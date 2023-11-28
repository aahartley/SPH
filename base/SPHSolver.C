#include "SPHSolver.h"

SPHSolver::SPHSolver(int num, float radius)
{
    this->num = num;
    dt = 0.016f;
    this->radius = radius;
    particles = std::vector<Particle>(num);
    viscosity = .01f;
    coef = 0.9f;
    GenerateParticles();

    //DFSPH
    for(int i = 0; i < num; i++)
    {
        CalcDensityField(i);
        CalcFactors(i);
    }
}

SPHSolver::~SPHSolver()
{

}

void SPHSolver::GenerateParticles()
{
    float xMin = 100.0f;
    float yMin = 100.0f;
    float xMax = 320.0f; // Adjust rectangle dimensions
    float yMax = 320.0f; // Adjust rectangle dimensions
    int numPerRow = 18;
    for (int i = 0; i < num; i++)
    {
        float xPos = xMin + (i % numPerRow) * (radius * 2 * SCALE); // Place particles in rows
        float yPos = yMin + (i / numPerRow) * (radius * 2 * SCALE); // Place particles in columns

        if (xPos <= xMax && xPos >= xMin && yPos <= yMax && yPos >= yMin)
        {
            Particle p(Vec2(0, 0), Vec2(0, 0), Vec2(xPos, yPos), radius, 1000, 1000);
            particles[i] = p;
        }
    }

}

//smoothing kernel
float SPHSolver::W( Vec2& distance, float h)
{
    //cubic spline from Smoothed Particle HydrodynamicsTechniques for the Physics Based Simulation of Fluids and Solids
    //Dan Koschier1, Jan Bender2, Barbara Solenthaler3, and Matthias Teschner4
    //normalization factors: sigma dimensions=2  40/7 *pi *h^2            d=3 8/pi*h^3
    //smmothing length == kernel support radius == h, h == 4 * particle radius, h = 0.3 m
    float q = (1.0f/h) * (distance.magnitude());
    float sigma_2 = 40.0f / (7.0f * PI * (h * h)); 
    if(q >= 0 && q <= 1.0f/2.0f) 
    {
        //std::cout << q << ' ' << sigma_2 * ((6 * (std::pow(q, 3) - std::pow(q, 2))) + 1) << '\n';
        return sigma_2 * ((6 * ((q * q * q) - (q * q))) + 1);
    }
    else if (q > 1.0f/2.0f && q <= 1) 
    {
        //std::cout << q << ' ' << sigma_2 * (2 * std::pow((1 - q), 3)) << '\n';
        return sigma_2 * ( 2 * ((1 - q) * (1 - q) * (1 - q)) );
    }
    else 
    {
        //std::cout << 0 << '\n';
        return 0;
    }

}

Vec2 SPHSolver::W_Gradient(Vec2& pos, Vec2& distance, float h)
{
    //∇W = (∂W/∂x, ∂W/∂y)
    Vec2 gradient;
    float q = (1.0f/h) * distance.magnitude();
    float sigma_2 = 40.0f / (7.0f * PI * (h * h)); 
    float dx, dy;
    if(q >= 0 && q <= 1.0f/2.0f) 
    {
        dx = (sigma_2 * ( (18 * (q*q)) - (12 * q) )  *  distance.X() ) / (h * distance.magnitude()) ;
        dy = (sigma_2 * ( (18 * (q*q)) - (12 * q) )  *  distance.Y() ) / (h * distance.magnitude()) ;
    }
    else if (q > 1.0f/2.0f && q <= 1)
    { 
        dx = (sigma_2 * (-6 * ((1 - q) * (1-q)) )  *  distance.X() ) / (h * distance.magnitude()) ;
        dy = (sigma_2 * (-6 * ((1 - q) * (1-q)) )  * distance.Y() ) / (h * distance.magnitude()) ;
    }
    else 
    {
        dx = 0;
        dy = 0;
    }
    gradient.X() = dx;
    gradient.Y() = dy;
    //std::cout << "gradient: " << gradient.X() << ' ' << gradient.Y() << '\n';
    //std::cout << gradient.magnitude() << '\n';
    return gradient;
}

void SPHSolver::CalcFactors(int i)
{

    Vec2 factor1;
    for(int j = 0; j < num; j++)
    {
        Vec2 dist = particles[i].GetPosN() - particles[j].GetPosN();
        if((&particles[i] != &particles[j]) && dist.magnitude() != 0)
        {
            //std::cout << (particles[j].GetMass() * W(dist, radius * 4))  << '\n';
            factor1 += (particles[j].GetMass() * W_Gradient(particles[i].GetPosN(), dist, radius * 4));
        }
    }
    float factor2 = 0;
    for(int j = 0; j < num; j++)
    {
        Vec2 dist = particles[i].GetPosN() - particles[j].GetPosN();
        if((&particles[i] != &particles[j]) && dist.magnitude() != 0)
        {
            //std::cout << (particles[j].GetMass() * W(dist, radius * 4))  << '\n';
            factor2 += ( (particles[j].GetMass() * W_Gradient(particles[i].GetPosN(), dist, radius * 4).magnitude())
                        * (particles[j].GetMass() * W_Gradient(particles[i].GetPosN(), dist, radius * 4).magnitude()) );
        }
    }
    //float factor = std::pow(particles[i].GetDens(),2) / ((factor1.magnitude() * factor1.magnitude()) + factor2);
    float factor = 1 / ((factor1.magnitude() * factor1.magnitude()) + factor2);
    if(factor <  1e-6) factor = 1e-6;
    particles[i].GetFactor() = factor;
}

void SPHSolver::CorrectDensityError()
{
    float densAvg = 0;
    float densRest = 1150;
    int iter = 0;
    float threshold = 100.f;
    #pragma omp parallel for
    for(int i = 0; i < num; i++)
    {
        densAvg += particles[i].GetDens();
    }
    densAvg /= num;
    //std::cout << densAvg << '\n';
    while((((densAvg - densRest) > threshold) || iter < 2)&& iter < 100)
    {
        //std::cout << iter << '\n';
        float densPAvg = 0;
        //predict dens
        #pragma omp parallel for
        for(int i = 0; i < num; i++)
        {
            CalcMaterialDens(i);
            particles[i].GetDensP() = particles[i].GetDens() + (dt * particles[i].GetDensM() );
            //if(particles[i].GetDensP() > densRest * 1.2) particles[i].GetDensP() = densRest * 1.2;
            //std::cout << particles[i].GetDensM() <<'\n';
            //std::cout <<  particles[i].GetDensP() << '\n';
        }
        #pragma omp parallel for
        for(int i = 0; i < num; i++)
        {
            Vec2 sum{0,0};
            float ki = ((particles[i].GetDensP() - densRest) / (dt * dt)) * particles[i].GetFactor();
            //std::cout << "dp-d0/dt^2: " << ((particles[i].GetDensP() - densRest) / std::pow(dt, 2))<<" factor: "<< particles[i].GetFactor() << " ki: "<<ki<< '\n';
            for(int j = 0; j < num; j++)
            {
                Vec2 dist = particles[i].GetPosN() - particles[j].GetPosN();
                if((&particles[i] != &particles[j]) && dist.magnitude() != 0)
                {
                    float kj = ((particles[j].GetDensP() - densRest) / (dt * dt)) * particles[j].GetFactor();
                    // sum += (particles[j].GetMass() * (( (ki / ((particles[i].GetDens())*(particles[i].GetDens()))) + (kj / ((particles[j].GetDens())*(particles[j].GetDens()))) ) *
                    //         W_Gradient(particles[i].GetPosN(), dist, radius*4)));    
                     sum += (particles[j].GetMass() * (( (ki / ((particles[i].GetDens()))) + (kj / ((particles[j].GetDens()))) ) *
                             W_Gradient(particles[i].GetPosN(), dist, radius*4)));
                    //std::cout << sum.X() << ' ' << sum.Y() << '\n';
                }
            }
            //std::cout << sum.X() << ' ' << sum.Y() << '\n';
            particles[i].GetVel() = particles[i].GetVel() -  ((dt * sum));
        }
        #pragma omp parallel for
        for(int i = 0; i < num; i++)
        {
            densPAvg += particles[i].GetDensP();
            particles[i].GetDensP() = 0;
            particles[i].GetDensM() = 0;
        }
        densPAvg /= num;
        densAvg = densPAvg;
        iter++;
      
    }
    //std::cout << "done\n";
}

void::SPHSolver::CorrectDivergenceError()
{
    float densMAvg = 0;
    int iter = 0;
    float threshold = 11.5f;
    for(int i = 0; i < num; i++)
    {
        densMAvg += particles[i].GetDensM();
    }
    densMAvg /= num;
    while((densMAvg > threshold || iter < 1) && iter < 100)
    {
        float dMAvg = 0;
        #pragma omp parallel for
        for(int i = 0; i < num; i++)
        {
            CalcMaterialDens(i);
        }
        #pragma omp parallel for
        for(int i = 0; i < num; i++)
        {
            Vec2 sum;
            float ki = (1.f/dt) * particles[i].GetDensM() * particles[i].GetFactor();
            for(int j = 0; j < num; j++)
            {
                Vec2 dist = particles[i].GetPosN() - particles[j].GetPosN();
                if((&particles[i] != &particles[j]) && dist.magnitude() != 0)
                {
                    float kj = (1.f/dt) * particles[j].GetDensM() * particles[j].GetFactor();
                    sum += ((particles[j].GetMass() * (( (ki / particles[i].GetDens())  + (kj / particles[j].GetDens()))) *
                            W_Gradient(particles[i].GetPosN(), dist, radius*4)));
                }
            }
            particles[i].GetVel() = particles[i].GetVel() -  ((dt * sum));
        }
        #pragma omp parallel for
        for(int i = 0; i < num; i++)
        {
            dMAvg += particles[i].GetDensM();
            particles[i].GetDensM() = 0;
        }
        dMAvg /= num;
        densMAvg = dMAvg;
        iter++;
    }
}

void SPHSolver::CalcMaterialDens(int i)
{
    for(int j = 0; j < num; j++)
    {
        Vec2 dist = particles[i].GetPosN() - particles[j].GetPosN();
        if((&particles[i] != &particles[j]) && dist.magnitude() != 0)
        {
            particles[i].GetDensM() += ((particles[j].GetMass() * (particles[i].GetVel() - particles[j].GetVel())).dot( W_Gradient(particles[i].GetPosN(), dist, radius * 4)));

        }
    }
    // for(int i = 0; i < num; i++)
    // {
    //     if(particles[i].GetDensM() < 0) particles[i].GetDensM() = 1;
    // }
  


}

void SPHSolver::CalcDensityField(int i)
{
    for(int j = 0; j < num; j++)
    {
        Vec2 dist = particles[i].GetPosN() - particles[j].GetPosN();
        if((&particles[i] != &particles[j]) && dist.magnitude() != 0)
        {
            particles[i].GetDens() += (particles[j].GetMass() * W(dist, radius * 4));
      

        }
    }
    for(int i = 0; i < num; i++)
    {
        //if(particles[i].GetDens() > 1500) particles[i].GetDens() = 1500;
    }
}

Vec2 SPHSolver::CalcPressureGradient(int i)
{
    Vec2 pres;

    for(int j = 0; j < num; j++)
    {
        Vec2 dist = particles[i].GetPosN() - particles[j].GetPosN();
        if((&particles[i] != &particles[j]) && dist.magnitude() != 0)
        {
            pres += (particles[j].GetMass() *
                ( (particles[i].GetPress() / (particles[i].GetDens() * particles[i].GetDens())) + (particles[j].GetPress() / (particles[j].GetDens() * particles[j].GetDens())) ) *
                W_Gradient(particles[i].GetPosN(), dist, radius * 4));
        }
    } 
    
    pres = pres * particles[i].GetDens();
    return pres;
}

Vec2 SPHSolver::LaplaceVel(int i)
{
    Vec2 vel;
    for(int j = 0; j < num; j++)
    {  
        Vec2 dist = particles[i].GetPosN() - particles[j].GetPosN();
        if((&particles[i] != &particles[j]) && dist.magnitude() != 0)
        {
            Vec2 gradient = W_Gradient(particles[i].GetPosN(), dist, radius * 4);
            //std::cout << gradient.magnitude() << '\n';
            vel += ( (particles[j].GetMass() / particles[j].GetDens()) *
                    (particles[i].GetVel() - particles[j].GetVel()) *
                    ( (2 * gradient.magnitude()) / dist.magnitude() ) );
        }

    }
    vel = vel * -1;
    //std::cout << vel.X() << ' ' << vel.Y() << '\n';
    return vel;
}



void SPHSolver::BoundaryCollisions()
{
    #pragma omp parallel for
    for(int i = 0; i < num; i++)
    {
        if(particles[i].GetPos().X() > 512-radius * SCALE)
        { 
            particles[i].GetPos().X() = (512-radius* SCALE) - 1;
            particles[i].GetPosN().X() = particles[i].GetPos().X() /SCALE;
            particles[i].GetVel().X() = particles[i].GetVel().X() * -coef;
        }
        if(particles[i].GetPos().X() < 0+radius* SCALE)
        {
            particles[i].GetPos().X() = (0+radius* SCALE) + 1;
            particles[i].GetPosN().X() = particles[i].GetPos().X() /SCALE;
            particles[i].GetVel().X() = particles[i].GetVel().X() * -coef;
        }
        if(particles[i].GetPos().Y() > 512-radius* SCALE)
        { 
            particles[i].GetPos().Y() = (512-radius* SCALE) - 1;
            particles[i].GetPosN().Y() = particles[i].GetPos().Y() /SCALE;
            particles[i].GetVel().Y() = particles[i].GetVel().Y() * -coef;
        }
        if(particles[i].GetPos().Y() < 0+radius* SCALE)
        {
            particles[i].GetPos().Y() = (0+radius* SCALE) + 1;
            particles[i].GetPosN().Y() = particles[i].GetPos().Y() /SCALE;
            particles[i].GetVel().Y() = particles[i].GetVel().Y() * -coef;
        }


    }
}

void SPHSolver::CFL()
{   
    float vel_m = particles[0].GetVel().magnitude();
    float lamda = 0.4f;
    #pragma omp parallel for
    for(int i = 1; i < num; i++)
    {
        float curr_vel = particles[i].GetVel().magnitude();
        if (curr_vel > vel_m) vel_m = curr_vel;
    }
    if(vel_m != 0) 
        dt = lamda * ((radius * 2.0f) / vel_m);
    if(dt > 0.016f || dt == 0) dt = 0.016f;
}

//DFSPH
void SPHSolver::PerformSimulation()
{
    CFL();
    //std::cout << dt << '\n';
    #pragma omp parallel for
    for(int i = 0; i < num; i++)
    {
        Vec2 f_ext(0, -9.8f);
        Vec2 f_viscosity = (particles[i].GetMass() * viscosity) * LaplaceVel(i);
        particles[i].GetVel() += ( (dt / particles[i].GetMass()) * (f_viscosity + f_ext + userF) );
        //std::cout << particles[i].GetVel().X() << ' ' << particles[i].GetVel().Y() << '\n';

    }

    //correct density error
    CorrectDensityError();
    #pragma omp parallel for
    for(int i = 0; i < num; i++)
    {
        particles[i].GetPosN() += (dt * particles[i].GetVel());
        particles[i].GetPos() = particles[i].GetPosN() * SCALE; //to draw
        //std::cout << particles[i].GetPos().X() << ' ' << particles[i].GetPos().Y() << '\n';

    }
    BoundaryCollisions();
    #pragma omp parallel for
    for(int i = 0; i < num; i++)
    {
        particles[i].GetDens() = 1000;
        particles[i].GetDensM() = 0;
        particles[i].GetDensP() = 0;
        particles[i].GetFactor() = 0;
        userF.X() = 0;

    }
    #pragma omp parallel for
    for(int i = 0; i < num; i++)
    {
        CalcDensityField(i);
        //std::cout << particles[i].GetDens() << '\n';
        CalcFactors(i);
    }
    CorrectDivergenceError();

}

















// //weakly compressible
// void SPHSolver::PerformSimulation()
// {
//     CFL();
//     //std::cout << dt << '\n';
//     for(int i = 0; i < num; i++)
//     {
//         CalcDensityField(i);
//         //std::cout << particles[i].GetDens() << '\n';
//     }

//     for(int i = 0; i < num; i++)
//     {
//         Vec2 f_ext(0, -9.8f);
//         Vec2 f_viscosity = (particles[i].GetMass() * viscosity) * LaplaceVel(i);
//         particles[i].GetVel() += ( (dt / particles[i].GetMass()) * (f_viscosity + f_ext) );
//         //std::cout << particles[i].GetVel().X() << ' ' << particles[i].GetVel().Y() << '\n';

//     }
//     for(int i = 0; i < num; i++)
//     {
//         particles[i].GetForces() = (-1 * (1 / particles[i].GetDens())) * CalcPressureGradient(i);
//         //std::cout << particles[i].GetForces().X() << ' ' << particles[i].GetForces().Y() << '\n';

//     }
//     for(int i = 0; i < num; i++)
//     {
//         particles[i].GetVel() += ((dt / particles[i].GetMass()) * particles[i].GetForces());
//         particles[i].GetPosN() += (dt * particles[i].GetVel());
//         particles[i].GetPos() = particles[i].GetPosN() * SCALE; //to draw
//         //std::cout << particles[i].GetPos().X() << ' ' << particles[i].GetPos().Y() << '\n';

//     }
//     BoundaryCollisions();
//     for(int i = 0; i < num; i++)
//     {
//         particles[i].GetDens() = 1000;
//         particles[i].GetForces().X() = 0;
//         particles[i].GetForces().Y() = 0;


//     }
// }
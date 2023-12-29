#include "SPHSolver.h"

SPHSolver::SPHSolver(int num, float radius)
{
    dt = 0.016f;
    this->radius = radius;
    this->num = num;
    numB = num + ((WIDTH/(radius * 2 * SCALE)) * 4);
    particles = std::vector<Particle>(numB);
    viscosity = .01f;
    coef = 0.9f;
    grid = GridSearch(WIDTH, radius * 4);
    GenerateParticles();
    GenerateBoundary();
    grid.InitGrid(particles);
    //DFSPH
    #pragma omp parallel for
    for(int i = 0; i < num; i++)
    {
        CalcDensityField(i);
        CalcFactors(i);

    }

}

SPHSolver::~SPHSolver()
{

}
void SPHSolver::CFL()
{   
    float vel_m = particles[0].GetVel().magnitude();
    float lamda = 0.4f;
    for(int i = 1; i < num; i++)
    {
        float curr_vel = particles[i].GetVel().magnitude();
        if (curr_vel > vel_m) vel_m = curr_vel;
    }
    if(vel_m != 0) 
        dt = lamda * ((radius * 2.0f) / vel_m);
    if(dt > 0.016f || dt == 0) dt = 0.016f;
}
void SPHSolver::GenerateParticles()
{
    float xMin = 100.0f;
    float yMin = 100.0f;
    float xMax = 320.0f; // Adjust rectangle dimensions
    float yMax = 520.0f; // Adjust rectangle dimensions
    int numPerRow = (xMax - xMin) / (float)(radius * 2 * SCALE);
    #pragma omp parallel for
    for (int i = 0; i < num; i++)
    {
        float xPos = xMin + (i % numPerRow) * (radius * 2 * SCALE); // Place particles in rows
        float yPos = yMin + (i / numPerRow) * (radius * 2 * SCALE); // Place particles in columns

        if (xPos <= xMax && xPos >= xMin && yPos <= yMax && yPos >= yMin)
        {
            Particle p(Vec2(0, 0), Vec2(xPos, yPos), radius, 1000, 1000, false);
            particles[i] = p;
        }
    }

}
//TODO
void SPHSolver::GenerateBoundary()
{
    int perEdge = float(WIDTH /(radius * 2 * SCALE));
    int index = num;
    for (int i = 0; i < perEdge; ++i)
    {
        //top, right, bot, left
        Particle p1(Vec2(0,0), Vec2(i * (radius * 2 * SCALE)+ radius*SCALE, WIDTH - (radius * SCALE)), radius, 1000, 1000, true);
        particles[index] = p1; index++;
        Particle p2(Vec2(0,0), Vec2(WIDTH-(radius * SCALE), i * (radius * 2 * SCALE)+ radius*SCALE), radius,1000,1000, true);
        particles[index] = p2; index++;
        Particle p3(Vec2(0,0), Vec2(i * (radius * 2 * SCALE) + radius*SCALE, (radius * SCALE)), radius,1000, 1000, true);
        particles[index] = p3; index++;
        Particle p4(Vec2(0,0), Vec2((radius * SCALE), i * (radius * 2 * SCALE)+ radius*SCALE), radius,1000,1000, true);
        particles[index] = p4; index++;

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

Vec2 SPHSolver::W_Gradient(Vec2& distance, float h)
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
    for(int j = 0; j < (int)grid.GetNeighbors()[i].size(); j++)
    {
        int index = grid.GetNeighbors()[i][j];
        Vec2 dist = particles[i].GetPosN() - particles[index].GetPosN();
        if((&particles[i] != &particles[index]) )
        {
            //std::cout << (particles[j].GetMass() * W(dist, radius * 4))  << '\n';
            factor1 += (particles[index].GetMass() * W_Gradient(dist, radius * 4));
        }
    }
    float factor2 = 0;
    for(int j = 0; j < (int)grid.GetNeighbors()[i].size(); j++)
    {
        int index = grid.GetNeighbors()[i][j];
        Vec2 dist = particles[i].GetPosN() - particles[index].GetPosN();
        if((&particles[i] != &particles[index]) )
        {
            //std::cout << (particles[j].GetMass() * W(dist, radius * 4))  << '\n';
            factor2 += ( (particles[index].GetMass() * W_Gradient(dist, radius * 4).magnitude())
                        * (particles[index].GetMass() * W_Gradient(dist, radius * 4).magnitude()) );
        }
    }
    float factor = 1 / ((factor1.magnitude() * factor1.magnitude()) + factor2);
    if(factor <  1e-6) factor = 1e-6;
    particles[i].GetFactor() = factor;
}
// void SPHSolver::CalcFactors(int i)
// {

//     Vec2 factor1;
//     for(int j = 0; j < num; j++)
//     {
//         Vec2 dist = particles[i].GetPosN() - particles[j].GetPosN();
//         if((&particles[i] != &particles[j]) )
//         {
//             //std::cout << (particles[j].GetMass() * W(dist, radius * 4))  << '\n';
//             factor1 += (particles[j].GetMass() * W_Gradient(dist, radius * 4));
//         }
//     }
//     float factor2 = 0;
//     for(int j = 0; j < num; j++)
//     {
//         Vec2 dist = particles[i].GetPosN() - particles[j].GetPosN();
//         if((&particles[i] != &particles[j]) )
//         {
//             //std::cout << (particles[j].GetMass() * W(dist, radius * 4))  << '\n';
//             factor2 += ( (particles[j].GetMass() * W_Gradient(dist, radius * 4).magnitude())
//                         * (particles[j].GetMass() * W_Gradient(dist, radius * 4).magnitude()) );
//         }
//     }
//     float factor = 1 / ((factor1.magnitude() * factor1.magnitude()) + factor2);
//     if(factor <  1e-6) factor = 1e-6;
//     particles[i].GetFactor() = factor;
// }
void SPHSolver::CorrectDensityError()
{
    float densAvg = 0;
    float densRest = 1200;
    int iter = 0;
    float threshold = 100.f;
    #pragma omp parallel for reduction (+:densAvg)
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
            for(int j = 0; j < (int)grid.GetNeighbors()[i].size(); j++)
            {
                int index = grid.GetNeighbors()[i][j];
                Vec2 dist = particles[i].GetPosN() - particles[index].GetPosN();
                if((&particles[i] != &particles[index]) )
                {
                    float kj = ((particles[index].GetDensP() - densRest) / (dt * dt)) * particles[index].GetFactor();  
                     sum += (particles[index].GetMass() * (( (ki / ((particles[i].GetDens()))) + (kj / ((particles[index].GetDens()))) ) *
                             W_Gradient(dist, radius*4)));
                    //std::cout << sum.X() << ' ' << sum.Y() << '\n';
                }
            }
            //std::cout << sum.X() << ' ' << sum.Y() << '\n';
            particles[i].GetVel() = particles[i].GetVel() -  ((dt * sum));

        }
        #pragma omp parallel for reduction (+:densPAvg)
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
// void SPHSolver::CorrectDensityError()
// {
//     float densAvg = 0;
//     float densRest = 1250;
//     int iter = 0;
//     float threshold = 100.f;
//     #pragma omp parallel for reduction (+:densAvg)
//     for(int i = 0; i < num; i++)
//     {
//         densAvg += particles[i].GetDens();
//     }
//     densAvg /= num;
//     //std::cout << densAvg << '\n';
//     while((((densAvg - densRest) > threshold) || iter < 2)&& iter < 100)
//     {
//         //std::cout << iter << '\n';
//         float densPAvg = 0;
//         //predict dens
//         #pragma omp parallel for
//         for(int i = 0; i < num; i++)
//         {
//             CalcMaterialDens(i);
//             particles[i].GetDensP() = particles[i].GetDens() + (dt * particles[i].GetDensM() );
//             //if(particles[i].GetDensP() > densRest * 1.2) particles[i].GetDensP() = densRest * 1.2;
//             //std::cout << particles[i].GetDensM() <<'\n';
//             //std::cout <<  particles[i].GetDensP() << '\n';
//         }
//         #pragma omp parallel for
//         for(int i = 0; i < num; i++)
//         {
//             Vec2 sum{0,0};
//             float ki = ((particles[i].GetDensP() - densRest) / (dt * dt)) * particles[i].GetFactor();
//             //std::cout << "dp-d0/dt^2: " << ((particles[i].GetDensP() - densRest) / std::pow(dt, 2))<<" factor: "<< particles[i].GetFactor() << " ki: "<<ki<< '\n';
//             for(int j = 0; j < num; j++)
//             {
//                 Vec2 dist = particles[i].GetPosN() - particles[j].GetPosN();
//                 if((&particles[i] != &particles[j]) )
//                 {
//                     float kj = ((particles[j].GetDensP() - densRest) / (dt * dt)) * particles[j].GetFactor();  
//                      sum += (particles[j].GetMass() * (( (ki / ((particles[i].GetDens()))) + (kj / ((particles[j].GetDens()))) ) *
//                              W_Gradient(dist, radius*4)));
//                     //std::cout << sum.X() << ' ' << sum.Y() << '\n';
//                 }
//             }
//             //std::cout << sum.X() << ' ' << sum.Y() << '\n';
//             particles[i].GetVel() = particles[i].GetVel() -  ((dt * sum));
//             //if(std::isnan(particles[i].GetVel().X()) || std::isnan(particles[i].GetVel().Y())) particles[i].GetVel() = Vec2(1,0);

//         }
//         #pragma omp parallel for reduction (+:densPAvg)
//         for(int i = 0; i < num; i++)
//         {
//             densPAvg += particles[i].GetDensP();
//             particles[i].GetDensP() = 0;
//             particles[i].GetDensM() = 0;
//         }
//         densPAvg /= num;
//         densAvg = densPAvg;
//         iter++;
      
//     }
//     //std::cout << "done\n";
// }
void::SPHSolver::CorrectDivergenceError()
{
    float densMAvg = 0;
    int iter = 0;
    float threshold = 12.0f;
    #pragma omp parallel for reduction (+:densMAvg)
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
            for(int j = 0; j < (int)grid.GetNeighbors()[i].size(); j++)
            {
                int index = grid.GetNeighbors()[i][j];
                Vec2 dist = particles[i].GetPosN() - particles[index].GetPosN();
                if((&particles[i] != &particles[index]))
                {
                    float kj = (1.f/dt) * particles[index].GetDensM() * particles[index].GetFactor();
                    sum += ((particles[index].GetMass() * (( (ki / particles[i].GetDens())  + (kj / particles[index].GetDens()))) *
                            W_Gradient(dist, radius*4)));
                }
            }
            particles[i].GetVel() = particles[i].GetVel() -  ((dt * sum));

        }
        #pragma omp parallel for reduction (+:dMAvg)
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
// void::SPHSolver::CorrectDivergenceError()
// {
//     float densMAvg = 0;
//     int iter = 0;
//     float threshold = 12.5f;
//     #pragma omp parallel for reduction (+:densMAvg)
//     for(int i = 0; i < num; i++)
//     {
//         densMAvg += particles[i].GetDensM();
//     }
//     densMAvg /= num;
//     while((densMAvg > threshold || iter < 1) && iter < 100)
//     {
//         float dMAvg = 0;
//         #pragma omp parallel for
//         for(int i = 0; i < num; i++)
//         {
//             CalcMaterialDens(i);
//         }
//         #pragma omp parallel for
//         for(int i = 0; i < num; i++)
//         {
//             Vec2 sum;
//             float ki = (1.f/dt) * particles[i].GetDensM() * particles[i].GetFactor();
//             for(int j = 0; j < num; j++)
//             {
//                 Vec2 dist = particles[i].GetPosN() - particles[j].GetPosN();
//                 if((&particles[i] != &particles[j]))
//                 {
//                     float kj = (1.f/dt) * particles[j].GetDensM() * particles[j].GetFactor();
//                     sum += ((particles[j].GetMass() * (( (ki / particles[i].GetDens())  + (kj / particles[j].GetDens()))) *
//                             W_Gradient(dist, radius*4)));
//                 }
//             }
//             particles[i].GetVel() = particles[i].GetVel() -  ((dt * sum));
//             //if(std::isnan(particles[i].GetVel().X()) || std::isnan(particles[i].GetVel().Y())) particles[i].GetVel() = Vec2(1,0);

//         }
//         #pragma omp parallel for reduction (+:dMAvg)
//         for(int i = 0; i < num; i++)
//         {
//             dMAvg += particles[i].GetDensM();
//             particles[i].GetDensM() = 0;
//         }
//         dMAvg /= num;
//         densMAvg = dMAvg;
//         iter++;
//     }
// }
void SPHSolver::CalcMaterialDens(int i)
{
    for(int j = 0; j < (int)grid.GetNeighbors()[i].size(); j++)
    {
        int index = grid.GetNeighbors()[i][j];
        Vec2 dist = particles[i].GetPosN() - particles[index].GetPosN();
        if((&particles[i] != &particles[index]))
        {
            particles[i].GetDensM() += ((particles[index].GetMass() * (particles[i].GetVel() - particles[index].GetVel())).dot(W_Gradient(dist, radius * 4)));

        }
    }


}
// void SPHSolver::CalcMaterialDens(int i)
// {
//     for(int j = 0; j < num; j++)
//     {
//         Vec2 dist = particles[i].GetPosN() - particles[j].GetPosN();
//         if((&particles[i] != &particles[j]))
//         {
//             particles[i].GetDensM() += ((particles[j].GetMass() * (particles[i].GetVel() - particles[j].GetVel())).dot( W_Gradient(dist, radius * 4)));
//             //if(particles[i].GetDensM() < 0) particles[i].GetDensM() = 0;
//         }
//     }

// }
void SPHSolver::CalcDensityField(int i)
{
    //std::cout << grid.GetNeighbors()[i].size() << '\n';
    for(int j = 0; j < (int)grid.GetNeighbors()[i].size(); j++)
    {
        int index = grid.GetNeighbors()[i][j];
        Vec2 dist = particles[i].GetPosN() - particles[index].GetPosN();
        if((&particles[i] != &particles[index]))
        {
            particles[i].GetDens() += (particles[index].GetMass() * W(dist, radius * 4));


        }
        //std::cout << particles[i].GetDens() << '\n';
    }
 
 }

// void SPHSolver::CalcDensityField(int i)
// {
//     for(int j = 0; j < num; j++)
//     {
//         Vec2 dist = particles[i].GetPosN() - particles[j].GetPosN();
//         if((&particles[i] != &particles[j]))
//         {
//             particles[i].GetDens() += (particles[j].GetMass() * W(dist, radius * 4));
      

//         }
//         //std::cout << particles[i].GetDens() << '\n';

//     }
 
// }

// for weakly compressible
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
                W_Gradient(dist, radius * 4));
        }
    } 
    
    pres = pres * particles[i].GetDens();
    return pres;
}

Vec2 SPHSolver::LaplaceVel(int i)
{
    Vec2 vel;
    for(int j = 0; j < (int)grid.GetNeighbors()[i].size(); j++)
    {
        int index = grid.GetNeighbors()[i][j];
        Vec2 dist = particles[i].GetPosN() - particles[index].GetPosN();
        if((&particles[i] != &particles[index]) && dist.magnitude() != 0)
        {
            Vec2 gradient = W_Gradient(dist, radius * 4);
            //std::cout << gradient.magnitude() << '\n';
            vel += ( (particles[index].GetMass() / particles[index].GetDens()) *
                    (particles[i].GetVel() - particles[index].GetVel()) *
                    ( (2 * gradient.magnitude()) / dist.magnitude() ) );
        }

    }
    vel = vel * -1;
    //if(std::isnan(vel.X()) || std::isnan(vel.Y())){ vel = Vec2(0,0);}
    //std::cout << vel.X() << ' ' << vel.Y() << '\n';
    return vel;
}

// Vec2 SPHSolver::LaplaceVel(int i)
// {
//     Vec2 vel;
//     for(int j = 0; j < num; j++)
//     {  
//         Vec2 dist = particles[i].GetPosN() - particles[j].GetPosN();
//         if((&particles[i] != &particles[j]) && dist.magnitude() != 0)
//         {
//             Vec2 gradient = W_Gradient(dist, radius * 4);
//             //std::cout << gradient.magnitude() << '\n';
//             vel += ( (particles[j].GetMass() / particles[j].GetDens()) *
//                     (particles[i].GetVel() - particles[j].GetVel()) *
//                     ( (2 * gradient.magnitude()) / dist.magnitude() ) );
//         }

//     }
//     vel = vel * -1;
//     //if(std::isnan(vel.X()) || std::isnan(vel.Y())) vel = Vec2(1,0);
//     //std::cout << vel.X() << ' ' << vel.Y() << '\n';
//     return vel;
// }



void SPHSolver::BoundaryCollisions()
{
    #pragma omp parallel for
    for(int i = 0; i < num; i++)
    {
        if(particles[i].GetPos().X() > WIDTH-radius * SCALE)
        { 
            particles[i].GetPos().X() = (WIDTH-radius* SCALE) - 1;
            particles[i].GetPosN().X() = particles[i].GetPos().X() /SCALE;
            particles[i].GetVel().X() = particles[i].GetVel().X() * -coef;
        }
        if(particles[i].GetPos().X() < 0+radius* SCALE)
        {
            particles[i].GetPos().X() = (0+radius* SCALE) + 1;
            particles[i].GetPosN().X() = particles[i].GetPos().X() /SCALE;
            particles[i].GetVel().X() = particles[i].GetVel().X() * -coef;
        }
        if(particles[i].GetPos().Y() > WIDTH-radius* SCALE)
        { 
            particles[i].GetPos().Y() = (WIDTH-radius* SCALE) - 1;
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
    CorrectDensityError();
    #pragma omp parallel for
    for(int i = 0; i < num; i++)
    {
        particles[i].GetPosN() += (dt * particles[i].GetVel());
        particles[i].GetPos() = particles[i].GetPosN() * SCALE; //to draw
        //std::cout << particles[i].GetPos().X() << ' ' << particles[i].GetPos().Y() << '\n';

    }
    //BoundaryCollisions();
    #pragma omp parallel for
    for(int i = 0; i < num; i++)
    {
        particles[i].GetDens() = 1000;
        particles[i].GetDensM() = 0;
        particles[i].GetDensP() = 0;
        particles[i].GetFactor() = 0;
        userF.X() = 0;

    }
    grid.InitGrid(particles);
    #pragma omp parallel for
    for(int i = 0; i < num; i++)
    {
        CalcDensityField(i);
        CalcFactors(i);
        //std::cout << particles[i].GetDens() << '\n';
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
//         particles[i].GetPForces() = (-1 * (1 / particles[i].GetDens())) * CalcPressureGradient(i);
//         //std::cout << particles[i].GetPForces().X() << ' ' << particles[i].GetForces().Y() << '\n';

//     }
//     for(int i = 0; i < num; i++)
//     {
//         particles[i].GetVel() += ((dt / particles[i].GetMass()) * particles[i].GetPForces());
//         particles[i].GetPosN() += (dt * particles[i].GetVel());
//         particles[i].GetPos() = particles[i].GetPosN() * SCALE; //to draw
//         //std::cout << particles[i].GetPos().X() << ' ' << particles[i].GetPos().Y() << '\n';

//     }
//     BoundaryCollisions();
//     for(int i = 0; i < num; i++)
//     {
//         particles[i].GetDens() = 1000;
//         particles[i].GetPForces().X() = 0;
//         particles[i].GetPForces().Y() = 0;


//     }
// }
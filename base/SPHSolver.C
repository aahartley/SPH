#include "SPHSolver.h"

SPHSolver::SPHSolver(int num, float radius)
{
    this->num = num;
    dt = 0.016;
    this->radius = radius;
    particles = std::vector<Particle>(num);
    GenerateParticles();
}

SPHSolver::~SPHSolver()
{

}

void SPHSolver::GenerateParticles()
{
    for(int i = 0; i < num; i++)
    {
        Particle p(Vec2(0,0), Vec2(0,0), Vec2(250 + i*(radius*2*80),250), radius, 1000, 0);
        particles[i] = p;
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
    float sigma_2 = (40.0f/7.0f) * 3.14f * std::pow(h, 2); 
    if(q >= 0 && q <= 1.0f/2.0f) 
    {
        //std::cout << q << ' ' << sigma_2 * ((6 * (std::pow(q, 3) - std::pow(q, 2))) + 1) << '\n';
        return sigma_2 * ((6 * (std::pow(q, 3) - std::pow(q, 2))) + 1);
    }
    else if (q > 1.0f/2.0f && q <= 1) 
    {
        //std::cout << q << ' ' << sigma_2 * (2 * std::pow((1 - q), 3)) << '\n';
        return sigma_2 * (2 * std::pow((1 - q), 3));
    }
    else 
    {
        //std::cout << 0 << '\n';
        return 0;
    }

}

float SPHSolver::W_Gradient(Vec2& pos, Vec2& distance, float h)
{
    //∇W = (∂W/∂x, ∂W/∂y)
    float gradient;
    float q = (1.0f/h) * distance.magnitude();
    float sigma_2 = (40.0f/7.0f) * 3.14f * std::pow(h, 2); 
    if(q >= 0 && q <= 1.0f/2.0f) 
    {
        float dx = (sigma_2 * ( (18 * std::pow(q, 2)) - (12 * q) ) ) * ( pos.X() / (h * distance.magnitude()) );
        float dy = (sigma_2 * ( (18 * std::pow(q, 2)) - (12 * q) ) ) * ( pos.Y() / (h * distance.magnitude()) );
        return dx + dy;
    }
    else if (q > 1.0f/2.0f && q <= 1)
    { 
        float dx = (sigma_2 * (-6 * std::pow((1 - q), 2) ) ) * ( pos.X() / (h * distance.magnitude()) );
        float dy = (sigma_2 * (-6 * std::pow((1 - q), 2) ) )  * ( pos.Y() / (h * distance.magnitude()) );
        return dx + dy;    
    }
    else 
        return 0;


    return gradient;
}

void SPHSolver::CalcDensityField(int i)
{
    for(int j = 0; j < num; j++)
    {
        if(&particles[i] != &particles[j])
        {
            Vec2 dist = particles[i].GetPosN() - particles[j].GetPosN();
            //std::cout << (particles[j].GetMass() * W(dist, radius * 4))  << '\n';
            particles[i].GetDens() += (particles[j].GetMass() * W(dist, radius * 4));
        }
    }
    
}

Vec2 SPHSolver::CalcGradient(float& A)
{
    //forward difference formula
    Vec2 grad;

    //for(int i = 0; i < num; i++)
    //{
        for(int j = 0; j < num; j++)
        {
            //A += particles[j].GetMass() * (particles[j].GetMass() / particles[j].GetDensity()) * W_Gradient();
        }
    //}

    return grad;
}

void SPHSolver::Integrate()
{
      for(int i = 0; i < num; i++)
    {
        particles[i].SetAcc(Vec2(0,-9.8f));
        particles[i].GetVel() += particles[i].GetAcc() * dt;
        particles[i].GetPos() += particles[i].GetVel() * dt;
    }

}

void SPHSolver::BoundaryCollisions()
{
    
}

void SPHSolver::PerformSimulation()
{
    for(int i = 0; i < num; i++)
    {
        CalcDensityField(i);
        //std::cout << particles[i].GetDens() << '\n';
    }

    Integrate();

}
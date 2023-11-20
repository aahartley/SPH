#ifndef SPHSOLVER_H
#define SPHSOLVER_H

#include <vector>
#include <iostream>

#include "Particle.h"

class SPHSolver
{
public:
    SPHSolver(int num, float radius);
    ~SPHSolver();

    void GenerateParticles();
    void Integrate();
    void BoundaryCollisions();
    //smoothing kernel
    float W(Vec2& distance, float h); //distance and smoothing radius
    float W_Gradient(Vec2& pos, Vec2& distance, float h);
    //calc rho field
    void CalcDensityField(int i);
    Vec2 CalcGradient(float& A);

    void PerformSimulation();

    int& GetNum(){return num;}
    std::vector<Particle> GetParticles(){return particles;}


private:
    int num;
    float dt, radius;
    std::vector<Particle> particles;

};



#endif
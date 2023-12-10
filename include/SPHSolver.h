#ifndef SPHSOLVER_H
#define SPHSOLVER_H

#include <vector>
#include <iostream>
#include <cmath>

#include "Particle.h"
#include "GridSearch.h"

class SPHSolver
{
public:
    SPHSolver(int num, float radius);
    ~SPHSolver();

    void GenerateParticles();
    void BoundaryCollisions();
    //smoothing kernel
    float W(Vec2& distance, float h); //distance and smoothing radius
    Vec2 W_Gradient(Vec2& distance, float h);
    //calc rho field
    void CalcDensityField(int i);
    Vec2 LaplaceVel(int index);
    Vec2 CalcPressureGradient(int i);
    void CFL();
    void CalcFactors(int i);
    void CalcMaterialDens(int i);
    void CorrectDensityError();
    void CorrectDivergenceError();

    void PerformSimulation();

    int& GetNum(){return num;}
    float& GetDT(){return dt;}
    Vec2& GetUserF(){return userF;}
    std::vector<Particle>& GetParticles(){return particles;}


private:
    int num;
    float dt, radius, viscosity, coef;
    std::vector<Particle> particles;
    Vec2 userF;
    GridSearch grid;

};



#endif
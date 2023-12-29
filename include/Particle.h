#ifndef PARTICLE_H 
#define PARTICLE_H

#include "Constants.h"
#include "Vec2.h"

class Particle
{
public:
    Particle(){}
    Particle(Vec2 v, Vec2 p, float r, float d, float pr, bool b):
             vel(v), pos(p), radius(r), rho(d), pressure(pr), boundary(b)
             {mass = 18.0f; posN.X()=pos.X()/SCALE; posN.Y()=pos.Y()/SCALE; factor=0; rho_m = 0; rho_p = 0;}

    ~Particle(){}

    Vec2& GetVel(){return vel;}
    Vec2& GetPos(){return pos;}
    Vec2& GetForces(){return forces;}
    Vec2& GetPosN(){return posN;} 
    float& GetRadius(){return radius;}
    float& GetDens(){return rho;}
    float& GetDensM() {return rho_m;}
    float& GetDensP() {return rho_p;}
    float& GetPress(){return pressure;}
    float& GetMass(){return mass;}
    float& GetFactor(){return factor;}
    bool& IsBoundary(){return boundary;}



private:
    Vec2 vel, pos, posN, forces;
    float radius, rho, rho_m, rho_p, pressure, mass, factor;
    bool boundary;
};





#endif
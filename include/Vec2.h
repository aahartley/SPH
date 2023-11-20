#ifndef VEC2_H
#define VEC2_H

#include <cmath>

class Vec2
{
public:
    Vec2(){xy[0] = xy[1] =0;}
    Vec2(float x, float y){xy[0] =x; xy[1]=y;}
    ~Vec2(){}

    float magnitude(){ return sqrtf(std::pow(xy[0],2) + std::pow(xy[1],2)); }

    //! Add two vectors together
    const Vec2 operator+(const Vec2& v) const 
    { 
       return Vec2(xy[0]+v.xy[0], xy[1]+v.xy[1]); 
    }
    //! Subtract two vectors together
    const Vec2 operator-(const Vec2& v) const 
    { 
       return Vec2(xy[0]-v.xy[0], xy[1]-v.xy[1]); 
    }
    Vec2& operator=       (const Vec2& v)
    { xy[0] = v.xy[0]; xy[1] = v.xy[1]; return *this; }
  
    Vec2& operator+=      (const Vec2& v)
    { xy[0] += v.xy[0]; xy[1] += v.xy[1]; return *this; }

    //! Multiplication of a constant with a vector
    friend const Vec2 operator* (const float w, const Vec2& v)
    { return v*w; }
	  
    //! Multiplication of a vector with a constant
    const Vec2 operator*        (const float v) const
    { return Vec2(xy[0]*v, xy[1]*v); }

    //! Division of a vector with a constant
    const Vec2 operator/        (const float v) const
    { return Vec2(xy[0]/v, xy[1]/v); }

    const float& operator[] (const int v) const { return xy[v]; }
    float& operator[] (const int v)       { return xy[v]; }
    const float& operator() (const int v) const { return xy[v]; }
    float& X()  { return xy[0]; }
    float& Y()  { return xy[1]; }

private:
    float xy[2];

};

#endif
#ifndef VEC_H
#define VEC_H

/*
 * Simple class for 2D/3D vectors
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
#include "dipole.h"
#include <iostream>

class Vec
{
    public:
        Vec();
        Vec(REAL x_, REAL y_);
        Vec(REAL x_, REAL y_, REAL z_);
        REAL GetX();
        REAL GetY();
        REAL GetZ();
        void SetX(REAL x_);
        void SetY(REAL y_);
        void SetZ(REAL z_);
        void operator+=(Vec& v);
        void operator-=(Vec& v);
        void operator=(Vec& v);
        Vec operator+(Vec& v);
        Vec operator-(Vec& v);
        void operator*=(REAL c);
        Vec operator*(REAL c);
        
        REAL Len();
        REAL LenSqr();
    
    private:
        REAL x,y,z;
};

std::ostream& operator<<(std::ostream& os, const Vec& ic);


#endif //VEC_H


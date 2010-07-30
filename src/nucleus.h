#ifndef NUCLEUS_H
#define NUCLEUS_H

/* 
 * Class to handle all nucleus related information, e.g. density
 * distribution etc.
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
#include "vector.h"
#include "dipole.h"
#include "gdist.h"
#include <iostream>
 
class Nucleus 
{
    public:
        Nucleus();
        Nucleus(int A_);
        int Intialize();
        int GetA();
        void SetA(int A_);
        REAL WS(REAL r);
        REAL GetWS_RA();
        REAL GetWS_delta();
        
        // Proton shape
        REAL Tp(REAL r);
        REAL Tp(Vec b);
        
        GDist* GetGDist();
        void SetGDist(GDist* gdist_);
        
        
    private:
        int A;      // Number of nucleons
        
        // Woods-Saxon parameters
        static const REAL WS_delta=0.54*FMGEV;
        REAL WS_N;
        REAL WS_RA;
        GDist* gdist;
        
        // Proton shape
        static const REAL B_G=4;

};

std::ostream& operator<<(std::ostream& os, Nucleus& ic);
 
#endif


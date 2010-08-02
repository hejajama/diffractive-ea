#ifndef NUCLEUS_H
#define NUCLEUS_H

/* 
 * Class to handle all nucleus related information, e.g. density
 * distribution etc.
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 

#include "dipole.h"
#include "gdist.h"
#include "vector.h"
#include <iostream>
#include <vector>
 
class Nucleus 
{
    public:
        Nucleus();
        Nucleus(int A_);
        ~Nucleus();
        int Intialize();
        int GetA();
        void SetA(int A_);
        
        // W-S distribution
        REAL WS(REAL r);
        REAL WS_unnorm(REAL r); // Unnormalized WS
        REAL T_WS(REAL b);
        REAL T_WS(Vec b);
        REAL T_WS_unnorm(REAL b); // Unnormalized T_WS
        
        // Nucleon configuration
        std::vector<Vec>& RandomNucleonConfiguration();
        //std::vector<Vec>& RandomNucleonConfiguration2d();
        
        int GenerateRandomNucleonConfigurations(int n,REAL mindist, REAL maxdist);
        
        REAL GetWS_RA();
        REAL GetWS_delta();
        
        // Proton shape
        REAL Tp(REAL r);
        REAL Tp(Vec b); // b = |b|
        
        GDist* GetGDist();
        void SetGDist(GDist* gdist_);
        
        REAL MaxR();    // Integrate up to this limit
        
        
    private:
        int A;      // Number of nucleons
        
        // Woods-Saxon parameters
        static const REAL WS_delta=0.54*FMGEV;
        REAL WS_RA;
        GDist* gdist;
    
        // Normalization
        REAL WS_N;
        REAL T_WS_N;
        REAL T_WS_0;    // T_WS(0)
        
        // Proton shape
        static const REAL B_G=4;
        
        // Generated random nucleon configurations
        std::vector< std::vector<Vec> > nucleon_configs;
        int next_config;  // id of the configuration to give out next


};

std::ostream& operator<<(std::ostream& os, Nucleus& ic);
 
#endif


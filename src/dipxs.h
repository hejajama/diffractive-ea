#ifndef DIPXS_H
#define DIPXS_H

/*
 * Parent class for dipole cross section
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
//#include "vector.h"
#include "nucleus.h"

class DipXS
{
    public:
        DipXS(Nucleus nucleus_);
        
        // Impact parameter dependent cross section averaged over 
        // initial nucleon configurations
        // virtual REAL DipXSection_b(REAL r, REAL xbj, Vec b ) = 0;
       // virtual REAL DipXSection_b(REAL r, REAL xbj, REAL b ) = 0;  // b = |b|
        REAL Alphas_r(REAL r);
        Nucleus& GetNucleus();
    protected:
        REAL Mu2(REAL rsqr);   // 4/r^2 + mu_0
        Nucleus nucleus;   
        REAL Sigmap(REAL r, REAL xbjork);  // Total dipole proton cross section

};

std::ostream& operator<<(std::ostream& os, DipXS& ic);


#endif

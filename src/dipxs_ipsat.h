#ifndef DIPXS_IPSAT_H
#define DIPXS_IPSAT_H

/*
 * Dipole cross section in IP sat model
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
#include "dipxs.h"
#include "nucleus.h"
#include <vector>
 
class DipXS_IPSat : public DipXS
{
    public:
        DipXS_IPSat(Nucleus nucleus_);
        REAL DipXSection_b_sqr(REAL rsqr, REAL r2sqr, Vec b, Vec b2, REAL xbjork ); // Impact parameter representation
        REAL DipXSection_b(REAL rsqr, REAL xbjork, REAL b ); // b=|b|
        
        REAL DipXSection_delta(REAL rsqr, REAL r2sqr, REAL xbjork, Vec delta);
        
        REAL DipXSection(REAL rsqr, REAL xbjork, Vec b, 
            std::vector<Vec> &nucleons); // Non-averaged dipole cross section

    
    private:
        static const REAL B_p=4.0; 

};


#endif  // DIPXS_IPSAT_H

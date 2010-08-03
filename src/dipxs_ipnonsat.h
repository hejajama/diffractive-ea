#ifndef DIPXS_IPNONSAT_H
#define DIPXS_IPNONSAT_H

/*
 * Dipole cross section in IP non sat model
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
#include "dipxs.h"
#include "nucleus.h"
#include <vector>
 
class DipXS_IPNonSat : public DipXS
{
    public:
        DipXS_IPNonSat(Nucleus nucleus_);
        REAL DipXSection_avg_sqr(REAL rsqr, REAL r2sqr, Vec b, Vec b2, REAL xbjork ); // Impact parameter representation
        REAL DipXSection(REAL rsqr, REAL xbjork, REAL b ); // b=|b|
        REAL DipXSection(REAL rsqr, REAL xbjork, Vec b, 
                std::vector<Vec>& nucleons); 

        
        REAL DipXSection_delta(REAL rsqr, REAL r2sqr, REAL xbjork, Vec delta);
        
        REAL DipXSection_b_nonaveraged(REAL rsqr, REAL xbjork, Vec b, 
            std::vector<Vec> &nucleons); // Non-averaged dipole cross section
    
    private:
        static const REAL B_p=4.0; 

};


#endif

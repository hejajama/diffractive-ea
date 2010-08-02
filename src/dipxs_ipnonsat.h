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
        REAL DipXSection_b_sqr(REAL r, REAL r2, Vec b, Vec b2, REAL xbjork ); // Impact parameter representation
        REAL DipXSection_b(REAL r, REAL xbjork, REAL b ); // b=|b|
        
        REAL DipXSection_delta(REAL r, REAL r2, REAL xbjork, Vec delta);
        
        REAL DipXSection_b_nonaveraged(REAL r, REAL xbjork, Vec b, 
            std::vector<Vec> &nucleons); // Non-averaged dipole cross section
    
    private:
        static const REAL B_p=4.0; 

};


#endif

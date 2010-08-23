#ifndef Dipxs_IPNONSAT_H
#define Dipxs_IPNONSAT_H

/*
 * Dipole cross section in IP non sat model
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
#include "dipxs.h"
#include "nucleus.h"
#include <vector>
 
class Dipxs_IPNonSat : public Dipxs
{
    public:
        Dipxs_IPNonSat(Nucleus &nucleus_);
        REAL Dipxsection_avg_sqr(REAL rsqr, REAL r2sqr, Vec b, Vec b2, REAL xbjork ); // Impact parameter representation
        REAL Dipxsection(REAL rsqr, REAL xbjork, REAL b ); // b=|b|
        REAL Dipxsection(REAL rsqr, REAL xbjork, Vec b, 
                std::vector<Vec>& nucleons); 

        
        REAL Dipxsection_delta(REAL rsqr, REAL r2sqr, REAL xbjork, Vec delta);
        
        REAL Dipxsection_b_nonaveraged(REAL rsqr, REAL xbjork, Vec b, 
            std::vector<Vec> &nucleons); // Non-averaged dipole cross section
    
    private:
        static const REAL B_p=4.0; 

};


#endif

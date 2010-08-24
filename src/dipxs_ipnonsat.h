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
        REAL Dipxsection_sqr_avg(REAL rsqr, REAL r2sqr, REAL xbjork, REAL delta);
        REAL Dipxsection(REAL rsqr, REAL xbjork, Vec b, 
                std::vector<Vec>& nucleons); 
       

    
    private:
        REAL Sigmap(REAL rsqr, REAL xbjork);
        static const REAL B_p=4.0; 

};


#endif

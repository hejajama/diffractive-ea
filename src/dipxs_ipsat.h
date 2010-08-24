#ifndef Dipxs_IPSAT_H
#define Dipxs_IPSAT_H

/*
 * Dipole cross section in IP sat model
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
#include "dipxs.h"
#include "nucleus.h"
#include <vector>
 
class Dipxs_IPSat : public Dipxs
{
    public:
        Dipxs_IPSat(Nucleus &nucleus_);
        // Amplitude squared averaged over nucleon configurations as a 
        // function of \Delta and r,r' (and x)
        // \int d^2 b_1 ... d^2 b_A T_A(b_1)...T_A(B_A) 
        //      * \int d^2 b d^2 b' e^(-i(b-b')*\Delta) 
        //      * (d\sigma^2 / d^2 b)(b,r) (d\sgima^2 / d^2b)(b',r')
        REAL Dipxsection_sqr_avg(REAL rsqr, REAL r2sqr, REAL xbj, 
                REAL delta);
        
        REAL Dipxsection(REAL rsqr, REAL xbjork, Vec b, 
            std::vector<Vec> &nucleons); // Non-averaged dipole cross section


        //REAL Dipxsection_b_avg_sqr(REAL rsqr, REAL r2sqr, Vec b, Vec b2, REAL xbjork ); // Impact parameter representation
        
        //REAL Dipxsection_delta(REAL rsqr, REAL r2sqr, REAL xbjork, Vec delta);
    
        REAL GetB_p();
        REAL FactorC(REAL rsqr, REAL xbjork);
    private:
        static const REAL B_p=4.0;       

};

static const int N_MAX=2;  // Upper limit for the sum in DipXSection_sqr_avg


#endif  // Dipxs_IPSAT_H

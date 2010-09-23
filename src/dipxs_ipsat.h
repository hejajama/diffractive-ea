#ifndef Dipxs_IPSAT_H
#define Dipxs_IPSAT_H

/*
 * Dipole cross section in IP sat model
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
#include "dipxs.h"
#include "nucleus.h"
#include "dipole.h"
#include <vector>
 
class Dipxs_IPSat : public Dipxs
{
    public:
        Dipxs_IPSat(Nucleus &nucleus_);
        Dipxs_IPSat(Nucleus &nucleus_, int mode_, REAL bp=DEFAULT_B_p);
        // Amplitude squared averaged over nucleon configurations as a 
        // function of \Delta and r,r' (and x)
        // \int d^2 b_1 ... d^2 b_A T_A(b_1)...T_A(B_A) 
        //      * \int d^2 b d^2 b' e^(-i(b-b')*\Delta) 
        //      * (d\sigma^2 / d^2 b)(b,r) (d\sgima^2 / d^2b)(b',r')
        REAL Dipxsection_sqr_avg(REAL rsqr, REAL r2sqr, REAL xbj, 
                REAL delta);
        
        REAL Dipxsection(REAL rsqr, REAL xbjork, Vec b, 
            std::vector<Vec> &nucleons); // Non-averaged dipole cross section

        // Amplitude squared for coherent scattering
        // |\int d^2 b_1...d^2 b_A T_A(b_1)...T_A(B_A)
        //      * \int d^2 b e^(-ib*\Delta)
        //      * (d\sigma^A/d^2 b)(b,r,x) |^2
        REAL CoherentDipxsection_avg(REAL rsqr, REAL xbj, REAL delta);

        // Dipole-proton amplitude
        REAL Dipxsection_proton(REAL rsqr, REAL xbj, REAL delta);
        REAL Dipxsection_proton(REAL rsqr, REAL xbj);


        //REAL Dipxsection_b_avg_sqr(REAL rsqr, REAL r2sqr, Vec b, Vec b2, REAL xbjork ); // Impact parameter representation
        
        //REAL Dipxsection_delta(REAL rsqr, REAL r2sqr, REAL xbjork, Vec delta);
    
        REAL GetB_p();
        REAL FactorC(REAL rsqr, REAL xbjork);
        void SetFactorize(bool f);
    private:
        REAL Sigmap(REAL rsqr, REAL xbjork);
        REAL B_p;   
        REAL mode;    
        bool factorize;

};

const int N_MAX=1;  // Upper limit for the sum in Dipxsection_sqr_avg

const int IPSAT_MODE_DEFAULT=1;
const int IPSAT_MODE_NONSAT_P=2;
/* In case of IPSAT_MODE_NONSAT_P the scattering amplitude for
 * dipole-proton scattering is not unitarized, but we however take into
 * account the possibility to scatter on multiple nucleons at the same time
 *
 * In IPSAT_MODE_DEFAULT the dipole-proton scattering amplitude is also
 * unitarized
 */



#endif  // Dipxs_IPSAT_H

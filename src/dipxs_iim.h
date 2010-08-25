#ifndef Dipxs_IIM_H
#define Dipxs_IIM_H

/*
 * Dipole cross section in heavy quark improved IIM model
 * taking into account the impact parameter dependence
 * Source: Article by Cyrillie Marquet
 * arXiv:0706.2682v1
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
#include "dipxs.h"
#include "nucleus.h"
#include <vector>
 
class Dipxs_IIM : public Dipxs
{
    public:
        Dipxs_IIM(Nucleus &nucleus_);
        // Amplitude squared averaged over nucleon configurations as a 
        // function of \Delta and r,r' (and x)
        // \int d^2 b_1 ... d^2 b_A T_A(b_1)...T_A(B_A) 
        //      * \int d^2 b d^2 b' e^(-i(b-b')*\Delta) 
        //      * (d\sigma^2 / d^2 b)(b,r) (d\sgima^2 / d^2b)(b',r')
        REAL Dipxsection_sqr_avg(REAL rsqr, REAL r2sqr, REAL xbj, 
                REAL delta);
        
        REAL Dipxsection(REAL rsqr, REAL xbjork, Vec b, 
            std::vector<Vec> &nucleons); // Non-averaged dipole cross section

    private:
        REAL Q_s(REAL x);       // Saturation scale
        REAL DipoleAmplitude(REAL r, REAL x);
        
    
        // Parameters for the model
        static const REAL alpha = 0.615065;
        static const REAL beta  = 1.00642;
        static const REAL x0    = 1.632e-5;
        static const REAL N0    = 0.7;
        static const REAL kappa = 9.9;
        static const REAL lambda= 0.2197;
        static const REAL gammac= 0.7376;
        static const REAL B_D   = 5.591;    // GeV^{-2}
        

};



#endif  // Dipxs_IPSAT_H

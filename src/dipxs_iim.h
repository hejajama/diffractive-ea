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

// How to generalize IIM model into qq-nucleus scattering
// IPNONSAT: amplitude(b) = \sum_i amplitude_qq(b-b_i),
// IPSAT: for S matrix element S^A(b) = \prod_i S(b-b_i)
const int IIM_IPSAT=1; const int IIM_IPNONSAT=2;
 
class Dipxs_IIM : public Dipxs
{
    public:
        Dipxs_IIM(Nucleus &nucleus_);
        Dipxs_IIM(Nucleus &nucleus_, std::string file);
        // Amplitude squared averaged over nucleon configurations as a 
        // function of \Delta and r,r' (and x)
        // \int d^2 b_1 ... d^2 b_A T_A(b_1)...T_A(B_A) 
        //      * \int d^2 b d^2 b' e^(-i(b-b')*\Delta) 
        //      * (d\sigma^2 / d^2 b)(b,r) (d\sgima^2 / d^2b)(b',r')
        REAL Dipxsection_sqr_avg(REAL rsqr, REAL r2sqr, REAL xbj, 
                REAL delta);
        REAL Dipxsection_proton(REAL rsqr, REAL xbj, REAL delta);
        REAL Dipxsection(REAL rsqr, REAL xbjork, Vec b, 
            std::vector<Vec> &nucleons); // Non-averaged dipole cross section
        REAL DipoleAmplitude(REAL r, REAL x);
        REAL Q_s(REAL x);       // Saturation scale
        REAL GetB_D();
    private:
        int ReadParameters(std::string file);
        REAL prevft, prevdelta;     // To optimize Dipxsection_sqr_avg
        
    
        // Parameters for the model
        REAL alpha;     // = 0.615065;
        REAL beta;      // = 1.00642;
        REAL x0;        // = 1.632e-5;
        REAL N0;        // = 0.7;
        REAL kappa;     // = 9.9;      // \eta in documentation
        REAL lambda;    // = 0.2197;
        REAL gammac;    // = 0.7376;
        REAL B_D;       // = 5.591 or 4.0;    // GeV^{-2}   
        static const int IIM_MODE=IIM_IPSAT;
};

static const int N_MAX_IIM=1;   // How many terms we take into account from
                                // the amplitude series

#endif  // Dipxs_IPSAT_H

#ifndef Dipxs_BK_H
#define Dipxs_BK_H

/*
 * Dipole cross section
 * Dipole amplitude from BK
 * B-dependence like in IIM
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2014
 */
 
#include "dipxs.h"
#include "nucleus.h"
#include <vector>
#include <gsl/gsl_integration.h>
#include <amplitudelib/amplitudelib.hpp>

 
class Dipxs_BK : public Dipxs
{
    public:
        Dipxs_BK(Nucleus &nucleus_);
        Dipxs_BK(Nucleus &nucleus_, std::string file);
        ~Dipxs_BK();
        // Amplitude squared averaged over nucleon configurations as a 
        // function of \Delta and r,r' (and x)
        // \int d^2 b_1 ... d^2 b_A T_A(b_1)...T_A(B_A) 
        //      * \int d^2 b d^2 b' e^(-i(b-b')*\Delta) 
        //      * (d\sigma^2 / d^2 b)(b,r) (d\sgima^2 / d^2b)(b',r')
        REAL DipoleAmplitude_sqr_avg(REAL rsqr, REAL r2sqr, REAL xbj, 
                REAL delta);
        REAL Dipxsection(REAL rsqr, REAL xbjork, Vec b, 
            std::vector<Vec> &nucleons); // Non-averaged dipole cross section
        REAL DipoleAmplitude(REAL r, REAL x);
        
        REAL DipoleAmplitude_proton(REAL rsqr, REAL xbj, REAL delta);
        
        // Scattering amplitude for coherent gamma^*A scattering
        REAL CoherentDipoleAmplitude_avg(REAL rsqr, REAL xbj, REAL delta);

        // Total dipole-proton cross section (integrated over d^2 b) in 1/Gev^2
        REAL TotalDipxsection_proton(REAL rsqr, REAL xbj);
        
        // 1/2*d\sigma/d^2b = q\barq-proton scattering amplitude
        REAL Qq_proton_amplitude(REAL rsqr, REAL xbj, REAL b);

        REAL Q_s(REAL x);       // Saturation scale
        REAL GetB_D();

        REAL Bp();
        REAL Sigma0();
    private:
        void Intialize();
        REAL prevft, prevdelta;     // To optimize Dipxsection_sqr_avg
        
        gsl_integration_workspace *ft_workspace_coh;
        
        AmplitudeLib* N;

        REAL B_D;       // = 5.591 or 4.0;    // GeV^{-2}
        REAL sigma0; 
};

static const int N_MAX_BK=1;   // How many terms we take into account from
                                // the amplitude series

#endif  // Dipxs_IPSAT_H

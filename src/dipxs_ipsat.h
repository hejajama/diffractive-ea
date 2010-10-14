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
#include <gsl/gsl_integration.h>
 
class Dipxs_IPSat : public Dipxs
{
    public:
        Dipxs_IPSat(Nucleus &nucleus_);
        Dipxs_IPSat(Nucleus &nucleus_, int mode_, REAL bp=DEFAULT_B_p);
        ~Dipxs_IPSat();
        // Amplitude squared averaged over nucleon configurations as a 
        // function of \Delta and r,r' (and x)
        // \int d^2 b_1 ... d^2 b_A T_A(b_1)...T_A(B_A) 
        //      * \int d^2 b d^2 b' e^(-i(b-b')*\Delta) 
        //      * (d\sigma^2 / d^2 b)(b,r) (d\sgima^2 / d^2b)(b',r')
        REAL DipoleAmplitude_sqr_avg(REAL rsqr, REAL r2sqr, REAL xbj, 
                REAL delta);

        // Averaged coherent scattering amplitude
        // \int d^2 b_1...d^2 b_A T_A(b_1)...T_A(B_A)
        //      * \int d^2 b e^(-ib*\Delta)
        //      * 1/2*(d\sigma^A / d^2 b) A(b,r,x) 
        REAL CoherentDipoleAmplitude_avg(REAL rsqr, REAL xbj, REAL delta);

        // Dipole-proton amplitude
        REAL DipoleAmplitude_proton(REAL rsqr, REAL xbj, REAL delta);
        
        // Total dipole-proton cross section (integrated over d^2 b) in 1/Gev^2
        REAL TotalDipxsection_proton(REAL rsqr, REAL xbj);
        
        // 1/2*d\sigma/d^2b = q\barq-proton scattering amplitude
        REAL Qq_proton_amplitude(REAL rsqr, REAL xbj, REAL b);
        
        

        //REAL Dipxsection_b_avg_sqr(REAL rsqr, REAL r2sqr, Vec b, Vec b2, REAL xbjork ); // Impact parameter representation
        
        //REAL Dipxsection_delta(REAL rsqr, REAL r2sqr, REAL xbjork, Vec delta);
    
        REAL GetB_p();
        REAL FactorC(REAL rsqr, REAL xbjork);
        void SetFactorize(bool f);
    private:
        void Intialize();
        REAL Sigmap(REAL rsqr, REAL xbjork);
        REAL B_p;   
        int mode;    
        bool factorize;
        gsl_integration_workspace *ft_workspace;
        gsl_integration_workspace *ft_workspace_proton;
        
        static const size_t MAXITER_FT=1000;  // Max number of interations when 
            // calculating the fourier transformation of the amplitude

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

#ifndef Dipxs_IPSAT2012_H
#define Dipxs_IPSAT2012_H

/*
 * Dipole class for diffractive calculations
 * Uses IPsat model from combined fit from 1212.2974
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>
 */
 
#include "dipxs.h"
#include "nucleus.h"
#include "dipole.h"
#include <vector>
#include "ipsat_mz/dipoleamplitude.hpp"
#include <gsl/gsl_integration.h>
 
class Dipxs_IPSat2012 : public Dipxs
{
    public:
        Dipxs_IPSat2012(Nucleus &n);
        ~Dipxs_IPSat2012();
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

        REAL Sigma0();
        REAL Bp();

       double DipoleAmplitude(double r, double xBj, double b=0, int param=2 ); 
		// Dipole amplitude where T_p is taken away,
		// C = 1 - exp(-Pi^2/(2*NC)*r^2*\alpha_s*xg/(2*Pi*Bp))
		// Stupid name, but use same as for the old ipsat
        REAL FactorC(REAL rsqr, REAL x);	 

        //REAL Dipxsection_b_avg_sqr(REAL rsqr, REAL r2sqr, Vec b, Vec b2, REAL xbjork ); // Impact parameter representation
        
        //REAL Dipxsection_delta(REAL rsqr, REAL r2sqr, REAL xbjork, Vec delta);
    
        REAL GetB_p();
		void SetB_p(double bp) { B_p = bp; }

		REAL SaturationScale(double x, int A, double b);

        void SetFactorize(bool f);
		bool GetFactorize();
    private:
        void Intialize();
        REAL B_p;   
        bool factorize;
		IPsat_MZ::DipoleAmplitude *ipsat_mz;
};

const int N_MAX2012=1;  // Upper limit for the sum in Dipxsection_sqr_avg

// Integration helper struct
// First the required integral helpers
struct inthelper_ipsatavg2012
{
    REAL rsqr,r2sqr;
    REAL xbj;
    REAL delta;
    
    Dipxs_IPSat2012* dip;
    Nucleus* nuke;
};

const int MAXITER_FT2012=1000;


#endif  // Dipxs_IPSAT_H

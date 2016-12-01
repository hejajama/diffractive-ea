#ifndef CALCULATOR_H
#define CALCULATOR_H

/*
 * Calculates different cross sections
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
#include "wave_function.h"
#include "dipxs.h"
#include "dipole.h"
#include <gsl/gsl_integration.h>

/*******************
 * \gamma^* N -> J/\Psi N cross section
 * Calculates:   d\sigma / dt = 1/(16*\pi)*|
 * \int d^2 r d^2 r' (jpsi)(r)*(jpsi)(r')*qqamplitude_sqr_avg(r,delta)|^2
 * And total cross section by integrating over t
 * Here (jpsi) is the inner product between \gamma^* and J/\Psi (or VM) wave 
 * functions integrated over z \in [0,1]
 *
 * Final result is multiplied by two correction factors, see arXiv:0712.260v2
 */

enum Diffraction
{
	COHERENT, INCOHERENT
};

class Calculator
{
    public:
        Calculator(Dipxs* amplitude_, WaveFunction* wavef_);
        REAL CrossSection_dt(REAL t, REAL Qsqr, REAL bjorkx);  // d\sigma / dt
        REAL CoherentCrossSection_dt(REAL t, REAL Qsqr, REAL bjorkx);
		REAL CoherentCrossSection_avgqsqr(REAL t, REAL minq2, REAL maxq2, REAL W, REAL M_v);
        REAL TotalCoherentCrossSection(REAL Qsqr, REAL bjorkx, 
            REAL mint, REAL maxt);
        REAL TotalCoherentCrossSection(REAL Qsqr, REAL bjorkx);
        REAL ProtonCrossSection_dt(REAL t, REAL Qsqr, REAL bjorkx);
        REAL ProtonCrossSection_dt_qsqravg(REAL dt, REAL minQsqr, REAL maxQsqr, REAL W, REAL M_v);
        
        REAL TotalCrossSection_noint(REAL Qsqr, REAL bjorkx);
        REAL TotalCrossSection(REAL Qsqr, REAL bjorkx); // incoherent
        REAL TotalCrossSection(REAL Qsqr, REAL bjorkx, REAL mint, REAL maxt); 
        REAL TotalProtonCrossSection(REAL Qsqr, REAL bjorkx);
        void SetPolarization(int pol);
        Dipxs* GetAmplitude();
        REAL RIntAmplitude(REAL t, REAL Qsqr, REAL bjorkx, 
            REAL(*helperf)(REAL x, void* p));
        void SetCorrections(bool c);
        void SetTAccuracy(REAL acc);
        
        // Flux of photons, ra is radius in GeV^(-1), z is charge 
        double NuclearPhotonFlux(double y, double sqrts, bool pa=false, int z=82); // if pa=true, we have pA collision -> min. impact param is R_A, not 2R_A
        double ProtonPhotonFlux(double y, double sqrts);
        double DiffractiveAAtoJpsi(double y, double sqrts, Diffraction d=COHERENT, bool pa=false, int z=82);	// d\\sigma/dy
        double DiffractiveAAtoJpsi_dt(double y, double sqrts, double t, Diffraction d=COHERENT, bool pa=false, int z=82);	// d\\sigma/dydt
        
        double CoherentIncoherent(double Qsqr, double bjorkx);	// find t where coherent=incoherent
        
    private:
        
        REAL Rg(REAL lambda);   // Skewedness effect
        REAL Beta(REAL lambda); // Ratio of the real part to im. part of
                                // the scattering amplitude
        REAL rgsqr_l;           // rgsqr and betasqr caches
        REAL rgsqr_t;           // t and l: transverse / longitudinal
        REAL betasqr_l;
        REAL betasqr_t;
        bool cached_corrections;    // true if we have cached corrections
        REAL cache_Q2;          // Value of Q^2 for which we cached corrections
        Dipxs* amplitude;
        WaveFunction* wavef;
        int polarization; // If polarization is set to VM_MODE_TOT, we must sum
                    // transversial and longitudinal polarization.
        bool corrections;       // Whether or not to use corrections
        REAL relaccuracy_t;     // Relative accuracy for integrations over t
};


// Integral helpers for outern r integral
struct inthelper_r
{
    Calculator* calculator;
    Dipxs* amplitude;
    WaveFunction* vm;  // Vector meson wave function
    REAL r;         // r for inner r' integral, not used in outern inthelperf_r1
    REAL Qsqr;
    REAL bjorkx;
    REAL delta;
    bool analytic_t;    // Perfrom t integral analytically
    bool coherent;  // Calculate coherent cross section, not quasi-elastic
    inthelper_r(){ analytic_t=false; coherent=false; }
}; 
 
// Integral helpers
REAL inthelperf_r1(REAL r, void* p);
REAL inthelperf_r2(REAL r, void* p);
REAL inthelper_totxs(REAL t, void* p);
REAL inthelper_totcohxs(REAL t, void* p);
REAL inthelperf_proton(REAL r, void* p);
REAL inthelperf_coherent(REAL r, void* p);

 



#endif


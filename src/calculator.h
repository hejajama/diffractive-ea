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

/*******************
 * \gamma^* N -> J/\Psi N cross section
 * Calculates:   d\sigma / dt = 1/(16*\pi)*
 * \int d^2 r d^2 r' (jpsi)(r)*(jpsi)(r')*qqamplitude_sqr_avg(r,delta)
 * Here (jpsi) is the inner product between \gamma^* and J/\Psi (or VM) wave 
 * functions integrated over z \in [0,1]
 */

class Calculator
{
    public:
        Calculator(Dipxs* amplitude_, WaveFunction* wavef_);
        REAL CrossSection_dt(REAL t, REAL Qsqr, REAL bjorkx);  // d\sigma / dt
        REAL ProtonCrossSection_dt(REAL t, REAL Qsqr, REAL Bjorkx);
        REAL TotalCrossSection(REAL Qsqr, REAL bjorkx); 
        REAL TotalProtonCrossSection(REAL Qsqr, REAL bjorkx);
        
        Dipxs* GetAmplitude();
        
    private:
        Dipxs* amplitude;
        WaveFunction* wavef;

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
    inthelper_r(){ analytic_t=false; }
}; 
 
// Integral helpers
REAL inthelperf_r1(REAL r, void* p);
REAL inthelperf_r2(REAL r, void* p);
REAL inthelper_totxs(REAL t, void* p);
REAL inthelperf_proton(REAL r, void* p);

// Integration settings
const REAL MAXR=4;
const REAL MINR=0.05;   // r=0 doesn't work, K_{0,1}(0)=inf
const REAL RINTACCURACY=0.0001;
const REAL TINTACCURACY=0.0001;
const REAL TOTXS_MAXT=2;  // Max |t| in GeV^2 when calculating total xs
 
 
#endif


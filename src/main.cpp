/*
 * Calculates eA cross section in dipole model
 * 
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
#include <iostream>
#include <ctime>
#include <vector>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_integration.h>

#include "dipole.h"
#include "vm_photon.h"
#include "vector.h"
#include "nucleus.h"
#include "dipxs_ipnonsat.h"
#include "dipxs_ipsat.h"
#include "gdist/gdist_dglap.h"
#include "mersenne/mersenne.h"

using namespace std;

// Integration settings
const REAL MAXR=4;
const REAL MINR=0.05;   // r=0 doesn't work, K_{0,1}(0)=inf
const REAL RINTACCURACY=0.0001;

// Integral helpers for outern r integral
struct inthelper_r
{
    Dipxs* amplitude;
    VM_Photon* vm;  // Vector meson wave function
    REAL r;         // r for inner r' integral, not used in outern inthelperf_r1
    REAL Qsqr;
    REAL bjorkx;
    REAL delta;
};

REAL inthelperf_r1(REAL r, void* p);

// For inner r' integral
REAL inthelperf_r2(REAL r2, void* p);

int main(int argc, char* argv[])
{    

    // Parameters
    REAL bjorkx=1e-4;
    REAL x=1e-4;
    REAL Qsqr=0;
    
    // Intialize random number generator
    seed_mersenne(std::time(NULL));
    
    // J/Psi wave function:  e_f, N_T, N_L, R_T, R_L, m_f, M_V, delta
    VM_Photon JPsi(2.0/3.0, 1.23, 0.83, sqrt(6.5), sqrt(3.0), 1.4, 3.097, 0);
    
    // Intialize Dipxs and Nucleus
    Nucleus nuke(197);
    //GDist *gdist = new DGLAPDist();
    GDist *gdist = new GDist_Toy();
    nuke.SetGDist(gdist);    
    //Dipxs* dsigmadb = new Dipxs_IPSat(nuke);
    Dipxs* dsigmadb = new Dipxs_IPNonSat(nuke);
    
    /******************
     * Dipole-nucleus cross section for a fixed r 
     */
    /*REAL rsqr=SQR(0.4);
    REAL maxt=2.0;
    int points=100;
    for (REAL delta=0; SQR(delta)<maxt; delta+=sqrt(maxt)/points)
    {
        REAL tmpxs=1.0/(16.0*M_PI)*dsigmadb->Dipxsection_sqr_avg(rsqr, rsqr, x, delta);
        cout << SQR(delta) << " " << tmpxs << endl;
    }*/
    
    /*******************
     * \gamma^* N -> J/\Psi N cross section
     * Calculates
     * \int d^2 r d^2 r' (jpsi)(r)*(jpsi)(r')*qqamplitude_sqr_avg(r,delta)
     * Here (jpsi) is the inner product between \gamma^* and J/\Psi wave 
     * functions integrated over z \in [0,1]
     */
    REAL maxt=0.3; REAL mint=0; int points=100;
    for (REAL delta=sqrt(mint); SQR(delta)<maxt; delta+=sqrt(maxt-mint)/points)
    {
        gsl_function fun;
        inthelper_r inthelp;
        inthelp.amplitude=dsigmadb;
        inthelp.vm=&JPsi; inthelp.bjorkx=bjorkx;
        inthelp.delta=delta; inthelp.Qsqr=Qsqr;
        fun.function=&inthelperf_r1;
        fun.params=&inthelp;
        
        REAL result,abserr; size_t eval;
        int status = gsl_integration_qng(&fun, MINR, MAXR, RINTACCURACY, RINTACCURACY, 
            &result, &abserr, &eval);
        result *= 1.0/(16.0*M_PI);
        cout << SQR(delta) << " " << result << endl;
    }


   
    return 0;
}


// Integral helpers

REAL inthelperf_r1(REAL r, void* p)
{
    inthelper_r* par = (inthelper_r*)p;
    par->r=r;
    gsl_function fun;
    fun.function=&inthelperf_r2;
    fun.params=par;
    
    REAL result,abserr; size_t eval;
    
    int status = gsl_integration_qng(&fun, MINR, MAXR, RINTACCURACY, RINTACCURACY, 
        &result, &abserr, &eval);
    
    return 2*M_PI*r*par->vm->PsiSqr_tot_intz(par->Qsqr, r)*result;
}

// Inner r' integral: \int dr' r' (jpsi)(r') * qqamplitude_sqr(r,r')
REAL inthelperf_r2(REAL r2, void* p)
{
    inthelper_r* par = (inthelper_r*)p;
    return 2*M_PI*r2 * par->vm->PsiSqr_tot_intz(par->Qsqr, r2) 
            * par->amplitude->Dipxsection_sqr_avg(SQR(par->r), SQR(r2), 
                    par->bjorkx, par->delta);

}

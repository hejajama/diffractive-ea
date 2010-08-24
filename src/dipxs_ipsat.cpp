/*
 * Dipole cross section in IP sat model
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
#include <gsl/gsl_math.h>
#include <iostream>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_gamma.h>
#include <vector>

#include "dipxs_ipsat.h"
#include "dipole.h"

const int AVERAGE_INTEGRAL_ITERATIONS=15;
const int MONTE_CARLO_ITERATIONS=1000;
const int MAX_ITERATIONS=80000;
const REAL MAXB=80;

const REAL AVGITACCURACY = 0.0001;

using std::cout; using std::endl; using std::cerr;

Dipxs_IPSat::Dipxs_IPSat(Nucleus &nucleus_) :
    Dipxs(nucleus_)
{
    // Nothing to do here
}

/* Amplitude squared averaged over nucleon configurations as a 
 * function of \Delta and r,r' (and x)
 * \int d^2 b_1 ... d^2 b_A T_A(b_1)...T_A(B_A) 
 *      * \int d^2 b d^2 b' e^(-i(b-b')*\Delta) 
 *      * (d\sigma^2 / d^2 b)(b,r) (d\sgima^2 / d^2b)(b',r')
 *
 * This quantity is calculated by using the following approximation which
 * works only for heavy nucleai if delta is big :
 * 16\pi B_p \int d^2 b \sum_{n=1}^A 1/n*exp(-B_p*\Delta^2/n)
 *    * exp(-2*A*\pi*B_p*T_A(b)*(C(r)+C(r'))) 
      * [ Pi*Bp*C(r)*C(r')*T_A(b)
 *    / (1 - 2*Pi*Bp*T_A(b)(C(r)-C(r'))) ]^n
 *
 * This sum converges rapidly, we can even assume that n=1
 * Upper limit is defined as N_MAX
 */
 
// First the required integral helpers
struct inthelper_ipsatavg
{
    REAL rsqr,r2sqr;
    REAL xbj;
    REAL delta;
    
    Dipxs_IPSat* dip;
    Nucleus* nuke;
};

REAL inthelperf_ipsatavg(REAL b, void* p)
{
    // Calculate the sum
    REAL sum=0;
    inthelper_ipsatavg* par=(inthelper_ipsatavg*)p;
    REAL A = par->nuke->GetA();
    REAL x = par->xbj;
    REAL B_p=par->dip->GetB_p();
    REAL twstmp = par->nuke->T_WS(b);
    for (int i=1; i<=N_MAX; i++)
    {
        sum+=gsl_sf_choose(A,i)/((REAL)i)*exp(-B_p*SQR(par->delta)/i)
             * exp(-2*A*M_PI*B_p*twstmp*(par->dip->FactorC(par->rsqr,x)
               + par->dip->FactorC(par->r2sqr,x) ) )
             * pow( 
                M_PI*B_p*par->dip->FactorC(par->rsqr,x)
                * par->dip->FactorC(par->r2sqr,x)*twstmp
                / (1-2*M_PI*B_p*twstmp
                *(par->dip->FactorC(par->rsqr,x)+par->dip->FactorC(par->r2sqr,x)
                ))  ,i);    
    }
    
    sum*=16*M_PI*B_p;
    
    // 2D integral -> 1D
    sum*=2*M_PI*b;
    return sum;
}

REAL Dipxs_IPSat::Dipxsection_sqr_avg(REAL rsqr, REAL r2sqr, REAL xbj, 
                REAL delta)
{

    inthelper_ipsatavg helper;
    helper.rsqr=rsqr; helper.r2sqr=r2sqr; helper.xbj=xbj; helper.delta=delta;
    helper.nuke=&nucleus; helper.dip=this;
    
    gsl_function int_helper;
    int_helper.function=&inthelperf_ipsatavg;
    int_helper.params=&helper;
    
    size_t eval;
    REAL result,abserr;

    int status = gsl_integration_qng(&int_helper, 0, MAXB, 
            AVGITACCURACY, AVGITACCURACY, &result, &abserr, &eval);
    
    return result;
}

/*
 * Non-averaged dipole cross section as a function of
 * nucleon transversial positions 
 */
REAL Dipxs_IPSat::Dipxsection(REAL rsqr, REAL xbjork, Vec b, 
            std::vector<Vec> &nucleons)
{
    if (nucleons.size() != nucleus.GetA())
    {
        std::cerr << "Dipxs_IPSat::Dipxsection: Got list of " 
            << nucleons.size() << " nucleons but "
            << "there should be " << nucleus.GetA() << " of them..." << std::endl;
            return -1;
    }
    REAL result=0; REAL ex=0;
    Vec tmp(0,0,0);
    for (int i=0; i<nucleons.size(); i++)
    {
        tmp.SetX(b.GetX()-nucleons[i].GetX());
        tmp.SetY(b.GetY()-nucleons[i].GetY());
        ex+=nucleus.Tp(tmp);        
    }
    result = 2.0*(1.0-exp(-2*nucleus.GetGDist()->Gluedist(xbjork,rsqr)/2.0*ex));
    return result;
}

/*
 * FactorC
 * Defined to simplify equations
 * C = 1 - exp(-Pi^2/NC*r^2*\alpha_s*xg/(4*Pi*Bp))
 *   = 1 - exp(-sigmap / (4*Pi*Bp) )
 *   = 1 - exp(-Gluedist()*r^2 / (2*Pi*Bp))
 */
REAL Dipxs_IPSat::FactorC(REAL rsqr, REAL xbjork)
{
    return 1-exp(-nucleus.GetGDist()->Gluedist(xbjork,rsqr)*rsqr/(2.0*M_PI*B_p));    
}


REAL Dipxs_IPSat::GetB_p()
{
    return B_p;
}


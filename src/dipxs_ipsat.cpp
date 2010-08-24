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
#include "cubature/cubature.h" // Multi dimensional integral

const int AVERAGE_INTEGRAL_ITERATIONS=15;
const int MONTE_CARLO_ITERATIONS=1000;
const int MAX_ITERATIONS=80000;
const REAL MAXB=80;
const int N_MAX=1;  // Upper limit for the sum in DipXSectionsqr_avg
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
 *    * exp(-2*A*\pi*B_p*T_A(b)*(C(r)+C(r'))) * Pi*Bp*C(r)*C(r')*T_A(b)
 *    / (1 - 2*Pi*Bp*T_A(b)(C(r)-C(r')))
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
             * M_PI*B_p*par->dip->FactorC(par->rsqr,x)
             * par->dip->FactorC(par->r2sqr,x)*twstmp
             / (1-2*M_PI*B_p*twstmp
              *(par->dip->FactorC(par->rsqr,x)+par->dip->FactorC(par->r2sqr,x)));    
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


//*******************************************
/* The following methods are not used here and these may not work,
 * I'm leaving them here for now as they might be useful when one 
 * tries to write a code for a dipole cross section model with which it is not
 * possible to average over nucleon configs analytically.
 *
 */
 
#ifdef DONT_HIDE_STUPID_CODE


struct inthelper_param_sat
{ 
    Dipxs_IPSat* Dipxs;
    REAL r1,r2,xbjork; 
};

REAL inthelperf_sat(REAL x[], size_t dim, void* p)
{
    inthelper_param_sat* fp = (inthelper_param_sat*)p;
    if (dim!=4) { cerr << "Dimension is not 4 " << endl; return 0; }
    
    Vec b1(x[0],x[1]); Vec b2(x[2],x[3]);
    REAL result = fp->Dipxs->Dipxsection_b_avg_sqr(fp->r1,fp->r2,b1,b2,fp->xbjork);
    return result;

}


void inthelperf_cuba_sat(unsigned int ndim, const double x[], void* fdata, 
    unsigned int fdim, REAL *fval)
{
    inthelper_param_sat* fp = (inthelper_param_sat*)fdata;
    if (ndim!=4) { cerr << "Dimension is not 4 " << endl; return; }
    
    Vec b1(x[0],x[1]); Vec b2(x[2],x[3]);
   fval[0] = fp->Dipxs->Dipxsection_b_avg_sqr(fp->r1,fp->r2,b1,b2,fp->xbjork);
}

struct inthelper_param2_sat
{
    Dipxs_IPSat* Dipxs;
    REAL r1sqr,r2sqr,xbjork;
    Vec* delta;
    std::vector<Vec>* nucleons;
};

void inthelperf2_cuba_sat(unsigned int ndim, const double x[], void* fdata, 
    unsigned int fdim, REAL *fval)
{
    // To calculate \int d^2 b d^2 b' dsgima/db(r,b_i) dsigma/db(r',b_i')
    inthelper_param2_sat* fp = (inthelper_param2_sat*)fdata;
    Vec b1(x[0],x[1]); Vec b2(x[2],x[3]); Vec& delta = *(fp->delta);
    fval[0]=fp->Dipxs->Dipxsection(fp->r1sqr, fp->xbjork, b1, *(fp->nucleons))
           *fp->Dipxs->Dipxsection(fp->r2sqr, fp->xbjork, b2, *(fp->nucleons));
    fval[0]*=cos(x[0]*delta.GetX())*cos(x[1]*delta.GetY())
            *cos(x[2]*delta.GetX())*cos(x[3]*delta.GetY());
    }

/* 
 * Dipole nucleus cross section as a function of delta
 * NOTE: Does not work if delta != 0!
 */
REAL Dipxs_IPSat::Dipxsection_delta(REAL r1sqr, REAL r2sqr, REAL xbjork, Vec delta)
{
    REAL xl[4]={-MAXB,-MAXB,-MAXB,-MAXB};
    REAL xu[4]={MAXB,MAXB,MAXB,MAXB};
    REAL sum,errorsum;
    
    // This is easy to parallerize as each iteration is independent
    #pragma omp parallel for
    for (int i=0; i<AVERAGE_INTEGRAL_ITERATIONS; i++)
    {
        REAL result,error;
        inthelper_param2_sat params;
        params.r1sqr=r1sqr; params.r2sqr=r2sqr; params.xbjork=xbjork;
        params.nucleons=&nucleus.RandomNucleonConfiguration();
        params.Dipxs=this;
        params.delta=&delta;
       
        adapt_integrate(1, inthelperf2_cuba_sat, &params, 4, xl, xu, MAX_ITERATIONS,
            1e-5,1e-5,&result,&error);
        cout << "Average " << i << " / " << AVERAGE_INTEGRAL_ITERATIONS << " gave " << result << endl;
        sum+=result; errorsum+=error;
    }
        
    cout << "Absolute error: " << errorsum/(AVERAGE_INTEGRAL_ITERATIONS*16.0*M_PI) << endl;
    return sum/(16.0*M_PI*AVERAGE_INTEGRAL_ITERATIONS);
    // Test by setting \delta=0 and calculating dipole-proton cross section
    /*
    Vec delta(0,0);
    
    gsl_monte_function F;
    inthelper_param params;
    params.r1=r1; params.r2=r2; params.xbjork=xbjork;
    params.Dipxs=this;
    
    F.f=&inthelperf; F.dim=4; F.params=&params;
    
    // Calculate by Monte Carlo (stupid?)
    size_t calls = MONTE_CARLO_ITERATIONS;
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    
    gsl_monte_plain_state *s  = gsl_monte_plain_alloc(4);
    REAL result,error;

    REAL xl[4]={-MAXB,-MAXB,-MAXB,-MAXB};
    REAL xu[4]={MAXB,MAXB,MAXB,MAXB};
    gsl_monte_plain_integrate(&F, xl, xu, 4, calls, r, s, &result, &error);
       
    */

    //return result/(16.0*M_PI);
}

/*
 * Averaged dipole nucleus cross section squared
 */
REAL Dipxs_IPSat::Dipxsection_b_avg_sqr(REAL rsqr, REAL r2sqr, Vec b, Vec b2, REAL xbjork )
{
    if (nucleus.GetA()==1) // Dipole-proton
        return Sigmap(rsqr,xbjork)*nucleus.Tp(b.Len());

    /* Dipole-nuclues, so we have to average over all possible nucleon
     * configurations. Monte Carlo method is used here as it 
     * should converge quite rapidly, all nuclei should be quite similar
     * when nucleon positions are generated from
     * Woods-Saxon distribution
     * Monte Carlo means here that we generate some random nucleai and average
     * the product of amplitudes over these nucleai. 
     * 
     * Calculates:
     * \int d^2 b_1 ... d^2 b_A T_A(b_1)...T_A(b_A) 
     *  * (d\sigma/d^2 b)(r,b_1,...,b_A) * (d\sigma/d^2 b')(r',b_1,...,b_A)
     */
     
    int A = nucleus.GetA(); // Dimension of the integration
    
    // Note: this would be easy to parallelize

    REAL result=0;
    for (int i=0; i<AVERAGE_INTEGRAL_ITERATIONS; i++)
    {
                // TODO: Define mindist and maxdist somewhere 
        result+=Dipxsection(rsqr,xbjork,b,
                nucleus.RandomNucleonConfiguration())
                * Dipxsection(r2sqr,xbjork,b2, 
                nucleus.RandomNucleonConfiguration());         
     }
     result*=1.0/AVERAGE_INTEGRAL_ITERATIONS;
     return result;
}

#endif


/*
 * Dipole cross section in IP non sat model
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
#include <gsl/gsl_math.h>
#include <iostream>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_exp.h>
#include <vector>

#include "dipxs_ipnonsat.h"
#include "dipole.h"
#include "cubature/cubature.h" // Multi dimensional integral

const int AVERAGE_INTEGRAL_ITERATIONS=130;
const int MONTE_CARLO_ITERATIONS=5000;
const REAL MAXB=30;

using std::cout; using std::endl; using std::cerr;

DipXS_IPNonSat::DipXS_IPNonSat(Nucleus nucleus_) :
    DipXS(nucleus_)
{
    // Nothing to do here
}
/*
struct inthelper
{
    DipXS_IPNonSat* dipxs;
    REAL b;
}

REAL inthelperf_b2(REAL r, void* params)
{
    // Inner b2 integral, 1 is given as a parameter
    REAL r = 0.3;   // Constant for testing purposes
    return ((*inthelper)params)->dipxs->DipXSection_b_sqr(r,r, 
    
}

REAL inthelperf_b1(REAL b, void* params)
{
    // Outern b1 integral
    gsl_integration_workspace* w2 = gsl_integration_workspace_alloc(1000);
    REAL result, abserr; 
    gsl_function F2;
    F2.function=&inthelperf_b2;

    // Give b1 as a parameter to inner integration
    inthelper params_2;
    params_2.dipxs=(inthelper*)params->dipxs;
    params_2.b=b;
    gsl_integration_qags(&F2, 0, MAXB, 1e-5,1e-7, 1000, w2, &result, &abserr);
    
    return result;
    
    gsl_integration_workspace_free(w2);

}
*/

struct inthelper_param 
{ 
    DipXS_IPNonSat* dipxs;
    REAL r1,r2,xbjork; 
};

REAL inthelperf(REAL x[], size_t dim, void* p)
{
    inthelper_param* fp = (inthelper_param*)p;
    if (dim!=4) { cerr << "Dimension is not 4 " << endl; return 0; }
    
    Vec b1(x[0],x[1]); Vec b2(x[2],x[3]);
    REAL result = fp->dipxs->DipXSection_b_sqr(fp->r1,fp->r2,b1,b2,fp->xbjork);
    return result;

}

void inthelperf_cuba(unsigned int ndim, const double x[], void* fdata, 
    unsigned int fdim, REAL *fval)
{
    inthelper_param* fp = (inthelper_param*)fdata;
    if (ndim!=4) { cerr << "Dimension is not 4 " << endl; return; }
    
    Vec b1(x[0],x[1]); Vec b2(x[2],x[3]);
   fval[0] = fp->dipxs->DipXSection_b_sqr(fp->r1,fp->r2,b1,b2,fp->xbjork);
}

/* 
 * Dipole nucleus cross section as a function of delta
 */
REAL DipXS_IPNonSat::DipXSection_delta(REAL r1, REAL r2, REAL xbjork, Vec delta)
{
    REAL xl[4]={-MAXB,-MAXB,-MAXB,-MAXB};
    REAL xu[4]={MAXB,MAXB,MAXB,MAXB};
    REAL result,error;
    
    inthelper_param params;
    params.r1=r1; params.r2=r2; params.xbjork=xbjork;
    params.dipxs=this;
    const int MAX_ITERATIONS=10000;
    adapt_integrate(1, inthelperf_cuba, &params, 4, xl, xu, MAX_ITERATIONS,
        1e-5,1e-5,&result,&error);
        
    cout << "Absolute error: " << error/(16.0*M_PI) << endl;

    // Test by setting \delta=0 and calculating dipole-proton cross section
    /*
    Vec delta(0,0);
    
    gsl_monte_function F;
    inthelper_param params;
    params.r1=r1; params.r2=r2; params.xbjork=xbjork;
    params.dipxs=this;
    
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

    return result/(16.0*M_PI);
}

/*
 * Dipole cross section in impact parameter reprsesentation
 * Averaged over all initial configurations
 */
REAL DipXS_IPNonSat::DipXSection_b(REAL r, REAL xbjork, REAL b )
{
    // Dipol-proton only for now
    return Sigmap(r,xbjork)*nucleus.Tp(b);
}

/*
 * Averaged dipole nucleus cross section squared
 */
REAL DipXS_IPNonSat::DipXSection_b_sqr(REAL r, REAL r2, Vec b, Vec b2, REAL xbjork )
{
    if (nucleus.GetA()==1) // Dipole-proton
        return Sigmap(r,xbjork)*nucleus.Tp(b.Len());

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

    REAL result=0;
    for (int i=0; i<AVERAGE_INTEGRAL_ITERATIONS; i++)
    {
                // TODO: Define mindist and maxdist somewhere 
        result+=DipXSection_b_nonaveraged(r,xbjork,b,
                nucleus.RandomNucleonConfiguration())
                * DipXSection_b_nonaveraged(r2,xbjork,b2, 
                nucleus.RandomNucleonConfiguration());         
     }
     result*=1.0/AVERAGE_INTEGRAL_ITERATIONS;
     return result;
}


/*
 * Non-averaged dipole cross section as a function of
 * nucleon transversial positions 
 */
REAL DipXS_IPNonSat::DipXSection_b_nonaveraged(REAL r, REAL xbjork, Vec b, 
            std::vector<Vec> &nucleons)
{
    if (nucleons.size() != nucleus.GetA())
    {
        std::cerr << "Got list of " << nucleons.size() << " nucleons but "
            << "there should be " << nucleus.GetA() << " of them...";
            return -1;
    }
    REAL result=0;
    for (int i=0; i<nucleons.size(); i++)
    {
        Vec tmp(b.GetX()-nucleons[i].GetX(),b.GetY()-nucleons[i].GetY());
//        cout << "tmp.LenSqr(): " << 
        result += exp(-tmp.LenSqr()/(2*B_p));
    }
    result *= Sigmap(r,xbjork)/(2*M_PI*B_p);
    return result;
}


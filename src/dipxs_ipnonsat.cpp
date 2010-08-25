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
#include <gsl/gsl_sf_bessel.h>
#include <vector>



#include "dipole.h"
#include "dipxs_ipnonsat.h"

const REAL MAXB=70;
const REAL FTINTACCURACY=0.0001;
const int INTPOINTS=1000;

using std::cout; using std::endl; using std::cerr;

Dipxs_IPNonSat::Dipxs_IPNonSat(Nucleus &nucleus_) :
    Dipxs(nucleus_)
{
    prevr=prevr2=prevft=-1;
}


/* Amplitude squared averaged over nucleon configurations as a 
 * function of \Delta and r,r' (and x)
 * \int d^2 b_1 ... d^2 b_A T_A(b_1)...T_A(B_A) 
 *      * \int d^2 b d^2 b' e^(-i(b-b')*\Delta) 
 *      * (d\sigma^2 / d^2 b)(b,r) (d\sgima^2 / d^2b)(b',r')
 *
 * The following result is derived by Caldwell&Kowalski, see arXiv:0909.1254v1:
 * |A|^2 = sigmap(r)*sigmap(r')*exp(-B_p*\Delta^2)*[ A + A(A-1)
 *            * |\int d^2 b e^(-ib*\delta) T_A(b)|^2
 * 
 * Here we neglect two-body correlations
 */
 
// We have to calculate the FT of T_A, so we need some integral helpers
struct inthelper_ipnonsat
{
    Nucleus* nuke;
    REAL delta;
};

REAL inthelperf_ipnonsat(REAL b, void* p)
{
    //TODO: Move this to Nucleus class?
    inthelper_ipnonsat* par = (inthelper_ipnonsat*)p;
    return 2*M_PI*b*gsl_sf_bessel_J0(b*par->delta)*par->nuke->T_WS(b);
}

REAL Dipxs_IPNonSat::Dipxsection_sqr_avg(REAL rsqr, REAL r2sqr, 
    REAL xbjork, REAL delta)
{
    size_t eval;
    REAL result,abserr;
        
    inthelper_ipnonsat helper;
    helper.delta=delta;
    helper.nuke=&nucleus;
    
    gsl_function int_helper;
    int_helper.function=&inthelperf_ipnonsat;
    int_helper.params=&helper;

    // If we just calculated this when this function was called previously
    // Yeah, this is quite an ugly hack, but this optimizes this quite much!
    if (rsqr==prevr and r2sqr==prevr2)
        result=prevft;
    else
        int status = gsl_integration_qng(&int_helper, 0, MAXB, 
                FTINTACCURACY, FTINTACCURACY, &result, &abserr, &eval);
    
    /*gsl_integration_workspace * w 
         = gsl_integration_workspace_alloc (1000);
    gsl_integration_qag(&int_helper, 0, MAXB, FTINTACCURACY, FTINTACCURACY,
        INTPOINTS, GSL_INTEG_GAUSS51, w, &result, &abserr); 
    gsl_integration_workspace_free(w);
    */
    REAL A = nucleus.GetA();
    
    //std::cerr << "Sigmap: " << Sigmap(rsqr,xbjork) << " , int: " << result << ", delta=" << delta << " err:" << abserr << endl;
    
    REAL bp=B_p;
    prevr=rsqr; prevr2=r2sqr; prevft=result;
    
    return Sigmap(rsqr,xbjork)*Sigmap(r2sqr,xbjork)*exp(-bp*SQR(delta))
        *A*(1.0+(A-1.0)*result*result);
}


/*
 * Dipole cross section in impact parameter reprsesentation as a 
 * function of nucleon transversal positions
 * _NOT_ Averaged over nucleon connfigurations
 */
REAL Dipxs_IPNonSat::Dipxsection(REAL rsqr, REAL xbjork, Vec b, 
                std::vector<Vec>& nucleons)
{
    REAL x=b.GetX(); REAL y=b.GetY();
    if (nucleons.size() != nucleus.GetA())
    {
        std::cerr << "Got list of " << nucleons.size() << " nucleons but "
            << "there should be " << nucleus.GetA() << " of them..." << endl;
            return -1;
    }
    REAL result=0;
    Vec tmp;
    for (int i=0; i<nucleons.size(); i++)
    {
        tmp.SetX(x-nucleons[i].GetX());
        tmp.SetY(y-nucleons[i].GetY());
        result += exp(-tmp.LenSqr()/(2*B_p));
    }
    result *= Sigmap(rsqr,xbjork)/(2*M_PI*B_p);
    return result;
} 

/*
 * Sigmap
 * Total dipole-proton cross section
 * \pi^2/Nc*r^2*alpha_s(mu(r)^2)*xg(x,mu(r)^2)
 * = 2*r^2*Gluedist(x,r)
 */
REAL Dipxs_IPNonSat::Sigmap(REAL rsqr, REAL xbjork)
{
    return 2*rsqr*nucleus.GetGDist()->Gluedist(xbjork,rsqr);
}



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
#include <gsl/gsl_sf_bessel.h>
#include <vector>

#include "dipxs_ipsat.h"
#include "dipole.h"

// Integral helpers
REAL inthelperf_satp(REAL b, void* p);  // Integrate over impact parameter 

const REAL MAXB=90;

const REAL AVGITACCURACY = 0.001;
const REAL FTINTACCURACY = 0.001; 

using std::cout; using std::endl; using std::cerr;

Dipxs_IPSat::Dipxs_IPSat(Nucleus &nucleus_) :
    Dipxs(nucleus_)
{
    mode = IPSAT_MODE_DEFAULT;
    B_p=DEFAULT_B_p;
    factorize=true;
    Intialize();
}

Dipxs_IPSat::Dipxs_IPSat(Nucleus &nucleus_, int mode_, REAL bp) : Dipxs(nucleus_)
{
    mode = mode_;
    B_p=bp;
    factorize=true;
    Intialize();
}

void Dipxs_IPSat::Intialize()
{
    factorize=true,
    ft_workspace = gsl_integration_workspace_alloc(MAXITER_FT);
    ft_workspace_proton = gsl_integration_workspace_alloc(MAXITER_FT);
}

Dipxs_IPSat::~Dipxs_IPSat()
{
    gsl_integration_workspace_free(ft_workspace);
    gsl_integration_workspace_free(ft_workspace_proton);
}

/* Amplitude squared averaged over nucleon configurations as a 
 * function of \Delta and r,r' (and x)
 * \int d^2 b_1 ... d^2 b_A T_A(b_1)...T_A(B_A) 
 *      * \int d^2 b d^2 b' e^(-i(b-b')*\Delta) 
 *      * (d\sigma^2 / d^2 b)(b,r) (d\sgima^2 / d^2b)(b',r')
 *
 * This quantity is calculated by using the following approximation which
 * works only for heavy nucleai if delta is large :
 * <|d\sigma/d^2 b|^2> = 
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
    int A = par->nuke->GetA();
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
                / (1.0-2.0*M_PI*B_p*twstmp
                *(par->dip->FactorC(par->rsqr,x)+par->dip->FactorC(par->r2sqr,x)
                ))  ,i);    
    }
    
    sum*=4*M_PI*B_p;   // Factor 4, not 16, as we return the amplitude, not
                        // d\sigma/d^2b
    
    // 2D integral -> 1D
    sum*=2*M_PI*b;
    return sum;
}

REAL Dipxs_IPSat::DipoleAmplitude_sqr_avg(REAL rsqr, REAL r2sqr, REAL xbj, 
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
            0, AVGITACCURACY, &result, &abserr, &eval);
    
    if (status) std::cerr << "Error " << status << " at " << __FILE__ << ":"
        << __LINE__ << ": Result " << result << ", abserror: " << abserr 
        << " (t=" << delta*delta <<")" << endl;
    return result;
}


/* 
 * Amplitude for coherent dipole-nucleus scattering
 * \int d^2 b_1...d^2 b_A T_A(b_1)...T_A(B_A)
 *      *\int d^2 b e^(-ib*\Delta)
 *      * 1/2*(d\sigma^A/d^2 b)(b,r,x)
 *
 * In IPSat model this can be derived to be
 * (when assuming smooth and heavy nucleus)
 *
 * A = \int d^2 e^(-i b*\Delta) (1 - exp(-A/2*T_A(b)*\sigma_dip^p(r,x)) )
 *
 * The nasty part here is that we have to perform the fourier transformation
 * numerically. 
 */

REAL inthelperf_ipsat_coherentavg(REAL b, void* p)
{
    inthelper_ipsatavg* par=(inthelper_ipsatavg*)p;
    int A = par->nuke->GetA();
    return 2*M_PI*b*gsl_sf_bessel_J0(b*par->delta)*
        (    1-exp(-A/2.0*par->nuke->T_WS(b)  
        * par->dip->TotalDipxsection_proton(par->rsqr, par->xbj) )  );
}

REAL Dipxs_IPSat::CoherentDipoleAmplitude_avg(REAL rsqr, REAL xbj, REAL delta)
{
    inthelper_ipsatavg helper;
    helper.rsqr=rsqr; helper.xbj=xbj; helper.delta=delta;
    helper.nuke=&nucleus; helper.dip=this;
    
    gsl_function int_helper;
    int_helper.function=&inthelperf_ipsat_coherentavg;
    int_helper.params=&helper;
    
    size_t eval; 
    REAL result,abserr;

    //int status = gsl_integration_qng(&int_helper, 0, MAXB, 
    //        0, OSCAVGITACCURACY, &result, &abserr, &eval);
    
    int status=gsl_integration_qag(&int_helper, 0, MAXB, 0, FTINTACCURACY, 
        MAXITER_FT, GSL_INTEG_GAUSS41, ft_workspace, &result, &abserr);

    if (status) std::cerr << "Error " << status << " at " << __FILE__ << ":"
        << __LINE__ << ": Result " << result << ", abserror: " << abserr 
        << " relerror: " << abserr/result << " (t=" << delta*delta <<")" << endl;
    return result;

}

/*
 * Total dipole-proton cross section
 * Calculated as \int d^2 b d\sigma/d^2 b
 * = 2 * 2*\pi*B_p*(1-exp(-r^2*F(x,r)))
 * = 4\pi B_p * (1-exp(-r^2 * Gluedist() / (2*\pi*B_p) ))
 *
 * Units: 1/GeV^2
 */
REAL Dipxs_IPSat::TotalDipxsection_proton(REAL rsqr, REAL xbj)
{
    if (!factorize)
    {
        std::cerr << "Nonfactorized Dipxs_IPSat::TotalDipxsecton_proton " 
            << " is not implemented, using factorized one. " << std::endl;
    }   
    return 4.0*M_PI*B_p*FactorC(rsqr, xbj);
}

/*
 * Dipole-proton amplitude as a function of \Delta
 * Integrated over impact parameter dependence
 * A = 2*Pi*B_p*C exp(-B_p \Delta^2 / 2)
 * 
 * Otherwise same as Dipxs_IPNonSAT::Dipxsection_proton but
 * we take into account the unitarity requirement.
 *
 * It is also approximated that 1-exp(-r^2T(b)...) = T(b)(1-exp(...))
 */
REAL Dipxs_IPSat::DipoleAmplitude_proton(REAL rsqr, REAL xbj, REAL delta)
{
    // Note: We cannot use FactorC() here as FactorC() depends on the mode
    // used: in NONSAT_P mode it doesn't saturate
    
    if (factorize) // Factorize T(b) dependency
    {
        /* If we want to allways use saturated dipole-proton amplitude
        return 2.0*M_PI*B_p*exp(-B_p*SQR(delta)/2.0)
            * ( 1.0-exp(-nucleus.GetGDist()->Gluedist(xbj,rsqr)
                *rsqr/(2.0*M_PI*B_p)) );
                */
        return 2.0*M_PI*B_p*exp(-B_p*SQR(delta)/2.0)*FactorC(rsqr, xbj);
    }
    // Else: Integrate numerically over impact parameter dependence, oscillatory
    // integral...
    inthelper_ipsatavg helper;
    helper.rsqr=rsqr; helper.xbj=xbj; helper.delta=delta;
    helper.nuke=&nucleus; helper.dip=this;
    
    gsl_function int_helper;
    int_helper.function=&inthelperf_satp;
    int_helper.params=&helper;
    
    size_t eval;
    REAL result,abserr;
    //int status = gsl_integration_qng(&int_helper, 0, MAXB, 
    //        0, FTINTACCURACY, &result, &abserr, &eval);
    int status=gsl_integration_qag(&int_helper, 0, MAXB, 0, FTINTACCURACY, 
        MAXITER_FT, GSL_INTEG_GAUSS41, ft_workspace_proton, &result, &abserr);
    
    if (status) std::cerr << "Error " << status << " at " << __FILE__ << ":"
        << __LINE__ << ": Result " << result << ", abserror: " << abserr <<
        ", t=" << SQR(delta) << " GeV^2 " << std::endl;
        
    return result;
}



REAL inthelperf_satp(REAL b, void* p)
{
    inthelper_ipsatavg* par = (inthelper_ipsatavg*)p;
    return 2.0*M_PI*b*gsl_sf_bessel_J0(b*par->delta)*(1-exp(-par->rsqr
         * par->nuke->GetGDist()->Gluedist(par->xbj, par->rsqr)
         * exp(-SQR(b)/(2*par->dip->GetB_p()))/(2.0*M_PI*par->dip->GetB_p()) ));
}


/*
 * Non-averaged dipole cross section as a function of
 * nucleon transversial positions 
 */
/*REAL Dipxs_IPSat::Dipxsection(REAL rsqr, REAL xbjork, Vec b, 
            std::vector<Vec> &nucleons)
{
    std::cerr << "Dipxs_IPSat::Dipxsection() may not work, test before using!"
        << std::endl;
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
    result = 1.0*(1.0-exp(-2*nucleus.GetGDist()->Gluedist(xbjork,rsqr)/2.0*ex));
    return result;
}*/

/*
 * Dipole-proton scattering amplitude
 */
REAL Dipxs_IPSat::Qq_proton_amplitude(REAL rsqr, REAL xbj, REAL b)
{
    return 1.0-exp(-rsqr*nucleus.GetGDist()->Gluedist(xbj, rsqr)
     *exp(-SQR(b)/(2.0*B_p))/(2*M_PI*B_p));
}

/*
 * FactorC
 * Defined to simplify equations
 * C = 1 - exp(-Pi^2/(2*NC)*r^2*\alpha_s*xg/(2*Pi*Bp))
 *   = 1 - exp(-F(x,r^2)*r^2 )
 *   = 1 - exp(-sigmap/2 * / (2*Pi*Bp) )
 *   = 1 - exp(-Gluedist()*r^2 / (2Pi*Bp))
 *
 * In IPSAT_MODE_NONSAT_P case we expand this to the first order in r^2
 */
REAL Dipxs_IPSat::FactorC(REAL rsqr, REAL xbjork)
{
    if (mode==IPSAT_MODE_DEFAULT)
        return 1.0-exp(-nucleus.GetGDist()->Gluedist(xbjork,rsqr)
            *rsqr/(2.0*M_PI*B_p));    
    else if (mode==IPSAT_MODE_NONSAT_P) 
        return nucleus.GetGDist()->Gluedist(xbjork,rsqr)*rsqr / (2.0*M_PI*B_p);
    else
        std::cerr << "Error: mode not set for Dipxs_IPSat. " << std::endl;
}

/*
 * Sigmap
 * \pi^2/Nc*r^2*alpha_s(mu(r)^2)*xg(x,mu(r)^2)
 * = 2*r^2*Gluedist(x,r)
 */
REAL Dipxs_IPSat::Sigmap(REAL rsqr, REAL xbjork)
{
    return 2*rsqr*nucleus.GetGDist()->Gluedist(xbjork,rsqr);
}


REAL Dipxs_IPSat::GetB_p()
{
    return B_p;
}

void Dipxs_IPSat::SetFactorize(bool f)
{  
    factorize=f;
}


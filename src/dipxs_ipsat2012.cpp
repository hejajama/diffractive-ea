/*
 * Dipole class for diffractive calculations
 * Uses IPsat model from combined fit from 1212.2974
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>
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

#include "dipxs_ipsat2012.h"
#include "dipole.h"

// Integral helpers
REAL inthelperf_satp2012(REAL b, void* p);  // Integrate over impact parameter 
REAL inthelperf_totxs_satp2012(REAL b, void* p);

const REAL MAXB=100;
const int MAXITER_BINT=1000;

const REAL AVGINTACCURACY = 0.001;
const REAL FTINTACCURACY = 0.0001; 

using std::cout; using std::endl; using std::cerr;

// Actual routine to compute dipole amplitude
extern "C" {
  double dipole_amplitude_(double* xBj, double* r, double* b, int* param);
};
double DipoleAmplitude(double r, double xBj, double b=0, int param=2)
{
		// param 2: m_c=1.4
	return 0.5*dipole_amplitude_(&xBj, &r, &b, &param);
}


Dipxs_IPSat2012::Dipxs_IPSat2012(Nucleus &n) : Dipxs(n)
{
    B_p=4;
    Intialize();
}

void Dipxs_IPSat2012::Intialize()
{
    factorize=false;
}

Dipxs_IPSat2012::~Dipxs_IPSat2012()
{

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
 *
 * NB! Uses factorized dipole amplitude even if mode is set to 
 * IPSAT_MODE_NOFACTOR
 */
 
REAL inthelperf_ipsatavg2012(REAL b, void* p)
{
    // Calculate the sum
    REAL sum=0;
    inthelper_ipsatavg2012* par=(inthelper_ipsatavg2012*)p;
    int A = par->nuke->GetA();
    REAL x = par->xbj;
    REAL r = std::sqrt(par->rsqr);
    REAL r2 = std::sqrt(par->r2sqr);
    REAL B_p=par->dip->GetB_p();
    REAL twstmp = par->nuke->T_WS(b);
    for (int i=1; i<=N_MAX2012; i++)
    {
        sum+=gsl_sf_choose(A,i)/((REAL)i)*exp(-B_p*SQR(par->delta)/i)
             * exp(-2.0*A*M_PI*B_p*twstmp*(DipoleAmplitude(r,x)
               + DipoleAmplitude(r2, x) ) )
             * pow( 
                M_PI*B_p*DipoleAmplitude(r, x)
                * DipoleAmplitude(r2, x)*twstmp
                / (1.0-2.0*M_PI*B_p*twstmp
                *(DipoleAmplitude(r, x)+DipoleAmplitude(r2, x)
                ))  ,i);    
    }
    
    sum*=4.0*M_PI*B_p;   // Factor 4, not 16, as we return 1/2 d\sigma/d^2b
    
    // 2D integral -> 1D
    sum*=2.0*M_PI*b;
    return sum;
}

REAL Dipxs_IPSat2012::DipoleAmplitude_sqr_avg(REAL rsqr, REAL r2sqr, REAL xbj, 
                REAL delta)
{

    inthelper_ipsatavg2012 helper;
    helper.rsqr=rsqr; helper.r2sqr=r2sqr; helper.xbj=xbj; helper.delta=delta;
    helper.nuke=&nucleus; helper.dip=this;
    
    gsl_function int_helper;
    int_helper.function=&inthelperf_ipsatavg2012;
    int_helper.params=&helper;
    
    size_t eval;
    REAL result,abserr;

    //int status = gsl_integration_qng(&int_helper, 0, MAXB, 
    //        0, AVGINTACCURACY, &result, &abserr, &eval);
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(MAXITER_BINT);
     // maxb: MAXB
    int status=gsl_integration_qag(&int_helper, 0, MAXB, 0, AVGINTACCURACY, 
        MAXITER_BINT, GSL_INTEG_GAUSS61, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    
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
 *
 * NB! Does not take into account whether the mode is default or NONSAT_P
 */

REAL inthelperf_ipsat_coherentavg2012(REAL b, void* p)
{
		cerr << "TODO: b dependence for proton correctly!" << endl;
    inthelper_ipsatavg2012* par=(inthelper_ipsatavg2012*)p;
    int A = par->nuke->GetA();
    return 2.0*M_PI*b*gsl_sf_bessel_J0(b*par->delta)*
        (    1.0-exp(-A/2.0*par->nuke->T_WS(b)  
        * par->dip->TotalDipxsection_proton(par->rsqr, par->xbj) )  );
}

REAL Dipxs_IPSat2012::CoherentDipoleAmplitude_avg(REAL rsqr, REAL xbj, REAL delta)
{
    inthelper_ipsatavg2012 helper;
    helper.rsqr=rsqr; helper.xbj=xbj; helper.delta=delta;
    helper.nuke=&nucleus; helper.dip=this;
    
    gsl_function int_helper;
    int_helper.function=&inthelperf_ipsat_coherentavg2012;
    int_helper.params=&helper;
    
    size_t eval; 
    REAL result,abserr;

    //int status = gsl_integration_qng(&int_helper, 0, MAXB, 
    //        0, OSCAVGINTACCURACY, &result, &abserr, &eval);
    gsl_integration_workspace *ft_workspace 
     = gsl_integration_workspace_alloc(MAXITER_FT2012);
    int status=gsl_integration_qag(&int_helper, 0, MAXB, 0, FTINTACCURACY, 
        MAXITER_FT2012, GSL_INTEG_GAUSS61, ft_workspace, &result, &abserr);
    gsl_integration_workspace_free(ft_workspace);

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
 * (in factorized approximation)
 *
 * Non-factorized form: \int d^2 b 2(1-exp(-r^2F(x,r)))
 *
 * Units: 1/GeV^2
 */
REAL Dipxs_IPSat2012::TotalDipxsection_proton(REAL rsqr, REAL xbj)
{
    if (!factorize)
    {
        inthelper_ipsatavg2012 helper;
        helper.rsqr=rsqr; helper.xbj=xbj; helper.delta=-1;
        helper.nuke=&nucleus; helper.dip=this;
        
        gsl_function int_helper;
        int_helper.function=&inthelperf_totxs_satp2012;
        int_helper.params=&helper;
        
        size_t eval;
        REAL result,abserr;

        gsl_integration_workspace *ft_workspace 
         = gsl_integration_workspace_alloc(MAXITER_FT2012);
        int status=gsl_integration_qag(&int_helper, 0, MAXB, 0, FTINTACCURACY, 
            MAXITER_FT2012, GSL_INTEG_GAUSS61, ft_workspace, &result, &abserr);
        gsl_integration_workspace_free(ft_workspace);
        if (status) std::cerr << "Error " << status << " at " << __FILE__ << ":"
            << __LINE__ << ": Result " << result << ", abserror: " << abserr <<
            ", (tot qqp-xs)" << std::endl;
        return result;
    }   
	cerr << "Factorization is not implemented in TotalDipxsection_proton" << endl;
}


REAL inthelperf_totxs_satp2012(REAL b, void* p)
{
    inthelper_ipsatavg2012* par = (inthelper_ipsatavg2012*)p;
    return 2.0*M_PI*b*2.0
        *par->dip->Qq_proton_amplitude(par->rsqr, par->xbj, b);
}

/*
 * Diffractive dipole-proton amplitude as a function of \Delta
 * Integrated over impact parameter dependence
 * 1/2 * \int d^2 b exp(-ib*\Delta) d\sigma/d^2b
 * (factorized version) = 2*Pi*B_p*C exp(-B_p \Delta^2 / 2)
 * 
 * Otherwise same as Dipxs_IPNonSAT::Dipxsection_proton but
 * we take into account the unitarity requirement.
 *
 * It is also approximated that 1-exp(-r^2T(b)...) = T(b)(1-exp(...))
 */
REAL Dipxs_IPSat2012::DipoleAmplitude_proton(REAL rsqr, REAL xbj, REAL delta)
{
    // Note: We cannot use FactorC() here as FactorC() depends on the mode
    // used: in NONSAT_P mode it doesn't saturate
    
    if (factorize) // Factorize T(b) dependency
    {
			cerr << "Factorized ipsat2012 is todo!" << endl;
    }
    //cerr << "Check non-factorized version! dipxs_ipsat2012.cpp DipoleAmplitude_proton" << endl;
    // Else: Integrate numerically over impact parameter dependence, oscillatory
    // integral...
    inthelper_ipsatavg2012 helper;
    helper.rsqr=rsqr; helper.xbj=xbj; helper.delta=delta;
    helper.nuke=&nucleus; helper.dip=this;
    
    gsl_function int_helper;
    int_helper.function=&inthelperf_satp2012;
    int_helper.params=&helper;
    
    size_t eval;
    REAL result,abserr;
    //int status = gsl_integration_qng(&int_helper, 0, MAXB, 
    //        0, FTINTACCURACY, &result, &abserr, &eval);
    gsl_integration_workspace *ft_workspace 
     = gsl_integration_workspace_alloc(MAXITER_FT2012);
    int status=gsl_integration_qag(&int_helper, 0, MAXB, 0, FTINTACCURACY, 
        MAXITER_FT2012, GSL_INTEG_GAUSS61, ft_workspace, &result, &abserr);
    gsl_integration_workspace_free(ft_workspace);
    if (status) std::cerr << "Error " << status << " at " << __FILE__ << ":"
        << __LINE__ << ": Result " << result << ", abserror: " << abserr <<
        ", t=" << SQR(delta) << " GeV^2 " << std::endl;
        
    return result;
}


// Inthelper for case without factorization
REAL inthelperf_satp2012(REAL b, void* p)
{
    inthelper_ipsatavg2012* par = (inthelper_ipsatavg2012*)p;
    return 2.0*M_PI*b*gsl_sf_bessel_J0(b*par->delta)
        *par->dip->Qq_proton_amplitude(par->rsqr, par->xbj, b);
}

/*
 * Dipole-proton scattering amplitude
 */
REAL Dipxs_IPSat2012::Qq_proton_amplitude(REAL rsqr, REAL xbj, REAL b)
{
    if (factorize)
	cerr << "Factorized ipsat2012 is not implemented" << endl;	
    return DipoleAmplitude(std::sqrt(rsqr), xbj, b);
}



REAL Dipxs_IPSat2012::GetB_p()
{
    return B_p;
}

void Dipxs_IPSat2012::SetFactorize(bool f)
{  
    factorize=f;
}

REAL Dipxs_IPSat2012::Bp()
{
    return B_p;
}

REAL Dipxs_IPSat2012::Sigma0()
{
    return 4.0*M_PI*Bp();
}

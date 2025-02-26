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
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_bessel.h>
#include <vector>
#include "ipsat_mz/dipoleamplitude.hpp"
#include "dipxs_ipsat2012.h"
#include "dipole.h"

// Integral helpers
REAL inthelperf_satp2012(REAL b, void* p);  // Integrate over impact parameter 
REAL inthelperf_totxs_satp2012(REAL b, void* p);

const REAL MAXB=120;
const int MAXITER_BINT=1000;

const REAL AVGINTACCURACY = 0.001;
const REAL FTINTACCURACY = 0.0001; 

using std::cout; using std::endl; using std::cerr;

// Actual routine to compute dipole amplitude
extern "C" {
  double dipole_amplitude_(double* xBj, double* r, double* b, int* param);
};
//double Dipxs_IPSat2012::DipoleAmplitude(double r, double xBj, double b=0, int param=2, double Bp=4)
double Dipxs_IPSat2012::DipoleAmplitude(double r, double xBj, double b, int param )
{
	//	if (r==0) return 0;
	//	return ipsat_mz->N(r, xBj, b);
	//	return 0.5*dipole_amplitude_(&xBj, &r, &b, &param);
		// param 2: m_c=1.4
	double n = 0.5*dipole_amplitude_(&xBj, &r, &b, &param);
	if (n ==1 or n==0) return n;
	// Change Bp
	double lns = log(1.0-n);	
	double newlns = lns / (1.0 / (2.0*M_PI*4.0) * exp(-b*b/(2.0*4.0))) * 1.0 / (2.0*M_PI*Bp()) * exp(-b*b/(2.0*Bp()) );

	return 1.0 - exp(newlns);
}


Dipxs_IPSat2012::Dipxs_IPSat2012(Nucleus &n) : Dipxs(n)
{
    B_p=4;
    Intialize();
}

void Dipxs_IPSat2012::Intialize()
{
    factorize=false;
	ipsat_mz = new IPsat_MZ::DipoleAmplitude(2.289363553168, std::sqrt(1.1), 0.08289088639946, 2.195310911936, 1.35277437092);
	ipsat_mz->SetSaturation(true);
//	ipsat_mz =  new IPsat_MZ::DipoleAmplitude(4.297444629517, std::sqrt(1.1), -0.006657294973805, 3.039134356321, 1.350367375905);
//    ipsat_mz->SetSaturation(false);
//	factorize=true;
	//factorize=true;
}

Dipxs_IPSat2012::~Dipxs_IPSat2012()
{
	delete ipsat_mz;
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
             * exp(-2.0*A*M_PI*B_p*twstmp*(par->dip->FactorC(r*r,x)
               + par->dip->FactorC(r2*r2, x) ) )
             * pow( 
                M_PI*B_p*par->dip->FactorC(r*r, x)
                * par->dip->FactorC(r2*r2, x)*twstmp
                / (1.0-2.0*M_PI*B_p*twstmp
                *(par->dip->FactorC(r*r, x)+par->dip->FactorC(r2*r2, x)
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

	// Factorized approximation: factorize e^(-b62/(2B)) in front
	double n = DipoleAmplitude(std::sqrt(rsqr), xbj, 0);

	return 4.0*M_PI*GetB_p()*n;
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
    
    if (factorize) // Factorize T(b) dependency
    {
		// Factorized approximation: factorize e^(-b62/(2B)) in front
		double n = DipoleAmplitude(std::sqrt(rsqr), xbj, 0);

		return 2.0*M_PI*GetB_p()*std::exp(-GetB_p()*SQR(delta)/2.0)*n;
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
	{
		// Factorized approximation: factorize e^(-b^2/(2B)) in front
		double n = DipoleAmplitude(std::sqrt(rsqr), xbj, b);
		double tp = std::exp(-SQR(b)/(2.0*GetB_p()));
		// Remove tmp from exponent
		double expon = std::log(1.0 - n);
		expon /= tp;
		return tp * (1.0 - std::exp( expon ) ); 
	}
    return DipoleAmplitude(std::sqrt(rsqr), xbj, b, 2 );
}

// Dipole amplitude without T_p, but with 1/(2pi B)
REAL Dipxs_IPSat2012::FactorC(REAL rsqr, REAL x)
{
		// Factorized approximation: factorize e^(-b^2/(2B)) in front
		double b=0;
		double n = DipoleAmplitude(std::sqrt(rsqr), x, b);
		double tp = std::exp(-SQR(b)/(2.0*GetB_p()));
		// Remove tmp from exponent
		double expon = std::log(1.0 - n);
		expon /= tp;
		return 1.0 - std::exp(expon);	
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

struct satscalehelper_ipsat2012
{
	Dipxs_IPSat2012* dipole;
	int A;
	double b;
	double xbj;
	Nucleus *nuke;
};

double satscalehelperf_ipsat2012(double r, void* p)
{
	satscalehelper_ipsat2012* par= (satscalehelper_ipsat2012*) p;
	if (par->A == 1)
		return par->dipole->Qq_proton_amplitude(r*r, par->xbj, par->b) - (1.0 - std::exp(-0.5));
	

// Factorized approximation: factorize e^(-b62/(2B)) in front
	 
     double B_p=4.0;
     double n = par->dipole->DipoleAmplitude(r, par->xbj, 0);
     double tp = 1.0/(2.0*B_p*M_PI)*std::exp(-SQR(0)/(2.0*B_p));
     // Remove tp from exponent
     double expon = std::log(1.0 - n);
     expon /= tp;
	 return  (1.0 - std::exp( expon  * par->A * par->nuke->T_WS(par->b) ) ) - (1.0-std::exp(-0.5));
	
}


REAL Dipxs_IPSat2012::SaturationScale(double x, int A, double b)
{
	 const gsl_root_fsolver_type *T = gsl_root_fsolver_bisection;
     gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
	 gsl_function f;
	 f.function = satscalehelperf_ipsat2012;
	 satscalehelper_ipsat2012 par;
	 par.dipole=this;
	 par.xbj=x;
	 par.b=b;
	 par.A=A;
	 Nucleus nuke(A);
	 nuke.Intialize();
	 par.nuke=&nuke;
	 f.params=&par;
	 double maxr = 50;
	 double ACC = 0.01;
	 int MAXITER = 100;
	 gsl_root_fsolver_set(s, &f, 0, maxr);
     int iter=0; int status; double min,max;
     do
     {
         iter++;
         gsl_root_fsolver_iterate(s);
         min = gsl_root_fsolver_x_lower(s);
         max = gsl_root_fsolver_x_upper(s);
         status = gsl_root_test_interval(min, max, 0, ACC);
     } while (iter <= MAXITER and status == GSL_CONTINUE);
																			     
	  double satr = gsl_root_fsolver_root(s);
											     
	  double qs = std::sqrt( 2.0) / satr;
																							     
	  gsl_root_fsolver_free(s);
	return qs*qs;
}


bool Dipxs_IPSat2012::GetFactorize()
{
	return factorize;
}

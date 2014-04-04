/*
 * Dipole cross section
 * Dipole amplitude from BK
 * B-dependence like in IIM
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2014
 */
 
#include "dipxs_bk.h"
#include "nucleus.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_bessel.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <amplitudelib/amplitudelib.hpp>


const REAL AVGITACCURACY = 0.001;
const REAL FTINTACCURACY = 0.001; 
const REAL MAXB=90;
const size_t MAXITER_FT=1000;  // Max number of interations when 
            // calculating the fourier transformation of the amplitude

Dipxs_BK::Dipxs_BK(Nucleus& nuke) : Dipxs(nuke)
{
    Intialize();
}

Dipxs_BK::Dipxs_BK(Nucleus& nuke, std::string file) : Dipxs(nuke)
{

    Intialize();
}

void Dipxs_BK::Intialize()
{
    prevdelta=prevft=-1;

    N = new AmplitudeLib("/nashome2/hejajama/rbk/data/bestfits/mve/bestfit.dat");
    cout <<"# Amplitude read from file /nashome2/hejajama/rbk/data/bestfits/mve/bestfit.dat, info: " << N->GetString() << endl;
    B_D=4.0;
    sigma0=83.985; // GeV^(-2)
    
    ft_workspace_coh = gsl_integration_workspace_alloc(MAXITER_FT);
}

Dipxs_BK::~Dipxs_BK()
{
    gsl_integration_workspace_free(ft_workspace_coh);
    delete N;
}

/* Amplitude squared averaged over nucleon configurations as a 
 * function of \Delta and r,r' (and x)
 * \int d^2 b_1 ... d^2 b_A T_A(b_1)...T_A(B_A) 
 *      * \int d^2 b d^2 b' e^(-i(b-b')*\Delta) 
 *      * 1/2*(d\sigma^2 / d^2 b)(b,r) 1/2*(d\sgima^2 / d^2b)(b',r')
 *
 * IIM model we use:
 * A = S(b)N(r,x), where
 * S(b) = exp(-b^2/(2*B_D))   impact parameter dependence
 * N(r,x) is the dipole scattering amplitude
 * eq. (13) in paper arXiv 0706.2682v1
 */
 

struct inthelper_bkavg
{
    REAL r1q,r2q;
    REAL xbj;
    REAL delta;
    
    Dipxs_BK* dip;
    Nucleus* nuke;
};

REAL inthelperf_bkavg(REAL b, void* p)
{
    // Calculate the sum
    REAL sum=0;
    inthelper_bkavg* par=(inthelper_bkavg*)p;
    int A = par->nuke->GetA();
    REAL x = par->xbj;
    REAL B_D=par->dip->GetB_D();
    REAL twstmp = par->nuke->T_WS(b);

    if (N_MAX_BK > 1)
        cerr << "Incoherent BK only supports the first term!" << endl;
    

    // Assume i=1 so that we can write 16\pi^2 B^2 -> 
    //for (int i=1; i<=N_MAX_BK; i++)
    //{
        sum=A*exp(-B_D*SQR(par->delta))
             * exp(-par->dip->Sigma0()/2.0*A*twstmp*(par->dip->DipoleAmplitude(par->r1q,x)
               + par->dip->DipoleAmplitude(par->r2q,x) ) )
             *  
                par->dip->DipoleAmplitude(par->r1q,x)
                * par->dip->DipoleAmplitude(par->r2q,x)*twstmp
                / (1.0-par->dip->Sigma0()/2.0*twstmp
                *(par->dip->DipoleAmplitude(par->r1q,x)
                    + par->dip->DipoleAmplitude(par->r2q,x) )
                 ) ;    
    //}
    
    sum*=1.0/4.0*SQR(par->dip->Sigma0());
    
    // 2D integral -> 1D
    sum*=2.0*M_PI*b;
    return sum;
}
 
REAL Dipxs_BK::DipoleAmplitude_sqr_avg(REAL rsqr, REAL r2sqr, REAL xbj, 
                REAL delta)
{
    REAL r1 = sqrt(rsqr); REAL r2 = sqrt(r2sqr);    // Stupid..
    int A = nucleus.GetA();
    
    REAL bd=B_D;
    
 
        // Generalize as a product of S matrix elements
        inthelper_bkavg helper;
        helper.r1q=r1; helper.r2q=r2;
        helper.xbj=xbj; helper.delta=delta;
        helper.dip=this;
        helper.nuke=&nucleus;
        
        gsl_function fun;
        fun.function=&inthelperf_bkavg;
        fun.params=&helper;
        
        size_t eval;
        REAL result,abserr;

        int status = gsl_integration_qng(&fun, 0, MAXB, 
                0, AVGITACCURACY, &result, &abserr, &eval);
        
        if (status) std::cerr << "Error " << status << " at " << __FILE__ << ":"
            << __LINE__ << ": Result " << result << ", abserror: " << abserr 
            << " (t=" << delta*delta <<")" << std::endl;
        
        return result;
    

}

/* 
 * Amplitude for coherent dipole-nucleus scattering
 * \int d^2 b_1...d^2 b_A T_A(b_1)...T_A(B_A)
 *      *\int d^2 b e^(-ib*\Delta)
 *      * 1/2*(d\sigma^A / d^2b )(b,r,x) 
 *
 * In bk model this can be derived to be
 * \int d^2b e^(-b*\Delta) (1 - exp(-A/2*T_A(b)*\sigma_dip^p(r,x)) )
 *
 * The nasty part here is that we have to perform the fourier transformation
 * numerically. 
 */

struct inthelper_coherent
{
    REAL rsqr, delta, xbj;
    Dipxs* dip;
    Nucleus* nuke;
};

REAL inthelperf_bk_coherentavg(REAL b, void* p)
{
    inthelper_coherent* par=(inthelper_coherent*)p;
    int A = par->nuke->GetA();
    return 2.0*M_PI*b*gsl_sf_bessel_J0(b*par->delta)*
         (1.0-exp(-A/2.0*par->nuke->T_WS(b)  
        * par->dip->TotalDipxsection_proton(par->rsqr, par->xbj) )  );
}

REAL Dipxs_BK::CoherentDipoleAmplitude_avg(REAL rsqr, REAL xbj, REAL delta)
{
    inthelper_coherent helper;
    helper.rsqr=rsqr; helper.xbj=xbj; helper.delta=delta;
    helper.nuke=&nucleus; helper.dip=this;
    
    gsl_function int_helper;
    int_helper.function=&inthelperf_bk_coherentavg;
    int_helper.params=&helper;
    
    size_t eval; 
    REAL result,abserr;
    
    int status=gsl_integration_qag(&int_helper, 0, MAXB, 0, FTINTACCURACY, 
        MAXITER_FT, GSL_INTEG_GAUSS41, ft_workspace_coh, &result, &abserr);

    if (status) std::cerr << "Error " << status << " at " << __FILE__ << ":"
        << __LINE__ << ": Result " << result << ", abserror: " << abserr 
        << " relerror: " << abserr/result << " (t=" << delta*delta <<")" << std::endl;
    return result;

}

/*
 * Dipole-proton amplitude as a function of \Delta
 * Integrated over impact parameter dependence
 * 
 * 2 \pi B_D exp(-B_D*\Delta^2/2) * DipoleAmplitude(rq,x)
 *
 */
 ///TODO: What is the impact parameter profile??? Factor 2 correctly?
REAL Dipxs_BK::DipoleAmplitude_proton(REAL rsqr, REAL xbj, REAL delta)
{
    REAL rq = sqrt(rsqr);
    //if (delta > 0)
    //    cerr << "Dipxs_BK::DipoleAmplitude_proton does not support delta dependence " << LINEINFO << endl;
        
    return Sigma0()/2.0*DipoleAmplitude(rq,xbj);
}

/*
 * Total dipole-proton cross section
 * Calculated as \int d^2 b d\sigma/d^2 b
 *
 * Units: 1/GeV^2
 */
REAL Dipxs_BK::TotalDipxsection_proton(REAL rsqr, REAL xbj)
{
    REAL rq =sqrt(rsqr);
    return Sigma0()*DipoleAmplitude(rq,xbj);    // todo: sigma0
}

/* 
 * Dipole scattering amplitude N(r,x)
 */
REAL Dipxs_BK::DipoleAmplitude(REAL r, REAL x)
{
    if (x > N->X0())
        return 0;
    return N->N(r, std::log(N->X0()/x));   

}
 

REAL Dipxs_BK::Dipxsection(REAL rsqr, REAL xbjork, Vec b, 
            std::vector<Vec> &nucleons)
{
    // Not implemented
    return 0;
}

/*
 * Dipole-proton scattering amplitude as a function of b and r
 * returns 1/2*d\sigma/d^2b
 */
REAL Dipxs_BK::Qq_proton_amplitude(REAL rsqr, REAL xbj, REAL b)
{
    return exp(-SQR(b)/(2.0*B_D))*DipoleAmplitude(std::sqrt(rsqr),xbj);
}
/*
 * Saturation scale
 */
REAL Dipxs_BK::Q_s(REAL x)
{
    cerr << "Qs not implemented for bk! " << LINEINFO << endl;
}


REAL Dipxs_BK::GetB_D()
{
    return B_D;
}

REAL Dipxs_BK::Bp()
{
    return B_D;
}

REAL Dipxs_BK::Sigma0()
{
    return 2.0/3.0*sigma0;
}

/*
 * Dipole cross section in IIM model
 * Source: Article by Cyrillie Marquet
 * arXiv:0706.2682v1
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
#include "dipxs_iim.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>

const REAL AVGITACCURACY = 0.001;
const REAL MAXB=90;

Dipxs_IIM::Dipxs_IIM(Nucleus& nuke) : Dipxs(nuke)
{
    prevdelta=prevft=-1;
    if(ReadParameters("iim.dat"))
    {
        std::cerr << "Could not read parameters from file iim.dat" 
            << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

Dipxs_IIM::Dipxs_IIM(Nucleus& nuke, std::string file) : Dipxs(nuke)
{
    prevdelta=prevft=-1;
    if(ReadParameters(file))
    {
        std::cerr << "Could not read parameters from file " << file
            << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

/* Amplitude squared averaged over nucleon configurations as a 
 * function of \Delta and r,r' (and x)
 * \int d^2 b_1 ... d^2 b_A T_A(b_1)...T_A(B_A) 
 *      * \int d^2 b d^2 b' e^(-i(b-b')*\Delta) 
 *      * (d\sigma^2 / d^2 b)(b,r) (d\sgima^2 / d^2b)(b',r')
 *
 * IIM model we use:
 * A = S(b)N(r,x), where
 * S(b) = exp(-b^2/(2*B_D))   impact parameter dependence
 * N(r,x) is the dipole scattering amplitude
 * eq. (13) in paper arXiv 0706.2682v1
 * 
 * This is generalized to a electron-nucleus scattering by two different ways:
 * 1) IIM_MODE=IIM_IPNONSAT:
 * * S(b) = \sum_i S(b-b_i), as it is done with IPNonSat model
 * 2) IIM_MODE=IIM_IPSAT
 * * S-matrix is the product of dipole-proton S matrix elements
 */
 
// Some integral helpers for IIM_IPSAT case 
// First the required integral helpers
struct inthelper_iimavg
{
    REAL r1q,r2q;
    REAL xbj;
    REAL delta;
    
    Dipxs_IIM* dip;
    Nucleus* nuke;
};

REAL inthelperf_iimavg(REAL b, void* p)
{
    // Calculate the sum
    REAL sum=0;
    inthelper_iimavg* par=(inthelper_iimavg*)p;
    int A = par->nuke->GetA();
    REAL x = par->xbj;
    REAL B_D=par->dip->GetB_D();
    REAL twstmp = par->nuke->T_WS(b);
    for (int i=1; i<=N_MAX_IIM; i++)
    {
        sum+=gsl_sf_choose(A,i)/((REAL)i)*exp(-B_D*SQR(par->delta)/i)
             * exp(-2*A*M_PI*B_D*twstmp*(par->dip->DipoleAmplitude(par->r1q,x)
               + par->dip->DipoleAmplitude(par->r2q,x) ) )
             * pow( 
                M_PI*B_D*par->dip->DipoleAmplitude(par->r1q,x)
                * par->dip->DipoleAmplitude(par->r2q,x)*twstmp
                / (1-2*M_PI*B_D*twstmp
                *(par->dip->DipoleAmplitude(par->r1q,x)
                    + par->dip->DipoleAmplitude(par->r2q,x) )
                 )  ,i);    
    }
    
    sum*=16*M_PI*B_D;
    
    // 2D integral -> 1D
    sum*=2*M_PI*b;
    return sum;
}
 
REAL Dipxs_IIM::Dipxsection_sqr_avg(REAL rsqr, REAL r2sqr, REAL xbj, 
                REAL delta)
{
    REAL r1 = sqrt(rsqr); REAL r2 = sqrt(r2sqr);    // Stupid..
    int A = nucleus.GetA();
    
    REAL Q_sA=Q_s(xbj);
    REAL r1Q = Q_sA*r1; REAL r2Q=Q_sA*r2;
    REAL bd=B_D;
    
 
    if (IIM_MODE==IIM_IPNONSAT)
    {
        // Generalize to qq-nucleus scattering as it was done with IPNonSat
        // model: amplitude=\sum_i amplitude_qq(b-b_i)
    
        // If we just calculated this when this function was called previously
        // Yeah, this is quite an ugly hack, but this optimizes this quite much!
        REAL ft;
        if (prevdelta==delta)
            ft=prevft;
        else
            ft=nucleus.FT_T_WS(delta);

        prevdelta=delta; prevft=ft;
        return 4*SQR(M_PI)*SQR(bd)*exp(-bd*delta*delta)*DipoleAmplitude(r1Q,xbj)
            * DipoleAmplitude(r2Q,xbj)*A*( 1.0 + (A-1.0)*ft*ft );
     }
     else if (IIM_MODE==IIM_IPSAT)
     {
        // Generalize as a product of S matrix elements
        inthelper_iimavg helper;
        helper.r1q=r1Q; helper.r2q=r2Q;
        helper.xbj=xbj; helper.delta=delta;
        helper.dip=this;
        helper.nuke=&nucleus;
        
        gsl_function fun;
        fun.function=&inthelperf_iimavg;
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
     

}

/*
 * Dipole-proton amplitude as a function of \Delta
 * Integrated over impact parameter dependence
 * 
 * 2 \pi B_D exp(-B_D*\Delta^2/2) * DipoleAmplitude(rq,x) * 2
 * Last * 2 is due to the fact that d\sigma/d^2b = 2 * amplitude
 */
REAL Dipxs_IIM::Dipxsection_proton(REAL rsqr, REAL xbj, REAL delta)
{
    REAL rq = Q_s(xbj)*sqrt(rsqr);
    return 2*2*M_PI*B_D*exp(-B_D*SQR(delta)/2.0)*DipoleAmplitude(rq,xbj);
}

/* 
 * Dipole scattering amplitude N(rQ,x)
 */
REAL Dipxs_IIM::DipoleAmplitude(REAL rq, REAL x)
{
    if (rq <= 2.0)
        return N0*pow(rq/2.0,2.0*gammac)*exp(-2*pow(log(rq/2.0),2)
            /(kappa*lambda*log(1/x) ) );
    else
        return 1-exp(-4*alpha*pow(log(beta*rq),2) );    

}
 

REAL Dipxs_IIM::Dipxsection(REAL rsqr, REAL xbjork, Vec b, 
            std::vector<Vec> &nucleons)
{
    // Not implemented
    return 0;
}

/*
 * Saturation scale
 */
REAL Dipxs_IIM::Q_s(REAL x)
{
    return pow(x0/x,lambda/2);
}

/*
 * Read parameters from the given file
 * Syntax is: key:value
 * Lines starting with # are ignored
 * Dimension of values is GeV^n
 * Keys: alpha, beta, x_0, N_0, kappa, lambda, gammac, B_D
 * 
 * Returns -1 in case of error, 0 otherwise
 */
int Dipxs_IIM::ReadParameters(std::string file)
{
    std::ifstream f(file.c_str());
  
    if (!f.is_open())
    {
        return -1;
    }
    std::string line;
    
    while(!f.eof() )
    {
        std::getline(f, line);
        if (line[0]=='#')   // Comment line
            continue;
        if (line.substr(0,5)=="alpha") 
            alpha=StrToReal(line.substr(6));
        if (line.substr(0,4)=="beta") 
            beta=StrToReal(line.substr(5));
        if (line.substr(0,3)=="x_0") 
            x0=StrToReal(line.substr(4));
        if (line.substr(0,3)=="N_0") 
            N0=StrToReal(line.substr(4));
        if (line.substr(0,5)=="kappa") 
            kappa=StrToReal(line.substr(6));
        if (line.substr(0,6)=="lambda") 
            lambda=StrToReal(line.substr(7));
        if (line.substr(0,6)=="gammac") 
            gammac=StrToReal(line.substr(7));
        if (line.substr(0,3)=="B_D")
            B_D=StrToReal(line.substr(4));
    }
    f.close();

    return 0;

}

REAL Dipxs_IIM::GetB_D()
{
    return B_D;
}



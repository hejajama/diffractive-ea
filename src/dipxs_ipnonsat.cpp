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
    prevdelta=prevft=-1;
    B_p=DEFAULT_B_p;
}

Dipxs_IPNonSat::Dipxs_IPNonSat(Nucleus &nucleus_, REAL bp) :
     Dipxs(nucleus_)
{
    prevdelta=prevft=-1;
    B_p=bp;
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
 *
 * Case A=1 (dipole-proton collision) is handled separately:
 * |A|^2 = 1/4*sigmap(r)*sigmap(r')*exp(-B_p*\Delta^2)
 * Factor 1/4 as this is amplitude, not Fourier transformed d\sigma/d^2b
 */

REAL Dipxs_IPNonSat::DipoleAmplitude_sqr_avg(REAL rsqr, REAL r2sqr, 
    REAL xbjork, REAL delta)
{
    REAL result;
    // In case dipole-proton use just the transverse profile of proton and
    // do not sum over nucleus states / average over nucleons
    if (nucleus.GetA()==1)
    {
        result=Sigmap(rsqr, xbjork)/2.0*Sigmap(r2sqr, xbjork)/2.0*exp(-B_p*SQR(delta));    
        return result;
    }    
    // Dipole-nucleus
    
    // If we just calculated this when this function was called previously
    // Yeah, this is quite an ugly hack, but this optimizes this quite much!
    if (prevdelta==delta)
        result=prevft;
    else
    {   
        result=nucleus.FT_T_WS(delta);
    }
    
    REAL A = nucleus.GetA();
    
    prevdelta=delta; prevft=result;
    
    return Sigmap(rsqr,xbjork)/2.0*Sigmap(r2sqr,xbjork)/2.0*exp(-B_p*SQR(delta))
        *A*(1.0+(A-1.0)*result*result);
}


/* 
 * Amplitude for coherent dipole-nucleus scattering
 * \int d^2 b_1...d^2 b_A T_A(b_1)...T_A(B_A)
 *      *\int d^2 b e^(-ib*\Delta)
 *      * 1/2*(d\sigma^A / d^2 B)(b,r,x)
 *
 * Calculated by Caldwell and Kowalski, arXiv:0909:1254v1
 *
 * A = 1/2*A*Sigmap*Exp(-B_p*\Delta^2/2)*\int d^2b exp(-ib.\Delta) T_A(b)
 */
REAL Dipxs_IPNonSat::CoherentDipoleAmplitude_avg(REAL rsqr, REAL xbj, 
                REAL delta)
{                

    return 1.0/2.0*nucleus.GetA()*Sigmap(rsqr, xbj)*exp(-B_p*SQR(delta)/2.0)
        *nucleus.FT_T_WS(delta);
}

/*
 * Dipole-proton amplitude as a function of \Delta
 * Integrated over impact parameter dependence
 * A = sigmap(r)/2 * exp(-B_p \Delta^2/2)
 */
REAL Dipxs_IPNonSat::DipoleAmplitude_proton(REAL rsqr, REAL xbj, REAL delta)
{
    return Sigmap(rsqr, xbj)/2.0*exp(-B_p*SQR(delta)/2.0);
}

/*
 * Total dipole-proton cross section
 * Calculated as \int d^2 b d\sigma/d^2 b
 * = \pi^2/3*r^2*\alpha_s(\mu^2)*xg(x,\mu^2) * \int d^2 b T_p(b)/(2*\pi*B_p)
 * = 2 * Gdist() * r^2
 *
 * Units: 1/GeV^2
 */
REAL Dipxs_IPNonSat::TotalDipxsection_proton(REAL rsqr, REAL xbj)
{
    return 2.0*nucleus.GetGDist()->Gluedist(xbj, rsqr)*rsqr;
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

/*
 * Dipole-proton scattering amplitude as a function of b and r
 * returns 1/2*d\sigma/d^2b
 */
REAL Dipxs_IPNonSat::Qq_proton_amplitude(REAL rsqr, REAL xbj, REAL b)
{
    return Sigmap(rsqr, xbj)*exp(-SQR(b)/(2.0*B_p))/(2.0*M_PI*B_p);
}

/*
 * Parent class for dipole cross section
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
 
#include "dipxs.h"
#include <iostream>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

const REAL TOTXS_MAXT=2;  // Max |t| in GeV^2 when calculating total xs
const REAL TINTACCURACY=0.0001;

struct inthelper_dipxs
{
    REAL xbj;
    Dipxs* amplitude;
    REAL rsqr, r2sqr;
};

REAL inthelperf_dipxs(REAL t, void* p);
REAL inthelperf_protxs(REAL t, void*p);

Dipxs::Dipxs()
{
    //InitGsl();
}

Dipxs::Dipxs(Nucleus &nucleus_)
{
    nucleus=nucleus_;
    //InitGsl();
}

Dipxs::~Dipxs()
{
    //FreeGsl();
}

/*
 * Dipole cross section integrated over |t| as a function of r and r'
 * This is virtual method, so if this can be calculated analytically 
 * in a given model then this can be implemented in the class impelmenting
 * the given model
 *
 * This integrates numerically Dipxsection_sqr_avg which is purely virtual
 * method
 */
REAL Dipxs::Dipxsection_sqr_avg(REAL rsqr, REAL r2sqr, REAL xbj)
{
    gsl_function fun;   
    inthelper_dipxs inthelp;
    inthelp.amplitude=this;
    inthelp.xbj=xbj;
    inthelp.rsqr=rsqr; inthelp.r2sqr=r2sqr;
    fun.function=&inthelperf_dipxs;
    fun.params=&inthelp; 
        
    REAL result,abserr; size_t eval;
    int status = gsl_integration_qng(&fun, 0, TOTXS_MAXT, 0, TINTACCURACY, 
        &result, &abserr, &eval);
    if (status)
        std::cerr << "Total cross section integral failed to reach tolerance: "
        << "Result: " << result << ", abserr: " << abserr << std::endl;
    
    std::cerr << "Dipxs::Dipxsection_sqr_avg may not work, test it before use!"
        << std::endl;
    
    return result;

}

/*
 * Set new nucleus
 */
void Dipxs::SetNucleus(Nucleus& nuke)
{
    nucleus=nuke;
}

/*
 * Coherent dipole-nucleus scattering amplitude averaged over nucleon
 * configurations
 */
 REAL Dipxs::CoherentDipoleAmplitude_avg(REAL rsqr, REAL xbj, 
                REAL delta)
{
    std::cerr << "CoherentDipoleAmplitude_avg is not yet implemented" << std::endl;
    return 0;
}

REAL Dipxs::CoherentDipxsection_avg(REAL rsqr, REAL xbj, 
                REAL delta)
{
    std::cerr << "Dipxs::CoherentDipxsection_avg is deprecated." << std::endl;
    return 0;
}

/*
 * Strong coupling constant as a function of r
 */
REAL Dipxs::Alphas_r(REAL rsqr)
{
    return Alpha_s(Mu2(rsqr));    // Alpha_s(Q^2) defined in dipole.h
}

Nucleus& Dipxs::GetNucleus()
{
    return nucleus;
}

REAL Dipxs::MaxB()
{
    return nucleus.MaxR();
}

REAL Dipxs::Dipxsection(REAL rsqr, REAL xbjork, REAL x, REAL y,
                std::vector<Vec>& nucleons)
{
    Vec tmpvec(x,y);
    return Dipxsection(rsqr, xbjork, tmpvec, nucleons);

}

std::ostream& operator<<(std::ostream& os, Dipxs& ic)
{
    return os << " Dipxs object ";

}



REAL inthelperf_dipxs(REAL t, void* p)
{
    inthelper_dipxs* par = (inthelper_dipxs*)p;
    return par->amplitude->DipoleAmplitude_sqr_avg(par->rsqr, par->r2sqr, par->xbj, sqrt(t));
}

// Integral helper for dipole-proton integration over t
REAL inthelperf_protxs(REAL t, void* p)
{
    inthelper_dipxs* par = (inthelper_dipxs*)p;
    return par->amplitude->DipoleAmplitude_proton(par->rsqr, par->xbj, sqrt(t));
}

/*
 * Not yet implemented errors
 */

REAL Dipxs::Dipxsection(REAL rsqr, REAL xbj, Vec b, 
                std::vector<Vec>& nucleons)
{
    std::cerr << "Dipxsection as a function of nucleon positions is not " <<
       " implemented" << std::endl;
    return 0;
}

// Total dipole-proton cross section
REAL Dipxs::TotalDipxsection_proton(REAL rsqr, REAL xbj)
{
    std::cerr << "TotalDipxsection_proton is not implemented " << std::endl;
    return 0;
}

// Deprecated
REAL Dipxs::Dipxsecton_sqr_avg(REAL r, REAL r2, REAL x, REAL d)
{
    std::cerr << "Dipxs::Dipxsection_sqr_avg is deprecated! Use DipoleAmplitude_sqr_avg"
        << " instead!" << std::endl;
}

// Qq-proton/nucleus amplitude
REAL Dipxs::Qq_proton_amplitude(REAL rsqr, REAL xbj, REAL b)
{
    std::cerr << "Dipxs::Qq_amplitude is not impelmented" << std::endl;
}

REAL Dipxs::Sigma0()
{
    std::cerr <<"Dipxs::Sigma0 is not implemented" << std::endl;
}


REAL Dipxs::Bp()
{
    std::cerr <<"Dipxs::Bp is not implemented" << std::endl;
}

REAL Dipxs::SaturationScale(double x, int A, double b)
{
	std::cerr <<"Dipxs::SaturationScale is nto implemented!" << std::endl;
	return 0;
}

/*
 * Class to handle all nucleus related information, e.g. density
 * distribution etc.
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */

#include "nucleus.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_exp.h>
#include <iostream>

const REAL WSINTACCURACY=0.00001;

// Integral helpers 

struct wsinthelper{
    Nucleus *nuke;
};

REAL  wsinthelperf(REAL r, void * p){
  return (4*M_PI*SQR(r)*((wsinthelper*)p)->nuke->WS(r));
}

Nucleus::Nucleus()
{
    A=0; gdist=0;
}

Nucleus::Nucleus(int A_)
{
    A=A_; gdist=0;
    Intialize();
}

void Nucleus::SetA(int A_)
{
    A=A_; gdist=0;
    Intialize();
}


/*
 * Intializes class: normalizes WS distribution and intializes gluon 
 * distribution
 * Returns <0 in case of error, 0 otherwise
 */
int Nucleus::Intialize()
{

    // Gluon distribution
    if (gdist) delete gdist;    // If we already have allocated gdist
    gdist = new GDist_Toy();    // Use toy model by default
    
    
    // Normalize WS distribution    

    WS_RA=(1.12*pow(A,1.0/3.0) - 0.86*pow(A,-1.0/3.0))*FMGEV;
    WS_N=1;
    
    // Normalize WS distribution to 1
    REAL result,abserr;
    size_t eval;
    
    wsinthelper param;
    param.nuke=this;
    
    gsl_function int_helper;
    int_helper.function=&wsinthelperf;
    int_helper.params=&param;
    
    int status = gsl_integration_qng(&int_helper, 0,90, WSINTACCURACY, WSINTACCURACY, 
        &result, &abserr, &eval);
    if (status) { 
        std::cerr << "WS integral in Nucleus constructor failed with code " 
            << gsl_strerror(status) << std::endl;
        return -1;
    }
        
    WS_N = 1/result;
    return 0;

}

void Nucleus::SetGDist(GDist* gdist_)
{
    if (gdist) delete gdist;
    gdist=gdist_;
}

GDist* Nucleus::GetGDist()
{
    return gdist;
}   
    

/* 
 * Woods-Saxon density distribution
 */
REAL Nucleus::WS(REAL r)
{
    return WS_N/(1+gsl_sf_exp((r-WS_RA)/WS_delta));

}

/*
 * Proton shape
 */
REAL Nucleus::Tp(REAL r)
{
    return gsl_sf_exp(-SQR(r)/(2.0*B_G))/(2*M_PI*B_G);
}

REAL Nucleus::Tp(Vec b)
{
    return Tp(b.Len());
}



int Nucleus::GetA() { return A; }
REAL Nucleus::GetWS_RA() { return WS_RA; }
REAL Nucleus::GetWS_delta() { return WS_delta; }

std::ostream& operator<<(std::ostream& os, Nucleus& ic)
{
  os << "Nucleus, A = " << ic.GetA() << ", WS parameters: R_A = " << 
    ic.GetWS_RA() << ", delta = " << ic.GetWS_delta() ;

  return os;

}



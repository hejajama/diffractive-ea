/*
 * Class to handle all nucleus related information, e.g. density
 * distribution etc.
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_exp.h>
#include <iostream>

#include "nucleus.h"
#include "mersenne/mersenne.h"

using std::cout; using std::endl; using std::cerr;

const REAL WSINTACCURACY=0.000001;

// Integral helpers 

struct wsinthelper{
    Nucleus *nuke;
    REAL b;    // For z integral in T_WS where b is constant
};

// 3 dimensional integral of WS
REAL  wsinthelperf3d(REAL r, void * p){
  return (4*M_PI*SQR(r)*((wsinthelper*)p)->nuke->WS(r));
}

// 1 dimensional integral of WS
REAL wsinthelperf1d(REAL z, void* p){
    return ((wsinthelper*)p)->nuke->WS(sqrt(SQR(z)
        + SQR( ((wsinthelper*)p)->b) ) )   ;
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

Nucleus::~Nucleus()
{
/*    for (int i=0; i<nucleon_configs.size(); i++) { 
        (*nucleon_configs)[i].clear(); 
    }
    delete nucleon_configs;*/
}

void Nucleus::SetA(int A_)
{
    A=A_; 
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
    int_helper.function=&wsinthelperf3d;
    int_helper.params=&param;
    
    int status = gsl_integration_qng(&int_helper, 0,MaxR(), WSINTACCURACY, WSINTACCURACY, 
        &result, &abserr, &eval);
    if (status) { 
        std::cerr << "WS integral in Nucleus constructor failed with code " 
            << gsl_strerror(status) << std::endl;
        return -1;
    }
        
    WS_N = 1/result;
    //T_WS_0=T_WS(0);
    
    next_config=-1;    
    
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
    return WS_N/(1+exp((r-WS_RA)/WS_delta));
}

/*
 * Unnormalized distributions
 * Returns the propability to find a nucleon with radius r
 */
REAL Nucleus::WS_unnorm(REAL r)
{
    return WS(r)/WS_N;
}
/*
REAL Nucleus::T_WS_unnorm(REAL r)
{
    return T_WS(r)/T_WS_0;
}*/

/*
 * Integrated density distribution
 */
REAL Nucleus::T_WS(REAL b)
{
    REAL result,abserr;
    size_t eval;
    
    wsinthelper param;
    param.nuke=this;
    param.b=b;
    
    gsl_function int_helper;
    int_helper.function=&wsinthelperf1d;
    int_helper.params=&param;
    
    int status = gsl_integration_qng(&int_helper, 0, MaxR(), 
            WSINTACCURACY, WSINTACCURACY, &result, &abserr, &eval);
    if (status) { 
        std::cerr << "WS integral in T_WS failed with code " 
            << gsl_strerror(status) << std::endl;
        return -1;
    }
    return result*2.0;  // Multiplied by 2.0 as integration is done from 0,
                        // not form -MaxR()

}

REAL Nucleus::T_WS(Vec b)
{
    // NOTE: There is one unecessary sqrt() as b.Len() is squared in
    // T_WS, optimize?
    return T_WS(b.Len());
}



/*
 * Proton shape
 */
REAL Nucleus::Tp(REAL rsqr)
{
    return exp(-rsqr/(2.0*B_G))/(2*M_PI*B_G);
}

REAL Nucleus::Tp(Vec b)
{
    return Tp(b.LenSqr());
}

/*
 * Nucleon configuration
 * Nucleon positions are generated from W-S distributio
 * mindist and maxdist are smallest and largest distance allowed 
 * between nucleons
 * 
 * Configuration is saved in private vector nucleon_configs
 * n is the number of configurations to generate
 *
 * Returns -1 in case of an error, 0 otherwise
 */

int Nucleus::GenerateRandomNucleonConfigurations(int n, REAL mindist, REAL maxdist)
{
    if (n<1) return -1;
    
    next_config=0;
    for (int i=0; i<nucleon_configs.size(); i++) { 
        nucleon_configs[i].clear(); 
    }
   
    nucleon_configs.clear();
    
    for (int i=0; i<n; i++)
    {
        std::vector<Vec> nucleons; 
        REAL smallestdist=999; REAL largestdist=0;
        REAL tmpdist;
        do
        {
            nucleons.clear(); smallestdist=999; largestdist=0;
            for (int i=0; i<A; i++)
            {
                Vec tmp;
                do {
                    Vec tmpvec (2.0*(mersenne()-0.5)*MaxR(),
                                2.0*(mersenne()-0.5)*MaxR(),
                                2.0*(mersenne()-0.5)*MaxR());
                    tmp=tmpvec;
                } while (mersenne() > WS_unnorm(tmp.Len())); // WS distribution!
                nucleons.push_back(tmp);
            }
             
                
            // Check smallestdist and largest dist
            for (int i=0; i<A; i++)
                for (int j=0; j<i; j++)
                {
                   tmpdist = sqrt( SQR(nucleons[j].GetX() - nucleons[i].GetX())
                        + SQR(nucleons[j].GetY() - nucleons[i].GetY()) 
                        + SQR(nucleons[j].GetZ()-nucleons[i].GetZ())); 
                   if (tmpdist < smallestdist) smallestdist=tmpdist;
                   if (tmpdist > largestdist) largestdist=tmpdist;    
            }
            
            // Try again if smallestdist/largestdist does not match the limits
            if (smallestdist < mindist or largestdist > maxdist)
                  cout << "Again... Smallest dist: " << smallestdist << " , largest: " << largestdist << endl;
        } while (smallestdist < mindist or largestdist > maxdist);
        
        nucleon_configs.push_back(nucleons);    
    }

    return 0;

}
 
std::vector<Vec>& Nucleus::RandomNucleonConfiguration()
{
    if (next_config<0) {
       cerr << "Random nucleon configuration asked before anything was " 
            << "generated, returning zero config" << endl;
            std::vector<Vec> *zero = new std::vector<Vec>; (*zero).push_back(Vec()); 
            return *zero;
    }
    next_config++;
    if (next_config >= nucleon_configs.size())
        next_config=0;

    return nucleon_configs[next_config];
}
/*
std::vector<Vec>& Nucleus::RandomNucleonConfiguration2d()
{
    std::vector<Vec> nucleons = RandomNucleonConfiguration();
    // 2d = integrated over z, so just move all nucleons to z=0
    for (int i=0; i<nucleons.size(); i++) nucleons[i].SetZ(0);
    return nucleons;

}*/

/*
 * Maximum radius
 * This is used e.g. as a upper limit of the integrals
 */
REAL Nucleus::MaxR()
{
    return WS_RA+4.0*WS_delta;
}

/* Some trivial functions */

int Nucleus::GetA() { return A; }
REAL Nucleus::GetWS_RA() { return WS_RA; }
REAL Nucleus::GetWS_delta() { return WS_delta; }

std::ostream& operator<<(std::ostream& os, Nucleus& ic)
{
  os << "Nucleus, A = " << ic.GetA() << ", WS parameters: R_A = " << 
    ic.GetWS_RA() << ", delta = " << ic.GetWS_delta() ;

  return os;

}



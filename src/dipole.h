#ifndef DIPOLE_H
#define DIPOLE_H

/*
 * Some constants and definitions
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */


#include <gsl/gsl_math.h>
#include <complex>
#include <string>
#include <iostream>
#include <tools/config.hpp>

// *** Datatypes ***
// REAL = double
//#define REAL double
typedef  double  REAL;
typedef  std::complex<REAL> COMPLEX;
//#define COMPLEX std::complex<REAL>

// *** Macros ***
//#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define ABS(val) ((val) < 0 ? -(val) : (val))


// *** Constants ***
const int NC=3;
const int Nf=3;
const REAL ALPHA_e = 1.0/137.035999679; 
const REAL e = sqrt(4*M_PI*ALPHA_e);
const REAL mu2_0 = 1.1699;  // \mu_0^2 = 1.1699 (Lappi) or 0.8 GeV^2 (KT) 
const REAL LAMBDAQCD2=0.21416*0.21416;   // GeV^2
const REAL FOURPI = (4.0*M_PI);
const REAL DEFAULT_B_p=4.0;	// Default proton shape for ipsat and ipnonsat
const REAL NBGEVSQR = 1.0/SQR(5.068)*1.0e7; // 1/Gev^4 -> nb/Gev^2
const REAL FMGEV = 5.068;       // fm * GeV
#define LINEINFO __FILE__ << ":" << __LINE__


// *** Functions ***
// eps(z,Q,m)^2 = z(1-z)*Q^2 + m^2
REAL epsfun(REAL z, REAL Qsqr, REAL msqr);
REAL epsfunsqr(REAL z, REAL Qsqr, REAL msqr); // eps^2
REAL Alpha_s(REAL Qsqr); // Strong coupling constant as a function of Q^2
REAL Mu2(REAL rsqr);   // 4/r^2 + mu_0

REAL StrToReal(std::string str);
int StrToInt(std::string str);

// GSL error handler
void ErrHandler(const char * reason,
                        const char * file,
                        int line,
                        int gsl_errno);
                        
// Wrapper for gsl_sf_exp to handle underflows
// Or just use exp() from gsl/gsl_math.h
REAL exp_wrap(REAL x);   

// *** Classes ***


#endif   // DIPOLE_H

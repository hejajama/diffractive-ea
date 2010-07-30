#ifndef DIPOLE_H
#define DIPOLE_H

/*
 * Some constants and definitions
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */


#include <gsl/gsl_math.h>


// *** Datatypes ***
// REAL = double
#define REAL double

// *** Macros ***
#define SQR(x) ((x)*(x))

// *** Constants ***
const int NC=3;
const REAL ALPHA_e = 1.0/137.0; 
const REAL e = sqrt(4*M_PI*ALPHA_e);
const REAL mu2_0 = 0.8; // \mu_0^2 = 0.8 GeV^2
const REAL LAMBDAQCD2=0.04;   // 0.2 GeV^2
#define FMGEV 5.068       // fm * GeV



// *** Functions ***
// eps(z,Q,m)^2 = z(1-z)*Q^2 + m^2
REAL epsfun(REAL z, REAL Qsqr, REAL msqr);
REAL epsfunsqr(REAL z, REAL Qsqr, REAL msqr); // eps^2
REAL Alpha_s(REAL Qsqr); // Strong coupling constant as a function of Q^2

// *** Classes ***


#endif   // DIPOLE_H

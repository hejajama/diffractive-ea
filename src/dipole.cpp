/*
 * Some commonly used functions
 * 
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
#include "dipole.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_exp.h>

// eps(z,Q,r) = sqrt(z(1-z)*Q^2 - m^2)
REAL epsfunsqr(REAL z, REAL Qsqr, REAL msqr)
{
    return z*(1-z)*Qsqr + msqr;
}

REAL epsfun(REAL z, REAL Qsqr, REAL msqr)
{
    return sqrt(epsfunsqr(z,Qsqr,msqr));
}


/* 
 * Q^2 dependent strong coupling constant
 * Takes into account only u,d ands s quarks
 */
REAL Alpha_s(REAL Qsqr)
{
    return 12.0*M_PI/((33.0-2.0*3.0)*log(Qsqr/LAMBDAQCD2));
}

/*
 * Wrapper for gsl_sf_exp
 * If |x| < 1e-9 then returns 1
 * If x < -100 returns 0
 * gsl_sf_exp would return underflow and abort the program
 */
REAL exp_wrap(REAL x)
{
    if (ABS(x)<1e-9)
        return 1;
    if (x<-100)
        return 0;
    return gsl_sf_exp(x);
}




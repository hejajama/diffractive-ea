/*
 * Some commonly used functions
 * 
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
#include "dipole.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_exp.h>
#include <string>
#include <sstream>
#include <cstdlib>

// eps(z,Q,r) = sqrt(z(1-z)*Q^2 + m^2)
REAL epsfunsqr(REAL z, REAL Qsqr, REAL msqr)
{
    return z*(1.0-z)*Qsqr + msqr;
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
    return 12.0*M_PI/( (33.0-2.0*Nf)*log(Qsqr/LAMBDAQCD2) );
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

/*
 * \mu^2 = 4/r^2 + mu_0
 */
REAL Mu2(REAL rsqr)
{
    return 4.0/rsqr + mu2_0;
}

/*
 * Str to REAL/int
 */
 /*
REAL StrToReal(std::string str)
{
    std::stringstream buff(str);
    REAL tmp;
    buff >> tmp;
    return tmp;
}

int StrToInt(std::string str)
{
    std::stringstream buff(str);
    int tmp;
    buff >> tmp;
    return tmp;
}
*/
// GSL Error handler
/*
int errors;
void ErrHandler(const char * reason,
                        const char * file,
                        int line,
                        int gsl_errno)
{
    errors++;
    if (gsl_errno !=14  )  // 14 = failed to reach tolerance, handle when gsl_int
                  // is called
        std::cerr << file << ":"<< line <<": Error " << errors << ": " <<reason 
            << " (code " << gsl_errno << ")." << std::endl;
    
}
*/


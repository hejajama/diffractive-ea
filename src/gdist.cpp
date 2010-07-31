/*
 * Gluon distribution
 *
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_math.h>
#include "gdist.h"
 
REAL GDist_Toy::gluedist(REAL bjorkx, REAL rsqr)
{
    return 13.7681*gsl_sf_exp(-0.351381*sqrt(rsqr));
}

std::ostream& operator<<(std::ostream& os, GDist& ic)
{
    return os << "GDist object ";
}


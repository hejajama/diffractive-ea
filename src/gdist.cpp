/*
 * Gluon distribution
 *
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_math.h>
#include "gdist.h"
#include <iostream>
 
REAL GDist_Toy::Gluedist(REAL bjorkx, REAL rsqr)
{
    return 15.2426629*gsl_sf_exp(-0.36264237*sqrt(rsqr));
}

std::ostream& operator<<(std::ostream& os, GDist& ic)
{
    return os << "GDist object ";
}


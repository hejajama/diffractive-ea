/*
 * Gluon distribution
 *
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_math.h>
#include "gdist.h"
#include "dipole.h"
#include <iostream>

// Pi^2 / (2*Nc) * alphas(mu(r)^2) xg(x, mu(r))
REAL GDist_Toy::Gluedist(REAL bjorkx, REAL rsqr)
{
    std::cerr << "You are using the toy model of GDist, you really shoudn't!!" << std::endl;
    return SQR(M_PI)/(2*NC)*Alpha_s(Mu2(rsqr))*15.2426629*gsl_sf_exp(-0.36264237*sqrt(rsqr));
}

std::ostream& operator<<(std::ostream& os, GDist& ic)
{
    return os << "GDist object ";
}


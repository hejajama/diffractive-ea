#ifndef GDIST_H
#define GDIST_H

/*
 * Gluon distribution
 *
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */

#include "dipole.h"
#include <iostream>
 
 
class GDist
{
    public:
        // Pi^2 / (2*Nc) * alphas(mu(r)^2) xg(x, mu(r))
        virtual REAL Gluedist(REAL bjorkx, REAL rsqr)=0; 
    private:
};

std::ostream& operator<<(std::ostream& os, GDist& ic);


#endif


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
        virtual REAL Gluedist(REAL bjorkx, REAL rsqr)=0; // = xg(x,µ^2)
    
    private:


};


/*
 * Toy model for gluon distribution
 * Assumes x=1e-4 and gluon density is got by an exponential
 * fit to few xg(r) data points */
 
class GDist_Toy : public GDist
{
    public:
        REAL Gluedist(REAL bjorkx, REAL rsqr);

};

std::ostream& operator<<(std::ostream& os, GDist& ic);


#endif


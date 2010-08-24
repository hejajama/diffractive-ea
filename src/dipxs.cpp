/*
 * Parent class for dipole cross section
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
 
#include "dipxs.h"
#include <iostream>

Dipxs::Dipxs()
{
    InitGsl();
}

Dipxs::Dipxs(Nucleus &nucleus_)
{
    nucleus=nucleus_;
    InitGsl();
}

Dipxs::~Dipxs()
{
    FreeGsl();
}

/*
 * Strong coupling constant as a function of r
 */
REAL Dipxs::Alphas_r(REAL rsqr)
{
    return Alpha_s(Mu2(rsqr));    // Alpha_s(Q^2) defined in dipole.h
}

Nucleus& Dipxs::GetNucleus()
{
    return nucleus;
}

REAL Dipxs::MaxB()
{
    return nucleus.MaxR();
}

REAL Dipxs::Dipxsection(REAL rsqr, REAL xbjork, REAL x, REAL y,
                std::vector<Vec>& nucleons)
{
    Vec tmpvec(x,y);
    return Dipxsection(rsqr, xbjork, tmpvec, nucleons);

}

std::ostream& operator<<(std::ostream& os, Dipxs& ic)
{
    return os << " Dipxs object ";

}




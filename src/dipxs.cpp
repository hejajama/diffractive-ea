/*
 * Parent class for dipole cross section
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
 
#include "dipxs.h"

DipXS::DipXS(Nucleus nucleus_)
{
    nucleus=nucleus_;
}

/*
 * Strong coupling constant as a function of r
 */
REAL DipXS::Alphas_r(REAL r)
{
    return Alpha_s(Mu2(SQR(r)));    // Alpha_s(Q^2) defined in dipole.h
}

REAL DipXS::Mu2(REAL rsqr)
{
    return 4/rsqr + mu2_0;
}

REAL DipXS::Sigmap(REAL r, REAL xbjork)
{
    return SQR(M_PI)/NC*SQR(r)*Alphas_r(r)
        *nucleus.GetGDist()->gluedist(xbjork,SQR(r));
}

std::ostream& operator<<(std::ostream& os, DipXS& ic)
{
    return os << " DipXS object ";

}



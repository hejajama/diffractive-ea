/*
 * Parent class for dipole cross section
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
 
#include "dipxs.h"

DipXS::DipXS()
{
    InitGsl();
}

DipXS::DipXS(Nucleus nucleus_)
{
    nucleus=nucleus_;
    InitGsl();
}

DipXS::~DipXS()
{
    FreeGsl();
}

/*
 * Strong coupling constant as a function of r
 */
REAL DipXS::Alphas_r(REAL rsqr)
{
    return Alpha_s(Mu2(rsqr));    // Alpha_s(Q^2) defined in dipole.h
}

REAL DipXS::Mu2(REAL rsqr)
{
    return 4/rsqr + mu2_0;
}

REAL DipXS::Sigmap(REAL rsqr, REAL xbjork)
{
    return SQR(M_PI)/NC*rsqr*Alphas_r(rsqr)
        *nucleus.GetGDist()->gluedist(xbjork,rsqr);
}

Nucleus& DipXS::GetNucleus()
{
    return nucleus;
}

REAL DipXS::MaxB()
{
    return nucleus.MaxR();
}

REAL DipXS::DipXSection(REAL rsqr, REAL xbjork, REAL x, REAL y,
                std::vector<Vec>& nucleons)
{
    Vec tmpvec(x,y);
    return DipXSection(rsqr, xbjork, tmpvec, nucleons);

}

std::ostream& operator<<(std::ostream& os, DipXS& ic)
{
    return os << " DipXS object ";

}




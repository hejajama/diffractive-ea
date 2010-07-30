/*
 * Dipole cross section in IP non sat model
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
#include <gsl/gsl_math.h>
#include <iostream>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>

#include "dipxs_ipnonsat.h"

DipXS_IPNonSat::DipXS_IPNonSat(Nucleus nucleus_) :
    DipXS(nucleus_)
{
    // Nothing to do here
}

/*
 * Dipole cross section in impact parameter reprsesentation
 * Averaged over all initial configurations
 * NOTE: Averaging could be done analytically, but it is here done
 * numerically to test the algorithm
 */
REAL DipXS_IPNonSat::DipXSection_b(REAL r, REAL xbjork, Vec b )
{
    // TODO: Only for dipole-proton currently
    return DipXSection_b(r,xbjork,b.Len());
    //return SQR(M_PI)/NC*SQR(r)*Alphas_r(r)*nucleus.GetGDist()->gluedist(xbjork,SQR(r))*nucleus.Tp(b);
}

REAL DipXS_IPNonSat::DipXSection_b (REAL r, REAL xbjork, REAL b )
{
    if (nucleus.GetA()==1) // Dipole-proton
        return SQR(M_PI)/NC*SQR(r)*Alphas_r(r)*nucleus.GetGDist()->gluedist(xbjork,SQR(r))*nucleus.Tp(b);

    /* Dipole-nuclues, so we have to average over all possible nucleon configurations
     * Monte Carlo method is used here as it should converge quite rapidly, all
     * nucleons should be quite similar when nucleon positions are generated from
     * Woods-Saxon distribution
     */
     
     int A = nucleus.GetA(); // Dimension of the integration
     
        

}


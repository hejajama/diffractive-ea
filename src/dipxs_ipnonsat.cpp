/*
 * Dipole cross section in IP non sat model
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
#include <gsl/gsl_math.h>
#include <iostream>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_sf_exp.h>
#include <vector>

#include "dipxs_ipnonsat.h"

using std::cout; using std::endl;

DipXS_IPNonSat::DipXS_IPNonSat(Nucleus nucleus_) :
    DipXS(nucleus_)
{
    // Nothing to do here
}

/*
 * Dipole cross section in impact parameter reprsesentation
 * Averaged over all initial configurations
 */
REAL DipXS_IPNonSat::DipXSection_b(REAL r, REAL xbjork, REAL b )
{
    // Dipol-proton only for now
    return DipXSection_b(r,xbjork,b);
}

REAL DipXS_IPNonSat::DipXSection_b (REAL r, REAL xbjork, Vec b )
{
    if (nucleus.GetA()==1) // Dipole-proton
        return Sigmap(r,xbjork)*nucleus.Tp(b.Len());

    /* Dipole-nuclues, so we have to average over all possible nucleon
     * configurations. Monte Carlo method is used here as it 
     * should converge quite rapidly, all nuclei should be quite similar
     * when nucleon positions are generated from
     * Woods-Saxon distribution
     */
     
     int A = nucleus.GetA(); // Dimension of the integration
     
     
}

/*
 * Non-averaged dipole cross section as a function of
 * nucleon transversial positions 
 */
REAL DipXS_IPNonSat::DipXSectoin_b_nonaveraged(REAL r, REAL xbjork, Vec b, 
            std::vector<Vec> nucleons)
{
    if (nucleons.size() != nucleus.GetA())
    {
        std::cerr << "Got list of " << nucleons.size() << " nucleons but "
            << "there should be " << nucleus.GetA() << " of them...";
            return -1;
    }
    REAL result=0;
    for (int i=0; i<nucleons.size(); i++)
    {
        Vec tmp = b-nucleons[i];
//        cout << "tmp.LenSqr(): " << 
        result += gsl_sf_exp(-tmp.LenSqr()/(2*B_p));
    }
    result *= Sigmap(r,xbjork)/(2*M_PI*B_p);
    return result;
}


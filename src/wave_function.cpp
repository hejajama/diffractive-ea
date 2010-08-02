/*
 * General class for wave functions
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
 
#include "wave_function.h"

REAL WaveFunction::PsiSqr_tot(REAL Qsqr, REAL r, REAL z)
{
    return PsiSqr_T(Qsqr,r,z)+PsiSqr_L(Qsqr,r,z);
}

REAL WaveFunction::PsiSqr_tot_intz(REAL Qsqr, REAL r)
{
    return PsiSqr_T_intz(Qsqr,r)+PsiSqr_L_intz(Qsqr,r);
}
 

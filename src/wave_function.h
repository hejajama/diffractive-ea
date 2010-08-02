#ifndef WAVE_FUNCTION_H
#define WAVE_FUNCTION_H

/*
 * General class for wave functions
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */

#include "dipole.h"

class WaveFunction{
    public:
        virtual REAL PsiSqr_T(REAL Qsqr, REAL r, REAL z) = 0;
        virtual REAL PsiSqr_L(REAL Qsqr, REAL r, REAL z) = 0;
        virtual REAL PsiSqr_T_intz(REAL Qsqr, REAL r) = 0;
        virtual REAL PsiSqr_L_intz(REAL Qsqr, REAL r) = 0;
        REAL PsiSqr_tot(REAL Qsqr, REAL r, REAL z);
        REAL PsiSqr_tot_intz(REAL Qsqr, REAL r);
    private:

};

#endif  // WAVE_FUNCTION_H

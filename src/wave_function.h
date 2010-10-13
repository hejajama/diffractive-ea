#ifndef WAVE_FUNCTION_H
#define WAVE_FUNCTION_H

/*
 * General class for wave functions
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */

#include "dipole.h"
#include "string"

const int VM_MODE_TOT=1; const int VM_MODE_L=2; const int VM_MODE_T=3;

class WaveFunction{
    public:
        WaveFunction();
        virtual REAL PsiSqr_T(REAL Qsqr, REAL r, REAL z) = 0;
        virtual REAL PsiSqr_L(REAL Qsqr, REAL r, REAL z) = 0;
        virtual REAL PsiSqr_T_intz(REAL Qsqr, REAL r) = 0;
        virtual REAL PsiSqr_L_intz(REAL Qsqr, REAL r) = 0;
        virtual std::string GetParamString()=0;
        REAL PsiSqr_tot(REAL Qsqr, REAL r, REAL z);
        REAL PsiSqr_tot_intz(REAL Qsqr, REAL r);
        REAL PsiSqr_intz(REAL Qsqr, REAL r);
        void SetMode(int m);
    private:
        int mode;   // What to return when PsiSqr_intz
};

#endif  // WAVE_FUNCTION_H

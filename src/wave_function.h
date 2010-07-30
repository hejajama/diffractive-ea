#ifndef WAVE_FUNCTION_H
#define WAVE_FUNCTION_H

/*
 * General class for wave functions
 */

class WaveFunction{
    public:
        virtual REAL PsiSqr_T(REAL Qsqr, REAL r, REAL z) = 0;
        virtual REAL PsiSqr_L(REAL Qsqr, REAL r, REAL z) = 0;
        virtual REAL PsiSqr_T_intz(REAL Qsqr, REAL r) = 0;
        virtual REAL PsiSqr_L_intz(REAL Qsqr, REAL r) = 0;
    private:

};

#endif  // WAVE_FUNCTION_H

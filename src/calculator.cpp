/*
 * Calculates different cross sections
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
#include "calculator.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <iostream>


Calculator::Calculator(Dipxs* amplitude_, WaveFunction* wavef_)
{
    amplitude=amplitude_;
    wavef=wavef_;
}

/*
 * Differential cross section d\sigma / dt
 */

REAL Calculator::CrossSection_dt(REAL t, REAL Qsqr, REAL bjorkx)
{
    gsl_function fun;
    inthelper_r inthelp;
    inthelp.amplitude=amplitude;
    inthelp.vm=wavef; inthelp.bjorkx=bjorkx;
    inthelp.delta=sqrt(t); inthelp.Qsqr=Qsqr;
    fun.function=&inthelperf_r1;
    fun.params=&inthelp;
        
    REAL result,abserr; size_t eval;
    
    int status = gsl_integration_qng(&fun, MINR, MAXR, RINTACCURACY, RINTACCURACY, 
        &result, &abserr, &eval);
    if (status) std::cerr << "Error " << status << " at " << __FILE__ << ":"
        << __LINE__ << ": Result " << result << ", abserror: " << abserr 
        << " (t=" << t <<")" << std::endl;
        
    result *= 1.0/(16.0*M_PI);
    return result;
}

/*
 * Total cross section
 * Calculated as \int dt d\sigma/dt
 */
REAL Calculator::TotalCrossSection(REAL Qsqr, REAL bjorkx)
{
    gsl_function fun;   
    inthelper_r inthelp;
    inthelp.amplitude=amplitude;
    inthelp.vm=wavef; inthelp.bjorkx=bjorkx;
    inthelp.delta=-1; inthelp.Qsqr=Qsqr;
    inthelp.calculator=this;
    fun.function=&inthelper_totxs;
    fun.params=&inthelp; 
        
    REAL result,abserr; size_t eval;
    int status = gsl_integration_qng(&fun, 0, TOTXS_MAXT, 0, RINTACCURACY, 
        &result, &abserr, &eval);
    if (status)
        std::cerr << "Total cross section integral failed to reach tolerance: "
        << "Result: " << result << ", abserr: " << abserr << std::endl;
    
    return result;
    
}


/*
 * Integral helpers
 */
 
REAL inthelper_totxs(REAL t, void* p)
{
    inthelper_r* par = (inthelper_r*)p;
    return par->calculator->CrossSection_dt(t, par->Qsqr, par->bjorkx);

}

REAL inthelperf_r1(REAL r, void* p)
{
    inthelper_r* par = (inthelper_r*)p;
    par->r=r;
    gsl_function fun;
    fun.function=&inthelperf_r2;
    fun.params=par;
    
    REAL result,abserr; size_t eval;
    
    int status = gsl_integration_qng(&fun, MINR, MAXR, 0, RINTACCURACY, 
        &result, &abserr, &eval);
    
    if (status) std::cerr << "Error " << status << " at " << __FILE__ << ":"
        << __LINE__ << ": Result " << result << ", abserror: " << abserr 
        << " (t=" << par->delta*par->delta << ")" << std::endl;
    
    return 2*M_PI*r*par->vm->PsiSqr_tot_intz(par->Qsqr, r)*result;
}

// Inner r' integral: \int dr' r' (jpsi)(r') * qqamplitude_avg_sqr(r,r')
REAL inthelperf_r2(REAL r2, void* p)
{
    inthelper_r* par = (inthelper_r*)p;
    return 2*M_PI*r2 * par->vm->PsiSqr_tot_intz(par->Qsqr, r2) 
            * par->amplitude->Dipxsection_sqr_avg(SQR(par->r), SQR(r2), 
                    par->bjorkx, par->delta);

}
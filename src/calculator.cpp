/*
 * Calculates different cross sections
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */

#include "dipole.h"
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
    if (status and result>0.000001) 
    std::cerr << "Error " << status << " at " << __FILE__ << ":"
        << __LINE__ << ": Result " << result << ", abserror: " << abserr 
        << " (t=" << t <<")" << std::endl;
    
    // Factor 4 as we integrate amplitude, not d\sigma/d^2b = 2A
    result *= 4.0/(16.0*M_PI);
    return result;
}

/*
 * Differential cross section for coherent \gamma^*A scattering
 */
 
// Only one r integral for cohernet dipole-nucleus scattering
REAL inthelperf_coherent(REAL r, void* p)
{
    inthelper_r* par = (inthelper_r*)p;
    return 2*M_PI*r*par->vm->PsiSqr_intz(par->Qsqr, r)
        * par->amplitude->CoherentDipoleAmplitude_avg(SQR(r), par->bjorkx,
         par->delta);
}

REAL Calculator::CoherentCrossSection_dt(REAL t, REAL Qsqr, REAL bjorkx)
{
    gsl_function fun;
    inthelper_r inthelp;
    inthelp.amplitude=amplitude;
    inthelp.vm=wavef; inthelp.bjorkx=bjorkx;
    inthelp.delta=sqrt(t); inthelp.Qsqr=Qsqr;
    fun.function=&inthelperf_coherent;
    fun.params=&inthelp;
        
    REAL result,abserr; size_t eval;
    
    int status = gsl_integration_qng(&fun, MINR, MAXR, 0, RINTACCURACY, 
        &result, &abserr, &eval);
    if (status and result>0.00001) 
    std::cerr << "Error " << status << " at " << __FILE__ << ":"
        << __LINE__ << ": Result " << result << ", abserror: " << abserr 
        << " (t=" << t <<")" << std::endl;
        
    // Multiply by 4 as CoherentDipxsection_avg returns an amplitude which
    // must be multiplied by 2 and we square the result at the end
    result *= result;
    result *= 4.0/(16.0*M_PI);  
    return result;

}

/*
 * Differential cross section for dipole-proton scattering
 */
REAL Calculator::ProtonCrossSection_dt(REAL t, REAL Qsqr, REAL bjorkx)
{
    gsl_function fun;   
    inthelper_r inthelp;
    inthelp.amplitude=amplitude;
    inthelp.vm=wavef; inthelp.bjorkx=bjorkx;
    inthelp.delta=sqrt(t); inthelp.Qsqr=Qsqr;
    inthelp.calculator=this; inthelp.analytic_t=false;
    fun.function=&inthelperf_proton;
    fun.params=&inthelp; 
        
    REAL result,abserr; size_t eval;
    int status = gsl_integration_qng(&fun, MINR, MAXR, 0, RINTACCURACY, 
        &result, &abserr, &eval);
    if (status and result>0.000001) 
    std::cerr << "Error " << status << " at " << __FILE__ << ":"
        << __LINE__ << ": Result " << result << ", abserror: " << abserr 
        << " (t=" << t <<", Q^2=" << Qsqr<<", x=" << bjorkx <<")" << std::endl;
    
    // Multiply by 4 as we used amplitude, not d\sigma/d^2b = 2A
    return 4*result*result/(16.0*M_PI);

}

// Only one r integral for dipole-proton scattering
REAL inthelperf_proton(REAL r, void* p)
{
    inthelper_r* par = (inthelper_r*)p;
    if (par->analytic_t)
    {   // We are calculating total xs and performing t integral analytically
        //return 2.0*M_PI*r*par->vm->PsiSqr_intz(par->Qsqr, r)
        // * par->amplitude->Dipxsection_proton(SQR(r), par->bjorkx);
        std::cerr << "Analytic integration over t is not implemented " 
            << " (inthelperf_proton)" << std::endl;
    
    }
    return 2*M_PI*r*par->vm->PsiSqr_intz(par->Qsqr, r)
        * par->amplitude->DipoleAmplitude_proton(SQR(r), par->bjorkx,
         par->delta);
       
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
    int status = gsl_integration_qng(&fun, 0, TOTXS_MAXT, 0, TINTACCURACY, 
        &result, &abserr, &eval);
    if (status and result>0.00001)
        std::cerr << "Total cross section integral failed to reach tolerance: "
        << "Result: " << result << ", abserr: " << abserr << std::endl;
    
    return result;
    
}

// Integration over t 
REAL inthelper_totxs(REAL t, void* p)
{
    inthelper_r* par = (inthelper_r*)p;
    if (par->calculator->GetAmplitude()->GetNucleus().GetA()==1) // Proton
        return par->calculator->ProtonCrossSection_dt(t, par->Qsqr, par->bjorkx); 
    return par->calculator->CrossSection_dt(t, par->Qsqr, par->bjorkx);

}

/*
 * Total cross section for dipole-proton scattering
 * One can also use TotalCrossSection with A=1, but this method performs the 
 * integration over |t| analytically.
 */
REAL Calculator::TotalProtonCrossSection(REAL Qsqr, REAL bjorkx)
{
    std::cerr << "Proton Cross Section calculation may not work!" << std::endl;
    /*gsl_function fun;
    inthelper_r inthelp;
    inthelp.amplitude=amplitude;
    inthelp.vm=wavef; inthelp.bjorkx=bjorkx;
    inthelp.delta=-1; inthelp.Qsqr=Qsqr;
    inthelp.analytic_t=true;
    fun.function=&inthelperf_proton;
    fun.params=&inthelp;
        
    REAL result,abserr; size_t eval;
    
    int status = gsl_integration_qng(&fun, MINR, MAXR, RINTACCURACY, RINTACCURACY, 
        &result, &abserr, &eval);
    if (status) std::cerr << "Error " << status << " at " << __FILE__ << ":"
        << __LINE__ << ": Result " << result << ", abserror: " << abserr 
        << " (total proton cross section)" << std::endl;
        
    return result*result/(16.0*M_PI);    
    */
    return 0;
    

}
/*
 * Integral helpers
 */



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
    
    if (status and result>0.000001) 
        std::cerr << "Error " << status << " at " << __FILE__ << ":"
        << __LINE__ << ": Result " << result << ", abserror: " << abserr 
        << " (t=" << par->delta*par->delta << ")" << std::endl;
    
    return 2*M_PI*r*par->vm->PsiSqr_intz(par->Qsqr, r)*result;
}

// Inner r' integral: \int dr' r' (jpsi)(r') * qqamplitude_avg_sqr(r,r')
REAL inthelperf_r2(REAL r2, void* p)
{
    inthelper_r* par = (inthelper_r*)p;
    return 2*M_PI*r2 * par->vm->PsiSqr_intz(par->Qsqr, r2) 
            * par->amplitude->DipoleAmplitude_sqr_avg(SQR(par->r), SQR(r2), 
                    par->bjorkx, par->delta);

}

// Return pointer to dipole amplitude
Dipxs* Calculator::GetAmplitude()
{
    return amplitude;
}


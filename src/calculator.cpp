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
#include <gsl/gsl_sf_gamma.h>
#include <iostream>

// Integration settings
const REAL MAXR=20;
const REAL MINR=0.0001;   // r=0 doesn't work, K_{0,1}(0)=inf
const REAL RINTACCURACY=0.003;
const REAL TINTACCURACY=0.001;
const REAL TOTXS_MAXT=2;  // Max |t| in GeV^2 when calculating total xs
const REAL epsfact = 4.0;   // eps = x_pom/epsfact

Calculator::Calculator(Dipxs* amplitude_, WaveFunction* wavef_)
{
    amplitude=amplitude_;
    wavef=wavef_;
    cached_corrections=false; cache_Q2=-1;
    corrections=true;
}

/*
 * Differential cross section d\sigma / dt
 *
 * Here we assume that the corrections due to the real part of the
 * scattering amplitude and the skewedness effect are exactly same
 * as in case of dipole-proton scattering.
 *
 * As we use only factorized IPsat and IIM models here, these correction
 * factors do not depend on r -> we cache them
 *
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
    result=0;
    REAL eps = bjorkx/epsfact;
    if (polarization == VM_MODE_TOT)
    {
        wavef->SetMode(VM_MODE_T);
        int status = gsl_integration_qng(&fun, MINR, MAXR, 0, RINTACCURACY, 
            &result, &abserr, &eval);
        if (status and result>0.000001) 
        std::cerr << "Error " << status << " at " << __FILE__ << ":"
            << __LINE__ << ": Result " << result << ", abserror: " << abserr 
            << " (t=" << t <<")" << std::endl;
        
        // Calculate real part and skewedness corrections or use cahced ones
        //if (!cached_corrections or Qsqr != cache_Q2)
        //{
            polarization=VM_MODE_T;
            REAL xs = RIntAmplitude(t, Qsqr, bjorkx, &inthelperf_proton);
            REAL xseps = RIntAmplitude(t, Qsqr, bjorkx+eps, &inthelperf_proton);
            REAL lambda = log(xs/xseps)*(bjorkx/eps);
            betasqr_t = SQR(Beta(lambda));
            rgsqr_t = SQR(Rg(lambda));
            
            //cache_Q2=Qsqr; cached_corrections=true;
        //}
        result = result*(1.0+betasqr_t)*rgsqr_t;
        
        REAL tmpres;    
        wavef->SetMode(VM_MODE_L);
        status = gsl_integration_qng(&fun, MINR, MAXR, 0, RINTACCURACY, 
            &tmpres, &abserr, &eval);
        
        if (status and tmpres>0.000001) 
        std::cerr << "Error " << status << " at " << __FILE__ << ":"
            << __LINE__ << ": Result " << tmpres << ", abserror: " << abserr 
            << " (t=" << t <<")" << std::endl;
        
        // Calculate real part and skewedness corrections or use cahced ones
        //if (!cached_corrections or Qsqr != cache_Q2)
        //{
            polarization=VM_MODE_L;
            xs = RIntAmplitude(t, Qsqr, bjorkx, &inthelperf_proton);
            xseps = RIntAmplitude(t, Qsqr, bjorkx+eps, &inthelperf_proton);
            lambda = log(xs/xseps)*(bjorkx/eps);
            betasqr_l = SQR(Beta(lambda));
            rgsqr_l = SQR(Rg(lambda));
            
            //cache_Q2=Qsqr; cached_corrections=true;
        //}
        tmpres = tmpres*(1.0+betasqr_l)*rgsqr_l;
        result = result + tmpres;

        polarization=VM_MODE_TOT;
        
    }
    else    // Correct polarization is already set
    {
        int status = gsl_integration_qng(&fun, MINR, MAXR, 0, RINTACCURACY, 
            &result, &abserr, &eval);
        if (status and result>0.000001) 
        std::cerr << "Error " << status << " at " << __FILE__ << ":"
            << __LINE__ << ": Result " << result << ", abserror: " << abserr 
            << " (t=" << t <<")" << std::endl;
            
        // Calculate real part and skewedness corrections or use cahced ones
        // We use betasqr_l and rgsqr_l even though we don't know wheter we are
        // calculating diff xs for transverse or longitudinal scattering
        //if (!cached_corrections or Qsqr != cache_Q2)
        //{
            REAL xs = RIntAmplitude(t, Qsqr, bjorkx, &inthelperf_proton);
            REAL xseps = RIntAmplitude(t, Qsqr, bjorkx+eps, &inthelperf_proton);
            REAL lambda = log(xs/xseps)*(bjorkx/eps);
            betasqr_l = SQR(Beta(lambda));
            rgsqr_l = SQR(Rg(lambda));
            cache_Q2=Qsqr; cached_corrections=true;
        //}
        result *= (1.0 + betasqr_l)*rgsqr_l;
    }
    

    
    
    // Factor 4 as we integrate amplitude, not d\sigma/d^2b = 2A
    result *= 4.0/(16.0*M_PI);
    return result;
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
    
    return 2.0*M_PI*r*par->vm->PsiSqr_intz(par->Qsqr, r)*result;
}

// Inner r' integral: \int dr' r' (jpsi)(r') * qqamplitude_avg_sqr(r,r')
REAL inthelperf_r2(REAL r2, void* p)
{
    inthelper_r* par = (inthelper_r*)p;
    return 2.0*M_PI*r2 * par->vm->PsiSqr_intz(par->Qsqr, r2) 
            * par->amplitude->DipoleAmplitude_sqr_avg(SQR(par->r), SQR(r2), 
                    par->bjorkx, par->delta);

}

/*
 * Differential cross section for coherent \gamma^*A scattering
 * We calculate correction factors from dipole-proton amplitude!
 * See ProtonCrossSection_dt()
 */
 
REAL Calculator::CoherentCrossSection_dt(REAL t, REAL Qsqr, REAL bjorkx)
{      
    REAL result,abserr; size_t eval;
    REAL eps = bjorkx/epsfact;
    result=0;
    if (polarization == VM_MODE_TOT)
    {
        wavef->SetMode(VM_MODE_T);
        polarization=VM_MODE_T;
        REAL tmpres = RIntAmplitude(t, Qsqr, bjorkx, &inthelperf_coherent);
        REAL xs = RIntAmplitude(t, Qsqr, bjorkx, &inthelperf_proton);
        REAL xseps = RIntAmplitude(t, Qsqr, bjorkx+eps, &inthelperf_proton);
        REAL lambda = log(xs/xseps)*(bjorkx/eps);
        result = tmpres*tmpres*(1.0+SQR(Beta(lambda)))*SQR(Rg(lambda));
            
        wavef->SetMode(VM_MODE_L);
        polarization=VM_MODE_L;
        tmpres = RIntAmplitude(t, Qsqr, bjorkx, &inthelperf_coherent);
        xs = RIntAmplitude(t, Qsqr, bjorkx, &inthelperf_proton);
        xseps = RIntAmplitude(t, Qsqr, bjorkx+eps, &inthelperf_proton);
        lambda = log(xs/xseps)*(bjorkx/eps);
        result += tmpres*tmpres*(1.0+SQR(Beta(lambda)))*SQR(Rg(lambda)); 
        polarization=VM_MODE_TOT;
         // Sum d\sigma/dt_L + d\sgima/dt_T, not amplitudes!
    }
    else    // Correct polarization is already set
    {
        result = RIntAmplitude(t, Qsqr, bjorkx, &inthelperf_coherent);
        REAL xs = RIntAmplitude(t, Qsqr, bjorkx, &inthelperf_proton);
        REAL xseps = RIntAmplitude(t, Qsqr, bjorkx+eps, &inthelperf_proton);
        REAL lambda = log(xs/xseps)*(bjorkx/eps);
        result=result*result*(1.0+SQR(Beta(lambda)))*SQR(Rg(lambda));
    }
        
    result *= 1.0/(16.0*M_PI);  
    return result;
}

// Only one r integral for cohernet dipole-nucleus scattering
REAL inthelperf_coherent(REAL r, void* p)
{
    inthelper_r* par = (inthelper_r*)p;
    return 2*M_PI*r*par->vm->PsiSqr_intz(par->Qsqr, r)
        * par->amplitude->CoherentDipoleAmplitude_avg(SQR(r), par->bjorkx,
         par->delta);
}

/*
 * Differential cross section for dipole-proton scattering
 * Result is easy to calculate, just integrate amplitude->DipoleAmplitudeProton
 * over r. However, as we have neglected the real part of the scattering
 * amplitude and because there is the skewedness effect, we have to multiply
 * the result by two factors Rg(lambda) and Beta(lambda), see
 * arXiv:0712.2670
 */
REAL Calculator::ProtonCrossSection_dt(REAL t, REAL Qsqr, REAL bjorkx)
{
    //NOTE: This could be optimized as lambda does not depend
    // on t in factorized IPsat or in IIM model, but in case of 
    // nonfactorized IPsat there is quite significant t-dependence
    REAL eps = bjorkx/epsfact;
    REAL result=0;
        
    if (polarization == VM_MODE_TOT)
    {
        wavef->SetMode(VM_MODE_T);
        REAL xs = RIntAmplitude(t, Qsqr, bjorkx, &inthelperf_proton);
        REAL xseps = RIntAmplitude(t, Qsqr, bjorkx+eps, &inthelperf_proton);
        REAL lambda = log(xs/xseps)*(bjorkx/eps);
        
        result=xs*xs*(1+SQR(Beta(lambda)) ) * SQR(Rg(lambda));
        
        wavef->SetMode(VM_MODE_L);
        xs = RIntAmplitude(t, Qsqr, bjorkx, &inthelperf_proton);
        xseps = RIntAmplitude(t, Qsqr, bjorkx+eps, &inthelperf_proton);
        lambda = log(xs/xseps)*(bjorkx/eps);

        result+= xs*xs*(1+SQR(Beta(lambda)) ) * SQR(Rg(lambda));
        
    }
    else     // Correct polarization is already set
    {
        REAL xs = RIntAmplitude(t, Qsqr, bjorkx, &inthelperf_proton);
        REAL xseps = RIntAmplitude(t, Qsqr, bjorkx+eps, &inthelperf_proton);        
        REAL lambda = log(xs/xseps)*(bjorkx/eps);
        if (Beta(lambda)>1)
        {
            std::cerr << "beta=" << Beta(lambda) << " at t=" << t 
                << std::endl;
        
        }
        result=xs*xs*(1+SQR(Beta(lambda)) ) * SQR(Rg(lambda));
	    //std::cout << Qsqr << " " << xs/xseps << std::endl;
	    //bjorkx << " " << eps << " " << lambda << std::endl;
        //std::cout << Qsqr << " " << (1+SQR(Beta(lambda)) ) * SQR(Rg(lambda)) << std::endl;
   }
    
     
    return result/(16.0*M_PI);
}


// Only one r integral for dipole-proton scattering
REAL inthelperf_proton(REAL r, void* p)
{
    inthelper_r* par = (inthelper_r*)p;
    return 2*M_PI*r*par->vm->PsiSqr_intz(par->Qsqr, r)
        * par->amplitude->DipoleAmplitude_proton(SQR(r), par->bjorkx,
         par->delta);
       
}


/*
 * Integrate dipole amplitude over r
 * To be used with dipole-proton and coherent dipole-nucleus amplitudes,
 * as in these cases we only have to perform one r integral and square the
 * result at the end.
 * inthelperf* is a pointer to the integral helper function, which should return
 * 2 pi r \int dz/(4 pi) \Psi^*\Psi(Q,r) * 1/2 d\sigma/d^2 b
 * Note: there is factor 1/2, and we multiply the result by 2 at the end
 * of this method
 */ 
REAL Calculator::RIntAmplitude(REAL t, REAL Qsqr, REAL bjorkx, 
    REAL(*inthelperf)(REAL r, void* p) )
{
    gsl_function fun;   
    inthelper_r inthelp;
    inthelp.amplitude=amplitude;
    inthelp.vm=wavef; inthelp.bjorkx=bjorkx;
    inthelp.delta=sqrt(t); inthelp.Qsqr=Qsqr;
    inthelp.calculator=this;
    fun.function=inthelperf;
    fun.params=&inthelp;
        
    REAL result,abserr; size_t eval;
    result=0;
    int status = gsl_integration_qng(&fun, MINR, MAXR, 0, RINTACCURACY, 
        &result, &abserr, &eval);
    if (status and result>0.000001) 
        std::cerr << "Error " << status << " at " << __FILE__ << ":"
            << __LINE__ << ": Result " << result << ", abserror: " << abserr 
            << " (t=" << t <<")" << std::endl;
    return result*2.0;  // Multiply by 2 as DipoleAmplitude returns 
                        // 1/2 * d\sigma/d^2 b 
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
 * Skewedness factor
 * = 2^(2*\lambda+3)/sqrt(pi) * Gamma(lambda+5/2)/Gamma(lambda+4)
 * see e.g. arXiv:0712.2670v2
 */
REAL Calculator::Rg(REAL lambda)
{
    if (!corrections)
        return 1.0;
    REAL result=pow(2.0,2.0*lambda+3.0)/sqrt(M_PI);
    result *= gsl_sf_gamma(lambda+5.0/2.0)/gsl_sf_gamma(lambda+4.0);
    return result;
}

REAL Calculator::Beta(REAL lambda)
{
    if (!corrections)
        return 0;
    return tan(M_PI*lambda/2.0);
}

// Return pointer to dipole amplitude
Dipxs* Calculator::GetAmplitude()
{
    return amplitude;
}

void Calculator::SetPolarization(int pol)
{
    polarization=pol;
    wavef->SetMode(pol);
}

void Calculator::SetCorrections(bool c)
{
    corrections=c;
}



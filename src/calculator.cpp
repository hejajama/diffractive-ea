/*
 * Calculates different cross sections
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010-2014
 */

#include "dipole.h"
#include "calculator.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_bessel.h>
#include <iostream>

// Integration settings
const REAL MAXR=20;
const REAL MINR=0.0001;   // r=0 doesn't work, K_{0,1}(0)=inf
const REAL RINTACCURACY=0.01;
const REAL TINTACCURACY=0.001;
const REAL DEFAULTACCURACY=0.001;
const REAL TOTXS_MAXT=2;  // Max |t| in GeV^2 when calculating total xs
const REAL epsfact = 50.0;   // eps = x_pom/epsfact
const REAL TOTXS_COHERENT_MINT  = 0.0;  // Min t for total cohernet xs integral
const REAL TOTXS_COHERENT_MAXT = 0.1;

const double eps_y = 0.07  ;	// when computing corrections amplitude
							// is evaluated at y and y+eps_y

using std::cout; using std::cerr; using std::endl;

Calculator::Calculator(Dipxs* amplitude_, WaveFunction* wavef_)
{
    amplitude=amplitude_;
    wavef=wavef_;
    cached_corrections=false; cache_Q2=-1;
    corrections=true;
    relaccuracy_t=DEFAULTACCURACY;
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
		cerr <<"# Please improve numerics here... (derivative computation)" << endl;
        wavef->SetMode(VM_MODE_T);
        int status = gsl_integration_qng(&fun, MINR, MAXR, 0, RINTACCURACY, 
            &result, &abserr, &eval);
        if (status and result>0.000001) 
        std::cerr << "Error " << status << " at " << __FILE__ << ":"
            << __LINE__ << ": Result " << result << ", abserror: " << abserr 
            << " (t=" << t <<")" << std::endl;
        
        // Calculate real part and skewedness corrections
        REAL xs, xseps, lambda;
        betasqr_t=0; rgsqr_t=1;
        if (corrections)
        {
            xs = RIntAmplitude(t, Qsqr, bjorkx, &inthelperf_proton);
            xseps = RIntAmplitude(t, Qsqr, bjorkx+eps, &inthelperf_proton);
            lambda = log(xs/xseps)*(bjorkx/eps);
            betasqr_t = SQR(Beta(lambda));
            rgsqr_t = SQR(Rg(lambda));
        }
        result = result*(1.0+betasqr_t)*rgsqr_t;
        
        REAL tmpres;    
        wavef->SetMode(VM_MODE_L);
        status = gsl_integration_qng(&fun, MINR, MAXR, 0, RINTACCURACY, 
            &tmpres, &abserr, &eval);
        
        if (status and tmpres>0.000001) 
        std::cerr << "Error " << status << " at " << __FILE__ << ":"
            << __LINE__ << ": Result " << tmpres << ", abserror: " << abserr 
            << " (t=" << t <<")" << std::endl;
        
        // Calculate real part and skewedness corrections
        betasqr_l=0; rgsqr_l=1.0;
        if (corrections)
        {
            xs = RIntAmplitude(t, Qsqr, bjorkx, &inthelperf_proton);
            xseps = RIntAmplitude(t, Qsqr, bjorkx+eps, &inthelperf_proton);
            lambda = log(xs/xseps)*(bjorkx/eps);
            betasqr_l = SQR(Beta(lambda));
            rgsqr_l = SQR(Rg(lambda));
            
            //cache_Q2=Qsqr; cached_corrections=true;
        }
        tmpres = tmpres*(1.0+betasqr_l)*rgsqr_l;
        result = result + tmpres;
        
        wavef->SetMode(VM_MODE_TOT);
        
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
        rgsqr_l=1.0; betasqr_l=0;
        if (corrections)
        {
            /*REAL xs_old = RIntAmplitude(t, Qsqr, bjorkx, &inthelperf_proton);
            REAL xseps_old = RIntAmplitude(t, Qsqr, bjorkx+eps, &inthelperf_proton);
            REAL lambda_old = log(xs_old/xseps_old)*(bjorkx/eps);*/
            
            double y0 = log(1.0/bjorkx);

            double y1 = y0-eps_y;
            double y2 = y0 + eps_y;
            double x1 = exp(-y1);
            double x2 = exp(-y2);
			// NOTE: there could be problem if xs1 and xs2 have different sign which is 
			// possible in principle
            double xs1 = RIntAmplitude(t, Qsqr, x1, &inthelperf_proton);
            double xs2 = RIntAmplitude(t, Qsqr, x2, &inthelperf_proton);
            double lambda = log(xs2 / xs1) / (2.0*eps_y);
            rgsqr_l = SQR(Rg(lambda));
            betasqr_l=SQR(Beta(lambda));
            cache_Q2=Qsqr; cached_corrections=true;
            //cout <<"# x=" << bjorkx << " realpartcorr " << betasqr_l <<" skewedness " << rgsqr_l << endl;
        }
        result *= (1.0 + betasqr_l)*rgsqr_l;
    }
    

    
    
    // Factor 4 as we integrate 1/2*d\sigma/d2b
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
	if (amplitude->GetNucleus().GetA()==1)
	{
		cerr <<"Coherent cross section is not defined for the proton! " << endl;
		return 0;
	}
	
    REAL result,abserr, xs, xseps, lambda; size_t eval;
    REAL eps = std::min(bjorkx/epsfact,0.0001);
    result=0;
    
    // Notice: CoherentDipoleAmplitude returns 1/2*d\sigma/dt, but it is 
    // multiplied by 2 in RIntAmplitude()
    
    if (polarization == VM_MODE_TOT)
    {
        wavef->SetMode(VM_MODE_T);
        REAL tmpres = RIntAmplitude(t, Qsqr, bjorkx, &inthelperf_coherent);
        if (corrections)
        {
            xs = RIntAmplitude(t, Qsqr, bjorkx, &inthelperf_proton);
            xseps = RIntAmplitude(t, Qsqr, bjorkx+eps, &inthelperf_proton);
            lambda = log(xs/xseps)*(bjorkx/eps);
            result = tmpres*tmpres*(1.0+SQR(Beta(lambda)))*SQR(Rg(lambda));
        }
        else result=tmpres*tmpres;
        
        wavef->SetMode(VM_MODE_L);
        tmpres = RIntAmplitude(t, Qsqr, bjorkx, &inthelperf_coherent);
        if (corrections)
        {
            xs = RIntAmplitude(t, Qsqr, bjorkx, &inthelperf_proton);
            xseps = RIntAmplitude(t, Qsqr, bjorkx+eps, &inthelperf_proton);
            lambda = log(xs/xseps)*(bjorkx/eps);
            result += tmpres*tmpres*(1.0+SQR(Beta(lambda)))*SQR(Rg(lambda)); 
        } 
        else result +=tmpres*tmpres;
       
       wavef->SetMode(VM_MODE_TOT);
         // Sum d\sigma/dt_L + d\sgima/dt_T, not amplitudes!
    }
    else    // Correct polarization is already set
    {
        result = RIntAmplitude(t, Qsqr, bjorkx, &inthelperf_coherent);
        if (corrections)
        {
			
            /*REAL xs_old = RIntAmplitude(t, Qsqr, bjorkx, &inthelperf_proton);
            REAL xseps_old = RIntAmplitude(t, Qsqr, bjorkx+eps, &inthelperf_proton);
            REAL lambda_old = log(xs_old/xseps_old)*(bjorkx/eps);*/
            
            
            double y0 = log(1.0/bjorkx);

            double y1 = y0-eps_y;
            double y2 = y0 + eps_y;
            double x1 = exp(-y1);
            double x2 = exp(-y2);
            double xs1 = RIntAmplitude(t, Qsqr, x1, &inthelperf_proton);
            double xs2 = RIntAmplitude(t, Qsqr, x2, &inthelperf_proton);
            REAL lambda = log(xs2 / xs1) / (2.0*eps_y);
            //cout <<"x " << bjorkx << " realpart " << 1.0+SQR(Beta(lambda)) << " skew " << SQR(Rg(lambda)) << " lambda " << lambda << endl;
            //cout <<"# x=" << bjorkx <<", lambda: " << lambda << " oldlambda " << lambda_old << " rg " << SQR(Rg(lambda)) << " realpart " << SQR(Beta(lambda)) << endl;
            result=result*result*(1.0+SQR(Beta(lambda)))*SQR(Rg(lambda));
        }
        else result=result*result;
    }
        
    result *= 1.0/(16.0*M_PI);  
    return result;
}

// Only one r integral for cohernet dipole-nucleus scattering
REAL inthelperf_coherent(REAL r, void* p)
{
    inthelper_r* par = (inthelper_r*)p;
    REAL result = 2.0*M_PI*r*par->vm->PsiSqr_intz(par->Qsqr, r)
        * par->amplitude->CoherentDipoleAmplitude_avg(SQR(r), par->bjorkx,
         par->delta);
    //cout << r << " " << result << endl;
    return result;
}

/*
 * Total coherent cross section, integrate CoherentCrossSection_dt 
 * from mint to maxt
 */
 
REAL Calculator::TotalCoherentCrossSection(REAL Qsqr, REAL bjorkx)
{
    return TotalCoherentCrossSection(Qsqr, bjorkx, TOTXS_COHERENT_MINT,
        TOTXS_COHERENT_MAXT);
}

REAL Calculator::TotalCoherentCrossSection(REAL Qsqr, REAL bjorkx, 
    REAL mint, REAL maxt)
{
    
    gsl_function fun;   
    inthelper_r inthelp;
    inthelp.amplitude=amplitude;
    inthelp.vm=wavef; inthelp.bjorkx=bjorkx;
    inthelp.delta=-1; inthelp.Qsqr=Qsqr;
    inthelp.calculator=this;
    fun.function=&inthelper_totcohxs;
    fun.params=&inthelp; 
        
    REAL result,abserr; size_t eval;
    int status = gsl_integration_qng(&fun, mint, maxt, 0, relaccuracy_t ,
        &result, &abserr, &eval);
    if (status and result>0.00001)
        std::cerr << "Total cross section integral failed to reach tolerance: "
        << "Result: " << result << ", abserr: " << abserr << std::endl;

    return result;
}

// Integration over t 
REAL inthelper_totcohxs(REAL t, void* p)
{
    inthelper_r* par = (inthelper_r*)p; 
    return par->calculator->CoherentCrossSection_dt(t, par->Qsqr, par->bjorkx);
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
        result=xs*xs*(1.0+SQR(Beta(lambda)) ) * SQR(Rg(lambda));

        wavef->SetMode(VM_MODE_L);
        xs = RIntAmplitude(t, Qsqr, bjorkx, &inthelperf_proton);
        xseps = RIntAmplitude(t, Qsqr, bjorkx+eps, &inthelperf_proton);
        lambda = log(xs/xseps)*(bjorkx/eps);

        result+= xs*xs*(1.0+SQR(Beta(lambda)) ) * SQR(Rg(lambda));
        
        wavef->SetMode(VM_MODE_TOT);
        
    }
    else     // Correct polarization is already set
    {
        REAL xs = RIntAmplitude(t, Qsqr, bjorkx, &inthelperf_proton);
        /*REAL xseps = RIntAmplitude(t, Qsqr, bjorkx+eps, &inthelperf_proton);        
        REAL lambda = log(xs/xseps)*(bjorkx/eps);
        if (Beta(lambda)>1)
        {
            std::cerr << "beta=" << Beta(lambda) << " at t=" << t 
                << std::endl;
        
        }*/
        double y0 = log(1.0/bjorkx);
        //eps=0.001;
        eps=0.07;
        double y1 = y0-eps;
        double y2 = y0 + eps;
        double x1 = exp(-y1);
        double x2 = exp(-y2);
        double xs1 = RIntAmplitude(t, Qsqr, x1, &inthelperf_proton);
        double xs2 = RIntAmplitude(t, Qsqr, x2, &inthelperf_proton);
        double lambda = log(xs2 / xs1) / (2.0*eps);

        result=xs*xs*(1.0+SQR(Beta(lambda)) ) * SQR(Rg(lambda));
   }
    

    return result/(16.0*M_PI);
}


// Only one r integral for dipole-proton scattering
REAL inthelperf_proton(REAL r, void* p)
{
    inthelper_r* par = (inthelper_r*)p;
    REAL result =  2.0*M_PI*r*par->vm->PsiSqr_intz(par->Qsqr, r)
        * par->amplitude->DipoleAmplitude_proton(SQR(r), par->bjorkx,
         par->delta);

    if (isnan(result) or isinf(result))
        cerr << result << "  result at " << LINEINFO << ", r=" << r << " delta " << par->delta << " dipole amplitude " << par->amplitude->DipoleAmplitude_proton(SQR(r), par->bjorkx,
         par->delta)
            << " vm overlap " << par->vm->PsiSqr_intz(par->Qsqr, r) << endl;

    return result;
       
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
        
    REAL result,abserr; 
    result=0;
    //size_t eval;
    //int status = gsl_integration_qng(&fun, MINR, MAXR, 0, RINTACCURACY, 
    //    &result, &abserr, &eval);
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(100);
	int status = gsl_integration_qag(&fun, MINR, MAXR, 0, RINTACCURACY,
		100, GSL_INTEG_GAUSS51, ws, &result, &abserr);
	gsl_integration_workspace_free(ws);
        
    if (status and result>0.000001) 
        std::cerr << "Error " << status << " at " << __FILE__ << ":"
            << __LINE__ << ": Result " << result << ", relerror: " << std::abs(abserr/result)
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
    REAL result = TotalCrossSection(Qsqr, bjorkx, 0, TOTXS_MAXT);
    return result;
}

/*
 * Total cross section \int dt d\sigma/dt
 * Calculated in a given t range
 */
REAL Calculator::TotalCrossSection(REAL Qsqr, REAL bjorkx, REAL mint, REAL maxt)
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
    int status = gsl_integration_qng(&fun, mint, maxt, 0, relaccuracy_t, 
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
    {
        return par->calculator->ProtonCrossSection_dt(t, par->Qsqr, par->bjorkx); 
    }
    return par->calculator->CrossSection_dt(t, par->Qsqr, par->bjorkx);

}

/*
 * Calculate total gammap -> jpsi p cross section directly without integrating
 * over t
 * Analytically one can show that
 * sigma = \int d^2 b T(b)^2 | \int d^2 r dz/4pi psi^*psi N(r) |^2
 *
 * We use an approximation \int d^2 b T(b)^2 = \pi B_p = \sigma_0/4
 *  = 1/4 \sigma_0^2/(4pi B_p)  = 1/(16pi) sigma_0^2 / (B_p)
 */
///TODO: For some weird reason this does not give the same result as TotalCrossSection,
///and the difference is even more than a constant factor????????
REAL inthelperf_totxs_notint(double r, void *p)
{
    inthelper_r* par = (inthelper_r*)p;
    REAL result =  2.0*M_PI*r*par->vm->PsiSqr_intz(par->Qsqr, r)
        * par->amplitude->Qq_proton_amplitude(SQR(r), par->bjorkx, 0);

    return result/2.0;  // divide by 2.0 as RIntAmplitude multiplies by 2
}

REAL Calculator::TotalCrossSection_noint(REAL Qsqr, REAL bjorkx)
{
    gsl_function fun;   
    inthelper_r inthelp;
    inthelp.amplitude=amplitude;
    inthelp.vm=wavef; inthelp.bjorkx=bjorkx;
    inthelp.Qsqr=Qsqr;
    inthelp.calculator=this;
    inthelp.delta=0;
    fun.function=inthelperf_totxs_notint;
    fun.params=&inthelp;

    
        
    REAL result,abserr; 
    result=0;
    //size_t eval;
    //int status = gsl_integration_qng(&fun, MINR, MAXR, 0, RINTACCURACY, 
    //    &result, &abserr, &eval);
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(100);
	int status = gsl_integration_qag(&fun, MINR, MAXR, 0, RINTACCURACY,
		100, GSL_INTEG_GAUSS51, ws, &result, &abserr);
	gsl_integration_workspace_free(ws);

    result = result * result * 1.0/(16.0*M_PI) * SQR(amplitude->Sigma0())/amplitude->Bp();

    return result;

}

/*
 * Find t where coherent and incoherent cross sections are equal
 * Returns that t (GeV^2)
 */
struct CohIncohHelper
{
	Calculator* calc;
	double qsqr, bjorkx;
};
double CohIncohHelperf(double t, void* p)
{
	CohIncohHelper* par = (CohIncohHelper*)p;
	return par->calc->CoherentCrossSection_dt(t, par->qsqr, par->bjorkx) - par->calc->CrossSection_dt(t, par->qsqr, par->bjorkx);
}
double Calculator::CoherentIncoherent(double qsqr, double bjorkx)
{
	CohIncohHelper helper;
	helper.calc=this; helper.qsqr=qsqr; helper.bjorkx=bjorkx;
	
	const int MAX_ITER = 30;
    const double ROOTFINDACCURACY = 0.01;
    gsl_function f;
    f.params = &helper;
        
    f.function = &CohIncohHelperf;

    const gsl_root_fsolver_type *T = gsl_root_fsolver_bisection;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
        
    gsl_root_fsolver_set(s, &f, 0, 0.1);
    int iter=0; int status; double min,max;
    do
    {
        iter++;
        gsl_root_fsolver_iterate(s);
        min = gsl_root_fsolver_x_lower(s);
        max = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(min, max, 0, ROOTFINDACCURACY);    
    } while (status == GSL_CONTINUE and iter < MAX_ITER);

    if (iter>=MAX_ITER)
        cerr << "Solving failed at x=" << bjorkx << " " << LINEINFO << endl;


    double res = gsl_root_fsolver_root(s);

    gsl_root_fsolver_free(s);
	
	return res;
}


/*
 * d\\sigma/dy for AA -> AA+J/\Psi
 */
 
//const double M_v=3.097;	// J/\Psi mass in GeV
// Cross section integrated over t
double Calculator::DiffractiveAAtoJpsi(double y, double sqrts, Diffraction d, bool pa, int z)
{
    double M_v = wavef->MesonMass();
	double wsqr = sqrts * M_v * std::exp(y);
	double wsqr2 = sqrts * M_v * std::exp(-y);
	double xbj=	SQR(M_v)/wsqr;	// y
	double xbj2= SQR(M_v)/wsqr2;	// -y
	SetPolarization(VM_MODE_T);
	if (xbj > 0.02 or xbj2>0.02) return 0;
	
	double switcht1=0, switcht2=0;
	
	if (d==INCOHERENT and !pa) // for dipole-proton integrate from t=0
	{
		switcht1 = CoherentIncoherent(0, xbj);
		switcht2 = CoherentIncoherent(0, xbj2);
		cout << "# at x=" << xbj << ", wsqr=" << wsqr <<" switch at t=" << switcht1 <<", at x=" << xbj2 <<", wsqr=" << wsqr2 <<", switch at t=" << switcht2 << endl;
	}
	
	const double incoh_maxt=1.5;
	double res=0;
	double res1=0,res2=0;
	if (pa)
	{
		if (d==INCOHERENT)
			cerr <<"pA is implemented only for coherent scattering!!!!" << endl;
		
		double flux1 = NuclearPhotonFlux(y, sqrts, pa, z);
		// change target to proton
		int tmpa = amplitude->GetNucleus().GetA();
		amplitude->GetNucleus().SetA(1);
		double xs1= TotalCrossSection(0, xbj, 0, incoh_maxt);
		// change target back to the original nuke
		amplitude->GetNucleus().SetA(tmpa);
		double flux2= ProtonPhotonFlux(-y, sqrts);
		double xs2= TotalCoherentCrossSection(0, xbj2,0,0.1);
		res1=flux1*xs1;
		res2=flux2*xs2;
		cout <<"#x1=" << xbj <<" flux1=" << flux1 <<" xs1=" << xs1 <<" res1=" << res1 << " x2=" << xbj2 << " flux2=" << flux2 <<" xs2=" << xs2<<" res2=" << res2<< endl;
		return res1+res2;
	}
	
	
	
	double res1_tmp;
	if (xbj<1)
	{
		if (d==INCOHERENT)
		{
			res1 = NuclearPhotonFlux(y, sqrts, pa, z) * TotalCrossSection(0.0, xbj,switcht1,incoh_maxt);
		}
		else
		{
			res1 = NuclearPhotonFlux(y, sqrts, pa, z) * TotalCoherentCrossSection(0.0, xbj,0,0.1);
		}
		//cout <<"x=" << xbj <<" nuclear flux=" << NuclearPhotonFlux(y, sqrts, z) << " photon flux=" << ProtonPhotonFlux(y, sqrts) << endl;
	}
	if (xbj2<1)
	{
		if (d==INCOHERENT)
			res2 = NuclearPhotonFlux(-y, sqrts, pa, z) * TotalCrossSection(0.0, xbj2,switcht2,incoh_maxt);
		else
			res2=NuclearPhotonFlux(-y, sqrts, pa, z) * TotalCoherentCrossSection(0.0, xbj2,0,0.1);
	}
		
	
	res=res1+res2;
	return res;
}

double Calculator::DiffractiveAAtoJpsi_dt(double y, double sqrts, double t, Diffraction d, bool pa, int z)
{
    double M_v = wavef->MesonMass();
	double wsqr = sqrts * M_v * std::exp(y);
	double wsqr2 = sqrts * M_v * std::exp(-y);
	double xbj=	SQR(M_v)/wsqr;	// y
	double xbj2= SQR(M_v)/wsqr2;	// -y

	SetPolarization(VM_MODE_T);
	if (xbj > 0.02 or xbj2>0.02) return 0;

	double res1=0,res2=0;
	
	if (pa)
	{
		cerr << "pa dt is not implemented yet...." << endl;
		return 0;
	}
	
	if (xbj<1)
	{
		if (d==INCOHERENT)
			res1 = NuclearPhotonFlux(y, sqrts, pa, z) * CrossSection_dt(t, 0, xbj);
		else
			res1 = NuclearPhotonFlux(y, sqrts, pa, z) * CoherentCrossSection_dt(t, 0, xbj);
	}
	if (xbj2<1)
	{
		if (d==INCOHERENT)
			res2 = NuclearPhotonFlux(-y, sqrts, pa, z) * CrossSection_dt(t, 0, xbj2);
		else
			res2=NuclearPhotonFlux(-y, sqrts, pa, z) * CoherentCrossSection_dt(t, 0, xbj2);
	}
	cout << "# x1 " << xbj << " W1^2 " << wsqr << " dsigma/dt " << res1 << " x2 " <<  xbj2 << " W2^2 " << wsqr2 << " dsigma/dt " <<  res2 << endl;
	return res1+res2;
	
}


/*
 * Photon flux for diffractive AA events
 * y: rapidity of produced meson, z: charge
 * if pa=true, we have pA collision -> min. impact param is R_A, not 2R_A
 */
double Calculator::NuclearPhotonFlux(double y,  double sqrts,  bool pa, int z)
{
    double M_v = wavef->MesonMass();
	///FIXME: argument z is useless

    double mass=0;
    z=0;
    if (sqrts==2760 or sqrts==5020)
    {
        // LHC
        mass=193.729;
        z=82;
    }
    else if (sqrts==200)
    {
        mass=183.47;
        z=79;
    }
    else
    {
        cerr << "Sqrts " << sqrts << " is unknown, don't know nuclear mass and z" << endl;
        exit(1);
    }
    
    double ma=mass;
    //const double ma=193.729; // Mass of Pb
	//const double ma = 
	
	
	double omega = M_v/2.0*std::exp(y);
	double wsqr = 2.0*omega*sqrts;
	double xp = M_v*M_v/wsqr;
	int a = amplitude->GetNucleus().GetA();

	if (a != 208 or z != 82)
	{
	//	std::cerr << "PhotonFlux uses hardcoded nucleus mass, so only A=208 (Pb) is allowed!" << std::endl;
		//return 0;
	}
	
	double ra = (1.12*std::pow(a, 1.0/3.0) - 0.86*std::pow(a, -1.0/3.0)) * FMGEV;
	//double ra = 1.13*std::pow(a, 1.0/3.0) * FMGEV;
	if (pa)
	{
		ra/=2.0;	// in pA computation now we integrate from R_A to infty, not from 2R_A.
	}
	double gamma = a*sqrts/(2.0*ma);
	
	double xsi = 2.0*omega*ra /gamma;
	
	double flux = xsi * gsl_sf_bessel_K0(xsi) * gsl_sf_bessel_K1(xsi)
		- xsi*xsi/2.0 * ( SQR( gsl_sf_bessel_K1(xsi)) - SQR( gsl_sf_bessel_K0(xsi)) ) ;

	return 2.0 * ALPHA_e * z*z / M_PI * flux;

}

double Calculator::ProtonPhotonFlux(double y, double sqrts)
{
    double M_v = wavef->MesonMass();
	double omega = M_v/2*std::exp(y);
	double wsqr = 2.0*omega*sqrts;
	double xp = M_v*M_v/wsqr;
	double mp = 0.938;
	double gamma = sqrts/(2.0*mp);
	
	double minqsqr = SQR(omega)/(SQR(gamma) * (1.0 - 2.0*omega/sqrts) );
	double xsi=1.0+(0.71/minqsqr);
	
	return ALPHA_e / (2.0*M_PI) * ( 1.0  + SQR(1.0 - 2.0*omega/sqrts) )
		* (log(xsi) - 11.0/6.0 + 3.0/xsi - 3.0/(2.0*SQR(xsi)) + 1.0/(3.0*xsi*xsi*xsi) );
	
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

void Calculator::SetTAccuracy(REAL acc)
{
    relaccuracy_t=acc;
}

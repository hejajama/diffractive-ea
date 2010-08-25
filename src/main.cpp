/*
 * Calculates eA cross section in dipole model
 * 
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
#include <iostream>
#include <ctime>
#include <vector>
#include <sstream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_integration.h>

#include "dipole.h"
#include "vm_photon.h"
#include "vector.h"
#include "nucleus.h"
#include "dipxs_ipnonsat.h"
#include "dipxs_ipsat.h"
#include "gdist/gdist_dglap.h"
#include "mersenne/mersenne.h"

using namespace std;

// Integration settings
const REAL MAXR=4;
const REAL MINR=0.05;   // r=0 doesn't work, K_{0,1}(0)=inf
const REAL RINTACCURACY=0.0001;

const int MODEL_IPSAT=1; const int MODEL_IPNONSAT=2;
const int GDIST_DGLAP=1; const int GDIST_TOY=2;

// Integral helpers for outern r integral
struct inthelper_r
{
    Dipxs* amplitude;
    VM_Photon* vm;  // Vector meson wave function
    REAL r;         // r for inner r' integral, not used in outern inthelperf_r1
    REAL Qsqr;
    REAL bjorkx;
    REAL delta;
};

REAL inthelperf_r1(REAL r, void* p);

// For inner r' integral
REAL inthelperf_r2(REAL r2, void* p);

int main(int argc, char* argv[])
{    

    // Parameters
    REAL bjorkx=1e-4;
    REAL Qsqr=0;
    int A=197;
    int model=MODEL_IPSAT;
    int gdist_model=GDIST_DGLAP;
    int points=100;
    REAL maxt=0.3; REAL mint=0; 
            
    // Parse parameters
    if (argc>1)
    {
        if (string(argv[1])=="--help")
        {
            cout << "Usage: -x bjorkx -Q2 Q^2" << endl;
            cout << "-dipole {ipsat,ipnonsat}" << endl;
            cout << "-gdist {dglap,toy}" << endl;
            cout << "-A number_of_nucleai" << endl;
            cout << "-N number_of_data_points" << endl;
            cout << "-mint t_value, -maxt t_value" << endl;
            cout << "Default values: x="<<bjorkx <<", Q^2="<<Qsqr 
                << " A="<<A<<", N="<<points<<", mint="<<mint<<", maxt="<<maxt<< endl;
            return 0;
        }
        for (int i=1; i<argc; i++)
        {
            if (string(argv[i])=="-x") 
                bjorkx=StrToReal(argv[i+1]); 
            if (string(argv[i])=="-Q2")
                Qsqr=StrToReal(argv[i+1]); 
            if (string(argv[i])=="-A")
                A=StrToInt(argv[i+1]);
            if (string(argv[i])=="-N")
                points=StrToInt(argv[i+1]);
            if (string(argv[i])=="-mint")
                mint=StrToReal(argv[i+1]);
            if (string(argv[i])=="-maxt")
                maxt=StrToReal(argv[i+1]);
            if (string(argv[i])=="-dipole")
            {
                if (string(argv[i+1])=="ipsat")
                    model=MODEL_IPSAT;
                else if (string(argv[i+1])=="ipnonsat")
                    model=MODEL_IPNONSAT;
                else
                    cerr << "Model " << argv[i+1] << " is not valid" << endl;
            }
            if (string(argv[i])=="-gdist")
            {
                if (string(argv[i+1])=="dglap")
                    gdist_model=GDIST_DGLAP;
                else if (string(argv[i+1])=="toy")
                    gdist_model=GDIST_TOY;
                else
                    cerr << "Gluon distribution " << argv[i+1] 
                        << "is not valid " << endl;
             }
                
        }
    }
    
    // Print values
    cout << "# x=" << bjorkx << ", Q^2=" << Qsqr << " A=" << A << endl;
    cout << "# GDist=" << gdist_model << ",  dipole model=" << model << endl;
    
    // Intialize random number generator
    seed_mersenne(std::time(NULL));
    
    // J/Psi wave function:  e_f, N_T, N_L, R_T, R_L, m_f, M_V, delta
    //VM_Photon JPsi(2.0/3.0, 1.23, 0.83, sqrt(6.5), sqrt(3.0), 1.4, 3.097, 1);
    VM_Photon JPsi("jpsi.dat");
    
    cout << JPsi.PsiSqr_tot_intz(Qsqr,1) << endl; return 0;

    
    // Intialize Dipxs and Nucleus
    Nucleus nuke(A);
    //GDist *gdist = new DGLAPDist();
    GDist *gdist;
    if (gdist_model==GDIST_DGLAP)
        gdist = new DGLAPDist();
    else if (gdist_model==GDIST_TOY)
        gdist = new GDist_Toy();
    nuke.SetGDist(gdist);    
    Dipxs *dsigmadb;
    if (model==MODEL_IPSAT)
        dsigmadb = new Dipxs_IPSat(nuke);
    else if (model==MODEL_IPNONSAT)
        dsigmadb = new Dipxs_IPNonSat(nuke);
    
    /*******************
     * \gamma^* N -> J/\Psi N cross section
     * Calculates
     * \int d^2 r d^2 r' (jpsi)(r)*(jpsi)(r')*qqamplitude_sqr_avg(r,delta)
     * Here (jpsi) is the inner product between \gamma^* and J/\Psi wave 
     * functions integrated over z \in [0,1]
     *
     * All iterations are independent, so this is straightforward to parallerize
     */

       
    #pragma omp parallel for
    for (int i=0; i<=points; i++)
    {
        REAL tmpt = (maxt-mint)/points*i;
        REAL delta = sqrt(tmpt);
        gsl_function fun;
        inthelper_r inthelp;
        inthelp.amplitude=dsigmadb;
        inthelp.vm=&JPsi; inthelp.bjorkx=bjorkx;
        inthelp.delta=delta; inthelp.Qsqr=Qsqr;
        fun.function=&inthelperf_r1;
        fun.params=&inthelp;
        
        REAL result,abserr; size_t eval;
        int status = gsl_integration_qng(&fun, MINR, MAXR, RINTACCURACY, RINTACCURACY, 
            &result, &abserr, &eval);
        result *= 1.0/(16.0*M_PI);
        cout << SQR(delta) << " " << result << endl;
    }



    /******************
     * Dipole-nucleus cross section for a fixed r 
     */
    /*REAL rsqr=SQR(0.4);
    REAL maxt=2.0;
    int points=100;
    for (REAL delta=0; SQR(delta)<maxt; delta+=sqrt(maxt)/points)
    {
        REAL tmpxs=1.0/(16.0*M_PI)*dsigmadb->Dipxsection_sqr_avg(rsqr, rsqr, x, delta);
        cout << SQR(delta) << " " << tmpxs << endl;
    }*/

   
    return 0;
}


// Integral helpers

REAL inthelperf_r1(REAL r, void* p)
{
    inthelper_r* par = (inthelper_r*)p;
    par->r=r;
    gsl_function fun;
    fun.function=&inthelperf_r2;
    fun.params=par;
    
    REAL result,abserr; size_t eval;
    
    int status = gsl_integration_qng(&fun, MINR, MAXR, RINTACCURACY, RINTACCURACY, 
        &result, &abserr, &eval);
    
    return 2*M_PI*r*par->vm->PsiSqr_tot_intz(par->Qsqr, r)*result;
}

// Inner r' integral: \int dr' r' (jpsi)(r') * qqamplitude_sqr(r,r')
REAL inthelperf_r2(REAL r2, void* p)
{
    inthelper_r* par = (inthelper_r*)p;
    return 2*M_PI*r2 * par->vm->PsiSqr_tot_intz(par->Qsqr, r2) 
            * par->amplitude->Dipxsection_sqr_avg(SQR(par->r), SQR(r2), 
                    par->bjorkx, par->delta);

}

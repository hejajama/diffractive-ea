/*
 * Calculates eA cross section in dipole model
 * 
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
#include <iostream>
#include <iomanip>  // Set cout precision
#include <ctime>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include "dipole.h"
#include "vm_photon.h"
#include "vector.h"
#include "nucleus.h"
#include "dipxs_ipnonsat.h"
#include "dipxs_ipsat.h"
#include "dipxs_iim.h"
#include "gdist/gdist_dglap.h"
#include "mersenne/mersenne.h"
#include "calculator.h"

using namespace std;

// Integration settings
/*const REAL MAXR=4;
const REAL MINR=0.05;   // r=0 doesn't work, K_{0,1}(0)=inf
const REAL RINTACCURACY=0.002;
const REAL TOTXS_MAXT=2;    // Max |t| in GeV^2 when calculating total xs
*/
const int MODEL_IPSAT=1; const int MODEL_IPNONSAT=2; const int MODEL_IIM=3;
const int MODEL_IPSAT_NONSATP=4;
const int GDIST_DGLAP=1; const int GDIST_TOY=2;

void Cleanup(); 

int main(int argc, char* argv[])
{    
    gsl_set_error_handler(&ErrHandler);   // We really don't want to use this!
    
    // Parameters
    REAL bjorkx=1e-4;
    REAL Qsqr=0;
    int A=197;
    int model=MODEL_IPSAT;
    int gdist_model=GDIST_DGLAP;
    int points=100;
    REAL maxt=0.3; REAL mint=0; 
    bool totxs=false;   // Calculate total cross section
    bool Ap=false;      // \sigma^A / A*\sigma_p as a function of Q^2
    REAL t=0.5;
    REAL maxQsqr=30;
    
    string iim_file="iim.dat";  // Read parameters for IIM model from this file
            
    // Parse parameters
    if (argc>1)
    {
        if (string(argv[1])=="--help")
        {
            cout << "Usage: -x bjorkx -Q2 Q^2" << endl;
            cout << "-dipole {ipsat,ipnonsat,iim,ipsat_nonsatp}" << endl;
            cout << "-gdist {dglap,toy}" << endl;
            cout << "-A number_of_nucleai" << endl;
            cout << "-N number_of_data_points" << endl;
            cout << "-mint t_value, -maxt t_value" << endl;
            cout << "-iimfile filename (parameters for the IIM model)" << endl;
            cout << "-totxs (calculates total cross section" << endl;
            cout << "-A/p (nucleus cross section / A* "
                    <<"proton cross section as a function of Q^2)" << endl;
            cout << "-t t, -maxQ2 maxq2 [GeV] (value of t and maximum of Q^2 in case of -A/p)" << endl;
            cout << "Default values: x="<<bjorkx <<", Q^2="<<Qsqr 
                << " A="<<A<<", N="<<points<<", mint="<<mint<<", maxt="<<maxt<< endl;
            cout << "                dipxs=false, A/p=false, iimfile=" << iim_file << endl;
            cout << "                t="<<t << endl;
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
            if (string(argv[i])=="-totxs")
                totxs=true;
            if (string(argv[i])=="-A/p")
                Ap=true;
            if (string(argv[i])=="-t")
                t=StrToReal(argv[i+1]);
            if (string(argv[i])=="-maxQ2")
                maxQsqr=StrToReal(argv[i+1]);
            if (string(argv[i])=="-iimfile")
                iim_file=string(argv[i+1]);
            if (string(argv[i])=="-dipole")
            {
                if (string(argv[i+1])=="ipsat")
                    model=MODEL_IPSAT;
                else if (string(argv[i+1])=="ipnonsat")
                    model=MODEL_IPNONSAT;
                else if (string(argv[i+1])=="iim")
                    model=MODEL_IIM;
                else if (string(argv[i+1])=="ipsat_nonastp")
                    model=MODEL_IPSAT_NONSATP;
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
    cout << "# x=" << bjorkx << ", Q^2=" << Qsqr << " A=" << A;
    if (A==1) cout << " (dipole-proton)";
    cout << endl;
    cout << "# GDist=" << gdist_model << ",  dipole model=" << model << endl;
 
    
    // Intialize random number generator
    seed_mersenne(std::time(NULL));
    
    // J/Psi wave function:  e_f, N_T, N_L, R_T, R_L, m_f, M_V, delta
    //VM_Photon JPsi(2.0/3.0, 1.23, 0.83, sqrt(6.5), sqrt(3.0), 1.4, 3.097, 1);
    WaveFunction *JPsi = new VM_Photon("jpsi.dat");
    
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
    else if (model==MODEL_IIM)
        dsigmadb = new Dipxs_IIM(nuke, iim_file);
    else if (model==MODEL_IPSAT_NONSATP)
        dsigmadb = new Dipxs_IPSat(nuke, IPSAT_MODE_NONSAT_P);

    Calculator calculator(dsigmadb, JPsi);

    /*******************
     * \gamma^* N -> J/\Psi N cross section
     * Calculates:   d\sigma / dt = 1/(16*\pi)*
     * \int d^2 r d^2 r' (jpsi)(r)*(jpsi)(r')*qqamplitude_sqr_avg(r,delta)
     * Here (jpsi) is the inner product between \gamma^* and J/\Psi (or VM) wave 
     * functions integrated over z \in [0,1]
     */

    if (totxs)  // Calculate total cross section
    {
        REAL result = calculator.TotalCrossSection(Qsqr, bjorkx);
        cout << "Total cross section: " << result*400.0*1000.0 << " nb" << endl;
    }
    else if (Ap)    // Calculate d\sigma^A/dt / A*d\sigma_p/dt as a function Q^2 at t
    {         
        cout << "# t=" << t << endl;
        cout << "# Q^2   nucleus_xs / A*proton_xs" << endl;
        
        #pragma omp parallel for
        for (int i=0; i<=points; i++)
        {
            REAL tmpqsqr=maxQsqr/points*i;
            REAL protonxs = calculator.ProtonCrossSection_dt(t, tmpqsqr, bjorkx);
            REAL nukexs = calculator.CrossSection_dt(t, tmpqsqr, bjorkx);
            #pragma omp critical
            {
                cout.precision(5);
                cout << fixed << tmpqsqr;
                cout.precision(8);
                cout << " " << nukexs / (A*protonxs) << endl;
            }
        }
    
    
    }
    else
    {
        // All iterations are independent, so this is straightforward to parallerize   
        #pragma omp parallel for
        for (int i=0; i<=points; i++)
        {
            REAL tmpt = (maxt-mint)/points*i;
            REAL delta = sqrt(tmpt);
            REAL result;
            if (A==1)
                result = calculator.ProtonCrossSection_dt(tmpt, Qsqr, bjorkx);
            else
                result = calculator.CrossSection_dt(tmpt, Qsqr, bjorkx);
            #pragma omp critical
            {
                cout.precision(5);
                cout << fixed << tmpt;
                cout.precision(8);
                cout << " " << result << endl;
            }
        }
    }
    
    delete dsigmadb;
    delete gdist;
    delete JPsi;
   
    return 0;
}



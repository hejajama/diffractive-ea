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

const int MODE_TOTXS=1; const int MODE_DIFFXS=2; const int MODE_Ap=3;
const int MODE_Ap_X=4; const int MODE_TOTXS_Q=5;

void Cleanup(); 

int main(int argc, char* argv[])
{    
    gsl_set_error_handler(&ErrHandler);   // We really don't want to use this!
    
    // Parameters
    REAL bjorkx=1e-4;
    REAL Qsqr=1;
    int A=197;
    int model=MODEL_IPSAT;
    int gdist_model=GDIST_DGLAP;
    int points=100;
    REAL maxt=0.3; REAL mint=0; 
    REAL t=0.5;
    REAL maxQsqr=30;
    int mode=MODE_DIFFXS;   // What to do
    REAL minx=1e-6;
    REAL maxx=1e-2;
    string xgfile="xg.dat";
    REAL W=-1;  // W^2 = (P+q)^2 = invariant mass^2 of the \gamma^*N system
    bool x_set=false;
    bool w_set=false;
    bool xgval=false;   // Print the value of xg(x,µ) and quit.
    REAL r=-1;
    REAL M_v=3.097; // Mass of the produced vector meson, 3.097 GeV = J/\Psi
    
    string iim_file="iim.dat";  // Read parameters for IIM model from this file
            
    // Parse parameters
    if (argc>1)
    {
        if (string(argv[1])=="--help")
        {
            cout << "Usage: -x bjorkx -Q2 Q^2 -W W (specify only x or W)" << endl;
            cout << "-dipole {ipsat,ipnonsat,iim,ipsat_nonsatp}" << endl;
            cout << "-gdist {dglap,toy} -xgfile file" << endl;
            cout << "-A number_of_nucleai" << endl;
            cout << "-N number_of_data_points" << endl;
            cout << "-mint t_value, -maxt t_value" << endl;
            cout << "-iimfile filename (parameters for the IIM model)" << endl;
            cout << "-totxs (calculates total cross section)" << endl;
            cout << "-totxs_q2 (total cross section as a function of Q^2)" << endl;
            cout << "-A/p (nucleus cross section / A* "
                    <<"proton cross section as a function of Q^2)" << endl;
            cout << "-t t, -maxQ2 maxq2 [GeV] (value of t and maximum of Q^2 in case of -A/p)" << endl;
            cout << "-A/p_x (same as -A/p but as a function of x), -minx minx, -maxx maxx" << endl;
            cout << "-xg rval (print the value of xg(x,r) and quit) " << endl;
            cout << "-Mv mass (mass of the produced vector meson) " << endl;
            cout << endl;
            cout << "Default values: x="<<bjorkx <<", Q^2="<<Qsqr 
                << " A="<<A<<", N="<<points<<", mint="<<mint<<", maxt="<<maxt<< endl;
            cout << "                dipxs=false, A/p=false, iimfile=" << iim_file << endl;
            cout << "                t="<<t << ", minx=" << minx << ", maxx=" << maxx << endl;
            cout << "                Mv=" << M_v << endl;

            return 0;
        }
        for (int i=1; i<argc; i++)
        {
            if (string(argv[i])=="-x") 
            {
                bjorkx=StrToReal(argv[i+1]); 
                x_set=true;
            }
            if (string(argv[i])=="-W")
            {
                W = StrToReal(argv[i+1]);                
                w_set=true;
            }
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
                mode=MODE_TOTXS;
            if (string(argv[i])=="-A/p")
                mode=MODE_Ap;
            if (string(argv[i])=="-A/p_x")
                mode=MODE_Ap_X;
            if (string(argv[i])=="-totxs_q2")
                mode=MODE_TOTXS_Q;
            if (string(argv[i])=="-t")
                t=StrToReal(argv[i+1]);
            if (string(argv[i])=="-maxQ2")
                maxQsqr=StrToReal(argv[i+1]);
            if (string(argv[i])=="-iimfile")
                iim_file=string(argv[i+1]);
            if (string(argv[i])=="-xgfile")
                xgfile=string(argv[i+1]);
            if (string(argv[i])=="-Mv")
                M_v=StrToReal(argv[i+1]);
            if (string(argv[i])=="-xg") {
                xgval=true;
                r=StrToReal(argv[i+1]);
            }
            if (string(argv[i])=="-dipole")
            {
                if (string(argv[i+1])=="ipsat")
                    model=MODEL_IPSAT;
                else if (string(argv[i+1])=="ipnonsat")
                    model=MODEL_IPNONSAT;
                else if (string(argv[i+1])=="iim")
                    model=MODEL_IIM;
                else if (string(argv[i+1])=="ipsat_nonsatp")
                    model=MODEL_IPSAT_NONSATP;
                else
                {
                    cerr << "Model " << argv[i+1] << " is not valid" << endl;
                    return -1;
                }
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
    
    if (x_set==true and w_set==true)
    {
        cerr << "Both x and W set, don't know what to do. Exiting..." << endl;
        return -1;
    }
    if (w_set)  // TODO: Check
    {
        //bjorkx = Qsqr/(Qsqr + SQR(W));
        bjorkx = (Qsqr + SQR(M_v))/SQR(W);
        bjorkx = Qsqr/(Qsqr+SQR(W))*(1+SQR(M_v)/Qsqr);
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
        gdist = new DGLAPDist(xgfile);
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

    if (xgval)  // Print the value of xg(x,µ) and quit
    {
        REAL tmp = gdist->Gluedist(bjorkx, r*r);
        tmp *= 2*NC/(M_PI*M_PI*Alpha_s(Mu2(SQR(r))) );
        cout << tmp << endl;
        return 0;
    }


    Calculator calculator(dsigmadb, JPsi);

    /*******************
     * \gamma^* N -> J/\Psi N cross section
     * Calculates:   d\sigma / dt = 1/(16*\pi)*
     * \int d^2 r d^2 r' (jpsi)(r)*(jpsi)(r')*qqamplitude_sqr_avg(r,delta)
     * Here (jpsi) is the inner product between \gamma^* and J/\Psi (or VM) wave 
     * functions integrated over z \in [0,1]
     */

    if (mode==MODE_TOTXS)  // Calculate total cross section
    {
        REAL result = calculator.TotalCrossSection(Qsqr, bjorkx);
        cout << "Total cross section: " << result*400.0*1000.0 << " nb" << endl;
    }
     
    else if (mode==MODE_TOTXS_Q)    // Total cross section as a function of Q^2
    {
        cout << "# Total cross section, W=" << W << endl;
        #pragma omp parallel for
        for (int i=1; i<=points; i++)
        {
            REAL tmpqsqr = maxQsqr/points*i;
            bjorkx = tmpqsqr/(tmpqsqr+SQR(W))*(1+SQR(M_v)/tmpqsqr);
            bjorkx = (tmpqsqr + SQR(M_v))/SQR(W);
            REAL xs = calculator.TotalCrossSection(tmpqsqr, bjorkx);
            xs *= 400.0 * 1000.0;   // Gev^{-4} => nb/GeV^2   
            #pragma omp critical
            {
                cout.precision(5);
                cout << fixed << tmpqsqr;
                cout.precision(8);
                cout << " " << xs << endl;
            }
        
        }
    
    
    
    }
    else if (mode==MODE_Ap)    // Calculate d\sigma^A/dt / A*d\sigma_p/dt as a function Q^2 at t
    {         
        cout << "# t=" << t << ", x = " << bjorkx <<  endl;
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
    
    else if (mode==MODE_Ap_X)   // d\sigma^A/dt / A*d\sigma_p/dt as a function of x
    {
        cout << "# t=" << t << ", Q^2=" << Qsqr << endl;
        cout << "# x-dependence of nucleus_xs / A*proton_xs" << endl;
        
        #pragma omp parallel for
        for (int i=0; i<points; i++)
        {
            REAL tmpx = minx + (maxx-minx)/points*i;
            REAL protonxs = calculator.ProtonCrossSection_dt(t, Qsqr, tmpx);
            REAL nukexs = calculator.CrossSection_dt(t, Qsqr, tmpx);
            #pragma omp critical
            {
                cout.precision(6);
                cout << fixed << tmpx;
                cout.precision(8);
                cout << " " << nukexs / (A*protonxs) << endl;
            }
        
        }
    
    }

    
    else    // dsigma/dt as a function of t
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



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
#include "gaus_lc.h"
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
const int MODE_Ap_X=4; const int MODE_TOTXS_Q=5; const int MODE_COHERENT_DT=6;
const int MODE_VM_INTZ=7; const int MODE_TOTXS_RATIO_Q=8;
const int MODE_DSIGMA_D2B_R=9; const int MODE_TOTXS_W=10;
const int MODE_DSIGMA_DT_QQ=11; const int MODE_DSIGMA_D2B_B=12;
const int MODE_TOTAL_DIPXS=13; const int MODE_XG=14; const int MODE_GLUEDIST=15;

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
    REAL minQsqr=0;
    REAL minW=0;
    REAL maxW=200;
    int mode=MODE_DIFFXS;   // What to do
    REAL minx=1e-6;
    REAL maxx=1e-2;
    string xgfile="xg.dat";
    REAL W=-1;  // W^2 = (P+q)^2 = invariant mass^2 of the \gamma^*N system
    bool x_set=false;
    bool w_set=false;
    bool scalex=false;  // Scale x so that x=x_bj*(1+M_V^2/Q^2)
    bool output_fm=false;   // Use fm's when printing the value of r
    bool factor_ipsat=true; // Factroize IPSat amplitude to T(b)(1-exp(-r^2...))
    REAL r=-1;
    REAL M_v=3.097; // Mass of the produced vector meson, 3.097 GeV = J/\Psi
    REAL M_n=0;     // Mass of nucleus/nucleon
    REAL bp=DEFAULT_B_p;
    int polarization = VM_MODE_TOT;
    REAL b=0;   // Impact parameter
    
    string iim_file="iim.dat";  // Read parameters for IIM model from this file
            
    // Parse parameters
    if (argc>1)
    {
        if (string(argv[1])=="--help")
        {
            cout << "Usage: -x x_pomeron -Q2 Q^2 -W W (specify only x or W)" << endl;
            cout << "-scalex (scale x by factor 1+M_V^2/Q^2) (can't be used with -W)" << endl;
            cout << "-dipole {ipsat,ipnonsat,iim,ipsat_nonsatp}" << endl;
            cout << "-gdist {dglap} -xgfile file" << endl;
            cout << "-A number_of_nucleai -Mn nucleus_mass" << endl;
            cout << "-N number_of_data_points" << endl;
            cout << "-coherent_dt (calculate coherent d\\sigma/dt)" << endl;
            cout << "-Bp proton_shape (for ipsat and ipnonsat)" << endl;
            cout << "-nofactor (do not factorize amplitude in IPSat model)" << endl;
            cout << "-mint t_value, -maxt t_value" << endl;
            cout << "-iimfile filename (parameters for the IIM model)" << endl;
            cout << "-totxs (calculates total cross section)" << endl;
            cout << "-totxs_q2 (total cross section as a function of Q^2)" << endl;
            cout << "-totxs_w (total cross section as a function of W)" << endl;
            cout << "-totxs_q2_l/t (longitudinal total xs / transversial total xs)" << endl;
            cout << "-A/p (nucleus cross section / A* "
                    <<"proton cross section as a function of Q^2)" << endl;
            cout << "-t t, -maxQ2 maxq2 -minQ2 minQ2[GeV] (value of t and max/min of Q^2)" << endl;
            cout << "-minW minW, -maxW maxW [GeV] (max/min of W)" << endl;
            cout << "-A/p_x (same as -A/p but as a function of x), -minx minx, -maxx maxx" << endl;
            cout << "-xg (print the value of xg(x,r) and quit) " << endl;
            cout << "-gdistval (print ONLY gluedist and quit)" << endl;
            cout << "-Mv mass (mass of the produced vector meson) " << endl;
            cout << "-vm_intz (print \\int d^z/(4\\pi) r Psi^*Psi) -fm (print r in fm)" << endl;
            cout << "-pol {l, t, sum} (polarization of the VM, sum = l + t is default)" << endl;
            cout << "-dsigma/d2b_r impact_par (print d\\sgima_{q\bar q}/d^2 b as a function of r)" << endl;
            cout << "-dsigma/d2b_b (print d\\sigma_{q\\bar q}/d^2 b as a function of b)" << endl;
            cout << "-dsigma/dt_qq (print d\\sigma_{q\bar q}/dt as a function of t), -r r (in GeV)" << endl;
            cout << "-totdipxs (print total dipole-proton cross section as a function of r)" << endl;
            cout << endl;
            cout << "Default values: x="<<bjorkx <<", Q^2="<<Qsqr 
                << " A="<<A<<", N="<<points<<", mint="<<mint<<", maxt="<<maxt<< endl;
            cout << "                dipxs=false, A/p=false, iimfile=" << iim_file << endl;
            cout << "                t="<<t << ", minx=" << minx << ", maxx=" << maxx 
                 << " xgfile=" << xgfile << endl;
            cout << "                Mv=" << M_v << ", Mn=" << M_n << ", B_p=" << bp << endl;
            cout << "                minW=" << minW << ", maxW=" << maxW << endl;

            return 0;
        }
        for (int i=1; i<argc; i++)
        {
            if (string(argv[i])=="-x") 
            {
                bjorkx=StrToReal(argv[i+1]); 
                x_set=true;
            }
            else if (string(argv[i])=="-W")
            {
                W = StrToReal(argv[i+1]);                
                w_set=true;
            }
            else if (string(argv[i])=="-Q2")
                Qsqr=StrToReal(argv[i+1]); 
            else if (string(argv[i])=="-A")
                A=StrToInt(argv[i+1]);
            else if (string(argv[i])=="-N")
                points=StrToInt(argv[i+1]);
            else if (string(argv[i])=="-mint")
                mint=StrToReal(argv[i+1]);
            else if (string(argv[i])=="-maxt")
                maxt=StrToReal(argv[i+1]);
            else if (string(argv[i])=="-maxx")
                maxx=StrToReal(argv[i+1]);
            else if (string(argv[i])=="-minx")
                minx=StrToReal(argv[i+1]);
            else if (string(argv[i])=="-coherent_dt")
                mode=MODE_COHERENT_DT;
            else if (string(argv[i])=="-nofactor")
                factor_ipsat=false;
            else if (string(argv[i])=="-totxs")
                mode=MODE_TOTXS;
            else if (string(argv[i])=="-A/p")
                mode=MODE_Ap;
            else if (string(argv[i])=="-A/p_x")
                mode=MODE_Ap_X;
            else if (string(argv[i])=="-totxs_q2")
                mode=MODE_TOTXS_Q;
            else if (string(argv[i])=="-totxs_q2_l/t")
                mode=MODE_TOTXS_RATIO_Q;
            else if (string(argv[i])=="-totxs_w")
                mode=MODE_TOTXS_W;
            else if (string(argv[i])=="-t")
                t=StrToReal(argv[i+1]);
            else if (string(argv[i])=="-r")
                r=StrToReal(argv[i+1]);
            else if (string(argv[i])=="-maxQ2")
                maxQsqr=StrToReal(argv[i+1]);
            else if (string(argv[i])=="-minQ2")
                minQsqr=StrToReal(argv[i+1]);
             else if (string(argv[i])=="-maxW")
                maxW=StrToReal(argv[i+1]);
            else if (string(argv[i])=="-minW")
                minW=StrToReal(argv[i+1]);
            else if (string(argv[i])=="-iimfile")
                iim_file=string(argv[i+1]);
            else if (string(argv[i])=="-xgfile")
                xgfile=string(argv[i+1]);
            else if (string(argv[i])=="-Mv")
                M_v=StrToReal(argv[i+1]);
            else if (string(argv[i])=="-Mn")
                M_n=StrToReal(argv[i+1]);
            else if (string(argv[i])=="-Bp")
                bp=StrToReal(argv[i+1]);
            else if (string(argv[i])=="-xg") 
                mode=MODE_XG;
            else if (string(argv[i])=="-gdistval")
                mode=MODE_GLUEDIST;
            else if (string(argv[i])=="-scalex")
                scalex=true;
            else if (string(argv[i])=="-vm_intz")
                mode=MODE_VM_INTZ;
            else if (string(argv[i])=="-fm")
                output_fm=true;
            else if (string(argv[i])=="-dsigma/d2b_r")
            {
                mode=MODE_DSIGMA_D2B_R;
                b = StrToReal(argv[i+1]);
            }
            else if (string(argv[i])=="-dsigma/d2b_b")
                mode=MODE_DSIGMA_D2B_B;
            else if (string(argv[i])=="-dsigma/dt_qq")
                mode=MODE_DSIGMA_DT_QQ;
            else if (string(argv[i])=="-totdipxs")
                mode=MODE_TOTAL_DIPXS;
            else if (string(argv[i])=="-dipole")
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
            else if (string(argv[i])=="-gdist")
            {
                if (string(argv[i+1])=="dglap")
                    gdist_model=GDIST_DGLAP;
                else if (string(argv[i+1])=="toy")
                    gdist_model=GDIST_TOY;
                else
                    cerr << "Gluon distribution " << argv[i+1] 
                        << "is not valid " << endl;
            } 
            else if (string(argv[i])=="-pol")
            {
                if (string(argv[i+1])=="tot")
                    polarization = VM_MODE_TOT;
                else if (string(argv[i+1])=="t")
                    polarization = VM_MODE_T;
                else if (string(argv[i+1])=="l")
                    polarization = VM_MODE_L;
                else
                    cerr << "Unrecognized polarization parameter " << argv[i+1]
                        << endl;
            }
            else if (string(argv[i]).substr(0,1)=="-")
            {
                cerr << "Unrecognized parameter " << argv[i] << endl;
                return -1;
            }
                
        }
    }
    
    if (x_set==true and w_set==true)
    {
        cerr << "Both x and W set, don't know what to do. Exiting..." << endl;
        return -1;
    }
    if (w_set and mode!=MODE_TOTXS_Q and mode!=MODE_TOTXS_W)  // TODO: Check
    {
        //bjorkx = Qsqr/(Qsqr + SQR(W));
        // x = x_bj*(1+M^2/Q^2) = Q^2/(Q^2+W^2+M_n^2)(1+M^2/Q^2)
        //   = (Q^2 + M_v^2) / (Q^2 + W^2 + M_n^2)
        std::cout << "# x scaled from " << bjorkx << " to " ;
        bjorkx = (Qsqr + SQR(M_v))/( Qsqr + SQR(W) + SQR(M_n) );
        cout << bjorkx << endl;
    }

    if (!w_set and scalex) // Scale x->x*(1+M^2/Q^2)
    {
        if (Qsqr<0.00001){ cerr << "Q^2=0, can't scale x" << endl; return -1;}
        bjorkx = bjorkx * (1 + SQR(M_v)/Qsqr);
    }
    
    // Print values
    cout << "# A=" << A;
    if (A==1) cout << " (dipole-proton)";
    cout << endl;
    cout << "# GDist=" << gdist_model << ",  dipole model=" << model << endl;
    
    // Intialize random number generator
    seed_mersenne(std::time(NULL));
    
    // J/Psi wave function:  e_f, N_T, N_L, R_T, R_L, m_f, M_V, delta
    //VM_Photon JPsi(2.0/3.0, 1.23, 0.83, sqrt(6.5), sqrt(3.0), 1.4, 3.097, 1);
    WaveFunction *JPsi = new GausLC("jpsi.dat");
    JPsi->SetMode(polarization);
    
    // Intialize Dipxs and Nucleus
    Nucleus nuke(A);
    //GDist *gdist = new DGLAPDist();
    GDist *gdist;
    if (gdist_model==GDIST_DGLAP)
        gdist = new DGLAPDist(xgfile);
    else
        {cerr << "Unknown gdist model" << endl; return -1; }
    if (mode==MODE_GLUEDIST)
    {
        cout << gdist->Gluedist(bjorkx,r*r);
        return 0;
    }
    
    nuke.SetGDist(gdist);    
    Dipxs *amplitude;
    if (model==MODEL_IPSAT)
    {
        amplitude = new Dipxs_IPSat(nuke, IPSAT_MODE_DEFAULT, bp);
        ((Dipxs_IPSat*)amplitude)->SetFactorize(factor_ipsat);
    }
    else if (model==MODEL_IPNONSAT)
        amplitude = new Dipxs_IPNonSat(nuke, bp);
    else if (model==MODEL_IIM)
        amplitude = new Dipxs_IIM(nuke, iim_file);
    else if (model==MODEL_IPSAT_NONSATP)
    {
        amplitude = new Dipxs_IPSat(nuke, IPSAT_MODE_NONSAT_P,bp);
        ((Dipxs_IPSat*)amplitude)->SetFactorize(factor_ipsat);
    }

    if (mode==MODE_XG)  // Print the value of xg(x,µ) and quit
    {
        cout << "# r = " << r << " GeV^(-1), x = " << bjorkx << endl;
        REAL gd = gdist->Gluedist(bjorkx, r*r);
        REAL tmp = gd* 2*NC/(M_PI*M_PI*Alpha_s(Mu2(SQR(r))) );
        cout <<"xg = " << tmp << " (Gluedist = " << gd <<")"  << endl;
        return 0;
    }
    

    Calculator calculator(amplitude, JPsi);

    /*******************
     * \gamma^* N -> J/\Psi N cross section
     * Calculates:   d\sigma / dt = 1/(16*\pi)*
     * \int d^2 r d^2 r' (jpsi)(r)*(jpsi)(r')*4*qqamplitude_sqr_avg(r,delta)
     * Here (jpsi) is the inner product between \gamma^* and J/\Psi (or VM) wave 
     * functions integrated over z \in [0,1]
     */

    if (mode==MODE_TOTXS)  // Calculate total cross section
    {
        cout << "# x = " << bjorkx << ", Q^2=" << Qsqr << " Gev^2" << endl;
        REAL xs=0;
        if (polarization==VM_MODE_TOT) // Sum transversial and longitudinal
                                // polarization components
            {
                JPsi->SetMode(VM_MODE_L);
                xs = calculator.TotalCrossSection(Qsqr, bjorkx);
                JPsi->SetMode(VM_MODE_T);
                xs += calculator.TotalCrossSection(Qsqr, bjorkx);
            }
            else
                xs = calculator.TotalCrossSection(Qsqr, bjorkx);
        cout << "Total cross section: " << xs*NBGEVSQR << " nb" << endl;
    }
     
    else if (mode==MODE_TOTXS_Q)    // Total cross section as a function of Q^2
    {                               // W fixed
        cout << "# Total cross section [nb], W=" << W << " GeV" << endl;
        if (minQsqr==0) minQsqr=0.0001; // Qsqr=0 doesn't work
        REAL multiplier = pow(maxQsqr/minQsqr, 1.0/points);
        #pragma omp parallel for
        for (int i=0; i<=points; i++)
        {
            REAL tmpqsqr = minQsqr*pow(multiplier, i);
            //bjorkx = tmpqsqr/(tmpqsqr+SQR(W))*(1+SQR(M_v)/tmpqsqr);
            //bjorkx = (tmpqsqr + SQR(M_v))/SQR(W);
            bjorkx = (tmpqsqr + SQR(M_v))/(SQR(W)+tmpqsqr);
            REAL xs;
            if (polarization==VM_MODE_TOT) // Sum transversial and longitudinal
                                // polarization components
            {
                JPsi->SetMode(VM_MODE_L);
                xs = calculator.TotalCrossSection(tmpqsqr, bjorkx);
                JPsi->SetMode(VM_MODE_T);
                xs += calculator.TotalCrossSection(tmpqsqr, bjorkx);
            }
            else
                xs = calculator.TotalCrossSection(tmpqsqr, bjorkx);

            xs *= NBGEVSQR;     // 1/Gev^2 -> nb  
            #pragma omp critical
            {
                cout.precision(5);
                cout << fixed << tmpqsqr;
                cout.precision(8);
                cout << " " << xs << endl;
            }
        
        }    
    }
    
    else if (mode==MODE_TOTXS_W)    // Total cross section as a function of W
    {                               // Q^2 fixed
        cout << "# W [GeV]  total cross section [nb], Q^2=" << Qsqr << " GeV" << endl;
        if (minW<0.0001) minW=0.0001; // W=0 doesn't work
        REAL multiplier = pow(maxW/minW, 1.0/points);
        #pragma omp parallel for
        for (int i=0; i<=points; i++)
        {
            REAL tmpw = minW*pow(multiplier, i);
            //bjorkx = tmpqsqr/(tmpqsqr+SQR(W))*(1+SQR(M_v)/tmpqsqr);
            //bjorkx = (tmpqsqr + SQR(M_v))/SQR(W);
            bjorkx = (Qsqr + SQR(M_v))/(SQR(tmpw)+Qsqr);
            REAL xs=0;
            if (polarization==VM_MODE_TOT) // Sum transversial and longitudinal
                                // polarization components
            {
                JPsi->SetMode(VM_MODE_L);
                xs = calculator.TotalCrossSection(Qsqr, bjorkx);
                JPsi->SetMode(VM_MODE_T);
                xs += calculator.TotalCrossSection(Qsqr, bjorkx);
            }
            else
                xs = calculator.TotalCrossSection(Qsqr, bjorkx);

            xs *= NBGEVSQR;     // 1/Gev^2 -> nb  
            #pragma omp critical
            {
                cout.precision(5);
                cout << fixed << tmpw;
                cout.precision(8);
                cout << " " << xs << endl;
            }
        
        }    
    }
    
    // Longitudinal cross section / transversial cross section as a 
    // functio of Q^2 
    else if (mode==MODE_TOTXS_RATIO_Q) 
    {
        cout << "# Longitudnal xs / transversioal xs, W=" << W << " GeV" << endl;
        if (minQsqr==0) minQsqr=0.0001; // Qsqr=0 doesn't work
        REAL multiplier = pow(maxQsqr/minQsqr, 1.0/points);
        #pragma omp parallel for
        for (int i=0; i<=points; i++)
        {
            REAL tmpqsqr = minQsqr*pow(multiplier, i);
            bjorkx = (tmpqsqr + SQR(M_v))/(SQR(W)+tmpqsqr);
            
            JPsi->SetMode(VM_MODE_L);
            REAL xsl = calculator.TotalCrossSection(tmpqsqr, bjorkx);
            JPsi->SetMode(VM_MODE_T);
            REAL xst = calculator.TotalCrossSection(tmpqsqr, bjorkx);
           
            #pragma omp critical
            {
                cout.precision(5);
                cout << fixed << tmpqsqr;
                cout.precision(8);
                cout << " " << xsl/xst << endl;
            }
        
        }
    }
    
    else if (mode==MODE_Ap)    // Calculate d\sigma^A/dt / A*d\sigma_p/dt as a function Q^2 at t
    {         
        cout << "# t=" << t << ", x = " << bjorkx <<  endl;
        cout << "# Q^2   nucleus_xs / A*proton_xs" << endl;
        if (minQsqr==0) minQsqr=0.0001; // Qsqr=0 doesn't work
        REAL multiplier = pow(maxQsqr/minQsqr, 1.0/points);
        std::cerr<< "TODO: Fix polarization sums" << endl;
        #pragma omp parallel for
        for (int i=0; i<=points; i++)
        {
            //REAL tmpqsqr=maxQsqr/points*i;
            REAL tmpqsqr = minQsqr*pow(multiplier, i);
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
        std::cerr<< "TODO: Fix polarization sums" << endl;
        if (minx==0) minx=1e-6; // x=0 doesn't work
        REAL multiplier = pow(maxx/minx, 1.0/points);
            
        #pragma omp parallel for
        for (int i=0; i<points; i++)
        {
            //REAL tmpx = minx + (maxx-minx)/points*i;
            REAL tmpx = minx*pow(multiplier,i);
            REAL protonxs = calculator.ProtonCrossSection_dt(t, Qsqr, tmpx);
            REAL nukexs = calculator.CrossSection_dt(t, Qsqr, tmpx);
            #pragma omp critical
            {
                cout.precision(8);
                cout << fixed << tmpx;
                cout.precision(8);
                cout << " " << nukexs / (A*protonxs) << endl;
            }
        
        }
    
    }
    
    else if (mode==MODE_COHERENT_DT)   // d\sigma^A/dt for coherent scattering
    {
        if (A==1) { 
            std::cerr << "A=1, can't be coherent scattering." << std::endl;
            return -1;
        }
        cout << "# d\\sigma/dt [1/GeV^4], coherent scattering " << endl;
        cout << "# x_pomeron = " << bjorkx << ", Q^2 = " << Qsqr << endl;
        std::cerr<< "TODO: Fix polarization sums" << endl;
        // All iterations are independent, so this is straightforward to parallerize   
        #pragma omp parallel for
        for (int i=0; i<=points; i++)
        {
            REAL tmpt = (maxt-mint)/points*i;
            REAL delta = sqrt(tmpt);
            REAL result;
            result = calculator.CoherentCrossSection_dt(tmpt, Qsqr, bjorkx);
            #pragma omp critical
            {
                cout.precision(5);
                cout << fixed << tmpt;
                cout.precision(8);
                cout << " " << result << endl;
            }
        }
    
    
    }

    
    else if (mode==MODE_DIFFXS)   // dsigma/dt as a function of t
    {
        cout << "# d\\sigma/dt [nb/GeV^2] " << endl;
        cout << "# x_pomeron = " << bjorkx << ", Q^2 = " << Qsqr << endl;
        // All iterations are independent, so this is straightforward to parallerize   
        #pragma omp parallel for
        for (int i=0; i<=points; i++)
        {
            REAL tmpt = mint+(maxt-mint)/points*i;
            REAL delta = sqrt(tmpt);
            REAL result=0;
            if (A==1)
            {
                if (polarization==VM_MODE_TOT) 
                {
                    JPsi->SetMode(VM_MODE_L);
                    result = calculator.ProtonCrossSection_dt(tmpt, Qsqr, bjorkx);
                    JPsi->SetMode(VM_MODE_T);
                    result += calculator.ProtonCrossSection_dt(tmpt, Qsqr, bjorkx);
                }
                else
                    result = calculator.ProtonCrossSection_dt(tmpt, Qsqr, bjorkx); 
            }
            else
            {
                if (polarization==VM_MODE_TOT) 
                {
                    JPsi->SetMode(VM_MODE_L);
                    result = calculator.CrossSection_dt(tmpt, Qsqr, bjorkx);
                    JPsi->SetMode(VM_MODE_T);
                    result += calculator.CrossSection_dt(tmpt, Qsqr, bjorkx);
                }
                else
                    result = calculator.CrossSection_dt(tmpt, Qsqr, bjorkx); 
            }
            #pragma omp critical
            {
                cout.precision(5);
                cout << fixed << tmpt;
                cout.precision(8);
                cout << " " << result*NBGEVSQR << endl;
            }
        }
    }
    
    else if (mode==MODE_VM_INTZ)
    {
        cout << "# 2\\pi r * \\int dz/(4\\pi) r Psi^*Psi, Q^2 = " << Qsqr << endl;
        cout << "# " << *((GausLC*)JPsi) << endl;
        if (output_fm) cout << "# [r] = fm"; else cout << "# [r] = GeV^(-1)";
        cout << endl;
        REAL maxr=5;
        for (int i=1; i<=points; i++)
        {
            REAL tmpr = (maxr-MINR)/points*i;
            REAL val=tmpr*JPsi->PsiSqr_intz(Qsqr, tmpr);
            
            // if we want r axis to be in units of fm
            // Note: [\int d^2 r \int dz \Psi^*\Psi]=1, so we don't have to 
            // multiply val by FMGEV
            if (output_fm) { tmpr/=FMGEV; }
            
            cout << tmpr << " " << 2*M_PI*val << endl;
            
        
        }
    }
    
    else if (mode==MODE_DSIGMA_D2B_R) 
    {
        cout << "# d\\sigma/d^2b for dipole-proton scattering" << endl;
        cout << "# b = " << b << ", x = " << bjorkx << endl;
        REAL maxr = 10; REAL minr=0.01;
        REAL multiplier = pow(maxr/minr, 1.0/points);
        
        for (int i=1; i<=points; i++)
        {
            REAL tmpr = minr*pow(multiplier, i);
            cout << tmpr << " " << 2.0*amplitude->Qq_proton_amplitude(tmpr*tmpr, 
                                    bjorkx, b) << endl;
        }    
    }
    
    else if (mode==MODE_DSIGMA_D2B_B) 
    {
        cout << "# d\\sigma/d^2b for dipole-proton scattering" << endl;
        cout << "# r = " << r << " GeV^(-1), x = " << bjorkx << endl;
        REAL minb=0; REAL maxb=10;
        for (int i=0; i<=points; i++)
        {
            REAL tmpb = minb+(maxb-minb)/points*i;
            cout << tmpb << " " << 2.0*amplitude->Qq_proton_amplitude(SQR(r), 
                                    bjorkx, tmpb) << endl;
        }    
    }
    
    else if (mode==MODE_DSIGMA_DT_QQ)
    {
        cout << "#t [Gev^2]  d\\sigma_{q\\bar q}/dt  [Gev^(-2)]" << endl;
        cout << "#r=" << r << "GeV^(-1), x=" << bjorkx << endl;
    
        for (int i=0; i<points; i++)
        {
            REAL tmpt = mint+(maxt-mint)/points*i;
            cout << tmpt << " ";
            REAL xs= 2.0*amplitude->DipoleAmplitude_proton(SQR(r), bjorkx, 
                        sqrt(tmpt));
            xs = SQR(xs);///(16*M_PI);  // Without 16\pi so we can compare with
                                // KT article, PRD68, 114005 hep-ph/0304189
            cout << xs << endl;
        }
    
    }
    else if (mode==MODE_TOTAL_DIPXS)
    {
        cout << "#r [Gev^(-1)]  \\sigma [mb] (total dipole-proton cross section)" << endl;
        cout << "# x=" << bjorkx << endl;
        REAL minr=0.1; REAL maxr=15;
        for (int i=0; i<=points; i++)
        {
            REAL tmpr = minr+(maxr-minr)/points*i;
            cout << tmpr << " ";
            REAL xs= amplitude->TotalDipxsection_proton(SQR(tmpr), bjorkx);
            xs = xs*NBGEVSQR / 1000.0 / 1000.0;  // Gev^(-2) => mb
            cout << xs << endl;
        }
    
    }
    
    delete amplitude;
    delete gdist;
    delete JPsi;
   
    return 0;
}



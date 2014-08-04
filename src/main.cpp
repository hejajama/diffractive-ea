/*
 * Calculates eA cross section in dipole model
 * 
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010-2013
 */
 
#include <iostream>
#include <iomanip>  // Set cout precision
#include <ctime>
#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include "dipole.h"
#include "gaus_lc.h"
#include "gauss_boost.h"
#include "vector.h"
#include "nucleus.h"
#include "dipxs_ipnonsat.h"
#include "dipxs_ipsat.h"
#include "dipxs_ipsat2012.h"
#include "dipxs_iim.h"
#include "dipxs_bk.h"
#include "gdist/gdist_dglap.h"
#include "mersenne/mersenne.h"
#include "calculator.h"

using namespace std;

const std::string VERSION = "1.02-dev";
const std::string DATE = "2013-xx-xx";

const int MODEL_IPSAT=1; const int MODEL_IPNONSAT=2; const int MODEL_IIM=3;
const int MODEL_IPSAT_NONSATP=4; const int MODEL_IPSAT_NOFACTOR=5;
const int MODEL_BK=6; const int MODEL_IPSAT2012=7;
const int GDIST_DGLAP=1; 

const int WAVEF_GAUS_LC=1; const int WAVEF_BOOSTED_GAUSSIAN=2;

enum MODE
{
        MODE_TOTXS, MODE_DIFFXS, MODE_Ap,
        MODE_Ap_X, MODE_Ap_A, MODE_Ap_A_COH, MODE_TOTXS_Q, MODE_COHERENT_DT,
        MODE_VM_INTZ, MODE_TOTXS_RATIO_Q, MODE_DSIGMA_D2B_R,
        MODE_TOTXS_W, MODE_DSIGMA_DT_QQ, MODE_DSIGMA_D2B_B, MODE_TOTAL_DIPXS,
        MODE_XG, MODE_GLUEDIST, MODE_XG_X, MODE_XG_R, MODE_DSIGMA_DT_A,
        MODE_QUASIELASTIC_COHERENT_Q, MODE_QUASIELASTIC_COHERENT_X,
        MODE_QUASIELASTIC_COHERENT_A,
        MODE_COHERENT_AA, MODE_INCOHERENT_AA,
};


void Cleanup(); 

int main(int argc, char* argv[])
{    
    gsl_set_error_handler(&ErrHandler);   // We really don't want to use this!
    
    // Print the cmdline args
    cout << "# ";
    for (int i=0; i<argc; i++)
        cout << argv[i] << " " ;
    cout << endl;
    
    // Parameters
    REAL bjorkx=1e-4;
    REAL Qsqr=1;
    int A=197;
    int minA=10; int maxA=300;  // minA and maxA from A-dep. calculations
    int Astep=1;  
    int model=MODEL_IPSAT;
    int gdist_model=GDIST_DGLAP;
    int points=100;
    REAL maxt=0.3; REAL mint=0; 
    REAL t=0.5;
    REAL maxQsqr=-1;
    REAL minQsqr=-1;
    REAL minW=0;
    REAL maxW=200;
    MODE mode=MODE_DIFFXS;   // What to do
    REAL minx=1e-6;
    REAL maxx=1e-2;
    string xgfile="xg.dat";
    REAL W=-1;  // W^2 = (P+q)^2 = invariant mass^2 of the \gamma^*N system
    REAL M_v=3.097; // Mass of the produced vector meson, 3.097 GeV = J/\Psi
    bool x_set=false;
    bool w_set=false;
    bool scalex=false;  // Scale x so that x=x_bj*(1+M_V^2/Q^2)
    bool output_fm=false;   // Use fm's when printing the value of r
    bool corrections=true;
    REAL r=-1;
    bool pa=false;	// pA instead of AA collision
    REAL M_n=0;     // Mass of nucleus/nucleon
    REAL bp=DEFAULT_B_p;
    int polarization = VM_MODE_TOT;
    REAL b=0;   // Impact parameter
    int wavef=WAVEF_BOOSTED_GAUSSIAN;
    double sqrts=2760;	// AA and pA events
    std::string waveffile="";   // file where wavefunction parameters are read, empty=default
    
    string iim_file="iim.dat";  // Read parameters for IIM model from this file
    string bk_file="";
    REAL bk_sigma0=0;   // sigma0 in GeV^-2, required for BK dipole amplitude
            
    // Parse parameters
    if (argc>1)
    {
        if (string(argv[1])=="--help" or string(argv[1]) == "-help")
        {
            cout << "Usage: -x x_pomeron -Q2 Q^2 -W W (specify only x or W), -sqrts sqrts" << endl;
            cout << "-scalex (scale x by factor 1+M_V^2/Q^2) (can't be used with -W)" << endl;
            cout << "-dipole {ipsat,ipnonsat,iim,ipsat_nonsatp,ipsat-nofactor,bk,ipsat2012} [bkfilename bksigma0]" << endl;
            cout << "-gdist {dglap} -xgfile file" << endl;
            cout << "-wavef {gaus-lc, boosted-gaussian} (specify VM wave function)" << endl;
            cout << "-wavef_file filename: file where wavefunction parameters is read" << endl;
            cout << "-A number_of_nucleai -Mn nucleus_mass" << endl;
            cout << "-minA A, -maxA A, -Astep n (min and max mass of the nucleai, step in computations)" << endl;
            cout << "-N number_of_data_points" << endl;
            cout << "-coherent_dt (calculate coherent d\\sigma/dt)" << endl;
            cout << "-Bp proton_shape (for ipsat and ipnonsat)" << endl;
            cout << "-mint t_value, -maxt t_value" << endl;
            cout << "-iimfile filename (parameters for the IIM model)" << endl;
            cout << "-totxs (calculates total cross section)" << endl;
            cout << "-totxs_q2 (total cross section as a function of Q^2)" << endl;
            cout << "-totxs_w (total cross section as a function of W)" << endl;
            cout << "-totxs_q2_l/t (longitudinal total xs / transversial total xs)" << endl;
            cout << "-A/p, -A/p_x, -A/p_A, -A/p_A_coh (nucleus cross section / A* " << endl;
            cout << "     proton cross section as a function of Q^2, x or A, coh=coherent scatt.)" << endl;
            cout << "-t t, -maxQ2 maxq2 -minQ2 minQ2[GeV] (value of t and max/min of Q^2)" << endl;
            cout << "-minW minW, -maxW maxW [GeV] (max/min of W)" << endl;
            cout << "-xg (print the value of xg(x,r) and quit) " << endl;
            cout << "-xg_x (print xg as a function of x) " << endl;
			cout << "-xg_r (print xg as a function of r) " << endl;
            cout << "-gdistval (print ONLY gluedist and quit)" << endl;
            cout << "-vm_intz (print \\int d^z/(4\\pi) r Psi^*Psi) -fm (print r in fm)" << endl;
            cout << "-pol {l, t, sum} (polarization of the VM, sum = l + t is default)" << endl;
            cout << "-dsigma/d2b_r impact_par (print d\\sgima_{q\bar q}/d^2 b as a function of r)" << endl;
            cout << "-dsigma/d2b_b (print d\\sigma_{q\\bar q}/d^2 b as a function of b)" << endl;
            cout << "-dsigma/dt_qq (print d\\sigma_{q\bar q}/dt as a function of t), -r r (in GeV)" << endl;
            cout << "-totdipxs (print total dipole-proton cross section as a function of r)" << endl;
            cout << "-nocorrections (don't calculate real part and skewedness corrections)" << endl;
            cout << "-dsigma/dt-A  (quasielastic cross section as a function of A "
                << "[" << minA << " - " << maxA << "])" << endl;
            cout << "-quasielastic/coherent_x, -quasielastic/coherent_q2, -quasielastic/coherent_A" << endl
                << "   (\\int quasieal. from mint to maxt / total coherent) as a function of x, A or Q)" << endl;
            cout << "-coherent_AA, -incoherent_AA: compute d\\sigma/dy in coherent AA -> j/\\psi + AA" << endl;
            cout << "-pA: compute pA instead of AA, allways incoherent!" << endl;
            cout << "-v, -version (print version and quit)" << endl;
            cout << endl;
            cout << "Default values: x="<<bjorkx <<", Q^2="<<Qsqr 
                << " A="<<A<<", N="<<points<<", mint="<<mint<<", maxt="<<maxt<< endl;
            cout << "                A/p=false, iimfile=" << iim_file << endl;
            cout << "                minQ^2=" << minQsqr << ", maxQ^2=" << maxQsqr << endl;
            cout << "                t="<<t << ", minx=" << minx << ", maxx=" << maxx 
                 << " xgfile=" << xgfile << endl;
            cout << "                Mv=" << M_v << ", Mn=" << M_n << ", B_p=" << bp << endl;
            cout << "                minW=" << minW << ", maxW=" << maxW << endl;
            cout << "                wavef=gaus-lc" << endl;

            return 0;
        }
        for (int i=1; i<argc; i++)
        {
            if (string(argv[i])=="-version" or string(argv[i])=="-v")
            {
                std::cout << "Dipole v. " << VERSION << " (" << DATE << ")" << endl;
                return 0;
            }
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
            else if (string(argv[i])=="-sqrts")
				sqrts=StrToReal(argv[i+1]);
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
            else if (string(argv[i])=="-minA")
                minA=StrToInt(argv[i+1]);
            else if (string(argv[i])=="-maxA")
                maxA=StrToInt(argv[i+1]);
            else if (string(argv[i])=="-Astep")
                Astep=StrToInt(argv[i+1]);
            else if (string(argv[i])=="-coherent_dt")
                mode=MODE_COHERENT_DT;
            else if (string(argv[i])=="-nofactor")
                cerr << "Option -nofactor is deprecated , ignoring" << endl;
            else if (string(argv[i])=="-totxs")
                mode=MODE_TOTXS;
            else if (string(argv[i])=="-A/p")
                mode=MODE_Ap;
            else if (string(argv[i])=="-A/p_x")
                mode=MODE_Ap_X;
            else if (string(argv[i])=="-A/p_A")
                mode=MODE_Ap_A;
            else if (string(argv[i])=="-A/p_A_coh")
                mode=MODE_Ap_A_COH;
            else if (string(argv[i])=="-totxs_q2")
                mode=MODE_TOTXS_Q;
            else if (string(argv[i])=="-totxs_q2_l/t")
                mode=MODE_TOTXS_RATIO_Q;
            else if (string(argv[i])=="-totxs_w")
                mode=MODE_TOTXS_W;
            else if (string(argv[i])=="-coherent_AA")
				mode=MODE_COHERENT_AA;
			else if (string(argv[i])=="-incoherent_AA")
				mode=MODE_INCOHERENT_AA;
			else if (string(argv[i])=="-pA" or string(argv[i])=="-pa")
				pa=true;
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
            else if (string(argv[i])=="-Mn")
                M_n=StrToReal(argv[i+1]);
            else if (string(argv[i])=="-Bp")
                bp=StrToReal(argv[i+1]);
            else if (string(argv[i])=="-xg") 
                mode=MODE_XG;
            else if (string(argv[i])=="-xg_x")
                mode=MODE_XG_X;
			else if (string(argv[i])=="-xg_r")
				mode=MODE_XG_R;
            else if (string(argv[i])=="-gdistval")
                mode=MODE_GLUEDIST;
            else if (string(argv[i])=="-scalex")
                scalex=true;
            else if (string(argv[i])=="-vm_intz")
                mode=MODE_VM_INTZ;
            else if (string(argv[i])=="-dsigma/dt-A")
                mode=MODE_DSIGMA_DT_A;
            else if (string(argv[i])=="-fm")
                output_fm=true;
            else if (string(argv[i])=="-nocorrections")
                corrections=false;
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
            else if (string(argv[i])=="-quasielastic/coherent_q2")
                mode=MODE_QUASIELASTIC_COHERENT_Q;
            else if (string(argv[i])=="-quasielastic/coherent_x")
                mode=MODE_QUASIELASTIC_COHERENT_X;
            else if (string(argv[i])=="-quasielastic/coherent_A")
                mode=MODE_QUASIELASTIC_COHERENT_A;
            else if (string(argv[i])=="-dipole")
            {
                if (string(argv[i+1])=="ipsat" or string(argv[i+1])=="ipsat06")
                    model=MODEL_IPSAT;
                else if (string(argv[i+1])=="ipsat-nofactor")
                    model=MODEL_IPSAT_NOFACTOR;
                else if (string(argv[i+1])=="ipnonsat")
                    model=MODEL_IPNONSAT;
                else if (string(argv[i+1])=="ipsat2012" or string(argv[i+1])=="ipsat12")
                    model=MODEL_IPSAT2012;
                else if (string(argv[i+1])=="iim")
                    model=MODEL_IIM;
                else if (string(argv[i+1])=="ipsat_nonsatp")
                    model=MODEL_IPSAT_NONSATP;
                else if (string(argv[i+1])=="bk")
                {
                    model=MODEL_BK;
                    bk_file = string(argv[i+2]);
                    bk_sigma0 = StrToReal(argv[i+3]);
                }
                else
                {
                    cerr << "Model " << argv[i+1] << " is not valid" << endl;
                    return -1;
                }
            }
            else if (string(argv[i])=="-wavef_file")
            {
                waveffile = argv[i+1];
                    
            }
            else if (string(argv[i])=="-wavef")
            {
                if (string(argv[i+1])=="gaus-lc")
                    wavef=WAVEF_GAUS_LC;
                else if (string(argv[i+1])=="boosted-gaussian")
                    wavef=WAVEF_BOOSTED_GAUSSIAN;
                else
                {
                    cerr << "Wavef " << argv[i+1] << " is not valid" << endl;
                    return -1;
                }
            
            }
            else if (string(argv[i])=="-gdist")
            {
                if (string(argv[i+1])=="dglap")
                    gdist_model=GDIST_DGLAP;
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
        //std::cout << "# x scaled from " << bjorkx << " to " ;
        bjorkx = (Qsqr + SQR(M_v))/( Qsqr + SQR(W) + SQR(M_n) );
        //cout << bjorkx << endl;
    }

    if (!w_set and scalex) // Scale x->x*(1+M^2/Q^2)
    {
        cerr << "Check scalex if you need it!" << endl;
        //if (Qsqr<0.00001){ cerr << "Q^2=0, can't scale x" << endl; return -1;}
        //bjorkx = bjorkx * (1 + SQR(M_v)/Qsqr);
    }

    
    // Print values
    cout << "# A=" << A;
    if (A==1) cout << " (dipole-proton)";
    cout << endl;
    cout << "# GDist=" << gdist_model << ",  dipole model=" << model 
        << " wavef=" << wavef << " pol=" << polarization << endl;
    
    // Intialize random number generator
    seed_mersenne(std::time(NULL));
    
    // Wave function
    WaveFunction *JPsi;
    std::string fname;
    switch (wavef)
    {
        case WAVEF_GAUS_LC:
           // J/Psi wave function:  e_f, N_T, N_L, R_T, R_L, m_f, M_V, delta
           //VM_Photon JPsi(2.0/3.0, 1.23, 0.83, sqrt(6.5), sqrt(3.0), 1.4, 3.097, 1);
           fname="gaus-lc.dat";
           if (waveffile != "")
                fname=waveffile;
           JPsi = new GausLC(fname);
           cout << "# Wave function: " << *((GausLC*)JPsi) << endl;
           break;
        case WAVEF_BOOSTED_GAUSSIAN:
            fname = "gauss-boosted.dat";
            if (waveffile != "")
                fname=waveffile;
            JPsi = new BoostedGauss(fname);
            cout << "# Wave function: " << *((BoostedGauss*)JPsi) << endl;
            break;
        default:
            cerr << "Unknown wave function set! Quitting...." << std::endl;
            return -1;
    }

    M_v = JPsi->MesonMass();
    JPsi->SetMode(polarization);
    cout << "# Produced meson mass: " << M_v << endl;
    // Intialize Dipxs and Nucleus
    Nucleus nuke(A);
    //GDist *gdist = new DGLAPDist();
    GDist *gdist;
    if (gdist_model==GDIST_DGLAP)
        gdist = new DGLAPDist(xgfile);
    else
        {cerr << "Unknown gdist model" << endl; return -1; }
    
    nuke.SetGDist(gdist);    
    Dipxs *amplitude;
    switch (model)
    { 
    case MODEL_IPSAT:
        amplitude = new Dipxs_IPSat(nuke, IPSAT_MODE_DEFAULT, bp);
        ((Dipxs_IPSat*)amplitude)->SetFactorize(true);
        break;
    case MODEL_IPSAT_NOFACTOR:
        amplitude = new Dipxs_IPSat(nuke, IPSAT_MODE_DEFAULT, bp);
        ((Dipxs_IPSat*)amplitude)->SetFactorize(false);
        break;
    case MODEL_IPNONSAT:
        amplitude = new Dipxs_IPNonSat(nuke, bp);
        break;
    case MODEL_IPSAT2012:
        amplitude = new Dipxs_IPSat2012(nuke);
        ((Dipxs_IPSat2012*)amplitude)->SetFactorize(true);
        break;
    case MODEL_IIM:
        amplitude = new Dipxs_IIM(nuke, iim_file);
        break;
    case MODEL_IPSAT_NONSATP:
        amplitude = new Dipxs_IPSat(nuke, IPSAT_MODE_NONSAT_P,bp);
        ((Dipxs_IPSat*)amplitude)->SetFactorize(true);
        break;
    case MODEL_BK:
        amplitude = new Dipxs_BK(nuke, bk_file);
        ((Dipxs_BK*)amplitude)->SetSigma0(bk_sigma0);
        break;
    }

    if (Qsqr < 0.00001 and polarization != VM_MODE_T)
    {
        polarization = VM_MODE_T;
        std::cerr << "#Q^2=0 -> only transverse component" << endl;
    }

    Calculator calculator(amplitude, JPsi);
    calculator.SetPolarization(polarization);
    calculator.SetCorrections(corrections);
    if (corrections==false)
        cout << "# Ignoring corrections " << endl;
    
/////////////////////////////////////////////////////////////////////////////
////////////// Different modes //////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

    if (mode==MODE_XG)  // Print the value of xg(x,µ) and quit
    {
        cout << "# r = " << r << " GeV^(-1), x = " << bjorkx << endl;
        REAL gd = gdist->Gluedist(bjorkx, r*r);
        REAL tmp = gd* 2*NC/(M_PI*M_PI*Alpha_s(Mu2(SQR(r))) );
        cout <<"xg = " << tmp << " (Gluedist = " << gd <<")"  << endl;
    }
    
    else if (mode==MODE_XG_X) // xg as a function of x
    {
        if (minx<1e-7) minx=1e-7;
        REAL m = pow(maxx/minx, 1.0/points);
        cout << "# x   xg(x,r)    r=" << r << " GeV^(-1)" << endl;
        cout << "# mu_0^2=" << mu2_0 << " GeV^2, \\lambda_{QCD}^2=" 
            << sqrt(LAMBDAQCD2) << " GeV" << endl;
        for (int i=0; i<=points; i++)
        {
            REAL tmpx = minx*pow(m,i);
            REAL gd = gdist->Gluedist(tmpx, SQR(r));
            REAL res = gd*2.0*NC/(M_PI*M_PI*Alpha_s(Mu2(SQR(r))) );
            cout << tmpx << " " << res << endl;
        }   
    }
	else if (mode==MODE_XG_R)
	{
		REAL maxr=50; double minr=1e-8;
		REAL m = pow(maxr/minr, 1.0/points);
		cout <<"# r   xg(x,r)   x=" << bjorkx << endl;
		for (int i=0; i<points; i++)
		{
			REAL tmpr = minr*pow(m, i);
			REAL gd = gdist->Gluedist(bjorkx, SQR(tmpr));
			REAL res = gd*2.0*NC/(M_PI*M_PI*Alpha_s(Mu2(SQR(r))) );
			cout << tmpr << " " << res << endl;
		} 		
	}
    
    else if (mode==MODE_GLUEDIST)
    {
        cout << gdist->Gluedist(bjorkx,r*r);
    }   

    /*******************
     * \gamma^* N -> J/\Psi N cross section
     * Calculates:   d\sigma / dt = 1/(16*\pi)*
     * \int d^2 r d^2 r' (jpsi)(r)*(jpsi)(r')*4*qqamplitude_sqr_avg(r,delta)
     * Here (jpsi) is the inner product between \gamma^* and J/\Psi (or VM) wave 
     * functions integrated over z \in [0,1]
     */

    else if (mode==MODE_TOTXS)  // Calculate total cross section
    {
        cout << "# x = " << bjorkx << ", Q^2=" << Qsqr << " Gev^2" << endl;
        REAL xs=0;
        xs = calculator.TotalCrossSection(Qsqr, bjorkx);
        cout << "Total cross section: " << xs*NBGEVSQR << " nb" << endl;
    }
     
    else if (mode==MODE_TOTXS_Q)    // Total cross section as a function of Q^2
    {                               // W fixed
        cout << "#Q^2 [GeV^2]   Total cross section [nb], W=" << W << " GeV" << endl;
        if (minQsqr==0) minQsqr=0.0001; // Qsqr=0 doesn't work
        REAL multiplier = pow(maxQsqr/minQsqr, 1.0/points);
        ////#pragma omp parallel for
        for (int i=0; i<=points; i++)
        {
            REAL tmpqsqr = minQsqr*pow(multiplier, i);
            //bjorkx = tmpqsqr/(tmpqsqr+SQR(W))*(1+SQR(M_v)/tmpqsqr);
            //bjorkx = (tmpqsqr + SQR(M_v))/SQR(W);
            bjorkx = (tmpqsqr + SQR(M_v))/(SQR(W)+tmpqsqr);

            REAL xs = 0;
            if (A==1)
                xs = 1.0/calculator.GetAmplitude()->Bp() * calculator.ProtonCrossSection_dt(0, tmpqsqr, bjorkx);
            else
                xs =  calculator.TotalCrossSection(tmpqsqr, bjorkx);

            xs *= NBGEVSQR;     // 1/Gev^2 -> nb  
            //#pragma omp critical
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
        cout << "# W [GeV]  total cross section [nb]; " << endl;
        if (minQsqr>=0)
            cout << "# minQ^2: " << minQsqr << ", maxQ^2: " << maxQsqr << endl;
        else cout << "# Q^2=" << Qsqr << " GeV" << endl;
        if (minW<0.0001) minW=10; // small W -> huge x -> model doesn't work
        //REAL multiplier = pow(maxW/minW, 1.0/points);
        ////#pragma omp parallel for
        for (int i=0; i<=points; i++)
        {
            //REAL tmpw = minW*pow(multiplier, i);
            REAL tmpw = minW + (maxW-minW)/points*i;
            bjorkx = (Qsqr + SQR(M_v))/(SQR(tmpw)+Qsqr);
            REAL xs=0;
            if (A==1)
            {
                if (minQsqr>=0)  // average
                {
                    xs = 1.0/calculator.GetAmplitude()->Bp() * calculator.ProtonCrossSection_dt_qsqravg(0, minQsqr, maxQsqr, tmpw, M_v)*NBGEVSQR;
                }
                else
                    xs = 1.0/calculator.GetAmplitude()->Bp() * calculator.ProtonCrossSection_dt(0, Qsqr, bjorkx)*NBGEVSQR;
            }
            else
                xs = calculator.TotalCrossSection(Qsqr, bjorkx)*NBGEVSQR;

            //#pragma omp critical
            {
                cout.precision(5);
                cout << fixed << tmpw;
                cout.precision(8);
                cout << " " << xs <<  endl;
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
        //#pragma omp parallel for
        for (int i=0; i<=points; i++)
        {
            REAL tmpqsqr = minQsqr*pow(multiplier, i);
            bjorkx = (tmpqsqr + SQR(M_v))/(SQR(W)+tmpqsqr);
            
            calculator.SetPolarization(VM_MODE_L);
            REAL xsl = calculator.TotalCrossSection(tmpqsqr, bjorkx);
            calculator.SetPolarization(VM_MODE_T);
            REAL xst = calculator.TotalCrossSection(tmpqsqr, bjorkx);
           
            //#pragma omp critical
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
        cout << "# t=" << t << ", x = " << bjorkx << " pol = " << polarization
            << endl;
        cout << "# Q^2   nucleus_xs / A*proton_xs" << endl;
        calculator.SetCorrections(false);
        if (minQsqr==0) minQsqr=0.0001; // Qsqr=0 doesn't work
        REAL multiplier = pow(maxQsqr/minQsqr, 1.0/points);
        //#pragma omp parallel for
        for (int i=0; i<=points; i++)
        {
            //REAL tmpqsqr=maxQsqr/points*i;
            REAL tmpqsqr = minQsqr*pow(multiplier, i);
            REAL protonxs = calculator.ProtonCrossSection_dt(t, tmpqsqr, bjorkx);
            REAL nukexs = calculator.CrossSection_dt(t, tmpqsqr, bjorkx);
            //#pragma omp critical
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
        calculator.SetCorrections(false);
        if (minx==0) minx=1e-6; // x=0 doesn't work
        REAL multiplier = pow(maxx/minx, 1.0/points);
            
        //#pragma omp parallel for
        for (int i=0; i<points; i++)
        {
            //REAL tmpx = minx + (maxx-minx)/points*i;
            REAL tmpx = minx*pow(multiplier,i);
            REAL protonxs = calculator.ProtonCrossSection_dt(t, Qsqr, tmpx);
            REAL nukexs = calculator.CrossSection_dt(t, Qsqr, tmpx);
            //#pragma omp critical
            {
                cout.precision(8);
                cout << fixed << tmpx;
                cout.precision(8);
                cout << " " << nukexs / (A*protonxs) << endl;
            }
        
        }
    }
    
    else if (mode==MODE_Ap_A)   // d\sigma^A/dt / A*d\sigma_p/dt as a function of A
    {
        cout << "# t=" << t << ", Q^2=" << Qsqr << ", x=" << bjorkx << endl;
        cout << "# A-dependence of nucleus_xs / A*proton_xs" << endl;
        calculator.SetCorrections(false);
            
        for (int tmpA=minA; tmpA<=maxA; tmpA+=Astep)
        {
            Nucleus tmpnuke(tmpA);
            tmpnuke.SetGDist(gdist);
            amplitude->SetNucleus(tmpnuke);
            REAL protonxs = calculator.ProtonCrossSection_dt(t, Qsqr, bjorkx);
            REAL nukexs = calculator.CrossSection_dt(t, Qsqr, bjorkx);
            cout << tmpA;
            cout.precision(8);
            cout << " " << nukexs / (tmpA*protonxs) << endl;
        }
    }
    
    else if (mode==MODE_Ap_A_COH)   // d\sigma^A/dt / A^2*d\sigma_p/dt as a function of A
    {
        // Same as above, but now for coherent scattering
        cout << "# t=" << t << ", Q^2=" << Qsqr << ", x=" << bjorkx << endl;
        cout << "# A-dependence of coherent nucleus_xs / A^2*proton_xs" << endl;
        calculator.SetCorrections(false);
            
        for (int tmpA=minA; tmpA<=maxA; tmpA+=Astep)
        {
            Nucleus tmpnuke(tmpA);
            tmpnuke.SetGDist(gdist);
            amplitude->SetNucleus(tmpnuke);
            REAL protonxs = calculator.ProtonCrossSection_dt(t, Qsqr, bjorkx);
            REAL nukexs = calculator.CoherentCrossSection_dt(t, Qsqr, bjorkx);
            cout << tmpA;
            cout.precision(8);
            cout << " " << nukexs / (tmpA*tmpA*protonxs) << endl;
        }
    }
    
    
    
    
    // \\int quasielstic from minQsqr to maxQsqr / \\int coherent
    else if (mode==MODE_QUASIELASTIC_COHERENT_Q)
    {
        cout << "# \\int quasielstic from mint to maxt/ \\int coherent" 
            << endl;
        cout << "# mint=" << mint << ", maxt=" << maxt << ", x=" << bjorkx << endl;
       
        if (minQsqr<0.000001) minQsqr=0.01;
        REAL multiplier = pow(maxQsqr/minQsqr, 1.0/points);
        calculator.SetTAccuracy(0.01);
        
        if (mint < 0.1)
            cerr << "mint=" << mint << ", are you sure?" << endl;
        
        for (int i=0; i<=points; i++)
        {
            REAL tmpQsqr = minQsqr*pow(multiplier, i);
            REAL coherent = calculator.TotalCoherentCrossSection(tmpQsqr, 
                                                                bjorkx);
            REAL quasiel = calculator.TotalCrossSection(tmpQsqr, bjorkx, 
                                                            mint, maxt);
            cout << tmpQsqr << " " << quasiel/coherent << endl;
        }
    }
    
    // \\int quasielstic from minQsqr to maxQsqr / \\int coherent
    else if (mode==MODE_QUASIELASTIC_COHERENT_X)
    {
        cout << "# \\int quasielstic from mint to maxt / \\int coherent" 
            << endl;
        cout << "# Q^2=" << Qsqr << endl;
        
        if (minx<1e-8) minx=1e-8;
        REAL multiplier = pow(maxx/minx, 1.0/points);
        calculator.SetTAccuracy(0.01);
        
        if (mint < 0.1)
            cerr << "mint=" << mint << ", are you sure?" << endl;
        
        for (int i=0; i<=points; i++)
        {
            REAL tmpx = minx*pow(multiplier, i);
            REAL coherent = calculator.TotalCoherentCrossSection(Qsqr, 
                                                                tmpx);
            REAL quasiel = calculator.TotalCrossSection(Qsqr, tmpx, 
                                                            mint, maxt);
            cout << tmpx << " " << quasiel/coherent << endl;
        }
    }
    
    else if (mode==MODE_QUASIELASTIC_COHERENT_A)
    {
        cout << "# \\int quasielstic from mint to maxt / \\int coherent" 
            << endl;
        cout << "# Q^2=" << Qsqr << ", bjorkx="<< bjorkx << endl;
        

        calculator.SetTAccuracy(0.01);
        
        if (mint < 0.1)
            cerr << "mint=" << mint << ", are you sure?" << endl;
        
        for (int tmpA=minA; tmpA<=maxA; tmpA+=Astep)
        {
            Nucleus tmpnuke(tmpA);
            tmpnuke.SetGDist(gdist);
            amplitude->SetNucleus(tmpnuke);
            REAL coherent = calculator.TotalCoherentCrossSection(Qsqr, 
                                                                bjorkx);
            REAL quasiel = calculator.TotalCrossSection(Qsqr, bjorkx, 
                                                            mint, maxt);
            cout << tmpA << " " << quasiel/coherent << endl;
        }
    }
    
    else if (mode==MODE_COHERENT_DT)   // d\sigma^A/dt for coherent scattering
    {
        if (A==1) { 
            std::cerr << "A=1, can't be coherent scattering." << std::endl;
            return -1;
        }
        cout << "# t [Gev^2]  d\\sigma/dt [nb/GeV^2], coherent scattering " << endl;
        cout << "# x_pomeron = " << bjorkx << ", Q^2 = " << Qsqr << endl;
        // All iterations are independent, so this is straightforward to parallerize   
        //#pragma omp parallel for
        for (int i=0; i<=points; i++)
        {
            REAL tmpt = (maxt-mint)/points*i;
            REAL delta = sqrt(tmpt);
            REAL result = calculator.CoherentCrossSection_dt(tmpt, Qsqr, bjorkx);
            
            //#pragma omp critical
            {
                cout.precision(5);
                cout << fixed << tmpt;
                cout.precision(8);
                cout << " " << result*NBGEVSQR << endl;
            }
        }
    
    
    }    
  
    else if (mode==MODE_DIFFXS)   // dsigma/dt as a function of t
    {
        cout << "# t [GeV^2]  d\\sigma/dt [nb/GeV^2] " << endl;
        cout << "# x_pomeron = " << bjorkx << ", Q^2 = " << Qsqr << endl;
        // All iterations are independent, so this is straightforward to parallerize   
        //#pragma omp parallel for
        for (int i=0; i<=points; i++)
        {
            REAL tmpt = mint+(maxt-mint)/points*i;
            REAL delta = sqrt(tmpt);
            REAL result=0;
            if (A==1)
            {
                result = calculator.ProtonCrossSection_dt(tmpt, Qsqr, bjorkx); 
            }
            else
            {
                result = calculator.CrossSection_dt(tmpt, Qsqr, bjorkx); 
            }
            //#pragma omp critical
            {
                cout.precision(5);
                cout << fixed << tmpt;
                cout.precision(8);
                cout << " " << result*NBGEVSQR << endl;
            }
        }
    }
    
    else if (mode == MODE_DSIGMA_DT_A) // d\sigma/dt as a function of A [minA-maxA]
    {
        cout << "# A    d\\sigma/dt [nb/GeV^2] " << endl;
        cout << "# x_pomeron = " << bjorkx << ", Q^2 = " << Qsqr << endl;
        
        for (unsigned int tmpA=minA; tmpA<=maxA; tmpA+=Astep)
        {
            Nucleus tmpnuke(tmpA);
            tmpnuke.SetGDist(gdist);
            amplitude->SetNucleus(tmpnuke);
            
            REAL result = calculator.CrossSection_dt(t, Qsqr, bjorkx); 
            cout.precision(3);
            cout << fixed << tmpA;
            cout.precision(8);
            cout << " " << result*NBGEVSQR << endl;
            
        }
        
        
    
    }
    
    else if (mode == MODE_COHERENT_AA or mode==MODE_INCOHERENT_AA)
    {
		if (!pa)
			cout << "# d\\sigma/dy for AA -> AA + J/\\Psi" << endl;
		else 
			cout <<"# d\\sigma/dy for pA -> J/\\Psi pA" << endl;

		cout <<"# sqrts=" << sqrts << " GeV" << endl;
		double miny=-2.2;   // -3.5
		Diffraction d=COHERENT;
		if (mode==MODE_INCOHERENT_AA) d=INCOHERENT;
		if (d==COHERENT) cout <<"# Coherent" << endl;
		else cout <<"# Incoherent" << endl;
	    cout << "# y totxs dsigma/dt" << endl;	
		double maxy=0;
		if (pa)
			maxy=-miny;
		

        for(double y=miny; y<=maxy+0.00001; y+=0.1)
		{
			double res = calculator.DiffractiveAAtoJpsi(y, sqrts, d, pa);
			if (res>0)
				cout << y << " " << res  << " " << calculator.DiffractiveAAtoJpsi_dt(y, sqrts, 0, d) << endl;
            //else
            //    cout << y << " 0 0" << endl;
		}
	
        //cout << calculator.DiffractiveAAtoJpsi_dt(0, sqrts, 0, d) << endl;
		/*cout << "# t    d\\sigma/(dtdy)" << endl;
		for (double t=0; t<0.3; t+=0.001)
		//for (double t=0.1; t<=0.3001; t+=0.002)
		{
			cout << t << " " << calculator.DiffractiveAAtoJpsi_dt(0, sqrts, t, d) << endl;
		}*/

	}
    
    else if (mode==MODE_VM_INTZ)
    {
        cout << "# 2\\pi r * \\int dz/(4\\pi) r Psi^*Psi, Q^2 = " << Qsqr << endl;
        cout << "# " << *((GausLC*)JPsi) << endl;
        if (output_fm) cout << "# [r] = fm"; else cout << "# [r] = GeV^(-1)";
        cout << endl;
        REAL maxr=5; REAL minr=0.0001;
        for (int i=0; i<=points; i++)
        {
            REAL tmpr = minr + (maxr-minr)/points*i;
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
        cout << "# r [GeV^-1]   dsigma/d2b" << endl;
        REAL maxr = 10; REAL minr=0.00001;
        REAL multiplier = pow(maxr/minr, 1.0/points);
        
        for (int i=0; i<=points; i++)
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
        cout << "# b [GeV^-1]   dsigma/d2b" << endl;
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



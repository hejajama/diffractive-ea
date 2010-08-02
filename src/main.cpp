/*
 * Calculates eA cross section in dipole model
 * 
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
#include <iostream>
#include <ctime>
#include <vector>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_exp.h>

#include "dipole.h"
#include "vm_photon.h"
#include "vector.h"
#include "nucleus.h"
#include "dipxs_ipnonsat.h"
#include "mersenne/mersenne.h"

using namespace std;
 
int main(int argc, char* argv[])
{    

    // Parameters
    double xbjork=1e-4;
    REAL x;
    
    // Intialize random number generator
    seed_mersenne(std::time(NULL));
    cout << time(NULL) << endl;
    
    // J/Psi:      e_f, N_T, N_L, R_T, R_L, m_f, M_V, delta
    VM_Photon JPsi(2.0/3.0, 1.23, 0.83, sqrt(6.5), sqrt(3.0), 1.4, 3.097, 0);
    
                    //    Q^2,r,z
    REAL Qsqr=1,r=0.3,z=0.5;
    
    Nucleus nuke(197);
    //cout << nuke << endl;
    //DipXS* dsigmadb = new DipXS_IPNonSat(nuke);
    DipXS_IPNonSat dsigmadb(nuke);
    
    cout << "Generating 100 random nucleon configurations " << endl;
    dsigmadb.GetNucleus().GenerateRandomNucleonConfigurations(100,0,100);
    cout << "Done!"<< endl;
    
     Vec b(0,1); Vec b2(1,1);
     Vec c;
     c=b;
     cout << c << endl;

     //cout << c << endl;
     //cout << dsigmadb->DipXSection_b(0.3, bjorkx, b) << endl;
   // delete dsigmadb;
    
    /*cout << "Random nucleon: " << endl;
    vector<Vec> nucleons = dsigmadb.GetNucleus().RandomNucleonConfiguration();
    //cout << "A = " << nucleons.size() << endl;
    for (int i=0; i<nucleons.size(); i++)
    {
        cout << nucleons[i].GetX() << " " << nucleons[i].GetY() << endl;
    }*/
    
    // Speedtest
    REAL sum=0;
    for (int i=1; i<1e4; i++)
    {
        sum+=dsigmadb.DipXSection_b_sqr(r, r, Vec(i/100.0,0), Vec(0,i/100.0-40),xbjork);
    
    }
    cout << sum << endl;
    
    //cout << dsigmadb.DipXSection_b_sqr(r, r, b, b2, xbjork ) << endl;
    //cout << "delta=0 cross section with r=0.3 : " << dsigmadb.DipXSection_delta(0.3, 0.3, xbjork, Vec(0,0)) << endl;
    
    //for (int i=1; i<1e8; i*=10)
    //    cout << "exp(-" << i << ") = " << exp(-(REAL)i) << endl;
    
   

    return 0;
}


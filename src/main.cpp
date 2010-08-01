/*
 * Calculates eA cross section in dipole model
 * 
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
#include <iostream>
#include <ctime>
#include <vector>
#include <gsl/gsl_math.h>

#include "dipole.h"
#include "vm_photon.h"
#include "vector.h"
#include "nucleus.h"
#include "dipxs_ipnonsat.h"
#include "mersenne/mersenne.h"

using namespace std;
 
int main()
{

    // Parameters
    double xbjork=1e-4;
    REAL x;
    
    // Intialize random number generator
    seed_mersenne(std::time(NULL));
    
    // J/Psi:      e_f, N_T, N_L, R_T, R_L, m_f, M_V, delta
    VM_Photon JPsi(2.0/3.0, 1.23, 0.83, sqrt(6.5), sqrt(3.0), 1.4, 3.097, 0);
    
                    //    Q^2,r,z
    REAL Qsqr=1,r=0.3,z=0.5;
    cout << JPsi << endl;
     
    Nucleus nuke(197);
    //cout << nuke << endl;
    //DipXS* dsigmadb = new DipXS_IPNonSat(nuke);
    DipXS_IPNonSat dsigmadb(nuke);
     
     Vec b(0,1);
     //cout << dsigmadb->DipXSection_b(0.3, bjorkx, b) << endl;
   // delete dsigmadb;
    
    cout << "Random nucleon: " << endl;
    vector<Vec> nucleons = nuke.RandomNucleonConfiguration(0,100);
    for (int i=0; i<nucleons.size(); i++)
    {
        cout << nucleons[i].GetX() << " " << nucleons[i].GetY() 
            << dsigmadb.DipXSectoin_b_nonaveraged(r, xbjork, Vec(0,0), 
            nucleons) << endl;
    }
    

    return 0;
}

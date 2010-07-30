/*
 * Calculates eA cross section in dipole model
 * 
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
#include <iostream>
#include <gsl/gsl_math.h>

#include "dipole.h"
#include "vm_photon.h"
#include "vector.h"
#include "nucleus.h"
#include "dipxs_ipnonsat.h"


using std::cout; using std::endl;
 
int main()
{

    // Parameters
    double bjorkx=1e-4;
    REAL x;
    
    // J/Psi:      e_f, N_T, N_L, R_T, R_L, m_f, M_V, delta
    VM_Photon JPsi(2.0/3.0, 1.23, 0.83, sqrt(6.5), sqrt(3.0), 1.4, 3.097, 0);
    
                    //    Q^2,r,z
    REAL Qsqr=1,r=0.3,z=0.5;
     
    Nucleus nuke(197);
    cout << nuke << endl;
    DipXS* dsigmadb = new DipXS_IPNonSat(nuke);
     
     REAL b = 1;
     cout << dsigmadb->DipXSection_b(0.3, bjorkx, b) << endl;
     
    delete dsigmadb;
    
    cout << log(2.7) << endl;
     

    return 0;
}

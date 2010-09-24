#ifndef Dipxs_IPNONSAT_H
#define Dipxs_IPNONSAT_H

/*
 * Dipole cross section in IP non sat model
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
#include "dipxs.h"
#include "nucleus.h"
#include "dipole.h"
#include <vector>
 
class Dipxs_IPNonSat : public Dipxs
{
    public:
        Dipxs_IPNonSat(Nucleus &nucleus_);
        Dipxs_IPNonSat(Nucleus &nucleus_, REAL bp);
        REAL Dipxsection_sqr_avg(REAL rsqr, REAL r2sqr, REAL xbjork, REAL delta);
        REAL Dipxsection(REAL rsqr, REAL xbjork, Vec b, 
                std::vector<Vec>& nucleons); 
        
       // Amplitude squared for coherent scattering
       REAL CoherentDipxsection_avg(REAL rsqr, REAL xbj, 
                REAL delta);
                
       // Dipole-proton amplitude
       REAL Dipxsection_proton(REAL rsqr, REAL xbj, REAL delta);
       
       // Total dipole-proton cross section (integrated over d^2 b) in 1/Gev^2
       REAL TotalDipxsection_proton(REAL rsqr, REAL xbj);
    
    private:
        REAL Sigmap(REAL rsqr, REAL xbjork);
        REAL prevft, prevdelta; // To optimize Dipxsection_sqr_avg
        REAL B_p; 

};


#endif

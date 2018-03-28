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
#include "ipsat_mz/dipoleamplitude.hpp"
#include <vector>
 
class Dipxs_IPNonSat : public Dipxs
{
    public:
        Dipxs_IPNonSat(Nucleus &nucleus_);
        Dipxs_IPNonSat(Nucleus &nucleus_, REAL bp);
		~Dipxs_IPNonSat();
        REAL DipoleAmplitude_sqr_avg(REAL rsqr, REAL r2sqr, REAL xbjork, REAL delta);
        REAL Dipxsection(REAL rsqr, REAL xbjork, Vec b, 
                std::vector<Vec>& nucleons); 
        
       // Averaged amplitude for coherent scattering
       REAL CoherentDipoleAmplitude_avg(REAL rsqr, REAL xbj, 
                REAL delta);
                
       // Dipole-proton amplitude
       REAL DipoleAmplitude_proton(REAL rsqr, REAL xbj, REAL delta);
       
       // Total dipole-proton cross section (integrated over d^2 b) in 1/Gev^2
       REAL TotalDipxsection_proton(REAL rsqr, REAL xbj);
       
       // 1/2*d\sigma/d^2b = q\barq-proton scattering amplitude
       REAL Qq_proton_amplitude(REAL rsqr, REAL xbj, REAL b);
    
		REAL Bp();
	
		void Initialize();
    private:
        REAL Sigmap(REAL rsqr, REAL xbjork);
        REAL prevft, prevdelta; // To optimize Dipxsection_sqr_avg
        REAL B_p; 
		IPsat_MZ::DipoleAmplitude *ipsat_mz;

};


#endif

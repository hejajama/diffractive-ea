#ifndef DIPXS_H
#define DIPXS_H

/*
 * Parent class for dipole cross section
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
//#include "vector.h"
#include "nucleus.h"
#include "dipole.h"
#include <gsl/gsl_integration.h>

class DipXS
{
    public:
        DipXS();
        DipXS(Nucleus nucleus_);
        ~DipXS();
        
        // Cross section as a function of impact parameter and nucleon config.
        // These methods must be implemented by derived classes
        virtual REAL DipXSection(REAL rsqr, REAL xbj, Vec b, 
                std::vector<Vec>& nucleons) = 0;
        
        
        REAL DipXSection(REAL rsqr, REAL xbj, REAL x, REAL y,
                std::vector<Vec>& nucleons);


        REAL Alphas_r(REAL r);
        Nucleus& GetNucleus();
        
        // Fourier transformation of non-averaged cros sect. DipXSection_b
        COMPLEX FTDipXSection(REAL rsqr, REAL xbjorken,REAL deltax,REAL deltay,
            std::vector<Vec>& nucleons);
        COMPLEX FTDipXSection(REAL rsqr, REAL xbjorken, Vec delta, 
            std::vector<Vec>& nucleons);
            
        // Fourier transformation squared averaged over nucleon configurations as
        // a function of |delta| (averages over angular dependence)
        REAL FTDipXSection_sqr_avg(REAL r1sqr, REAL r2sqr, REAL xbjorken, 
            REAL delta);
        
    protected:
        REAL Mu2(REAL rsqr);   // 4/r^2 + mu_0
        Nucleus nucleus;   
        REAL Sigmap(REAL rsqr, REAL xbjork);  // Total dipole proton cross section
           
        REAL MaxB();        // Integrate b up to this limit, MaxB()=nuke.MaxR()
     
    private:
        // GSL workspace
        void InitGsl();
        void FreeGsl();
        gsl_integration_workspace   * gsl_xworkspace;
        gsl_integration_workspace * gsl_xcycle_workspace;
        gsl_integration_qawo_table * gsl_xqawotable;
        gsl_integration_workspace * gsl_yworkspace;
        gsl_integration_workspace * gsl_ycycle_workspace;
        gsl_integration_qawo_table * gsl_yqawotable;

};

std::ostream& operator<<(std::ostream& os, DipXS& ic);


#endif

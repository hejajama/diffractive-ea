#ifndef Dipxs_H
#define Dipxs_H

/*
 * Parent class for dipole cross section
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
//#include "vector.h"
#include "nucleus.h"
#include "dipole.h"
#include <gsl/gsl_integration.h>

class Dipxs
{
    public:
        Dipxs();
        Dipxs(Nucleus &nucleus_);
        ~Dipxs();
        
        // Some of these methods must be implemented in the derived classes
        
        // Cross section (amplitude) as a function of impact parameter
        // and nucleon config. (d \sigma^2 / d^2 b)_(b_1, ..., b_A)
        virtual REAL Dipxsection(REAL rsqr, REAL xbj, Vec b, 
                std::vector<Vec>& nucleons);
        // Amplitude squared averaged over nucleon configurations as a 
        // function of \Delta and r,r' (and x)
        // \int d^2 b_1 ... d^2 b_A T_A(b_1)...T_A(B_A) 
        //      * \int d^2 b d^2 b' e^(-i(b-b')*\Delta) 
        //      * (d\sigma^A / d^2 b)(b,r,x) (d\sigma^A / d^2b)(b',r',x)
        virtual REAL DipoleAmplitude_sqr_avg(REAL rsqr, REAL r2sqr, REAL xbj, 
                REAL delta)=0;
        // Deprecated:
        REAL Dipxsecton_sqr_avg(REAL r, REAL r2, REAL x, REAL d);
        // Cross section integrated over |t|
        virtual REAL Dipxsection_sqr_avg(REAL rsqr, REAL r2sqr, REAL xbj);
        
        // Coherent scattering amplitude
        // \int d^2 b_1...d^2 b_A T_A(b_1)...T_A(B_A)
        //      * \int d^2 b e^(-ib*\Delta)
        //      * A(b,r,x) 
        virtual REAL CoherentDipoleAmplitude_avg(REAL rsqr, REAL xbj, 
                REAL delta);   
        REAL CoherentDipxsection_avg(REAL rsqr, REAL xbj, 
                REAL delta); // Deprecated
        
        // Scattering amplitude for dipole-proton scattering
        virtual REAL DipoleAmplitude_proton(REAL rsqr, REAL xbj, REAL delta)=0;
        
        // Total dipole-proton cross section in 1/Gev^2
        // \int d^2 b d\sigma/d^2 b
        virtual REAL TotalDipxsection_proton(REAL rsqr, REAL xbj);
                
        // Wrapper, calls virtual Dipxsection with parameter Vec(x,y)
        REAL Dipxsection(REAL rsqr, REAL xbj, REAL x, REAL y,
                std::vector<Vec>& nucleons);
                
        // 1/2*d\sigma/d^2b q\bar q-proton amplitude as a function of r/x/b
        virtual REAL Qq_proton_amplitude(REAL rsqr, REAL xbj, REAL b);

        // Strong coupling constant as a function of dipole size
        REAL Alphas_r(REAL rsqr);
        
        Nucleus& GetNucleus();
        void SetNucleus(Nucleus& nuke);

        virtual REAL Sigma0();
        virtual REAL Bp();

		virtual REAL SaturationScale(double x, int A, double b);   // Saturation scale in proton or nucleus
        
        // Fourier transformation of non-averaged cross sect. Dipxsection_b
        //COMPLEX FTDipxsection(REAL rsqr, REAL xbjorken,REAL deltax,REAL deltay,
        //    std::vector<Vec>& nucleons);
        //COMPLEX FTDipxsection(REAL rsqr, REAL xbjorken, Vec delta, 
        //    std::vector<Vec>& nucleons);
            
        // Fourier transformation squared averaged over nucleon configurations as
        // a function of |delta| (averages over angular dependence)
        //REAL FTDipxsection_sqr_avg(REAL r1sqr, REAL r2sqr, REAL xbjorken, 
        //    REAL delta);
        
    protected:
        Nucleus nucleus;   
           
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

std::ostream& operator<<(std::ostream& os, Dipxs& ic);


#endif

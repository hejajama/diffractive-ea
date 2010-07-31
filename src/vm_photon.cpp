/*
 * Overlap between the photon and the vector meson wave functions
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
#include "vm_photon.h"
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <iostream>

const REAL ZINTACCURACY=0.00001;


VM_Photon::VM_Photon(REAL e_f_, REAL N_T_, REAL N_L_, REAL R_T_, REAL R_L_, 
    REAL m_f_, REAL M_V_, int delta_=0 )
{
    e_f=e_f_;
    m_f=m_f_; M_V=M_V_;
    N_T=N_T_; N_L=N_L_;
    R_T=R_T_; R_L=R_L_;
    delta=delta_; 
     
}


/*
 * Transversially polarized component of the overlap
 */

REAL VM_Photon::PsiSqr_T(REAL Qsqr, REAL r, REAL z)
{
    REAL result;
    REAL epstmp=epsfun(z,Qsqr,SQR(m_f));
    result = SQR(m_f)*gsl_sf_bessel_K0(epstmp*r)*Psi_T(r,z) 
        - (SQR(z)+SQR(1-z))*epstmp*gsl_sf_bessel_K1(epstmp*r)*Psi_T_DR(r,z);
    result *= e_f*e*NC/(M_PI*z*(1-z));
    
    return result;
}

/*
 * Longitudinally polarized component of the overlap
 */

REAL VM_Photon::PsiSqr_L(REAL Qsqr, REAL r, REAL z)
{
    REAL epstmp=epsfun(z,Qsqr,SQR(m_f));
    REAL result;
    result = M_V*Psi_L(r,z)+delta/(M_V*z*(1-z))*(SQR(m_f)*Psi_L(r,z) 
        - 1/r*Psi_L_DR(r,z) - Psi_L_D2R(r,z));
    result *= e_f*e*NC/M_PI*2*sqrt(Qsqr)*z*(1-z)*gsl_sf_bessel_K0(epstmp*r);
    return result;
}

/* 
 * Wave function overlap integrated over z=[0,1]
 * Normalization factor 1/(4*Pi) is included here!
 * PsiSqr_T/L is quite a smooth function so there is 
 * nothing tricky to do here, just use gsl_integration_qng
 */
 
/* As we have to integrate a member function of this class by GSL,
 * we need some helper structures
 */
struct zinthelper{
    VM_Photon *vm_p;
    double  Qsqr;
    double  r;
};

REAL  zhelperfuncT(REAL z, void * p){
  return ((zinthelper*)p)->vm_p->PsiSqr_T(
							((zinthelper*)p)->Qsqr,
							((zinthelper*)p)->r,
							z);
}

REAL  zhelperfuncL(REAL z, void * p){
  return ((zinthelper*)p)->vm_p->PsiSqr_L(
							((zinthelper*)p)->Qsqr,
							((zinthelper*)p)->r,
							z);
}

REAL VM_Photon::PsiSqr_T_intz(REAL Qsqr, REAL r)
{
    REAL result,abserr;
    size_t eval;
    struct zinthelper zintpar;
    zintpar.vm_p=this;
    zintpar.Qsqr=Qsqr;
    zintpar.r=r;
    gsl_function int_helper;
    int_helper.function=&zhelperfuncT;
    int_helper.params=&zintpar;
    
    int status = gsl_integration_qng(&int_helper, 0,1, ZINTACCURACY, ZINTACCURACY, 
        &result, &abserr, &eval);
    
    if(status){ std::cerr<< "z integral in Photon failed with code " 
        << gsl_strerror(status) << std::endl;}

    result*=1/(4.0*M_PI); // Normalization
    return result;
}

REAL VM_Photon::PsiSqr_L_intz(REAL Qsqr, REAL r)
{
    REAL result,abserr;
    size_t eval;
    struct zinthelper zintpar;
    zintpar.vm_p=this;
    zintpar.Qsqr=Qsqr;
    zintpar.r=r;
    gsl_function int_helper;
    int_helper.function=&zhelperfuncL;
    int_helper.params=&zintpar;
    
    int status = gsl_integration_qng(&int_helper, 0,1, ZINTACCURACY, ZINTACCURACY, 
        &result, &abserr, &eval);
    
    if(status){ std::cerr<< "z integral in Photon failed with code " 
        << gsl_strerror(status) << std::endl;}

    result*=1/(4.0*M_PI); // Normalization
    return result;
}

/*
 * Gaussian lightcone functions 
 */
REAL VM_Photon::Psi_T(REAL r, REAL z)
{
    return N_T*SQR(z*(1-z))*gsl_sf_exp(-SQR(r)/(2.0*SQR(R_T)));
}

REAL VM_Photon::Psi_L(REAL r, REAL z)
{
    return N_L*z*(1-z)*gsl_sf_exp(-SQR(r)/(2.0*SQR(R_L)));
}

// \partial_r Psi_T(r,z)
REAL VM_Photon::Psi_T_DR(REAL r, REAL z)
{
    return -1/(SQR(R_T))*N_T*r*SQR(z*(1-z))*gsl_sf_exp(-SQR(r)/(2.0*SQR(R_T)));
}

// \partial_r PSI_L(r,z)
REAL VM_Photon::Psi_L_DR(REAL r, REAL z)
{
    return -1/(SQR(R_T))*N_L*r*z*(1-z)*gsl_sf_exp(-SQR(r)/(2.0*SQR(R_T)));
}

// \partial^2_r PSI_L(r,z)
REAL VM_Photon::Psi_L_D2R(REAL r, REAL z)
{
    return (z-1)*z/gsl_pow_4(R_L) * (SQR(R_L)-SQR(r)) * N_L * gsl_sf_exp(-SQR(r)/(2.0*SQR(R_L)));
}

std::ostream& operator<<(std::ostream& os, const VM_Photon& ic)
{
    return os << " Vector meson and photon wave function overlap. ";
    /*e_f = "
    << e_f << ", m_f = " << m_f << ", M_V = " << M_V << ", N_T = "
    << N_T << ", N_L = " << N_L << ", R_T = " << R_T << ", R_L = " << R_L
    << ", delta = " << delta << " ";
    */
}


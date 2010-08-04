/*
 * Methods to calclulate 2d Fourier transform of dipole nucleus cross section
 * This calculates
 * \int d^2 b  (d\sigma / d^2 b)(r,b_1,...,b_A) * exp(-i*b*\delta)
 * 
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 *
 * Most of this code is written by Tuomas Lappi under non-specified license
 * Changes made to that code:
 * - Removed parts not needed here
 * - Ported to use DipXS class
 * - Takes into account and averages over nucleon configurations
 */
 
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include <iostream>
#include <complex>
#include "dipxs.h"

using std::cout; using std::cerr; using std::endl;
 
 

// helper functions for 2d Fourier transforms. To use GSL we have to 
// separate space into four quadrants
//                 y   
//           2     |    1
//                 |   
//         -----------------x
//                 |
//           3     |    4
// 
// Need four differently symmetrized helper functions,
// and versions of both for integration over x (fixed y) and y (x integrated)
// \int_{-\infty}^{\infty} dx dy e(i qx x + i qy y)
// = \int_{-\infty}^{\infty} dx dy
//       cos(qx x)cos(qy y) ( f(x,y) + f(x,-y) + f(-x,y) + f(-x,-y) ) <-- cc
//       -  sin(qx x)sin(qy y) ( f(x,y) - f(x,-y) - f(-x,y) + f(-x,-y) )  <-- ss
//       + i cos(qx x)sin(qy y) ( f(x,y) - f(x,-y) + f(-x,y) - f(-x,-y) ) <-- cs
//       + i sin(qx x)cos(qy y) ( f(x,y) + f(x,-y) - f(-x,y) - f(-x,-y) )  <-- sc




#define WORKSPACE 10000
#define QAWOSIZE 10000

#define ACCURACY 0.01 // For qx,qy integrals
#define BESACCURACY 0.0001 // For J_0(qb)-integrals
#define QSACC 0.00001

#define MAXBRATIO 0.000001

#define AVGDIPXSACCLIMIT 0.01 // swich tactics for accuracy better than this

#define MINQ 0.01 // q smaller than this treat as zero for Fourier integrals
#define MAXR (1.0/MINQ)     // TODO: This shouldn't be hardcoded?

#define INT_NPOINTS 10000
// #define gslL 100.0 // interval size; required by gsl, but not used

#define MAXR_IN_FT  // If this is defined, integrate only up to maximum bt in 2d Fourier integrals
// This option does not affect the 1dim (cylindrical Bessel) integals.

// Number of integrals to calculate to average over different nucleon configurations
const int NUCLEON_CONFIG_INTEGRALS = 1;

// Number of integrals to calculate to average over angular dependence
const int ANGULAR_INTEGRALS = 1;



// helper struct to hold pointers to Dipxs and y value
struct xintparastruct{
  DipXS * dxs;
  double rsqr;
  double xbj;
  double y;
  std::vector<Vec>* nucleons;
};

struct yintparastruct{
  DipXS * dxs;
  double rsqr;
  double xbj;
  double qx;
  double xmax;
  gsl_integration_workspace * xworkspace;
  gsl_integration_workspace * xcycle_workspace;
  gsl_integration_qawo_table * xwf;
  std::vector<Vec>* nucleons;
};

double ccxintegrand(double x, void *p){
 double yval =  ((struct xintparastruct *)p)->y;
 double rsqrval =  ((struct xintparastruct *)p)->rsqr;
 double xbjval =  ((struct xintparastruct *)p)->xbj;
 DipXS *dipole = ((struct xintparastruct *)p)->dxs;
 std::vector<Vec> *nucleons = ((struct xintparastruct *)p)->nucleons;

 return dipole->DipXSection(rsqrval,xbjval,yval,x, *nucleons) + 
   dipole->DipXSection(rsqrval,xbjval,-yval,x, *nucleons) + 
   dipole->DipXSection(rsqrval,xbjval,yval,-x, *nucleons) + 
   dipole->DipXSection(rsqrval,xbjval,-yval,-x, *nucleons) 
   ;

}
double ssxintegrand(double x, void *p){
 double yval =  ((struct xintparastruct *)p)->y;
 double rsqrval =  ((struct xintparastruct *)p)->rsqr;
 double xbjval =  ((struct xintparastruct *)p)->xbj;
 DipXS *dipole = ((struct xintparastruct *)p)->dxs;
 std::vector<Vec> *nucleons = ((struct xintparastruct *)p)->nucleons;

 return dipole->DipXSection(rsqrval,xbjval,yval,x, *nucleons)
   - dipole->DipXSection(rsqrval,xbjval,-yval,x, *nucleons)  
   - dipole->DipXSection(rsqrval,xbjval,yval,-x, *nucleons)  
   + dipole->DipXSection(rsqrval,xbjval,-yval,-x, *nucleons) 
   ;
}

double csxintegrand(double x, void *p){
 double yval =  ((struct xintparastruct *)p)->y;
 double rsqrval =  ((struct xintparastruct *)p)->rsqr;
 double xbjval =  ((struct xintparastruct *)p)->xbj;
 DipXS *dipole = ((struct xintparastruct *)p)->dxs;
 std::vector<Vec> *nucleons = ((struct xintparastruct *)p)->nucleons;
 
 return dipole->DipXSection(rsqrval,xbjval,yval,x,*nucleons)
   - dipole->DipXSection(rsqrval,xbjval,-yval,x, *nucleons)  
   + dipole->DipXSection(rsqrval,xbjval,yval,-x, *nucleons)  
   - dipole->DipXSection(rsqrval,xbjval,-yval,-x, *nucleons) 
   ;
}
double scxintegrand(double x, void *p){
 double yval =  ((struct xintparastruct *)p)->y;
 double rsqrval =  ((struct xintparastruct *)p)->rsqr;
 double xbjval =  ((struct xintparastruct *)p)->xbj;
 DipXS *dipole = ((struct xintparastruct *)p)->dxs;
 std::vector<Vec> *nucleons = ((struct xintparastruct *)p)->nucleons;

 return dipole->DipXSection(rsqrval,xbjval,yval,x, *nucleons)
   + dipole->DipXSection(rsqrval,xbjval,-yval,x, *nucleons)  
   - dipole->DipXSection(rsqrval,xbjval,yval,-x, *nucleons)  
   - dipole->DipXSection(rsqrval,xbjval,-yval,-x, *nucleons) 
   ;
}

double ccyintegrand(double y, void *p){
  
  struct xintparastruct xintpara;
  xintpara.dxs = ((struct yintparastruct *)p)->dxs;
  xintpara.rsqr = ((struct yintparastruct *)p)->rsqr;
  xintpara.xbj = ((struct yintparastruct *)p)->xbj;
  xintpara.y = y;
  xintpara.nucleons = ((struct yintparastruct *)p)->nucleons;
  
  double abserr,result;
  gsl_function int_helper;
  gsl_set_error_handler_off();
  
  int status;
  status = gsl_integration_qawo_table_set(((struct yintparastruct *)p)->xwf
					  ,((struct yintparastruct *)p) ->qx, 
					  ((struct yintparastruct *)p)->xmax,  GSL_INTEG_COSINE);
  int_helper.function = &ccxintegrand;
  int_helper.params = (void *)(&xintpara);
  
  if( ABS(((struct yintparastruct *)p) ->qx  ) > MINQ){
#ifdef MAXR_IN_FT
    status = gsl_integration_qawo(&int_helper, 0.0 ,ACCURACY,ACCURACY, WORKSPACE,
				  ((struct yintparastruct *)p)-> xworkspace, 
				  ((struct yintparastruct *)p)-> xwf,
				  &result, &abserr);
#else
    status = gsl_integration_qawf(&int_helper, 0.0 ,ACCURACY, WORKSPACE,
				  ((struct yintparastruct *)p)-> xworkspace, 
				  ((struct yintparastruct *)p)-> xcycle_workspace,
				  ((struct yintparastruct *)p)-> xwf,
				  &result, &abserr);
#endif
    
  } else { // qx = 0, do not treat as oscillatory
    double maxr = MAXR;
    status = gsl_integration_qag(&int_helper, 0.0 ,maxr,ACCURACY,ACCURACY, 
				 WORKSPACE,GSL_INTEG_GAUSS61,
				 ((struct yintparastruct *)p)-> xworkspace, 
				 &result, &abserr);
    //	cerr << "y " << y << " qx =0, ccint(x) " << result << endl;
  }
  
  if(status){cerr<< "cc integral over x failed with code " 
		 << gsl_strerror(status) << endl;}
  
  //cerr << "y " << y << " qx " << ((struct yintparastruct *)p) ->qx << " int dx " << result << endl;
  
  return result;
}

double ssyintegrand(double y, void *p){
  
  struct xintparastruct xintpara;
  xintpara.dxs = ((struct yintparastruct *)p)->dxs;
  xintpara.rsqr = ((struct yintparastruct *)p)->rsqr;
  xintpara.xbj = ((struct yintparastruct *)p)->xbj;
  xintpara.y = y;
  xintpara.nucleons = ((struct yintparastruct *)p)->nucleons;
  
  double abserr,result;
  gsl_function int_helper;
  
  
  gsl_set_error_handler_off();
  
  int status;
  status = gsl_integration_qawo_table_set(((struct yintparastruct *)p)->xwf
					  ,((struct yintparastruct *)p) ->qx, 
					  ((struct yintparastruct *)p)->xmax,  GSL_INTEG_SINE);
  int_helper.function = &ssxintegrand;
  int_helper.params = (void *)(&xintpara);
  
  if( ABS(((struct yintparastruct *)p) ->qx  ) > MINQ){
#ifdef MAXR_IN_FT
    status = gsl_integration_qawo(&int_helper, 0.0 ,ACCURACY,ACCURACY, WORKSPACE,
				  ((struct yintparastruct *)p)-> xworkspace, 
				  ((struct yintparastruct *)p)-> xwf,
				  &result, &abserr);
#else
    status = gsl_integration_qawf(&int_helper, 0.0 ,ACCURACY, WORKSPACE,
				  ((struct yintparastruct *)p)-> xworkspace, 
				  ((struct yintparastruct *)p)-> xcycle_workspace,
				  ((struct yintparastruct *)p)-> xwf,
				  &result, &abserr);
#endif
    
  } else { // qx = 0, do not treat as oscillatory
    status = 0;
    result = 0;
    abserr =  ABS(((struct yintparastruct *)p) ->qx  );
  }
  if(status){cerr<< "ss integral over x failed with code " 
		 << gsl_strerror(status) << endl;}
  
  return result;
}

double csyintegrand(double y, void *p){
  
  struct xintparastruct xintpara;
  xintpara.dxs = ((struct yintparastruct *)p)->dxs;
  xintpara.rsqr = ((struct yintparastruct *)p)->rsqr;
  xintpara.xbj = ((struct yintparastruct *)p)->xbj;
  xintpara.y = y;
  
  double abserr,result;
  gsl_function int_helper;
  
  
  gsl_set_error_handler_off();
  
  int status;
  status = gsl_integration_qawo_table_set(((struct yintparastruct *)p)->xwf
					  ,((struct yintparastruct *)p) ->qx, 
					  ((struct yintparastruct *)p)->xmax,  GSL_INTEG_SINE);
  int_helper.function = &csxintegrand;
  int_helper.params = (void *)(&xintpara);
  
  if( ABS(((struct yintparastruct *)p) ->qx  ) > MINQ){
#ifdef MAXR_IN_FT
    status = gsl_integration_qawo(&int_helper, 0.0 ,ACCURACY,ACCURACY, WORKSPACE,
				  ((struct yintparastruct *)p)-> xworkspace, 
				  ((struct yintparastruct *)p)-> xwf,
				  &result, &abserr);
#else
   status = gsl_integration_qawf(&int_helper, 0.0 ,ACCURACY, WORKSPACE,
				  ((struct yintparastruct *)p)-> xworkspace, 
				  ((struct yintparastruct *)p)-> xcycle_workspace,
				  ((struct yintparastruct *)p)-> xwf,
				  &result, &abserr);
#endif
  } else { // qx = 0, do not treat as oscillatory
    double maxr = MAXR;
    status = gsl_integration_qag(&int_helper, 0.0 ,maxr,ACCURACY,ACCURACY, 
				 WORKSPACE,GSL_INTEG_GAUSS61,
				 ((struct yintparastruct *)p)-> xworkspace, 
				 &result, &abserr);
  }
  
  if(status){cerr<< "cs integral over x failed with code " 
		 << gsl_strerror(status) << endl;}
  
  return result;
}
double scyintegrand(double y, void *p){
  
  struct xintparastruct xintpara;
  xintpara.dxs = ((struct yintparastruct *)p)->dxs;
  xintpara.rsqr = ((struct yintparastruct *)p)->rsqr;
  xintpara.xbj = ((struct yintparastruct *)p)->xbj;
  xintpara.y = y;
  xintpara.nucleons = ((struct yintparastruct *)p)->nucleons;
  
  double abserr,result;
  gsl_function int_helper;
  
  
  gsl_set_error_handler_off();
  
  int status;
  status = gsl_integration_qawo_table_set(((struct yintparastruct *)p)->xwf
					  ,((struct yintparastruct *)p) ->qx, 
					 ((struct yintparastruct *)p)->xmax,  GSL_INTEG_COSINE);
  int_helper.function = &scxintegrand;
  int_helper.params = (void *)(&xintpara);
  
  if( ABS(((struct yintparastruct *)p) ->qx  ) > MINQ){
#ifdef MAXR_IN_FT
    status = gsl_integration_qawo(&int_helper, 0.0 ,ACCURACY,ACCURACY, WORKSPACE,
				  ((struct yintparastruct *)p)-> xworkspace, 
				  ((struct yintparastruct *)p)-> xwf,
				  &result, &abserr);
#else
    status = gsl_integration_qawf(&int_helper, 0.0 ,ACCURACY, WORKSPACE,
				  ((struct yintparastruct *)p)-> xworkspace, 
				  ((struct yintparastruct *)p)-> xcycle_workspace,
				  ((struct yintparastruct *)p)-> xwf,
				  &result, &abserr);
#endif
  } else { // qx = 0, do not treat as oscillatory
    status = 0;
    result = 0;
    abserr =  ABS(((struct yintparastruct *)p) ->qx  );
  }
  if(status){cerr<< "sc integral over x failed with code " 
		 << gsl_strerror(status) << endl;}
  
  return result;
}

// ********************** Actual FT method ***************

COMPLEX DipXS::FTDipXSection(REAL rsqr,REAL xbjorken, Vec delta, 
    std::vector<Vec>& nucleons)
{
    return FTDipXSection(rsqr, xbjorken, delta.GetX(), delta.GetY(), 
            nucleons);
}

COMPLEX DipXS::FTDipXSection(REAL rsquare,REAL xbjorken, REAL qx, REAL qy, 
    std::vector<Vec>& nucleons) {
  double ccint(0),ssint(0),scint(0),csint(0);
  gsl_function int_helper;
  struct yintparastruct yintpara;
  yintpara.rsqr= rsquare;
  yintpara.xbj= xbjorken;
  yintpara.dxs= this;
  yintpara.qx= qx;
  yintpara.xworkspace = gsl_xworkspace;
  yintpara.xcycle_workspace = gsl_xcycle_workspace;
  yintpara.xwf=gsl_xqawotable;
  yintpara.xmax=MaxB();
  yintpara.nucleons=&nucleons;
  double abserr;
  
  gsl_set_error_handler_off();
  int status;
  // First ccintegral
  // Set weight function to cos(qy*x)
  gsl_integration_qawo_table_set(gsl_yqawotable, qy, MaxB(),  GSL_INTEG_COSINE);
  int_helper.function = &ccyintegrand;
  int_helper.params = (void *)(&yintpara);
  
  if(ABS(qy) > MINQ){
#ifdef MAXR_IN_FT
	status = gsl_integration_qawo(&int_helper, 0.0 ,ACCURACY,ACCURACY, WORKSPACE, 
				      gsl_yworkspace, 
				      gsl_yqawotable, 
				      &ccint, &abserr);
#else
	status = gsl_integration_qawf(&int_helper, 0.0 ,ACCURACY, WORKSPACE, 
				      gsl_yworkspace, 
				      gsl_ycycle_workspace, gsl_yqawotable, 
				      &ccint, &abserr);
#endif
  }else{
    double maxr = MAXR;
    status = gsl_integration_qag(&int_helper, 0.0 ,maxr,ACCURACY,ACCURACY, 
				 WORKSPACE,GSL_INTEG_GAUSS61,
				 gsl_yworkspace, &ccint, &abserr);	
    //	cerr << " qy=0, ccint " << ccint<<endl;
  }
  if(status){cerr<< "cc integral over y failed with code " 
		 << gsl_strerror(status) << endl;}
  
  // Then ssintegral
  gsl_integration_qawo_table_set(gsl_yqawotable, qy, MaxB(),  GSL_INTEG_SINE);
  int_helper.function = &ssyintegrand;
  int_helper.params = (void *)(&yintpara);
  if(ABS(qy) > MINQ){
#ifdef MAXR_IN_FT
    status = gsl_integration_qawo(&int_helper, 0.0 ,ACCURACY,ACCURACY, WORKSPACE, gsl_yworkspace, 
				  gsl_yqawotable, &ssint, &abserr);
#else
    status = gsl_integration_qawf(&int_helper, 0.0 ,ACCURACY, WORKSPACE, gsl_yworkspace, 
				  gsl_ycycle_workspace, gsl_yqawotable, &ssint, &abserr);
#endif
  }else{
    status=0;
    ssint=0;
    abserr = ABS(qy);
  }
  
  if(status){cerr<< "cc integral over y failed with code " 
		 << gsl_strerror(status) << endl;}
  // Then csintegral
  gsl_integration_qawo_table_set(gsl_yqawotable, qy, MaxB(),  GSL_INTEG_SINE);
  int_helper.function = &csyintegrand;
  int_helper.params = (void *)(&yintpara);
  if(ABS(qy) > MINQ){
 #ifdef MAXR_IN_FT
  status = gsl_integration_qawo(&int_helper, 0.0 ,ACCURACY,ACCURACY, WORKSPACE, gsl_yworkspace, 
				  gsl_yqawotable, &csint, &abserr);
 #else
  status = gsl_integration_qawf(&int_helper, 0.0 ,ACCURACY, WORKSPACE, gsl_yworkspace, 
				  gsl_ycycle_workspace, gsl_yqawotable, &csint, &abserr);
 #endif
 }else{
    status=0;
    csint=0;
    abserr = ABS(qy);
  }
  if(status){cerr<< "cs integral over y failed with code " 
		 << gsl_strerror(status) << endl;}
  
  // Then scintegral
  gsl_integration_qawo_table_set(gsl_yqawotable, qy, MaxB(),  GSL_INTEG_COSINE);
  int_helper.function = &scyintegrand;
  int_helper.params = (void *)(&yintpara);
  if(ABS(qy) > MINQ){
#ifdef MAXR_IN_FT
   status = gsl_integration_qawo(&int_helper, 0.0 ,ACCURACY,ACCURACY, WORKSPACE, gsl_yworkspace, 
				   gsl_yqawotable, &scint, &abserr);
#else
   status = gsl_integration_qawf(&int_helper, 0.0 ,ACCURACY, WORKSPACE, gsl_yworkspace, 
				  gsl_ycycle_workspace, gsl_yqawotable, &scint, &abserr);
#endif
  }else{
    double maxr = MAXR;
    status = gsl_integration_qag(&int_helper, 0.0 ,maxr,ACCURACY,ACCURACY,
				 WORKSPACE,GSL_INTEG_GAUSS61,
				 gsl_yworkspace, &scint, &abserr);	
  }
  
  if(status){cerr<< "sc integral over y failed with code " 
		 << gsl_strerror(status) << endl;}
  
  return COMPLEX(ccint-ssint,scint+csint);

}


/* 
 * Averaged Fourier transformation squared as a function of r1 and r2
 *
 * Averages the result of FTDipXSection over different
 * random nucleon configuratons and the angular dependence on \delta
 *
 * Monte Carlo method is used here as it 
 * should converge quite rapidly, all nuclei should be quite similar
 * when nucleon positions are generated from the Woods-Saxon distribution
 * Monte Carlo means here that we generate some random nucleai and average
 * the product of  
 * 
 * Calculates:
 * \int d^2 b_1 ... d^2 b_A T_A(b_1)...T_A(b_A) 
 *  * (d\sigma/d^2 b)(r,b_1,...,b_A) * F[ (d\sigma/d^2 b')(r',b_1,...,b_A)],
 * where F[] is the Fourier transformation averaged over the angular dependence
 * on \delta
 *
 */

REAL DipXS::FTDipXSection_sqr_avg(REAL r1sqr, REAL r2sqr, REAL xbjork, REAL delta)
{
    int npoints=ANGULAR_INTEGRALS;  // Number of points to use to average over angular depend.
    if (delta < MINQ) npoints=1;
    REAL result=0;
    REAL deltaphi = 2.0*M_PI/npoints;
    
    cout << "Averaging over " << NUCLEON_CONFIG_INTEGRALS << " nucleon "
        << "configurations and " << npoints << " angles " << endl;
    
    // Note: trivial to parallerize
    #pragma omp parallel for
    for (int i=0; i<NUCLEON_CONFIG_INTEGRALS; i++)
    {
        std::vector<Vec> &config = nucleus.RandomNucleonConfiguration();
        REAL angular_result=0;
        for (int j=0; j<npoints; j++)    // Angular part
        {
            REAL qx = delta*cos(j*deltaphi);
            REAL qy = delta*sin(j*deltaphi);
            COMPLEX xs(0,0);
            if (r1sqr==r2sqr)
                xs = std::norm(FTDipXSection(r1sqr,xbjork, qx, qy, config));
            else
                xs = FTDipXSection(r1sqr,xbjork, qx, qy, config)
                    *std::conj(FTDipXSection(r2sqr,xbjork, qx, qy, config));
            
            angular_result += std::abs(xs);
            
            REAL foo = 1.0/(16.0*M_PI);
            cout << "Average integral " << i << " angular integral (" <<qx << " , " << qy << ") = " 
                << xs << " => res " << std::abs(xs)*foo << endl; 
        }
        result+=angular_result;
        cout << "Average result after " << i+1 << " average integrals: " << result/((i+1)*npoints*16*M_PI) << endl;
    }
    //cout << "Sum: " << result << endl;
    result*=1.0/(npoints*NUCLEON_CONFIG_INTEGRALS);
    return result;
}


/*
REAL DipXS::avgftdipxsectionsqr(REAL rsqr,REAL xbj, REAL qt) const {
  int npoints=7;
  if(qt < MINQ){npoints = 1;}
  REAL deltaphi = 2*M_PI/npoints;
  REAL sum=0;
  COMPLEX xs;
  for(int i=0;i< npoints;i++){
    xs = ftDipXSection(rsqr,xbj,qt*cos(i*deltaphi),qt*sin(i*deltaphi));	
//	cerr << " q = " <<qt*cos(i*deltaphi) << "," << qt*sin(i*deltaphi) << " Xs " << xs << endl;
	sum += norm(xs);gsl_yqawotable
  }

//	cerr << "DipXS::avgftDipXSection() not implemented yet"<< endl;
  return sum/npoints;
}
*/
/*
struct pftstruct{
	const DipXS * dip;
	REAL rsqr;
	REAL xbj;
	REAL q;
};

double pfthelper(REAL bt, void *p){
  return bt*((struct pftstruct *)p)->dip->pDipXSection(
					((struct pftstruct *)p)->rsqr,
					((struct pftstruct *)p)->xbj,
					bt  
				) *
				gsl_sf_bessel_J0(
									(((struct pftstruct *)p)->q )
									*
									bt
				);
}*/
/*
REAL DipXS::ftpDipXSection(REAL rsquare,REAL xbjorken, REAL qt) const {
  struct pftstruct pfstr;
  pfstr.dip = this;
  pfstr.xbj = xbjorken;
  pfstr.rsqr = rsquare;
  pfstr.q = qt;
  gsl_function int_helper;
  int_helper.function = &pfthelper;
  int_helper.params = (void *)(&pfstr);

  int status;
  double result, abserr;
// Works very bad in the glauber nucleus case, so use finite maxbt also here, although
// for the proton this is not bad
//  status = gsl_integration_qagiu(&int_helper, 0.0 ,BESACCURACY,BESACCURACY, WORKSPACE, 
//				      gsl_yworkspace,  &result, &abserr);
  status = gsl_integration_qag(&int_helper, 0.0 ,maxbt(),BESACCURACY,BESACCURACY, 
				   WORKSPACE,GSL_INTEG_GAUSS61,
				   gsl_yworkspace, 
				   &result, &abserr);

  return 2*M_PI*result; // 2 pi from angular integral of bt
}

double glauberfthelper(REAL bt, void *p){
  return bt*((struct pftstruct *)p)->dip->avgdipxsectionglauber(
					((struct pftstruct *)p)->rsqr,
					((struct pftstruct *)p)->xbj,
					bt  
				) *
				gsl_sf_bessel_J0(
									(((struct pftstruct *)p)->q )
									*
									bt
				);
}

REAL DipXS::ftglauberDipXSection(REAL rsquare,REAL xbjorken, REAL qt) const {
  struct pftstruct pfstr;
  pfstr.dip = this;
  pfstr.xbj = xbjorken;
  pfstr.rsqr = rsquare;
  pfstr.q = qt;
  gsl_function int_helper;
  int_helper.function = &glauberfthelper;
  int_helper.params = (void *)(&pfstr);

  int status;
  double result, abserr;
// qagiu works very bad, integrate to maxbt only
//  status = gsl_integration_qagiu(&int_helper, 0.0 ,BESACCURACY,BESACCURACY, WORKSPACE, 
//				      gsl_yworkspace,  &result, &abserr);
  status = gsl_integration_qag(&int_helper, 0.0 ,maxbt(),BESACCURACY,BESACCURACY, 
				   WORKSPACE,GSL_INTEG_GAUSS61,
				   gsl_yworkspace, 
				   &result, &abserr);

  return 2*M_PI*result; // 2 pi from angular integral of bt
}
*/

 

void DipXS::InitGsl(){
/*
  if(!gsl_xworkspace){  gsl_xworkspace = gsl_integration_workspace_alloc(WORKSPACE);}
  if(!gsl_xcycle_workspace){  gsl_xcycle_workspace 
      = gsl_integration_workspace_alloc(WORKSPACE);}
  if(!gsl_xqawotable){ gsl_xqawotable = 
      gsl_integration_qawo_table_alloc(1.0, 1.0, GSL_INTEG_SINE,
				       QAWOSIZE);}
  if(!gsl_yworkspace){  gsl_yworkspace = gsl_integration_workspace_alloc(WORKSPACE);}
  if(!gsl_ycycle_workspace){  gsl_ycycle_workspace = 
      gsl_integration_workspace_alloc(WORKSPACE);}
  if(!gsl_yqawotable){  gsl_yqawotable = 
      gsl_integration_qawo_table_alloc(1.0, 1.0, GSL_INTEG_SINE,   QAWOSIZE);}
  */
  gsl_xworkspace = gsl_integration_workspace_alloc(WORKSPACE);
  gsl_xcycle_workspace = gsl_integration_workspace_alloc(WORKSPACE);
  gsl_xqawotable = gsl_integration_qawo_table_alloc(1.0, 1.0, GSL_INTEG_SINE,
				       QAWOSIZE);
  gsl_yworkspace = gsl_integration_workspace_alloc(WORKSPACE);
  gsl_ycycle_workspace =  gsl_integration_workspace_alloc(WORKSPACE);
  gsl_yqawotable = 
      gsl_integration_qawo_table_alloc(1.0, 1.0, GSL_INTEG_SINE,   QAWOSIZE);
  
}
void DipXS::FreeGsl(){
    gsl_integration_workspace_free(gsl_xworkspace);
    gsl_integration_workspace_free(gsl_xcycle_workspace);
    gsl_integration_qawo_table_free (gsl_xqawotable);
    gsl_integration_workspace_free(gsl_yworkspace);
    gsl_integration_workspace_free(gsl_ycycle_workspace);
    gsl_integration_qawo_table_free (gsl_yqawotable); 
}

#ifndef GDist_h
#define GDist_h
#include "../dipole.h"
#include <fstream>
#include <iostream>

#define SPLINE

#ifdef SPLINE
#include "gsl/gsl_interp.h"
#include "gsl/gsl_spline.h"
#endif

class GDist {
 public: 
  GDist(){}; 
  ~GDist(){};

  virtual REAL gluedist(REAL x,REAL rsqr) const = 0;

 private:

};


class DGLAPDist : public GDist {
 public: 
  DGLAPDist(); 
  ~DGLAPDist();

  REAL gluedist(REAL x,REAL rsqr) const;

 private:
  REAL minxbj; // these are enforced
  REAL maxxbj;
  REAL minrsqr; // this just determines whether we interpolate or extrapolate in rsqr
  REAL maxrsqr;
  REAL log2maxxbj,deltalog2xbj,log2minrsqr,deltalog2rsqr;
  int Nxbj;
  int Nrsqr;
  REAL C;
  REAL mu0sqr;
  REAL * xbjvals;
  REAL * rsqrvals;
  REAL * gdistdata;
  REAL& gdistat(int xbjind,int rsqrind){return gdistdata[Nrsqr*xbjind + rsqrind];}
  REAL& gdistat(int xbjind,int rsqrind) const {return gdistdata[Nrsqr*xbjind + rsqrind];}

#ifdef SPLINE
  gsl_interp ** interpdata;
  REAL * rsqrlogs;
  gsl_interp_accel * intaccel;
#endif

  // make these mutable to be able to cache values
  // rsqr value being asked for is between rsqrvals[rind-1] and rsqrvals[rind],
  // if outside range then rind = 0 or rind == Nrsqr
  // xbj value being asked for is between xbjvals[xind-1] and xbjvals[xind],
  mutable int rind;
  mutable int xind;

};


std::ostream& operator<<(std::ostream& os,const GDist& ic);

#endif

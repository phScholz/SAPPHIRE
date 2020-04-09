#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <vector>
#include <complex>
#include "ParticleTransmissionFunc.h"

class CoulFunc; 

class Potential : public ParticleTransmissionFunc {
 public:
  Potential(int,int,int,int,double,int,double,int,
	    double,int,double,double,double,double, 
	    TransmissionFunc*);
  ~Potential();
  int GetZ1Z2() const {return z1_*z2_;};
  double GetBoundaryRadius() const {return boundaryRadius_;};
  double GetRedMass() const {return redmass_;};
  double GetRMax() const {return boundaryRadius_+1.;};
  virtual std::complex<double> Calculate(double,int,double,double,double) const = 0;
  std::complex<double> CalcBeta(double,int,double,double,double) const;
  double CalcTransmission(double,int,double);
  void NormalizeInternally(std::vector<std::complex<double> >&,double) const;
  void NormalizeOverAllSpace(std::vector<std::complex<double> >&,double) const;
  std::vector<std::complex<double> > Solve(double,int,double,double,double) const;
 protected:
  CoulFunc* coulFunc_;
  double boundaryRadius_;
  double coulombRadius_;
};

#endif

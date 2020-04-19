#ifndef JLMPOTENTIAL_H
#define JLMPOTENTIAL_H
#pragma once
#include "Potential.h"

class TransmissionFunc;

class JLMPotential : public Potential {
 public:
  JLMPotential(int,int,int,int,double,int,
	       double,int,double,int,double,
	       double,double,double,
	       TransmissionFunc*); 
  double CalculateDensity(double,int) const ;
  std::complex<double> Calculate(double,int,double,double,double) const;
  static double GetA(int i, int j) {return a[i][j];};
 private:
  static constexpr double aRho_ = 0.54;
  static const double a[3][3];
  static const double b[3][3];
  static const double c[3][3];
  static const double dHighE[4][4];
  static const double fHighE[4][4];
  static const double dLowE[4][4];
  static const double fLowE[4][4];
  double cRho_;
  double rho0_[2];
};

#endif

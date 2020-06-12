#include "JLMPotential.h"
#include <assert.h>

JLMPotential::JLMPotential(int z1, int m1, int z2, int m2, 
			   double jInitial, int piInitial,
			   double jFinal, int piFinal,
			   double spin, int parity, double maxL,
			   double totalWidthForCorrection,
			   double uncorrTotalWidthForCorrection,
			   double uncorrTotalWidthSqrdForCorrection,
			   TransmissionFunc* previous) : 
  Potential(z1,m1,z2,m2,jInitial,piInitial,
	    jFinal,piFinal,spin,parity,maxL,
	    totalWidthForCorrection,
	    uncorrTotalWidthForCorrection,
	    uncorrTotalWidthSqrdForCorrection,
	    previous) {
  cRho_ = (0.978+0.0206*pow(double(m2_),1./3.))*pow(double(m2_),1./3.);
  double rhoFactor = 3./4./pi/cRho_/cRho_/cRho_/
    (1.+pi*pi*aRho_*aRho_/cRho_/cRho_);
  rho0_[0] = rhoFactor*(m2_-z2_);
  rho0_[1] = rhoFactor*z2_;
  boundaryRadius_=cRho_+7.*aRho_;
  coulombRadius_=1.3*pow(double(m2_),1./3.);
}

double JLMPotential::CalculateDensity(double r, int which) const {
  assert(which<2);
  return rho0_[which]/(1.+exp((r-cRho_)/aRho_));
}

std::complex<double> JLMPotential::Calculate(double r, int l, double s, double j, double E) const {
  double rhoN = CalculateDensity(r,0);
  double rhoP = CalculateDensity(r,1);
  double rho = rhoN+rhoP;
  double asymParam = (rhoN-rhoP)/(rhoN+rhoP);
  if(GetZ1Z2()!=0) asymParam*=-1.;

  double coulombPotential=0.0;
  if(r<=coulombRadius_) {
    coulombPotential=hbarc*fstruc*GetZ1Z2()/coulombRadius_*
      (1.0+0.5*(1-pow(r/coulombRadius_,2.0)));
  } else {
    coulombPotential=hbarc*fstruc*GetZ1Z2()/r;
  }
  E-=coulombPotential;

  double realScaler = 0.;
  double derivRealScaler = 0.;
  double realVector = 0.;
  double mRatio = 1.;
  double imaginaryScaler = 0.;
  double imaginaryVector = 0.;
  for(int i = 0; i<4;i++) {
    for(int j = 0;j<4;j++) {
      double factor =  pow(rho,double(i+1))*pow(E,double(j));
      if(i<3&&j<3) {
	realScaler+=a[j][i]*factor;
	derivRealScaler+=a[j][i]*factor*double(j)/E;
	realVector+=b[j][i]*factor;
	mRatio-=c[j][i]*factor;
      }
      if(E>15.) {
	imaginaryScaler+=dHighE[j][i]*factor;
	imaginaryVector+=fHighE[j][i]*factor;
      } else {	
	imaginaryScaler+=dLowE[j][i]*factor;
	imaginaryVector+=fLowE[j][i]*factor;
      }
    }
  }
  realVector*=mRatio;
  double realTotal = realScaler+asymParam*realVector;
  double fermiEnergy = (E>15.) ? rho*(-510.8+3222.*rho-6250.*rho*rho) :
    -22.-rho*(298.52-3760.23*rho+12435.82*rho*rho);
  double dValue = (E>15.) ? 600. : 100.;
  imaginaryScaler/=(1.+dValue/(E-fermiEnergy)/(E-fermiEnergy));
  imaginaryVector*=mRatio/(1.-derivRealScaler)/(1.+1./(E-fermiEnergy));
  double imaginaryTotal =  imaginaryScaler+asymParam*imaginaryVector;
  std::complex<double> potential(realTotal+coulombPotential,imaginaryTotal);
  return potential;
}

const double JLMPotential::a[3][3] = {{-0.9740e3 , 0.7097e4 ,-0.1953e5 },
				      { 0.1126e2 ,-0.1257e3 , 0.4180e3 },
				      {-0.4250e-1, 0.5853e0 ,-0.2054e1 }};
const double JLMPotential::b[3][3] = {{ 0.3601e3 ,-0.2691e4 , 0.7733e4 },
				      {-0.5224e1 , 0.5130e2 ,-0.1717e3 },
				      { 0.2051e-1,-0.2470e0 , 0.8846e0 }};
const double JLMPotential::c[3][3] = {{ 0.4557e1 ,-0.2051e1 ,-0.6509e2 },
				      {-0.5291e-2,-0.4906e0 , 0.3095e1 },
				      { 0.6108e-5, 0.1812e-2,-0.1190e-1}};
const double JLMPotential::dHighE[4][4] = {{-0.1483e4 , 0.2988e5 ,-0.2128e6 , 0.5125e6 },
					   { 0.3718e2 ,-0.9318e3 , 0.7209e4 ,-0.1796e5 },
					   {-0.3549e0 , 0.9591e1 ,-0.7752e2 , 0.1980e3 },
					   { 0.1119e-2,-0.3160e-1, 0.2611e0 ,-0.6753e0 }};
const double JLMPotential::fHighE[4][4] = {{ 0.5461e3 ,-0.8471e4 , 0.5172e5 ,-0.1140e6},
					   {-0.1120e2 , 0.2300e3 ,-0.1520e4 , 0.3543e4},
					   { 0.1065e0 ,-0.2439e1 , 0.1717e2 ,-0.4169e2},
					   {-0.3541e-3, 0.8544e-2,-0.6211e-1, 0.1537e0}};
const double JLMPotential::dLowE[4][4] = {{-0.5138e3 , 0.9078e4,-0.6192e5, 0.1516e6},
					  {-0.2985e2 , 0.5757e3,-0.4155e4, 0.1037e5},
					  { 0.1452e1 ,-0.3435e2, 0.2657e3,-0.6748e3},
					  { 0.9428e-1,-0.2310e1, 0.1882e2,-0.5014e2}};
const double JLMPotential::fLowE[4][4] = {{ 0.6597e3 ,-0.1263e5 , 0.9428e5 ,-0.2453e6 },
					  { 0.4509e1 ,-0.6572e2 , 0.5972e3 ,-0.1800e4 },
					  {-0.2383e1 , 0.5866e2 ,-0.4923e3 , 0.1358e4 },
					  {-0.4324e-1, 0.1348e1 ,-0.1295e2 , 0.3836e2 }};

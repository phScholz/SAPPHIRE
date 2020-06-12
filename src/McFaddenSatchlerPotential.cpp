#include "McFaddenSatchlerPotential.h"

McFaddenSatchlerPotential::McFaddenSatchlerPotential(int z1, int m1, int z2, int m2, 
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
  //Potential Parameters			   
  aR_ = 0.52;
  aI_ = 0.52;
  V_ = 185.;
  W_ = 25.;
  double r0R = 1.4;
  double r0I = 1.4;
  
  coulombRadius_=1.3*pow(m2_,0.33333333333333333333);
  rR_ = r0R*pow(m2_,0.33333333333333333333);
  rI_ = r0I*pow(m2_,0.33333333333333333333);  
  boundaryRadius_=rR_+7.*aR_;
}


std::complex<double> 
McFaddenSatchlerPotential::Calculate(double r, int l, double s, double j, double energy) const {
  std::complex<double> nuclearPotential= std::complex<double>(-V_/(1.+exp((r-rR_)/aR_)),-W_/(1.+exp((r-rI_)/aI_)));
  
  double coulombPotential=0.0;

  if(r<=coulombRadius_)
  {
      coulombPotential=hbarc*fstruc*GetZ1Z2()/coulombRadius_*(1.0+0.5*(1-pow(r/coulombRadius_,2.0)));
  } 
  else
  {
    coulombPotential=hbarc*fstruc*GetZ1Z2()/r;
  }

  return nuclearPotential+coulombPotential;
}
 

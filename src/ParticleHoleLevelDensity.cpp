#include "ParticleHoleLevelDensity.h"
#include "gsl/gsl_sf_gamma.h"
#include "math.h"
#include <algorithm>
#include <iostream>
#include "NuclearMass.h"

ParticleHoleLevelDensity::ParticleHoleLevelDensity(int Z, int A, double J,
						   int neutronNumber,
						   int neutronHoleNumber,
						   int protonNumber,
						   int protonHoleNumber) : 
  A_(A),Z_(Z),totalExcitonNumber_(0),preFactor_(0.),pauliCorrection_(0.),
  pairingEnergy_(0.),nCritical_(0.),spinFactor_(0.),gPi_(0.),gNu_(0.),
  invalid_(false) {
  
  if(neutronNumber<0||
     neutronHoleNumber<0||
     protonNumber<0||
     protonHoleNumber<0||
     (neutronNumber+neutronHoleNumber+protonNumber+protonHoleNumber)<=0) {
    invalid_=true;
    return;
  } 

  incident_=0;
  if(neutronHoleNumber==0&&protonHoleNumber==1) {
    if(protonNumber==2&&neutronNumber==0) incident_=1;
    else if(protonNumber==1&&neutronNumber==1) incident_=2;
  } else if(neutronHoleNumber==1&&protonHoleNumber==0) {
    if(protonNumber==1&&neutronNumber==1) incident_=1;
    else if(protonNumber==0&&neutronNumber==2) incident_=2;
  }

  if(incident_==1) {
    double qValue = 0.;
    NuclearMass::QValue(Z,A,Z-1,A-1,qValue);
    separationEnergy_ = -1.*qValue;
  } else if(incident_==2) {
    double qValue = 0.;
    NuclearMass::QValue(Z,A,Z,A-1,qValue);
    separationEnergy_ = -1.*qValue;
  }

  gPi_ = Z/15.;
  gNu_ = (A-Z)/15.;

  totalExcitonNumber_ = 
    protonNumber+protonHoleNumber+
    neutronNumber+neutronHoleNumber;
  totalHoleNumber_ = protonHoleNumber+neutronHoleNumber;

  double spinCutoff = sqrt(0.24*totalExcitonNumber_*pow(A,0.6666667));
  spinFactor_ = 0.5*pow(2.*J+1.,2.)/sqrt(3.14159)/pow(totalExcitonNumber_,1.5)/
    pow(spinCutoff,3.)*exp(-pow(J+0.5,2.)/totalExcitonNumber_/pow(spinCutoff,2.));

  preFactor_ = pow(gNu_,double(neutronNumber+neutronHoleNumber))*
    pow(gPi_,double(protonNumber+protonHoleNumber))/
    gsl_sf_fact(protonNumber)/
    gsl_sf_fact(protonHoleNumber)/
    gsl_sf_fact(neutronNumber)/
    gsl_sf_fact(neutronHoleNumber)/
    gsl_sf_fact(totalExcitonNumber_-1);

  pauliCorrection_ = PauliCorrection(neutronNumber,neutronHoleNumber,
				     protonNumber,protonHoleNumber);

  if((A-Z)%2==0&&Z%2==0)
    pairingEnergy_ = 24./sqrt(A);
  else if(((A-Z)%2!=0&&Z%2==0)||
	  ((A-Z)%2==0&&Z%2!=0))
    pairingEnergy_ = 12./sqrt(A);
  else
    pairingEnergy_ = 0.;

  nCritical_ = 4.*sqrt(pairingEnergy_*4.*(gPi_+gNu_))/3.5*log(2.);
}

double ParticleHoleLevelDensity::PairingCorrection(double energy) {
  if(invalid_||pairingEnergy_==0.) return 0.;
  if(energy/pairingEnergy_>=0.716+2.44*pow(totalExcitonNumber_/nCritical_,2.17))
    return pairingEnergy_*(1.-pow(0.996-1.76*pow(totalExcitonNumber_/nCritical_,1.6)/pow(energy/pairingEnergy_,0.68),2.));
  else return pairingEnergy_;
}

double ParticleHoleLevelDensity::operator()(double energy, bool correct, bool spin) {
  if(invalid_) return 0.;
  double U = (correct) ? energy-PairingCorrection(energy)-pauliCorrection_ : energy-pauliCorrection_;
  if(U<=0.) return 0.;
  if(spin) return preFactor_*pow(U,totalExcitonNumber_-1)*FiniteDepth(energy)*spinFactor_;
  else return preFactor_*pow(U,totalExcitonNumber_-1)*FiniteDepth(energy);
}

double ParticleHoleLevelDensity::PauliCorrection(int neutronNumber,
						 int neutronHoleNumber,
						 int protonNumber,
						 int protonHoleNumber) {
  if(invalid_) return 0.;
  return pow(std::max(protonNumber,protonHoleNumber),2.)/gPi_+
    pow(std::max(neutronNumber,neutronHoleNumber),2.)/gNu_-
    (pow(protonNumber,2.)+pow(protonHoleNumber,2.)
     +protonNumber+protonHoleNumber)/4./gPi_-
    (pow(neutronNumber,2.)+pow(neutronHoleNumber,2.)
     +neutronNumber+neutronHoleNumber)/4./gNu_;
}

double ParticleHoleLevelDensity::FiniteDepth(double energy) {
  double V=38.;
  if(energy>separationEnergy_&&incident_!=0) {
    double ep4 = pow(energy-separationEnergy_,4.);
    if(incident_==1&&energy>separationEnergy_) {
      V=22.+16.*ep4/(ep4+pow(450./pow(A_,0.33333333),4.));
    } else if(incident_==2&&energy>separationEnergy_) {
      V=22.+16.*ep4/(ep4+pow(245./pow(A_,0.33333333),4.));
    }
  }
  double fact=1.;
  double coeff=1.;
  for(int i = 1;i<=totalHoleNumber_;i++) {
    coeff*=double(totalHoleNumber_+1-i)/i;
    if(energy<i*V) continue;
    fact+=pow(-1.,i)*coeff*pow((energy-i*V)/energy,totalExcitonNumber_-1);
  }
  return fact;
}

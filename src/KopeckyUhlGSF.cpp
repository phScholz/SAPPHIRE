#include "KopeckyUhlGSF.h"
#include "LevelDensity.h"
#include <iostream>

KopeckyUhlGSF::KopeckyUhlGSF(int z2, int m2, double jInitial, int piInitial,
			     double jFinal, int piFinal, double maxL, 
			     double totalWidthForCorrection,
			     double uncorrTotalWidthForCorrection,
			     double uncorrTotalWidthSqrdForCorrection,
			     TransmissionFunc* previous,
			     LevelDensity* levelDensity, double compoundE) :
  GammaTransmissionFunc(z2,m2,jInitial,piInitial,jFinal,piFinal,maxL,
			totalWidthForCorrection,uncorrTotalWidthForCorrection,
			uncorrTotalWidthSqrdForCorrection,previous),
  levelDensity_(levelDensity), compoundE_(compoundE){
  double qValueNeutron;
  if(!NuclearMass::QValue(z2,m2,z2,m2-1,qValueNeutron)) {
    std::cout << "Unknown masses requested.  Aborting." << std::endl;
    exit(1);
  }
  levelDenParam_ = levelDensity->CalcDensityParam(-1.*qValueNeutron-levelDensity_->backshift_);
}

double KopeckyUhlGSF::CalcStrengthFunction(double energy) {
  double strength = 0.;
  for(int i = 0;i<2;i++) {
    strength+=gdrParameters_.kSigmaGamma_[i]*
      (energy*CalcEnergyDependentWidth(energy,i)/
       (pow(pow(energy,2.)-pow(gdrParameters_.E_[i],2.),2.)+
	pow(CalcEnergyDependentWidth(energy,i)*energy,2.))+
       0.7*gdrParameters_.W_[i]*4.*pow(pi,2.)*CalcKUTempSqrd(energy)/
       pow(gdrParameters_.E_[i],5.));
  }
  return strength;
}

double KopeckyUhlGSF::CalcEnergyDependentWidth(double energy, int which) {
  double width =  gdrParameters_.W_[which]*(pow(energy,2.0)+4.*pow(pi,2.)*CalcKUTempSqrd(energy))/pow(gdrParameters_.E_[which],2.);
  return width;
}

double KopeckyUhlGSF::CalcKUTempSqrd(double energy) {
  double e = compoundE_-levelDensity_->backshift_-energy;
  return (e>0) ? e/levelDenParam_ : 0.;
}

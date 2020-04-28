#include "LevelDensity/RauscherLevelDensity.h"

#include <math.h>

void RauscherLevelDensity::CalcBackShift() {
  double neutronPairingGap,protonPairingGap;
  if(!NuclearMass::NeutronPairingGap(Z_,A_,neutronPairingGap)||
     !NuclearMass::ProtonPairingGap(Z_,A_,protonPairingGap)) {
    //std::cout << "Mass not found in calculation of pairing gaps. Aborting." 
    //	      << std::endl;
    //exit(1);
    double chi;
    if(Z_%2==0&&(A_-Z_)%2==0) chi = 1.;
    else if(Z_%2==0||(A_-Z_)%2==0) chi = 0.;
    else chi = -1.;
    backshift_ = chi*12./sqrt(double(A_));
  } else backshift_ = (neutronPairingGap+protonPairingGap)/2.;
}

double RauscherLevelDensity::CalcDensityParam(double u) {
  double f = 1.-exp(-1.*gamma_*u);
  return aTilda_*(1.+microEnergyCorr_*f/u);
}

double RauscherLevelDensity::CalcNuclearTemp(double u) {
  double a = CalcDensityParam(u);
  double inverseNuclearTemp=aTilda_/sqrt(u*a)*
    (1.+gamma_*microEnergyCorr_*exp(-1.*gamma_*u))-3./2./u;
  return 1./inverseNuclearTemp;
}

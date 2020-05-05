#include "Potential.h"
#include "Constants.h"
#include "CoulFunc.h"
#include <iostream>
#include <float.h>


Potential::Potential(int z1, int m1, int z2, int m2, 
		     double jInitial, int piInitial,
		     double jFinal, int piFinal,
		     double spin, int parity, double maxL,
		     double totalWidthForCorrection,
		     double uncorrTotalWidthForCorrection,
		     double uncorrTotalWidthSqrdForCorrection,
		     TransmissionFunc* previous) :
  ParticleTransmissionFunc(z1,m1,z2,m2,jInitial,piInitial,
			   jFinal,piFinal,spin,parity,maxL,
			   totalWidthForCorrection,
			   uncorrTotalWidthForCorrection,
			   uncorrTotalWidthSqrdForCorrection,
			   previous) {
  coulFunc_ = new CoulFunc(z1,z2,GetRedMass(),true);
}

Potential::~Potential() {
  delete coulFunc_;
}

std::vector<std::complex<double> > Potential::Solve(double energy,int L,double S,double J,double rStep) const {
  std::vector<std::complex<double> > waveVector;
  std::vector<std::complex<double> > potentialVector;
  std::vector<double> radiusVector;
  
  radiusVector.push_back(0.0);
  radiusVector.push_back(rStep);
  potentialVector.push_back(Calculate(radiusVector[0],L,S,J,energy));
  potentialVector.push_back(Calculate(radiusVector[1],L,S,J,energy));
  if(L%2==0) {
    waveVector.push_back(std::complex<double>(0.1,0.1));
    waveVector.push_back(std::complex<double>(0.1,0.1));
  } else {
    waveVector.push_back(std::complex<double>(0.0,0.0));
    waveVector.push_back(std::complex<double>(rStep,rStep));
  }

  int i=2;
  while(i*rStep<=GetBoundaryRadius()+3.*rStep) {
    radiusVector.push_back(i*rStep);
    potentialVector.push_back(Calculate(radiusVector[i],L,S,J,energy));
    std::complex<double> tempWave=(2.0+pow(rStep,2.0)*2.0*uconv/pow(hbarc,2.0)*GetRedMass()*(potentialVector[i-1]-energy)
				   -2.0*rStep/radiusVector[i-1]+
				   pow(rStep,2.0)*(float)L*((float)L+1.0)/pow(radiusVector[i-1],2.0))*
      waveVector[i-1]+waveVector[i-2]*(2.0*rStep/radiusVector[i-1]-1.0);
    waveVector.push_back(tempWave);
    i++;
  }
  return waveVector;
}

void Potential::NormalizeInternally(std::vector<std::complex<double> > &waveFunction, double rStep) const {
  double sum=0.0;
  for(int i=0;i<=int(GetBoundaryRadius()/rStep);i++) {
    sum+=real(waveFunction[i]*conj(waveFunction[i]))*pow(i*rStep,2.0)*rStep;
  }
  for(int i=0;i<waveFunction.size();i++) {
    double normalizedRealValue = real(waveFunction[i])*i*rStep/sqrt(sum);
    double normalizedImagValue = imag(waveFunction[i])*i*rStep/sqrt(sum);
    waveFunction[i]=std::complex<double>(normalizedRealValue,normalizedImagValue);
  }
}

void Potential::NormalizeOverAllSpace(std::vector<std::complex<double> > &waveFunction,double rStep) const {
  double sum=0.0;
  for(int i=0;i<waveFunction.size();i++) {
    sum+=real(waveFunction[i]*conj(waveFunction[i]))*pow(i*rStep,2.0)*rStep;
  }
  for(int i=0;i<waveFunction.size();i++) {
    double normalizedRealValue = real(waveFunction[i])*i*rStep/sqrt(sum);
    double normalizedImagValue = imag(waveFunction[i])*i*rStep/sqrt(sum);
    waveFunction[i]=std::complex<double>(normalizedRealValue,normalizedImagValue);
  }
}

std::complex<double> Potential::CalcBeta(double energy,int l, double s, double j, double rStep) const {
  std::vector<std::complex<double> > waveFunction=Solve(energy,l,s,j,rStep);
  NormalizeInternally(waveFunction,rStep);

  int iAtBoundary=int(GetBoundaryRadius()/rStep);
  std::complex<double> beta=1.0/12.0/rStep*GetBoundaryRadius()/waveFunction[iAtBoundary]*
    (-waveFunction[iAtBoundary+2]+8.0*waveFunction[iAtBoundary+1]
     -8.0*waveFunction[iAtBoundary-1]+waveFunction[iAtBoundary-2]);
  
  return beta;
}

double Potential::CalcTransmission(double s, int l, double energy)  {
  double rStep=1./(3.+sqrt(2.*uconv*m1_*energy)/hbarc);
  std::complex<double> beta = CalcBeta( energy, l,  s,  jInitial_,rStep);
  CoulWaves coulWaves = coulFunc_->operator()(l,GetBoundaryRadius(),energy);
  
  double rho = sqrt(2.*uconv)/hbarc*GetBoundaryRadius()*sqrt(GetRedMass()*energy);
  std::complex<double> out = std::complex<double>(coulWaves.G,coulWaves.F);  
  std::complex<double> diffOut = std::complex<double>(coulWaves.dG,coulWaves.dF);  
  std::complex<double> lFunction = rho*diffOut/out;

  double pene = rho/(coulWaves.F*coulWaves.F+coulWaves.G*coulWaves.G);
  double transmission = (pene>0.) ? -1.*real(4.*pene*imag(beta)/(lFunction-beta)/conj(lFunction-beta)) : 0.;
  if(m1_==1 && z1_==1) return pNorm_*transmission;
  if(m1_==1 && z1_==0) return nNorm_*transmission;
  if(m1_==4 && z1_==2) return aNorm_*transmission;
}

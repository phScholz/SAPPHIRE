#include "PreEqTransitionRateFunc.h"
#include "ParticleTransmissionFunc.h"
#include <iostream>

PreEqTransitionRateFunc::PreEqTransitionRateFunc(int z1, int m1, int z2, int m2,
						 int initialNeutronNumber,
						 int initialNeutronHoleNumber,
						 int initialProtonNumber,
						 int initialProtonHoleNumber,
						 int finalNeutronNumber,
						 int finalNeutronHoleNumber,
						 int finalProtonNumber,
						 int finalProtonHoleNumber,
						 double jInitial, int piInitial, 
						 double jFinal, int piFinal,
						 double spin, int parity, double maxL,
						 double compoundE, double qValue) :
  z1_(z1),m1_(m1),z2_(z2),m2_(m2),
  initialNeutronNumber_(initialNeutronNumber), initialNeutronHoleNumber_(initialNeutronHoleNumber),
  initialProtonNumber_(initialProtonNumber), initialProtonHoleNumber_(initialProtonHoleNumber),
  finalNeutronNumber_(finalNeutronNumber), finalNeutronHoleNumber_(finalNeutronHoleNumber),
  finalProtonNumber_(finalProtonNumber), finalProtonHoleNumber_(finalProtonHoleNumber),
  piInitial_(piInitial), piFinal_(piFinal), parity_(parity), jInitial_(jInitial), jFinal_(jFinal),
  spin_(spin),maxL_(maxL),compoundE_(compoundE),qValue_(qValue) {

  int n = initialNeutronNumber_+initialNeutronHoleNumber_+
    initialProtonNumber_+initialProtonHoleNumber_;
  M2_ =2.*3.14159/6.58211814e-22*c1_/pow(m2_+m1_,3.)*
    (7.48*c2_+4.62e5/pow(compoundE_/n+10.7*c3_,3.));

  if((initialNeutronNumber==finalNeutronNumber+1&&
      initialNeutronHoleNumber==finalNeutronHoleNumber&&
      initialProtonNumber==finalProtonNumber&&
      initialProtonHoleNumber==finalProtonHoleNumber)||
     (initialNeutronNumber==finalNeutronNumber&&
      initialNeutronHoleNumber==finalNeutronHoleNumber&&
      initialProtonNumber==finalProtonNumber+1&&
      initialProtonHoleNumber==finalProtonHoleNumber)) 
    CalcParticleDecay();
  else if(initialNeutronNumber==finalNeutronNumber-1&&
	  initialNeutronHoleNumber==finalNeutronHoleNumber-1&&
	  initialProtonNumber==finalProtonNumber&&
	  initialProtonHoleNumber==finalProtonHoleNumber) 
    CalcNeutronPairProduction();
  else if(initialNeutronNumber==finalNeutronNumber&&
	  initialNeutronHoleNumber==finalNeutronHoleNumber&&
	  initialProtonNumber==finalProtonNumber-1&&
	  initialProtonHoleNumber==finalProtonHoleNumber-1) 
    CalcProtonPairProduction();
  else if(initialNeutronNumber==finalNeutronNumber+1&&
	  initialNeutronHoleNumber==finalNeutronHoleNumber+1&&
	  initialProtonNumber==finalProtonNumber-1&&
	  initialProtonHoleNumber==finalProtonHoleNumber-1) 
    CalcNeutronPairConversion();
  else if(initialNeutronNumber==finalNeutronNumber-1&&
	  initialNeutronHoleNumber==finalNeutronHoleNumber-1&&
	  initialProtonNumber==finalProtonNumber+1&&
	  initialProtonHoleNumber==finalProtonHoleNumber+1) 
    CalcProtonPairConversion();
}

void PreEqTransitionRateFunc::CalcParticleDecay() {
  ParticleHoleLevelDensity* initialLevelDensity = 
    new ParticleHoleLevelDensity(z2_+z1_,m2_+m1_,jInitial_,
				 initialNeutronNumber_,
				 initialNeutronHoleNumber_,
				 initialProtonNumber_,
				 initialProtonHoleNumber_);
  ParticleHoleLevelDensity* finalLevelDensity = 
    new ParticleHoleLevelDensity(z2_,m2_,jFinal_,
				 finalNeutronNumber_,
				 finalNeutronHoleNumber_,
				 finalProtonNumber_,
				 finalProtonHoleNumber_);
  ParticleTransmissionFunc* transmissionFunc = 
    ParticleTransmissionFunc::CreateParticleTransmissionFunc(z1_,m1_,z2_,m2_,
							     jInitial_,piInitial_,
							     jFinal_,piFinal_,
							     spin_,parity_,maxL_,
							     0., 0., 0., NULL);
  if(!transmissionFunc->IsValid()) {
    std::cout << "An problem occurred creating transmission function. Aborting."
	      << std::endl;
    delete transmissionFunc;
    exit(1);
  }
  double sum=0.;
  double dE=0.01;
  double lowEnergy =  0.001;
  double highEnergy = compoundE_+qValue_;
  double initialDensity = initialLevelDensity->operator()(compoundE_,true,true);
  for(double ep = lowEnergy;ep<highEnergy;ep+=dE) {
    double rate = (2.*jInitial_+1.)/(2.*jFinal_+1.)/2./3.14159/6.58211814e-22*
      finalLevelDensity->operator()(compoundE_+qValue_-ep,true,true)*
      transmissionFunc->operator()(ep)*dE/initialDensity;
    sum+=rate;
    function_.push_back(XYPair(ep,rate));
    cumulativeSum_.push_back(XYPair(ep,sum));
  }
  integral_=sum;
  delete initialLevelDensity;
  delete finalLevelDensity;
  delete transmissionFunc;
}

void PreEqTransitionRateFunc::CalcNeutronPairProduction() {
  double sum=0.;
  int nn = initialNeutronNumber_;
  int nhn = initialNeutronHoleNumber_;
  int pn = initialProtonNumber_;
  int phn = initialProtonHoleNumber_;

  ParticleHoleLevelDensity* levelDensity[9];
  double lowerLimit[4], upperLimit[4], dE[4];

  //neutron particle scattering
  levelDensity[0] = new ParticleHoleLevelDensity(z2_,m2_,jFinal_,nn,nhn,pn,phn);
  levelDensity[1] = new ParticleHoleLevelDensity(z2_,m2_,jFinal_,nn-1,nhn,pn,phn);
  levelDensity[2] = new ParticleHoleLevelDensity(z2_,m2_,jFinal_,2,1,0,0);
  levelDensity[3] = new ParticleHoleLevelDensity(z2_,m2_,jFinal_,nn,nhn-1,pn,phn);
  levelDensity[4] = new ParticleHoleLevelDensity(z2_,m2_,jFinal_,1,2,0,0);
  levelDensity[5] = new ParticleHoleLevelDensity(z2_,m2_,jFinal_,nn,nhn,pn-1,phn);
  levelDensity[6] = new ParticleHoleLevelDensity(z2_,m2_,jFinal_,1,1,1,0);
  levelDensity[7] = new ParticleHoleLevelDensity(z2_,m2_,jFinal_,nn,nhn,pn,phn-1);
  levelDensity[8] = new ParticleHoleLevelDensity(z2_,m2_,jFinal_,1,1,0,1);

  double gNu = levelDensity[0]->gNu();
  double gPi = levelDensity[0]->gPi();
  double U = compoundE_-levelDensity[0]->PairingCorrection(compoundE_);

  double finalStatePauliEnergy = levelDensity[0]->PauliCorrection(nn+1,nhn+1,pn,phn);
  lowerLimit[0] = finalStatePauliEnergy-levelDensity[0]->PauliCorrection(nn-1,nhn,pn,phn);
  lowerLimit[1] = finalStatePauliEnergy-levelDensity[0]->PauliCorrection(nn,nhn-1,pn,phn);
  lowerLimit[2] = finalStatePauliEnergy-levelDensity[0]->PauliCorrection(nn,nhn,pn-1,phn);
  lowerLimit[3] = finalStatePauliEnergy-levelDensity[0]->PauliCorrection(nn,nhn,pn,phn-1);

  upperLimit[0] = U-levelDensity[0]->PauliCorrection(nn-1,nhn,pn,phn);
  upperLimit[1] = U-levelDensity[0]->PauliCorrection(nn,nhn-1,pn,phn);
  upperLimit[2] = U-levelDensity[0]->PauliCorrection(nn,nhn,pn-1,phn);
  upperLimit[3] = U-levelDensity[0]->PauliCorrection(nn,nhn,pn,phn-1);

  if(nn==1&&nhn==0&&pn==0&&phn==0) {
    sum = levelDensity[2]->operator()(U,false)*M2_*1.5;
  } else if(nn==0&&nhn==0&&pn==1&&phn==0) {
    sum = levelDensity[6]->operator()(U,false)*M2_;
  } else {
    for(int i = 0;i<4;++i) 
      dE[i] = fabs(upperLimit[i]-lowerLimit[i])/integrationSteps_;
    for(int i = 1;i<=integrationSteps_;++i) {
      for(int j = 0;j<4;++j) {
        double E = lowerLimit[j]+(double(i)-0.5)*dE[j];
        double singleParticleDensity = (j<2) ? gNu : gPi;		 
        double M2 = (j<2) ? M2_*1.5 : M2_;
        if(E<upperLimit[j]) 
  	  sum+=levelDensity[1+j*2]->operator()(U-E,false)*singleParticleDensity*
	    levelDensity[2+j*2]->operator()(E,false)*M2*dE[j];
      }
    }
    sum/=levelDensity[0]->operator()(U,false);
  }
  
  for(int i = 0;i<9;i++) delete levelDensity[i];
  
  function_.push_back(XYPair(0.,sum));
  cumulativeSum_.push_back(XYPair(0.,sum));
  integral_=sum;
}

void PreEqTransitionRateFunc::CalcProtonPairProduction() {
  double sum = 0.;
  int nn = initialNeutronNumber_;
  int nhn = initialNeutronHoleNumber_;
  int pn = initialProtonNumber_;
  int phn = initialProtonHoleNumber_;

  ParticleHoleLevelDensity* levelDensity[9];
  double lowerLimit[4], upperLimit[4], dE[4];

  //neutron particle scattering
  levelDensity[0] = new ParticleHoleLevelDensity(z2_,m2_,jFinal_,nn,nhn,pn,phn);
  levelDensity[1] = new ParticleHoleLevelDensity(z2_,m2_,jFinal_,nn-1,nhn,pn,phn);
  levelDensity[2] = new ParticleHoleLevelDensity(z2_,m2_,jFinal_,1,0,1,1);
  levelDensity[3] = new ParticleHoleLevelDensity(z2_,m2_,jFinal_,nn,nhn-1,pn,phn);
  levelDensity[4] = new ParticleHoleLevelDensity(z2_,m2_,jFinal_,0,1,1,1);
  levelDensity[5] = new ParticleHoleLevelDensity(z2_,m2_,jFinal_,nn,nhn,pn-1,phn);
  levelDensity[6] = new ParticleHoleLevelDensity(z2_,m2_,jFinal_,0,0,2,1);
  levelDensity[7] = new ParticleHoleLevelDensity(z2_,m2_,jFinal_,nn,nhn,pn,phn-1);
  levelDensity[8] = new ParticleHoleLevelDensity(z2_,m2_,jFinal_,0,0,1,2);

  double gNu = levelDensity[0]->gNu();
  double gPi = levelDensity[0]->gPi();
  double U = compoundE_-levelDensity[0]->PairingCorrection(compoundE_);

  double finalStatePauliEnergy = levelDensity[0]->PauliCorrection(nn,nhn,pn+1,phn+1);
  lowerLimit[0] = finalStatePauliEnergy-levelDensity[0]->PauliCorrection(nn-1,nhn,pn,phn);
  lowerLimit[1] = finalStatePauliEnergy-levelDensity[0]->PauliCorrection(nn,nhn-1,pn,phn);
  lowerLimit[2] = finalStatePauliEnergy-levelDensity[0]->PauliCorrection(nn,nhn,pn-1,phn);
  lowerLimit[3] = finalStatePauliEnergy-levelDensity[0]->PauliCorrection(nn,nhn,pn,phn-1);

  upperLimit[0] = U-levelDensity[0]->PauliCorrection(nn-1,nhn,pn,phn);
  upperLimit[1] = U-levelDensity[0]->PauliCorrection(nn,nhn-1,pn,phn);
  upperLimit[2] = U-levelDensity[0]->PauliCorrection(nn,nhn,pn-1,phn);
  upperLimit[3] = U-levelDensity[0]->PauliCorrection(nn,nhn,pn,phn-1);

  if(nn==1&&nhn==0&&pn==0&&phn==0) {
    sum = levelDensity[2]->operator()(U,false)*M2_;
  } else if(nn==0&&nhn==0&&pn==1&&phn==0) {
    sum = levelDensity[6]->operator()(U,false)*M2_;
  } else {
    for(int i = 0;i<4;++i) 
      dE[i] = fabs(upperLimit[i]-lowerLimit[i])/integrationSteps_;
    for(int i = 1;i<=integrationSteps_;++i) {
      for(int j = 0;j<4;++j) {
        double E = lowerLimit[j]+(double(i)-0.5)*dE[j];
        double singleParticleDensity = (j<2) ? gNu : gPi;		 
        if(E<upperLimit[j]) 
  	  sum+=levelDensity[1+j*2]->operator()(U-E,false)*singleParticleDensity*
	    levelDensity[2+j*2]->operator()(E,false)*M2_*dE[j];
      }
    }
    sum/=levelDensity[0]->operator()(U,false);
  }
  
  for(int i = 0;i<9;i++) delete levelDensity[i];
  
  function_.push_back(XYPair(0.,sum));
  cumulativeSum_.push_back(XYPair(0.,sum));
  integral_=sum;
}

void PreEqTransitionRateFunc::CalcNeutronPairConversion() {
  double sum=0.;
  int nn = initialNeutronNumber_;
  int nhn = initialNeutronHoleNumber_;
  int pn = initialProtonNumber_;
  int phn = initialProtonHoleNumber_;

  ParticleHoleLevelDensity* levelDensity[4];
  levelDensity[0] = new ParticleHoleLevelDensity(z2_,m2_,jFinal_,nn,nhn,pn,phn);
  levelDensity[1] = new ParticleHoleLevelDensity(z2_,m2_,jFinal_,nn-1,nhn-1,pn,phn);
  levelDensity[2] = new ParticleHoleLevelDensity(z2_,m2_,jFinal_,1,1,0,0);
  levelDensity[3] = new ParticleHoleLevelDensity(z2_,m2_,jFinal_,0,0,1,1);
  
  double U = compoundE_-levelDensity[0]->PairingCorrection(compoundE_);
  
  double lowerLimit = levelDensity[0]->PauliCorrection(nn,nhn,pn,phn)-
    levelDensity[0]->PauliCorrection(nn-1,nhn-1,pn,phn);
  double upperLimit = U-levelDensity[0]->PauliCorrection(nn-1,nhn-1,pn,phn);
  double dE = fabs(upperLimit-lowerLimit)/integrationSteps_;

  for(int i = 1;i<=integrationSteps_;++i) {
    double E = lowerLimit+(double(i)-0.5)*dE;
    if(E<upperLimit) 
      sum+=levelDensity[1]->operator()(U-E,false)*levelDensity[2]->operator()(E,false)*
	levelDensity[3]->operator()(E,false)*M2_*dE;
  }
  sum/=levelDensity[0]->operator()(U,false);

  for(int i = 0;i<4;i++) delete levelDensity[i];
  
  function_.push_back(XYPair(0.,sum));
  cumulativeSum_.push_back(XYPair(0.,sum));
  integral_=sum;
}

void PreEqTransitionRateFunc::CalcProtonPairConversion() {
  double sum=0.;
  int nn = initialNeutronNumber_;
  int nhn = initialNeutronHoleNumber_;
  int pn = initialProtonNumber_;
  int phn = initialProtonHoleNumber_;

  ParticleHoleLevelDensity* levelDensity[4];
  levelDensity[0] = new ParticleHoleLevelDensity(z2_,m2_,jFinal_,nn,nhn,pn,phn);
  levelDensity[1] = new ParticleHoleLevelDensity(z2_,m2_,jFinal_,nn,nhn,pn-1,phn-1);
  levelDensity[2] = new ParticleHoleLevelDensity(z2_,m2_,jFinal_,0,0,1,1);
  levelDensity[3] = new ParticleHoleLevelDensity(z2_,m2_,jFinal_,1,1,0,0);
  
  double U = compoundE_-levelDensity[0]->PairingCorrection(compoundE_);
  
  double lowerLimit = levelDensity[0]->PauliCorrection(nn,nhn,pn,phn)-
    levelDensity[0]->PauliCorrection(nn,nhn,pn-1,phn-1);
  double upperLimit = U-levelDensity[0]->PauliCorrection(nn,nhn,pn-1,phn-1);
  double dE = fabs(upperLimit-lowerLimit)/integrationSteps_;

  for(int i = 1;i<=integrationSteps_;++i) {
    double E = lowerLimit+(double(i)-0.5)*dE;
    if(E<upperLimit) 
      sum+=levelDensity[1]->operator()(U-E,false)*levelDensity[2]->operator()(E,false)*
	levelDensity[3]->operator()(E,false)*M2_*dE;
  }
  sum/=levelDensity[0]->operator()(U,false);

  for(int i = 0;i<4;i++) delete levelDensity[i];
  
  function_.push_back(XYPair(0.,sum));
  cumulativeSum_.push_back(XYPair(0.,sum));
  integral_=sum;
}

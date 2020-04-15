#include "NuclearLevels.h"
#include "CrossSection.h"
#include "Decayer.h"
#include "TransitionRateFunc.h"
#include "RauscherLevelDensity.h"
#include "Constants.h"
#include <math.h>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <iostream>
#include <sstream>
#include <gsl/gsl_integration.h>
#include <TGraph.h>
#include <algorithm>

/* CrossSection::CrossSection_simple(int Z, int A, int pType, std::string energyFile, bool forRates, int entranceState, std::vector<int> exitStates)
{ 
  // FindLevels in the levels data table 
  std::vector<Level> knownLevels = NuclearLevels::FindLevels(Z,A);

  if(!knownLevels.size()>entranceState){
    std::cout << "Initial state not known." << std::endl;
  }



}
 */
CrossSection::CrossSection(int Z, int A, int pType, std::string energyFile, bool forRates,
			   int entranceState, std::vector<int> exitStates) : 
  Z_(Z), A_(A), pType_(pType), skipEnergy_(1000.), entranceState_(entranceState), exitStates_(exitStates) {

  groundStateJ_ =-1.;
  
  std::vector<Level> knownLevels = NuclearLevels::FindLevels(Z,A);
  
  if(knownLevels.size()>entranceState) {
    Level groundState = knownLevels[entranceState];
    std::string parity = (groundState.Pi_<0) ? "-" : "+";
      if(entranceState!=0) std::cout << "Initial state: E=" << groundState.energy_
				                          << " MeV, J=" << groundState.J_ << parity << std::endl;
    
      if(entranceState>0||groundState.energy_==0.) {
        groundStateJ_ = groundState.J_;
        groundStatePi_ = groundState.Pi_;
      }
  } else {
    std::cout << "Initial state not known." << std::endl;
  }
  
  if(groundStateJ_==-1.) {
    isValid_=false;
    return;
  }

  double qValue;
  
  if(pType_ == 0) {
    qValue=0.;
    compoundZ_ = Z;
    compoundA_ = A;
    double spinNorm=2.*(2.*groundStateJ_+1.); 
    preFactor_=hbarc*hbarc/100.*pi/spinNorm;
  } else if(pType_ == 1) {
    
    if(!NuclearMass::QValue(Z,A+1,Z,A,qValue)) {
      isValid_=false;
      return;
    }
    compoundZ_ = Z;
    compoundA_ = A+1;
    double reducedMass =double(A)/double(A+1);
    double spinNorm=2.*(2.*groundStateJ_+1.);
    preFactor_ = hbarc*hbarc/200.*pi/uconv/reducedMass/spinNorm;
  
  } else if(pType_ == 2) {
    
    if(!NuclearMass::QValue(Z+1,A+1,Z,A,qValue)) {
      isValid_=false;
      return;
    }
    compoundZ_ = Z+1;
    compoundA_ = A+1;
    double reducedMass=double(A)/double(A+1);
    double spinNorm=2.*(2.*groundStateJ_+1.);
    preFactor_ = hbarc*hbarc/200.*pi/uconv/reducedMass/spinNorm;
  
  } else if(pType_ ==3) {
    
    if(!NuclearMass::QValue(Z+2,A+4,Z,A,qValue)) {
      isValid_=false;
      return;
    }
  
    compoundZ_ = Z+2;
    compoundA_ = A+4;
    double reducedMass=double(A)*4./double(A+4);
    double spinNorm=(2.*groundStateJ_+1.);
    preFactor_ = hbarc*hbarc/200.*pi/uconv/reducedMass/spinNorm;
  }
  
  seperationEnergy_ = -qValue+knownLevels[entranceState].energy_;

  specifiedExitSepE_ = std::vector<double>(4,0.);
  specifiedExitJ_ = std::vector<double>(4,-1.);
  specifiedExitPi_ = std::vector<int>(4,0);
  
  for(int i =0;i<exitStates.size();i++) {
    if(exitStates[i]>=0) {
      std::vector<Level> levels;  
      double qVal;
      
      if(i==0) {
	      levels = NuclearLevels::FindLevels(compoundZ_,compoundA_);
	      qVal = 0.;
      } else if(i==1) {
	      levels = NuclearLevels::FindLevels(compoundZ_,compoundA_-1);
	      if(!NuclearMass::QValue(compoundZ_,compoundA_,compoundZ_,compoundA_-1,qVal)) {
	        exitStates_[i]=-1;
	        continue;
	      }
      } else if(i==2) {
	        levels = NuclearLevels::FindLevels(compoundZ_-1,compoundA_-1);
	        
          if(!NuclearMass::QValue(compoundZ_,compoundA_,compoundZ_-1,compoundA_-1,qVal)) {
	          exitStates_[i]=-1;
	          continue;
	        }
      } else if(i==3) {
	      levels = NuclearLevels::FindLevels(compoundZ_-2,compoundA_-4);
	        if(!NuclearMass::QValue(compoundZ_,compoundA_,compoundZ_-2,compoundA_-4,qVal)) {
	          exitStates_[i]=-1;
	          continue;
	        }
      }
      
      if(levels.size()<=exitStates[i]) {
	      if(i==0) std::cout << "Cannot find specified gamma residual state.  Calculating total." << std::endl;
	      else if(i==1) std::cout << "Cannot find specified neutron residual state.  Calculating total." << std::endl;
	      else if(i==2) std::cout << "Cannot find specified proton residual state.  Calculating total." << std::endl;
	      else if(i==3) std::cout << "Cannot find specified alpha residual state.  Calculating total." << std::endl;
	      exitStates_[i]=-1;
      } else {
	      specifiedExitJ_[i] = levels[exitStates[i]].J_;
	      specifiedExitPi_[i] = levels[exitStates[i]].Pi_;
	      specifiedExitSepE_[i] = -qVal+levels[exitStates[i]].energy_;
	      std::string parity = (levels[exitStates[i]].Pi_<0) ? "-" : "+";
	      
        if(i==0) std::cout << "Gamma residual state: E=" << levels[exitStates[i]].energy_ 
			                    << " MeV J=" << levels[exitStates[i]].J_ << parity << std::endl;
	      else if(i==1) std::cout << "Neutron residual state: E=" << levels[exitStates[i]].energy_ 
			                    << " MeV J=" << levels[exitStates[i]].J_ << parity << std::endl;
	      else if(i==2) std::cout << "Proton residual state: E=" << levels[exitStates[i]].energy_ 
			                    << " MeV J=" << levels[exitStates[i]].J_ << parity << std::endl;
	      else if(i==3) std::cout << "Alpha residual state: E=" << levels[exitStates[i]].energy_ 
			                    << " MeV J=" << levels[exitStates[i]].J_ << parity << std::endl;
      }
    }
  }

  if(forRates) {
    std::cout << "Calculating Partition Function..." << std::endl;
    CalcPartitionFunc();
  }
  
  if(!CalcAllowedJPi(forRates)) isValid_=false;
  else isValid_=true;

  energiesGiven_=false;
  if(energyFile.length()==0 || !FillEnergies(energyFile)) {
    if(pType_!=0) CalculateEnergyGrid();
    else {
      for(double E = minEnergy_+seperationEnergy_; E<=maxEnergy_+seperationEnergy_; E+=dE_) {
	      crossSections_.push_back(std::pair<double,CrossSectionValues>(E-seperationEnergy_,
								      CrossSectionValues(0.,0.,0.,0.,0.,0.,0.,0.)));
      }
    }
  }
}

bool CrossSection::FillEnergies(std::string energyFile) {
  std::ifstream in(energyFile.c_str());
  if(!in) return false;
  while(!in.eof()) {
    std::string line;
    std::getline(in,line);
    if(!in.eof()) {
      std::istringstream stm(line);
      double energy;
      if(stm>>energy) {
	crossSections_.push_back(std::pair<double,CrossSectionValues>(energy,
								      CrossSectionValues(0.,0.,0.,0.,0.,0.,0.,0.)));
      }
    }
  }
  in.close();
  energiesGiven_ = crossSections_.size()>0;
  return energiesGiven_;
}

bool CrossSection::CalcAllowedJPi(bool forRates) {
  if(forRates) {
    if((pType_==1 || pType_==2)&&A_%2==0) {
      for(int i =0;i<10;i++) {
	allowedJPi_.push_back(std::pair<double,int>(0.5+i,-1));
	allowedJPi_.push_back(std::pair<double,int>(0.5+i,1));
      }
    } else {
      for(int i =0;i<10;i++) {
	allowedJPi_.push_back(std::pair<double,int>(i,-1));
	allowedJPi_.push_back(std::pair<double,int>(i,1));
      }
    }
    return true;
  }
  if(pType_ == 1 || pType_ ==2) {
    for(double s = fabs(groundStateJ_-0.5);s<=groundStateJ_+0.5;s+=1.) {
      for(double l = 0.; l<=Decayer::maxL_;l+=1.) {
	int piCompound = (int(l)%2==0) ? groundStatePi_ : -1*groundStatePi_;
	for(double jCompound = fabs(s-l);jCompound<=s+l;jCompound+=1.) {
	  bool found = false;
	  for(int i = 0;i<allowedJPi_.size();i++) {
	    if(allowedJPi_[i].first==jCompound&&
	       allowedJPi_[i].second==piCompound) {
	      found=true;
	      break;
	    }
	  }
	  if(!found) 	      
	    allowedJPi_.push_back(std::pair<double,int>(jCompound,piCompound));
	}
      }
    }
  } else if(pType_==3) {
    for(double l = 0.;l<=Decayer::maxL_;l+=1.) {
      int piCompound = (int(l)%2==0) ? groundStatePi_ : -1*groundStatePi_;
      for(double jCompound = fabs(groundStateJ_-l);jCompound<=groundStateJ_+l;jCompound+=1.) {
	bool found = false;
	for(int i = 0;i<allowedJPi_.size();i++) {
	  if(allowedJPi_[i].first==jCompound&&
	     allowedJPi_[i].second==piCompound) {
	    found=true;
	    break;
	  }
	}
	if(!found) 	      
	  allowedJPi_.push_back(std::pair<double,int>(jCompound,piCompound));
      }
    }
  } else {
    for(double jCompound = fabs(groundStateJ_-1.); jCompound<=groundStateJ_+1.;jCompound+=1.) {
      bool found = false;
      for(int i = 0;i<allowedJPi_.size();i++) {
	if(allowedJPi_[i].first==jCompound&&
	   allowedJPi_[i].second==-1*groundStatePi_) {
	  found = true;
	  break;
	}
      }
      if(!found) 
	allowedJPi_.push_back(std::pair<double,int>(jCompound,-1*groundStatePi_));
      found = false;
      for(int i = 0;i<allowedJPi_.size();i++) {
	if(allowedJPi_[i].first==jCompound&&
	   allowedJPi_[i].second==groundStatePi_) {
	  found=true;
	  break;
	}
      }
      if(!found) 
      	allowedJPi_.push_back(std::pair<double,int>(jCompound,groundStatePi_));
    }
    for(double jCompound = fabs(groundStateJ_-2.); jCompound<=groundStateJ_+2.;jCompound+=1.) {
      bool found = false;
      for(int i = 0;i<allowedJPi_.size();i++) {
	if(allowedJPi_[i].first==jCompound&&
	   allowedJPi_[i].second==groundStatePi_) {
	  found=true;
	  break;
	}
      }
      if(!found) 
       	allowedJPi_.push_back(std::pair<double,int>(jCompound,groundStatePi_));
    }
  }
  for(int i = 0;i<allowedJPi_.size();i++) {
    entranceTrans_[i]=std::vector<double>();
    gExitTrans_[i]=std::vector<double>();
    nExitTrans_[i]=std::vector<double>();
    pExitTrans_[i]=std::vector<double>();
    aExitTrans_[i]=std::vector<double>();
  }
  if(allowedJPi_.size()==0) {
    return false;
  } else return true;
}

bool CrossSection::CalcDecayerVector(double E, DecayerVector& decayerVector, bool forAverageWidth) {
  for(int i = 0;i<allowedJPi_.size();i++) {
    Decayer* newDecayer = new Decayer(compoundZ_,compoundA_,allowedJPi_[i].first,
				      allowedJPi_[i].second,E);
    newDecayer->CorrectWidthFluctuations();
    std::vector<SpinRatePair*> entrancePairs;
    for(int i = 0;i<newDecayer->widthCorrectedDecayer_->spinRatePairs_.size();
	i++) {
      if(newDecayer->widthCorrectedDecayer_->spinRatePairs_[i].Z_==Z_&&
	 newDecayer->widthCorrectedDecayer_->spinRatePairs_[i].A_==A_&&
	 newDecayer->widthCorrectedDecayer_->spinRatePairs_[i].spin_==groundStateJ_&&
	 newDecayer->widthCorrectedDecayer_->spinRatePairs_[i].parity_==groundStatePi_)  {
	entrancePairs.push_back(&newDecayer->widthCorrectedDecayer_->spinRatePairs_[i]);
      }
    }
    /*
    if(entrancePairs.size()==0&&!forAverageWidth) {
      delete newDecayer;
      return false;
    }
    */
    decayerVector.push_back(std::pair<Decayer*,std::vector<SpinRatePair*> >(newDecayer,entrancePairs));
  }
  return true;
}

void CrossSection::Calculate() {
  if(calculateGammaCutoff_) {
    TransitionRateFunc::SetGammaCutoffEnergy(10000.);
    gammaCutoffSet_ = false;
  } else {
    gammaCutoffSet_ = true;
  }

  std::cout <<  "\t [                         ] 0%";std::cout.flush();
  
  int pointIndex=0;
  int numPoints=crossSections_.size();
  time_t startTime = time(NULL);
  int skipCounter = 0;


  for(int i=0;i<numPoints;i++) {
    if(i>0&&!energiesGiven_) {
      skipCounter++;
      
      if((crossSections_[i-1].second.neutron_==0.|| crossSections_[i-1].second.alpha_==0.|| crossSections_[i-1].second.proton_==0.)&& !skipped_[i-1]) skipCounter=0;
      
      if(crossSections_[i].first>skipEnergy_&&i!=numPoints-1 && (i+1)%3!=0 && skipCounter>8) {
        skipped_.push_back(true);
        continue;
      }
    }

    skipped_.push_back(false);
    ++pointIndex;
    
    if(difftime(time(NULL),startTime)>0.25) {
      startTime=time(NULL);
      std::string progress=" [";
      double percent=0.;
      
      for(int j = 1;j<=25;j++) {
	      if(pointIndex>=percent*numPoints&&percent<1.) {
	        percent+=0.04;
	        progress+='*';
	      } else progress+=' ';
      } 
      
      progress+="] ";
      
      std::cout << "\r\t" << progress << percent*100 << "%"; std::cout.flush();
    }
    
    double E = crossSections_[i].first+seperationEnergy_;
    double geometricCrossSection = (pType_!=0) ? preFactor_/(E-seperationEnergy_) : preFactor_/E/E;
    
    DecayerVector decayerVector;
    
    if(!CalcDecayerVector(E,decayerVector)) {
      for(int j = 0;j<decayerVector.size();j++) 
	      delete decayerVector[j].first;
      continue;
    }

    double gammaSum = 0.;
    double neutronSum = 0.;
    double protonSum = 0.;
    double alphaSum = 0.;
    double gammaStellarSum = 0.;
    double neutronStellarSum = 0.;
    double protonStellarSum = 0.;
    double alphaStellarSum = 0.;
    
    for(int j = 0;j<decayerVector.size();j++) {
      std::vector<SpinRatePair*> entrancePairs = decayerVector[j].second;
      Decayer* decayer = decayerVector[j].first->widthCorrectedDecayer_;
      double entranceTransmission=0.;
      
      for(int k = 0;k<entrancePairs.size();k++) {
	      entranceTransmission += entrancePairs[k]->rateFunc_->CalcTransmissionFunc(E-seperationEnergy_);
      }

      double gammaExitTransmission=0.;
      double neutronExitTransmission=0.;
      double protonExitTransmission=0.;
      double alphaExitTransmission=0.;

      for(int k = 0;k<decayer->spinRatePairs_.size();k++) {
	      if(decayer->spinRatePairs_[k].A_==compoundA_&&
	        decayer->spinRatePairs_[k].Z_==compoundZ_) {
	  
          if(exitStates_[0]<0) {
	          double branching = 
	            (residualGamma_) ? decayer->spinRatePairs_[k].rateFunc_->ExclusiveBranching() : 1.;
	    gammaExitTransmission+=decayer->spinRatePairs_[k].integral_*branching;
	  } else {
	    if(decayer->spinRatePairs_[k].spin_==specifiedExitJ_[0]&&
	       decayer->spinRatePairs_[k].parity_==specifiedExitPi_[0])
	      if(E-specifiedExitSepE_[0]>0.)
		gammaExitTransmission+=decayer->spinRatePairs_[k].rateFunc_->CalcTransmissionFunc(E-specifiedExitSepE_[0]);
	  }
	} else if(decayer->spinRatePairs_[k].A_==compoundA_-1&&
		  decayer->spinRatePairs_[k].Z_==compoundZ_) {
	  if(exitStates_[1]<0) {
	    double branching = 
	      (residualNeutron_) ? decayer->spinRatePairs_[k].rateFunc_->ExclusiveBranching() : 1.;
	    neutronExitTransmission+=decayer->spinRatePairs_[k].integral_*branching;
	  } else {
	    if(decayer->spinRatePairs_[k].spin_==specifiedExitJ_[1]&&
	       decayer->spinRatePairs_[k].parity_==specifiedExitPi_[1])
	      if(E-specifiedExitSepE_[1]>0.)
		neutronExitTransmission+=decayer->spinRatePairs_[k].rateFunc_->CalcTransmissionFunc(E-specifiedExitSepE_[1]);
	  }
	} else if(decayer->spinRatePairs_[k].A_==compoundA_-1&&
		  decayer->spinRatePairs_[k].Z_==compoundZ_-1) {
	  if(exitStates_[2]<0) {
	    double branching = 
	      (residualProton_) ? decayer->spinRatePairs_[k].rateFunc_->ExclusiveBranching() : 1.;
	    protonExitTransmission+=decayer->spinRatePairs_[k].integral_*branching;	
	  } else {
	    if(decayer->spinRatePairs_[k].spin_==specifiedExitJ_[2]&&
	       decayer->spinRatePairs_[k].parity_==specifiedExitPi_[2])
	      if(E-specifiedExitSepE_[2]>0.)
		protonExitTransmission+=decayer->spinRatePairs_[k].rateFunc_->CalcTransmissionFunc(E-specifiedExitSepE_[2]);
	  }
	} else if(decayer->spinRatePairs_[k].A_==compoundA_-4&&
		  decayer->spinRatePairs_[k].Z_==compoundZ_-2) {
	  if(exitStates_[3]<0) {
	    double branching = 
	      (residualAlpha_) ? decayer->spinRatePairs_[k].rateFunc_->ExclusiveBranching() : 1.;
	    alphaExitTransmission+=decayer->spinRatePairs_[k].integral_*branching;
	  } else {
	    if(decayer->spinRatePairs_[k].spin_==specifiedExitJ_[3]&&
	       decayer->spinRatePairs_[k].parity_==specifiedExitPi_[3])
	      if(E-specifiedExitSepE_[3]>0.)
		alphaExitTransmission+=decayer->spinRatePairs_[k].rateFunc_->CalcTransmissionFunc(E-specifiedExitSepE_[3]);
	  }	
	}
      }
      double compoundTransmission;
      if(pType_==0) {
	compoundTransmission=gammaExitTransmission;
	if(exitStates_[0]<0||entranceState_==exitStates_[0]) 
	  gammaExitTransmission-=entranceTransmission;
      } else if(pType_==1) {
	compoundTransmission=neutronExitTransmission;
	if(exitStates_[1]<0||entranceState_==exitStates_[1]) 
	  neutronExitTransmission-=entranceTransmission;
      } else if(pType_==2) {
	compoundTransmission=protonExitTransmission;
	if(exitStates_[2]<0||entranceState_==exitStates_[2]) 
	  protonExitTransmission-=entranceTransmission;
      } else if(pType_==3) {
	compoundTransmission=alphaExitTransmission;
	if(exitStates_[3]<0||entranceState_==exitStates_[3]) 
	  alphaExitTransmission-=entranceTransmission;
      }
      entranceTrans_[j].push_back(entranceTransmission);
      nExitTrans_[j].push_back(neutronExitTransmission);
      gExitTrans_[j].push_back(gammaExitTransmission);
      pExitTrans_[j].push_back(protonExitTransmission);
      aExitTrans_[j].push_back(alphaExitTransmission);
      if(decayer->totalIntegral_>0.) {
	gammaSum+=(2.*decayer->jInitial_+1.)*entranceTransmission*gammaExitTransmission/
	  decayer->totalIntegral_;
	neutronSum+=(2.*decayer->jInitial_+1.)*entranceTransmission*neutronExitTransmission/
	  decayer->totalIntegral_;
	protonSum+=(2.*decayer->jInitial_+1.)*entranceTransmission*protonExitTransmission/
	  decayer->totalIntegral_;
	alphaSum+=(2.*decayer->jInitial_+1.)*entranceTransmission*alphaExitTransmission/
	  decayer->totalIntegral_;
	gammaStellarSum+=(2.*decayer->jInitial_+1.)*compoundTransmission*gammaExitTransmission/
	  decayer->totalIntegral_;
	neutronStellarSum+=(2.*decayer->jInitial_+1.)*compoundTransmission*neutronExitTransmission/
	  decayer->totalIntegral_;
	protonStellarSum+=(2.*decayer->jInitial_+1.)*compoundTransmission*protonExitTransmission/
	  decayer->totalIntegral_;
	alphaStellarSum+=(2.*decayer->jInitial_+1.)*compoundTransmission*alphaExitTransmission/
	  decayer->totalIntegral_;
      }
      delete decayerVector[j].first;
    }
    double totalSum = gammaSum + neutronSum + protonSum + alphaSum;
    if(totalSum>=2.*gammaSum&&!gammaCutoffSet_&&i>0) {
      TransitionRateFunc::SetGammaCutoffEnergy(crossSections_[i-1].first+seperationEnergy_);
      gammaCutoffSet_=true;
    }
    crossSections_[i].second.gamma_=gammaSum*geometricCrossSection;
    crossSections_[i].second.proton_=protonSum*geometricCrossSection;
    crossSections_[i].second.neutron_=neutronSum*geometricCrossSection;
    crossSections_[i].second.alpha_=alphaSum*geometricCrossSection;
    crossSections_[i].second.gammaStellar_=gammaStellarSum*geometricCrossSection;
    crossSections_[i].second.protonStellar_=protonStellarSum*geometricCrossSection;
    crossSections_[i].second.neutronStellar_=neutronStellarSum*geometricCrossSection;
    crossSections_[i].second.alphaStellar_=alphaStellarSum*geometricCrossSection;
  }
  std::cout  << "\r\t" << " [*************************] 100%" << std::endl;
}

void CrossSection::PrintCrossSections() {
  char filename[256];
  if(pType_==0) {
    sprintf(filename,"Sapphire_%d%s+g.dat",
    	A_,
    	NuclearMass::FindElement(Z_).c_str());  
  } else if(pType_==1) {
    sprintf(filename,"Sapphire_%d%s+n.dat",
    	A_,
    	NuclearMass::FindElement(Z_).c_str());  
  } else if(pType_==2) {
    sprintf(filename,"Sapphire_%d%s+p.dat",
    	A_,
    	NuclearMass::FindElement(Z_).c_str());  
  } else {
    sprintf(filename,"Sapphire_%d%s+a.dat",
      	    A_,
    	    NuclearMass::FindElement(Z_).c_str());  
  }
  std::ofstream out(filename);
  out << "#" << std::setw(14) << "Energy [MeV]"
      << std::setw(15) << "gamma [b]"
      << std::setw(15) << "neutron [b]"
      << std::setw(15) << "proton [b]"
      << std::setw(15) << "alpha [b]"
      << std::setw(0)  << std::endl; 
  for(int i = 0;i<crossSections_.size();i++) {
    if(skipped_[i]) continue;
    out << std::scientific
	<< std::setw(15) << crossSections_[i].first
	<< std::setw(15) << crossSections_[i].second.gamma_
	<< std::setw(15) << crossSections_[i].second.neutron_
        << std::setw(15) << crossSections_[i].second.proton_ 
        << std::setw(15) << crossSections_[i].second.alpha_ 
        << std::setw(0)  << std::endl; 
  }
  out.flush();
  out.close();
}

void CrossSection::PrintTransmissionTerms() {
  for(int j = 0;j<allowedJPi_.size();j++) {
    char filename[256];
    if(allowedJPi_[j].second>0)
      sprintf(filename,"Transmission_%d%s_J=%.1f+.dat",
	      compoundA_,
	      NuclearMass::FindElement(compoundZ_).c_str(),
	      allowedJPi_[j].first);  
    else
      sprintf(filename,"Transmission_%d%s_J=%.1f-.dat",
	      compoundA_,
	      NuclearMass::FindElement(compoundZ_).c_str(),
	      allowedJPi_[j].first);  
    std::ofstream out(filename);
    out << "#" << std::setw(14) << "Energy [MeV]"
	<< std::setw(15) << "entrance"
	<< std::setw(15) << "gamma"
	<< std::setw(15) << "neutron"
	<< std::setw(15) << "proton"
	<< std::setw(15) << "alpha"
	<< std::setw(0)  << std::endl; 
    int k = 0;
    for(int i = 0;i<crossSections_.size();i++) {
      if(skipped_[i]) continue;
      out << std::scientific
	  << std::setw(15) << crossSections_[i].first
	  << std::setw(15) << entranceTrans_[j][k]
	  << std::setw(15) << gExitTrans_[j][k]
	  << std::setw(15) << nExitTrans_[j][k]
	  << std::setw(15) << pExitTrans_[j][k]
	  << std::setw(15) << aExitTrans_[j][k]
	  << std::setw(0)  << std::endl; 
      k++;
    }
    out.flush();
    out.close();
  }
}

std::pair<double,double> CrossSection::CalcAverageSWaveResWidth() {
  double energy = seperationEnergy_;
  DecayerVector decayerVector;
  if(!CalcDecayerVector(energy,decayerVector,true)) {
    for(int j = 0;j<decayerVector.size();j++) 
      delete decayerVector[j].first;
    return std::pair<double,double>(0.,0.);
  }
  double upSum=0.;
  double downSum=0.;
  for(int j = 0;j<decayerVector.size();j++) {
    Decayer* decayer = decayerVector[j].first->widthCorrectedDecayer_;
    if(decayer->jInitial_!=groundStateJ_+0.5&&
       decayer->jInitial_!=groundStateJ_-0.5) continue;
    if(decayer->piInitial_!=groundStatePi_) continue;
    for(int k=0;k<decayer->spinRatePairs_.size();k++) {
      if(decayer->spinRatePairs_[k].A_!=compoundA_||
	 decayer->spinRatePairs_[k].Z_!=compoundZ_) continue;
      if(decayer->jInitial_==groundStateJ_+0.5)
	upSum+=decayer->spinRatePairs_[k].integral_;
      else if(decayer->jInitial_==groundStateJ_-0.5)
	downSum+=decayer->spinRatePairs_[k].integral_;
    }
    delete decayerVector[j].first;
  }
  LevelDensity* levelDensityUp = new RauscherLevelDensity(compoundZ_,compoundA_,groundStateJ_+0.5);
  LevelDensity* levelDensityDown = new RauscherLevelDensity(compoundZ_,compoundA_,groundStateJ_-0.5);
  double ldValueUp=0.;
  double ldValueDown=0.;  
  if(levelDensityUp->operator()(energy)>0.) {
    ldValueUp=levelDensityUp->operator()(energy);
    upSum/=2.*pi*ldValueUp;
  }
  if(levelDensityDown->operator()(energy)>0.) {
    ldValueDown=levelDensityDown->operator()(energy);
    downSum/=2.*pi*ldValueDown;
  }
  delete levelDensityUp;
  delete levelDensityDown;
  double sum = upSum*(2.*groundStateJ_+2.);
  double totalWeight = (2.*groundStateJ_+2.);
  if(groundStateJ_>=0.5) {
    totalWeight+=2.*groundStateJ_;
    sum+=downSum*2.*groundStateJ_;
  }
  return std::pair<double,double>(sum/totalWeight,
				  1./(ldValueUp+ldValueDown));
}

std::pair<double,double> CrossSection::CalcAveragePWaveResWidth() {
  double energy = seperationEnergy_;
  DecayerVector decayerVector;
  if(!CalcDecayerVector(energy,decayerVector,true)) {
    for(int j = 0;j<decayerVector.size();j++) 
      delete decayerVector[j].first;
    return std::pair<double,double>(0.,0.);
  }
  double upSum15=0.;
  double upSum05=0.;
  double downSum05=0.;
  double downSum15=0.;
  for(int j = 0;j<decayerVector.size();j++) {
    Decayer* decayer = decayerVector[j].first->widthCorrectedDecayer_;
    if(decayer->jInitial_!=groundStateJ_+1.5&&
       decayer->jInitial_!=groundStateJ_+0.5&&
       decayer->jInitial_!=groundStateJ_-0.5&&
       decayer->jInitial_!=groundStateJ_-1.5) continue;
    if(decayer->piInitial_==groundStatePi_) continue;
    for(int k=0;k<decayer->spinRatePairs_.size();k++) {
      if(decayer->spinRatePairs_[k].A_!=compoundA_||
	 decayer->spinRatePairs_[k].Z_!=compoundZ_) continue;
      if(decayer->jInitial_==groundStateJ_+1.5)
	upSum15+=decayer->spinRatePairs_[k].integral_;
      else if(decayer->jInitial_==groundStateJ_+0.5)
	upSum05+=decayer->spinRatePairs_[k].integral_;
      else if(decayer->jInitial_==groundStateJ_-0.5)
	downSum05+=decayer->spinRatePairs_[k].integral_;
      else if(decayer->jInitial_==groundStateJ_-1.5)
	downSum15+=decayer->spinRatePairs_[k].integral_;
    }
    delete decayerVector[j].first;
  }
  LevelDensity* levelDensityUp15 = new RauscherLevelDensity(compoundZ_,compoundA_,groundStateJ_+1.5);
  LevelDensity* levelDensityUp05 = new RauscherLevelDensity(compoundZ_,compoundA_,groundStateJ_+0.5);
  LevelDensity* levelDensityDown05 = new RauscherLevelDensity(compoundZ_,compoundA_,groundStateJ_-0.5);
  LevelDensity* levelDensityDown15 = new RauscherLevelDensity(compoundZ_,compoundA_,groundStateJ_-1.5);
  double ldValueUp15 = 0.;
  double ldValueUp05 = 0.;
  double ldValueDown05 = 0.;
  double ldValueDown15 = 0.;
  if(levelDensityUp15->operator()(energy)>0.) {
    ldValueUp15=levelDensityUp15->operator()(energy);
    upSum15/=2.*pi*ldValueUp15;
  }
  if(levelDensityUp05->operator()(energy)>0.) {
    ldValueUp05=levelDensityUp05->operator()(energy);
    upSum05/=2.*pi*ldValueUp05;
  }
  if(levelDensityDown05->operator()(energy)>0.) {
    ldValueDown05=levelDensityDown05->operator()(energy);
    downSum05/=2.*pi*ldValueDown05;
  }
  if(levelDensityDown15->operator()(energy)>0.) {
    ldValueDown15=levelDensityDown15->operator()(energy);
    downSum15/=2.*pi*ldValueDown15;
  }
  delete levelDensityUp15;
  delete levelDensityUp05;
  delete levelDensityDown05;
  delete levelDensityDown15;
  double sum = upSum15*(2.*groundStateJ_+4.)+upSum05*(2.*groundStateJ_+2.);
  double totalWeight = 4.*groundStateJ_+6.;
  if(groundStateJ_>=0.5) {
    totalWeight+= 2.*groundStateJ_;
    sum+=downSum05*(2.*groundStateJ_);
  }
  if(groundStateJ_>=1.5) {
    totalWeight+= 2.*groundStateJ_-2.;
     sum+=downSum15*(2.*groundStateJ_-2.);
  }
  
  return std::pair<double,double>(sum/totalWeight,
				  1./(ldValueUp15+ldValueUp05+ldValueDown05+ldValueDown15));
}

std::pair<double,double> CrossSection::CalcAverageDWaveResWidth() {
  double energy = seperationEnergy_;
  DecayerVector decayerVector;
  if(!CalcDecayerVector(energy,decayerVector,true)) {
    for(int j = 0;j<decayerVector.size();j++) 
      delete decayerVector[j].first;
    return std::pair<double,double>(0.,0.);
  }
  double upSum25=0.;
  double upSum15=0.;
  double upSum05=0.;
  double downSum05=0.;
  double downSum15=0.;
  double downSum25=0.;
  for(int j = 0;j<decayerVector.size();j++) {
    Decayer* decayer = decayerVector[j].first->widthCorrectedDecayer_;
    if(decayer->jInitial_!=groundStateJ_+2.5&&
       decayer->jInitial_!=groundStateJ_+1.5&&
       decayer->jInitial_!=groundStateJ_+0.5&&
       decayer->jInitial_!=groundStateJ_-0.5&&
       decayer->jInitial_!=groundStateJ_-1.5&&
       decayer->jInitial_!=groundStateJ_-2.5) continue;
    if(decayer->piInitial_!=groundStatePi_) continue;
    for(int k=0;k<decayer->spinRatePairs_.size();k++) {
      if(decayer->spinRatePairs_[k].A_!=compoundA_||
	 decayer->spinRatePairs_[k].Z_!=compoundZ_) continue;
      if(decayer->jInitial_==groundStateJ_+2.5)
	upSum25+=decayer->spinRatePairs_[k].integral_;
      else if(decayer->jInitial_==groundStateJ_+1.5)
	upSum15+=decayer->spinRatePairs_[k].integral_;
      else if(decayer->jInitial_==groundStateJ_+0.5)
	upSum05+=decayer->spinRatePairs_[k].integral_;
      else if(decayer->jInitial_==groundStateJ_-0.5)
	downSum05+=decayer->spinRatePairs_[k].integral_;
      else if(decayer->jInitial_==groundStateJ_-1.5)
	downSum15+=decayer->spinRatePairs_[k].integral_;
      else if(decayer->jInitial_==groundStateJ_-2.5)
	downSum25+=decayer->spinRatePairs_[k].integral_;
    }
    delete decayerVector[j].first;
  }
  LevelDensity* levelDensityUp25 = new RauscherLevelDensity(compoundZ_,compoundA_,groundStateJ_+2.5);
  LevelDensity* levelDensityUp15 = new RauscherLevelDensity(compoundZ_,compoundA_,groundStateJ_+1.5);
  LevelDensity* levelDensityUp05 = new RauscherLevelDensity(compoundZ_,compoundA_,groundStateJ_+0.5);
  LevelDensity* levelDensityDown05 = new RauscherLevelDensity(compoundZ_,compoundA_,groundStateJ_-0.5);
  LevelDensity* levelDensityDown15 = new RauscherLevelDensity(compoundZ_,compoundA_,groundStateJ_-1.5);
  LevelDensity* levelDensityDown25 = new RauscherLevelDensity(compoundZ_,compoundA_,groundStateJ_-2.5);
  double ldValueUp25 = 0.;
  double ldValueUp15 = 0.;
  double ldValueUp05 = 0.;
  double ldValueDown05 = 0.;
  double ldValueDown15 = 0.;
  double ldValueDown25 = 0.;
  if(levelDensityUp25->operator()(energy)>0.) {
    ldValueUp25=levelDensityUp25->operator()(energy);
    upSum25/=2.*pi*ldValueUp25;
  }
  if(levelDensityUp15->operator()(energy)>0.) {
    ldValueUp15=levelDensityUp15->operator()(energy);
    upSum15/=2.*pi*ldValueUp15;
  }
  if(levelDensityUp05->operator()(energy)>0.) {
    ldValueUp05=levelDensityUp05->operator()(energy);
    upSum05/=2.*pi*ldValueUp05;
  }
  if(levelDensityDown05->operator()(energy)>0.) {
    ldValueDown05=levelDensityDown05->operator()(energy);
    downSum05/=2.*pi*ldValueDown05;
  }
  if(levelDensityDown15->operator()(energy)>0.) {
    ldValueDown15=levelDensityDown15->operator()(energy);
    downSum15/=2.*pi*ldValueDown15;
  }
  if(levelDensityDown25->operator()(energy)>0.) {
    ldValueDown25=levelDensityDown25->operator()(energy);
    downSum25/=2.*pi*ldValueDown25;
  }
  delete levelDensityUp25;
  delete levelDensityUp15;
  delete levelDensityUp05;
  delete levelDensityDown05;
  delete levelDensityDown15;
  delete levelDensityDown25;
  double totalWeight = 4.*groundStateJ_+10.;
  double sum = upSum25*(2.*groundStateJ_+6.)+upSum15*(2.*groundStateJ_+4.);//+upSum05*(2.*groundStateJ_+2.);
  
  if(groundStateJ_>=1.5) {
    totalWeight+= 2.*groundStateJ_-2.;
     sum+=downSum15*(2.*groundStateJ_-2.);
  }
  if(groundStateJ_>=2.5) {
    totalWeight+= 2.*groundStateJ_-4.;
     sum+=downSum25*(2.*groundStateJ_-4.);
  }
  return std::pair<double,double>(sum/totalWeight,
				  1./(ldValueUp25+ldValueUp15+downSum15+downSum25));
}

/* 
 * Creates SMOKER-like energy grid. A bit crazy.
 */
void CrossSection::CalculateEnergyGrid() {
  
  double neutronSepE;
  double protonSepE;
  double alphaSepE;

  if(!NuclearMass::QValue(compoundZ_,compoundA_,
			  compoundZ_,compoundA_-1,neutronSepE)) {
    isValid_=false;
    return;
  } else neutronSepE *= -1.;
  if(!NuclearMass::QValue(compoundZ_,compoundA_,
			  compoundZ_-1,compoundA_-1,protonSepE)) {
    isValid_=false;
    return;
  } else protonSepE *= -1.;
  if(!NuclearMass::QValue(compoundZ_,compoundA_,
			  compoundZ_-2,compoundA_-4,alphaSepE)) {
    isValid_=false;
    return;
  } else alphaSepE *= -1.;


  //Determine energy window for neutrons
  double minNeutronEnergy = 0.086*(0.5)*rateTemps_[0]-
    0.6*0.097*sqrt(0.5)*rateTemps_[0]+
    neutronSepE;
  double maxNeutronEnergy = 
    0.086*(Decayer::maxL_+0.5)*rateTemps_[rateTemps_.size()-1]+
    2.*0.097*sqrt(Decayer::maxL_+0.5)*rateTemps_[rateTemps_.size()-1]+
    neutronSepE;
  
  //Determine max energy window for protons and alphas

  double maxProtonEnergy=0.122*pow((compoundA_-1.)/compoundA_,0.33333333333333)*
    pow((compoundZ_-1.)*rateTemps_[rateTemps_.size()-1],0.6666666666667)+
    2.*0.237*pow((compoundA_-1.)/compoundA_,0.166666666667)*
    pow((compoundZ_-1.),0.333333333333)*pow(rateTemps_[rateTemps_.size()-1],0.83333333333)+
    protonSepE;
  
  double maxAlphaEnergy=0.122*pow(4.*(compoundA_-4.)/compoundA_,0.33333333333333)*
    pow(2.*(compoundZ_-2.)*rateTemps_[rateTemps_.size()-1],0.6666666666667)+
    2.*0.237*pow(4.*(compoundA_-4.)/compoundA_,0.166666666667)*
    pow(2.*(compoundZ_-2.),0.333333333333)*pow(rateTemps_[rateTemps_.size()-1],0.83333333333)+
    alphaSepE;

  //Find minimum energy for protons
  double thresholdProton = (pType_==2) ? 0. : protonSepE - seperationEnergy_;
  double minProtonEnergy;
  if(thresholdProton<0.) minProtonEnergy=0.;
  else {
    minProtonEnergy=0.122*pow((compoundA_-1.)/compoundA_,0.33333333333333)*
      pow((compoundZ_-1.)*rateTemps_[0],0.6666666666667)+thresholdProton;
 
    double gamowEnergyIn;
    if(pType_==1) gamowEnergyIn=0.;
    else if(pType_==2) gamowEnergyIn=6.283186/137.0*(compoundZ_-1.)*
			 sqrt((compoundA_-1.)/compoundA_*465.739);
    else if(pType_==3) gamowEnergyIn=6.283186/137.0*2.*(compoundZ_-2.)*
			 sqrt(4*(compoundA_-4.)/compoundA_*465.739);
    double gamowEnergyOut = 
      (pType_==2) ? 0. : 6.283186/137.0*(compoundZ_-1.)*sqrt((compoundA_-1.)/compoundA_*465.739);

    bool continueLoop = true;
    while(continueLoop) {
      double newMinProtonEnergy = minProtonEnergy
	-(-11.605/rateTemps_[0]+0.5*gamowEnergyIn*pow(minProtonEnergy,-1.5)
	  +0.5*gamowEnergyOut*pow(minProtonEnergy-thresholdProton,-1.5))/
	(-0.75*gamowEnergyIn*pow(minProtonEnergy,-2.5)-
	 0.75*gamowEnergyOut*pow(minProtonEnergy-thresholdProton,-2.5));
      if(fabs((minProtonEnergy-newMinProtonEnergy)/newMinProtonEnergy)>=1.e-5) {
	if(newMinProtonEnergy<=thresholdProton)
	  newMinProtonEnergy=thresholdProton+(minProtonEnergy-thresholdProton)*0.8;
      } else continueLoop=false;
      minProtonEnergy = newMinProtonEnergy;
    }
 
    double firstIterMinProtonEnergy = minProtonEnergy;
    minProtonEnergy-=0.2*0.237*pow((compoundA_-1.)/compoundA_,0.166666666667)*
      pow((compoundZ_-1.),0.333333333333)*pow(rateTemps_[0],0.83333333333);
    continueLoop = true;
    while(continueLoop) {
      double newMinProtonEnergy = minProtonEnergy
	-(log(2.)+firstIterMinProtonEnergy*11.605/rateTemps_[0]+
	  gamowEnergyIn*pow(firstIterMinProtonEnergy,-0.5)
	  +gamowEnergyOut*pow(firstIterMinProtonEnergy-thresholdProton,-0.5)
	  -11.605*minProtonEnergy/rateTemps_[0]-gamowEnergyIn*pow(minProtonEnergy,-0.5)
	  -gamowEnergyOut*pow(minProtonEnergy-thresholdProton,-0.5))/
	(-11.605/rateTemps_[0]+0.5*gamowEnergyIn*pow(minProtonEnergy,-1.5)
	 +0.5*gamowEnergyOut*pow(minProtonEnergy-thresholdProton,-1.5));
      if(fabs((minProtonEnergy-newMinProtonEnergy)/newMinProtonEnergy)>=1.e-5) {
	if(newMinProtonEnergy<=thresholdProton)
	  newMinProtonEnergy=thresholdProton+(minProtonEnergy-thresholdProton)*0.8;
      } else continueLoop=false;
      minProtonEnergy = newMinProtonEnergy;
    }
  }
  minProtonEnergy+=seperationEnergy_;

  //Find minimum energy for alphas
  double thresholdAlpha = (pType_==3) ? 0. : alphaSepE - seperationEnergy_;
  double minAlphaEnergy;
  if(thresholdAlpha<0.) minAlphaEnergy=0.;
  else {
    minAlphaEnergy=0.122*pow(4.*(compoundA_-4.)/compoundA_,0.33333333333333)*
      pow(2.*(compoundZ_-2.)*rateTemps_[0],0.6666666666667)+thresholdAlpha;
    
    double gamowEnergyIn;
    if(pType_==1) gamowEnergyIn=0.;
    else if(pType_==2) gamowEnergyIn=6.283186/137.0*(compoundZ_-1.)*
			 sqrt((compoundA_-1.)/compoundA_*465.739);
    else if(pType_==3) gamowEnergyIn=6.283186/137.0*2*(compoundZ_-2.)*
			 sqrt(4.*(compoundA_-4.)/compoundA_*465.739);
    double gamowEnergyOut = 
      (pType_==3) ? 0. : 6.283186/137.0*2.*(compoundZ_-2.)*sqrt(4.*(compoundA_-4.)/compoundA_*465.739);

    bool continueLoop = true;
    while(continueLoop) {
      double newMinAlphaEnergy = minAlphaEnergy
	-(-11.605/rateTemps_[0]+0.5*gamowEnergyIn*pow(minAlphaEnergy,-1.5)
	  +0.5*gamowEnergyOut*pow(minAlphaEnergy-thresholdAlpha,-1.5))/
	(-0.75*gamowEnergyIn*pow(minAlphaEnergy,-2.5)-
	 0.75*gamowEnergyOut*pow(minAlphaEnergy-thresholdAlpha,-2.5));
      if(fabs((minAlphaEnergy-newMinAlphaEnergy)/newMinAlphaEnergy)>=1.e-5) {
	if(newMinAlphaEnergy<=thresholdAlpha) 
	  newMinAlphaEnergy=thresholdAlpha+(minAlphaEnergy-thresholdAlpha)*0.8;
      } else continueLoop=false;
      minAlphaEnergy = newMinAlphaEnergy;
    }
  }
  minAlphaEnergy+=seperationEnergy_;
  
  //Setup and sort energy min/max vectors
  std::vector<int_double_pair> minEnergies;
  minEnergies.push_back(int_double_pair(1,minNeutronEnergy));
  minEnergies.push_back(int_double_pair(2,minProtonEnergy));
  minEnergies.push_back(int_double_pair(3,minAlphaEnergy));
  sort(minEnergies.begin(),minEnergies.end(),int_double_pair_compare());
  std::vector<int_double_pair> maxEnergies;
  maxEnergies.push_back(int_double_pair(1,maxNeutronEnergy));
  maxEnergies.push_back(int_double_pair(2,maxProtonEnergy));
  maxEnergies.push_back(int_double_pair(3,maxAlphaEnergy));
  sort(maxEnergies.begin(),maxEnergies.end(),int_double_pair_compare());
  reverse(maxEnergies.begin(),maxEnergies.end());
  int entrancePosition;
  for(int i = 0; i< minEnergies.size();++i) {
    if(minEnergies[i].first==pType_) {
      entrancePosition=i;
      break;
    }
  }
  if(entrancePosition>0) 
    minEnergies.erase(minEnergies.begin(),minEnergies.begin()+entrancePosition);

  //Determine energy windows and number of point in window
  double maxEnergy = maxEnergies[0].second;
  double energyRange = maxEnergy-minEnergies[0].second;
  std::vector<double> energyWindow;
  std::vector<int> numPointsInWindow;
  int pointsUsed = 0;
  int totalPoints = 96;
  for(int i = 0; i< minEnergies.size();++i) {
    if(i==0) energyWindow.push_back(maxEnergy-minEnergies[minEnergies.size()-1].second);
    else energyWindow.push_back(minEnergies[minEnergies.size()-i].second-
				minEnergies[minEnergies.size()-i-1].second);
    double weight = (i==minEnergies.size()-1) ? 5. : 1.;
    numPointsInWindow.push_back(int(energyWindow[i]*totalPoints*weight/energyRange));
    pointsUsed+=numPointsInWindow[i];
  }
  double stretch = double(totalPoints)/double(pointsUsed);
  totalPoints=0;
  for(int i = 0; i< minEnergies.size();++i) {
    numPointsInWindow[i]=int(stretch*numPointsInWindow[i]);
    if(numPointsInWindow[i]==0) numPointsInWindow[i]=1;
    totalPoints+=numPointsInWindow[i];
  }
  reverse(energyWindow.begin(),energyWindow.end());
  reverse(numPointsInWindow.begin(),numPointsInWindow.end());

  //Determine energy points for each window
  std::vector<double> energyGrid;
  for(int i = 0; i< minEnergies.size();++i) {
    double sepE,redMass,z1z2;
    if(minEnergies[minEnergies.size()-i-1].first==1) {
      sepE=neutronSepE;
      redMass=(compoundA_-1.)/compoundA_;
      z1z2=0.0;
    } else if(minEnergies[minEnergies.size()-i-1].first==2) {
      sepE=protonSepE;
      redMass=(compoundA_-1.)/compoundA_;
      z1z2=(compoundZ_-1.);
    } else if(minEnergies[minEnergies.size()-i-1].first==3) {
      sepE=alphaSepE;
      redMass=4.*(compoundA_-4.)/compoundA_;
      z1z2=2.*(compoundZ_-2.);
    }

    double eLow = minEnergies[minEnergies.size()-i-1].second-sepE;
    double eHigh = minEnergies[minEnergies.size()-i-1].second+
      energyWindow[minEnergies.size()-i-1]-sepE;
    int numPoints = numPointsInWindow[minEnergies.size()-i-1];
    std::vector<double> de;
    double sum = 0.;
    for(int j = 0;j<numPoints;j++) {
      if(minEnergies[minEnergies.size()-i-1].first==1) 
	de.push_back(0.097*sqrt(0.5)*
		     exp(log(rateTemps_[0])+
			 log(rateTemps_[rateTemps_.size()-1]/rateTemps_[0])*j/numPoints));
      else 
	de.push_back(0.237*pow(redMass*z1z2,0.3333333333)*
		     pow(exp(log(rateTemps_[0])+
			     log(rateTemps_[rateTemps_.size()-1]/rateTemps_[0])*j/numPoints),
			 0.8333333333333));
      sum+=de[j];
    }
    sum=(eHigh-eLow)/sum;
    energyGrid.push_back(eLow+sepE-seperationEnergy_);
    for(int j = 0;j<numPoints;j++) {
      if(i>0&&j==numPoints-1) continue;
      energyGrid.push_back(energyGrid[energyGrid.size()-1]+de[j]*sum);
    }
    if(i>0&&minEnergies[minEnergies.size()-i].first==1) {
      double lastPoint = energyGrid[energyGrid.size()-1];
      energyGrid.push_back(lastPoint+de[numPoints-1]*sum*0.88);
      energyGrid.push_back(lastPoint+de[numPoints-1]*sum*0.95);      
    }
  }
  sort(energyGrid.begin(),energyGrid.end());

  //Calculate skipping energy
  std::vector<Level> compoundLevels = NuclearLevels::FindLevels(compoundZ_,compoundA_);
  std::vector<Level> neutronLevels = NuclearLevels::FindLevels(compoundZ_,compoundA_-1);
  std::vector<Level> protonLevels = NuclearLevels::FindLevels(compoundZ_-1,compoundA_-1);
  std::vector<Level> alphaLevels = NuclearLevels::FindLevels(compoundZ_-2,compoundA_-4);
  double highestCompoundLevel = 
    (compoundLevels.size()) ? compoundLevels[compoundLevels.size()-1].energy_ : 0.;
  double neutronHighestLevel = 
    (neutronLevels.size()) ? neutronLevels[neutronLevels.size()-1].energy_ : 0.;
  double protonHighestLevel = 
    (protonLevels.size()) ? protonLevels[protonLevels.size()-1].energy_ : 0.;
  double alphaHighestLevel = 
    (alphaLevels.size()) ? alphaLevels[alphaLevels.size()-1].energy_ : 0.;
  std::vector<double> highestLevels;
  highestLevels.push_back(highestCompoundLevel);
  highestLevels.push_back(neutronHighestLevel+neutronSepE);
  highestLevels.push_back(protonHighestLevel+protonSepE);
  highestLevels.push_back(alphaHighestLevel+alphaSepE);
  sort(highestLevels.begin(),highestLevels.end());
  double lowestCompoundLevel = (compoundLevels.size()) ? compoundLevels[0].energy_ : 0.;
  double neutronLowestLevel = (neutronLevels.size()) ? neutronLevels[0].energy_ : 0.;
  double protonLowestLevel = (protonLevels.size()) ? protonLevels[0].energy_ : 0.;
  double alphaLowestLevel = (alphaLevels.size()) ? alphaLevels[0].energy_ : 0.;
  std::vector<double> lowestLevels;
  lowestLevels.push_back(lowestCompoundLevel);
  lowestLevels.push_back(neutronLowestLevel+neutronSepE);
  lowestLevels.push_back(protonLowestLevel+protonSepE);
  lowestLevels.push_back(alphaLowestLevel+alphaSepE);
  sort(lowestLevels.begin(),lowestLevels.end());
  skipEnergy_ = 
    (lowestLevels[lowestLevels.size()-1]==highestLevels[highestLevels.size()-1]) ?
    lowestLevels[lowestLevels.size()-1]+1. : 
    std::max(lowestLevels[lowestLevels.size()-1],highestLevels[highestLevels.size()-1]);
  skipEnergy_-=seperationEnergy_;

  //Fill cross section vector
  for(int i = 0; i<energyGrid.size();i++) {
    crossSections_.push_back(std::pair<double,CrossSectionValues>(energyGrid[i],
								  CrossSectionValues(0.,0.,0.,0.,0.,0.,0.,0.)));
  }
}

void CrossSection::CreateTempVector() {
  /*
  rateTemps_.push_back(0.0001);
  rateTemps_.push_back(0.0005);
  rateTemps_.push_back(0.001);
  rateTemps_.push_back(0.005);
  rateTemps_.push_back(0.01);
  rateTemps_.push_back(0.05);
  */
  rateTemps_.push_back(0.10);
  rateTemps_.push_back(0.15);
  rateTemps_.push_back(0.20);
  //rateTemps_.push_back(0.25);
  rateTemps_.push_back(0.30);
  rateTemps_.push_back(0.40);
  rateTemps_.push_back(0.50);
  rateTemps_.push_back(0.60);
  rateTemps_.push_back(0.70);
  rateTemps_.push_back(0.80);
  rateTemps_.push_back(0.90);
  rateTemps_.push_back(1.00);
  rateTemps_.push_back(1.50);
  rateTemps_.push_back(2.00);
  rateTemps_.push_back(2.50);
  rateTemps_.push_back(3.00);
  rateTemps_.push_back(3.50);
  rateTemps_.push_back(4.00);
  rateTemps_.push_back(4.50);
  rateTemps_.push_back(5.00);
  rateTemps_.push_back(6.00);
  rateTemps_.push_back(7.00);
  rateTemps_.push_back(8.00);
  rateTemps_.push_back(9.00);
  rateTemps_.push_back(10.00);
}

void CrossSection::CreateMACSEnergiesVector() {
  macsEnergies_.push_back(0.005);
  macsEnergies_.push_back(0.010);
  macsEnergies_.push_back(0.015);
  macsEnergies_.push_back(0.020);
  macsEnergies_.push_back(0.025);
  macsEnergies_.push_back(0.030);
  macsEnergies_.push_back(0.040);
  macsEnergies_.push_back(0.050);
  macsEnergies_.push_back(0.060);
  macsEnergies_.push_back(0.080);
  macsEnergies_.push_back(0.100);
}

struct gsl_reactionrate_params {
  double temperature;
  TGraph* graph;
  bool useSpline;
};

double gsl_reactionrate_integrand(double x, void * p) {
  struct gsl_reactionrate_params *params= (struct gsl_reactionrate_params *)p;
  double temperature = params->temperature;
  TGraph* graph = params->graph;
  bool useSpline = params->useSpline;

  double crossSection= (useSpline) ? graph->Eval(x,0,"S") : graph->Eval(x);
  return crossSection*x*exp(-x/temperature/boltzConst);
}


void CrossSection::CalculateReactionRates(bool macs) {
  if(macs) std::cout << "Calculating MACS..." << std::endl;
  else std::cout << "Calculating Reaction Rates..." << std::endl;
  double energies[crossSections_.size()];
  double gamma[crossSections_.size()];
  double neutron[crossSections_.size()];
  double proton[crossSections_.size()];
  double alpha[crossSections_.size()];
  double gammaStellar[crossSections_.size()];
  double neutronStellar[crossSections_.size()];
  double protonStellar[crossSections_.size()];
  double alphaStellar[crossSections_.size()];
  int j = 0;
  for(int i = 0;i<crossSections_.size();i++) {
    if(skipped_[i]) continue;
    energies[j]=crossSections_[i].first;
    gamma[j]=crossSections_[i].second.gamma_;
    neutron[j]=crossSections_[i].second.neutron_;
    proton[j]=crossSections_[i].second.proton_;
    alpha[j]=crossSections_[i].second.alpha_;
    gammaStellar[j]=crossSections_[i].second.gammaStellar_;
    neutronStellar[j]=crossSections_[i].second.neutronStellar_;
    protonStellar[j]=crossSections_[i].second.protonStellar_;
    alphaStellar[j]=crossSections_[i].second.alphaStellar_;
    j++;
  }
  TGraph* g_gamma = new TGraph(j,energies,gamma);
  TGraph* g_neutron = new TGraph(j,energies,neutron);
  TGraph* g_proton = new TGraph(j,energies,proton);
  TGraph* g_alpha = new TGraph(j,energies,alpha);
  TGraph* g_gammaStellar = new TGraph(j,energies,gammaStellar);
  TGraph* g_neutronStellar = new TGraph(j,energies,neutronStellar);
  TGraph* g_protonStellar = new TGraph(j,energies,protonStellar);
  TGraph* g_alphaStellar = new TGraph(j,energies,alphaStellar);
  reactionRates_.clear();
  if(macs) 
    for(std::vector<double>::const_iterator it = macsEnergies_.begin();
	it<macsEnergies_.end();++it) {
      struct gsl_reactionrate_params paramsGamma = {(*it)/boltzConst,g_gamma,true};
      struct gsl_reactionrate_params paramsGammaStellar = {(*it)/boltzConst,g_gammaStellar,true};

       gsl_integration_workspace * w 
	= gsl_integration_workspace_alloc (1000);
      
      gsl_function F;
      F.function = &gsl_reactionrate_integrand;
      
      double error,gammaRate,gammaStellarRate;
      
      F.params=&paramsGamma;
      gsl_integration_qag(&F,energies[0],energies[j-1],0.0,1e-4,1000,3,w,&gammaRate,&error);
      F.params=&paramsGammaStellar;
      gsl_integration_qag(&F,energies[0],energies[j-1],0.0,1e-4,1000,3,w,&gammaStellarRate,&error);
            
      gsl_integration_workspace_free (w);
           
      gammaRate*=2./sqrt(pi)/(*it)/(*it);
      gammaStellarRate*=2./sqrt(pi)/(*it)/(*it);
    
      reactionRates_.push_back(std::pair<double,CrossSectionValues>((*it),
								    CrossSectionValues(gammaRate,
										       0.,
										       0.,
										       0.,
										       gammaStellarRate,
										       0.,
										       0.,
										       0.)));
    }
  else
    for(std::vector<double>::const_iterator it = rateTemps_.begin();
	it<rateTemps_.end();++it) {
      struct gsl_reactionrate_params paramsGamma = {(*it),g_gamma,true};
      struct gsl_reactionrate_params paramsNeutron = {(*it),g_neutron,false};
      struct gsl_reactionrate_params paramsProton = {(*it),g_proton,false};
      struct gsl_reactionrate_params paramsAlpha = {(*it),g_alpha,false};
      struct gsl_reactionrate_params paramsGammaStellar = {(*it),g_gammaStellar,true};
      struct gsl_reactionrate_params paramsNeutronStellar = {(*it),g_neutronStellar,false};
      struct gsl_reactionrate_params paramsProtonStellar = {(*it),g_protonStellar,false};
      struct gsl_reactionrate_params paramsAlphaStellar = {(*it),g_alphaStellar,false};
      
      gsl_integration_workspace * w 
	= gsl_integration_workspace_alloc (1000);
      
      gsl_function F;
      F.function = &gsl_reactionrate_integrand;
      
      double error,gammaRate=0.,neutronRate=0.,protonRate=0.,alphaRate=0.,gammaStellarRate=0.,
	neutronStellarRate=0.,protonStellarRate=0.,alphaStellarRate=0.;
      
      F.params=&paramsGamma;
      gsl_integration_qag(&F,energies[0],energies[j-1],0.0,1e-3,1000,3,w,&gammaRate,&error);
      F.params=&paramsNeutron;
      gsl_integration_qag(&F,energies[0],energies[j-1],0.0,1e-3,1000,3,w,&neutronRate,&error);
      F.params=&paramsProton;
      gsl_integration_qag(&F,energies[0],energies[j-1],0.0,1e-3,1000,3,w,&protonRate,&error);
      F.params=&paramsAlpha;
      gsl_integration_qag(&F,energies[0],energies[j-1],0.0,1e-3,1000,3,w,&alphaRate,&error);
      F.params=&paramsGammaStellar;
      gsl_integration_qag(&F,energies[0],energies[j-1],0.0,1e-3,1000,3,w,&gammaStellarRate,&error);
      F.params=&paramsNeutronStellar;
      gsl_integration_qag(&F,energies[0],energies[j-1],0.0,1e-3,1000,3,w,&neutronStellarRate,&error);
      F.params=&paramsProtonStellar;
      gsl_integration_qag(&F,energies[0],energies[j-1],0.0,1e-3,1000,3,w,&protonStellarRate,&error);
      F.params=&paramsAlphaStellar;
      gsl_integration_qag(&F,energies[0],energies[j-1],0.0,1e-3,1000,3,w,&alphaStellarRate,&error);
      
      gsl_integration_workspace_free (w);
      
      double redMass = (pType_==1||pType_==2) ? double(compoundA_-1.)/double(compoundA_) :
	4.*double(compoundA_-4.)/double(compoundA_);
      gammaRate*=1e-24*avagadroNum*lightSpeedInCmPerS*sqrt(8.0/pi/redMass/uconv)/
	pow(boltzConst*(*it),1.5);
      neutronRate*=1e-24*avagadroNum*lightSpeedInCmPerS*sqrt(8.0/pi/redMass/uconv)/
	pow(boltzConst*(*it),1.5);
      protonRate*=1e-24*avagadroNum*lightSpeedInCmPerS*sqrt(8.0/pi/redMass/uconv)/
	pow(boltzConst*(*it),1.5);
      alphaRate*=1e-24*avagadroNum*lightSpeedInCmPerS*sqrt(8.0/pi/redMass/uconv)/
	pow(boltzConst*(*it),1.5);
      gammaStellarRate*=1e-24*avagadroNum*lightSpeedInCmPerS*sqrt(8.0/pi/redMass/uconv)/
	pow(boltzConst*(*it),1.5);
      neutronStellarRate*=1e-24*avagadroNum*lightSpeedInCmPerS*sqrt(8.0/pi/redMass/uconv)/
	pow(boltzConst*(*it),1.5);
      protonStellarRate*=1e-24*avagadroNum*lightSpeedInCmPerS*sqrt(8.0/pi/redMass/uconv)/
	pow(boltzConst*(*it),1.5);
      alphaStellarRate*=1e-24*avagadroNum*lightSpeedInCmPerS*sqrt(8.0/pi/redMass/uconv)/
	pow(boltzConst*(*it),1.5);
      
      reactionRates_.push_back(std::pair<double,CrossSectionValues>((*it),
								    CrossSectionValues(gammaRate,
										       neutronRate,
										       protonRate,
										       alphaRate,
										       gammaStellarRate,
										       neutronStellarRate,
										       protonStellarRate,
										       alphaStellarRate)));
    }
  delete g_gamma;
  delete g_neutron;
  delete g_proton;
  delete g_alpha;
  delete g_gammaStellar;
  delete g_neutronStellar;
  delete g_protonStellar;
  delete g_alphaStellar;
}

void CrossSection::PrintReactionRates(bool macs) {
  char filename[256];
  if(pType_==0) {
    sprintf(filename,"Sapphire_%d%s+g_rates.dat",
	    A_,
	    NuclearMass::FindElement(Z_).c_str());  
  } else if(pType_==1) {
    if(macs) sprintf(filename,"Sapphire_%d%s+n_macs.dat",
		     A_,
		     NuclearMass::FindElement(Z_).c_str());  
    else sprintf(filename,"Sapphire_%d%s+n_rates.dat",
		 A_,
		 NuclearMass::FindElement(Z_).c_str());  
  } else if(pType_==2) {
    sprintf(filename,"Sapphire_%d%s+p_rates.dat",
	    A_,
	    NuclearMass::FindElement(Z_).c_str());  
  } else {
    sprintf(filename,"Sapphire_%d%s+a_rates.dat",
	    A_,
	    NuclearMass::FindElement(Z_).c_str());  
  }
  std::ofstream out(filename);
  if(macs) out << "#"
	       << std::setw(14) << "kT [MeV]"
	       << std::setw(15) << "G"
	       << std::setw(15) << "MACS [b]"
	       << std::setw(15) << "MACS* [b]"
	       << std::setw(15) << "SEF"
	       << std::setw(0) << std::endl;
  else out << "#" 
	   << std::setw(14) << "T9"
	   << std::setw(15) << "G"
	   << std::setw(15) << "gamma"
	   << std::setw(15) << "gamma*"
	   << std::setw(15) << "gamma SEF"
           << std::setw(15) << "neutron" 
           << std::setw(15) << "neutron*" 
	   << std::setw(15) << "neutron SEF"
	   << std::setw(15) << "proton"
	   << std::setw(15) << "proton*"
	   << std::setw(15) << "proton SEF"
	   << std::setw(15) << "alpha"
	   << std::setw(15) << "alpha*"
	   << std::setw(15) << "alpha SEF"
	   << std::setw(0)  << std::endl; 
  for(int i = 0;i<reactionRates_.size();i++) {
    double partFunc = (macs) ? partFuncMACS_[i].second : partFunc_[i].second;
    out << std::scientific
	<< std::setw(15) << reactionRates_[i].first
	<< std::setw(15) << partFunc
	<< std::setw(15) << reactionRates_[i].second.gamma_ 
	<< std::setw(15) << reactionRates_[i].second.gammaStellar_/partFunc
	<< std::setw(15) << reactionRates_[i].second.gammaStellar_/partFunc/reactionRates_[i].second.gamma_;
    if(macs) out << std::setw(0) << std::endl;
    else out << std::setw(15) << reactionRates_[i].second.neutron_
	     << std::setw(15) << reactionRates_[i].second.neutronStellar_/partFunc
	     << std::setw(15) << reactionRates_[i].second.neutronStellar_/partFunc/reactionRates_[i].second.neutron_
	     << std::setw(15) << reactionRates_[i].second.proton_ 
	     << std::setw(15) << reactionRates_[i].second.protonStellar_/partFunc
	     << std::setw(15) << reactionRates_[i].second.protonStellar_/partFunc/reactionRates_[i].second.protonStellar_
	     << std::setw(15) << reactionRates_[i].second.alpha_ 
	     << std::setw(15) << reactionRates_[i].second.alphaStellar_/partFunc
	     << std::setw(15) << reactionRates_[i].second.alphaStellar_/partFunc/reactionRates_[i].second.alpha_
	     << std::setw(0)  << std::endl; 
  }
  out.flush();
  out.close();
}

struct gsl_partfunc_params {
  double temperature;
  LevelDensity* density;
};

double gsl_partfunc_integrand(double x, void * p) {
  struct gsl_partfunc_params* params = (struct gsl_partfunc_params*) p;
  double temperature = params->temperature;
  LevelDensity* density = params->density;
  return (x>1000.||x/boltzConst/temperature>100.) ? 
    0. : exp(-x/boltzConst/temperature)*density->operator()(x);
}

void CrossSection::CalcPartitionFunc() {
  std::vector<Level> knownLevels = NuclearLevels::FindLevels(Z_,A_);
  double highestLevel = (knownLevels.size()) ? knownLevels[knownLevels.size()-1].energy_ : 0.;
  double spinOffset = (A_%2==0) ? 0. : 0.5;
  for(int i=0;i<2;++i) {
    std::vector<double>::const_iterator end = (i==0) ? rateTemps_.end() : macsEnergies_.end();
    for(std::vector<double>::const_iterator it = (i==0) ? rateTemps_.begin() : macsEnergies_.begin();it<end;++it) {
      double sum = 0.;
      double temperature = (i==0) ? (*it) : (*it)/boltzConst;
      for(int k=0; k<knownLevels.size();k++) {
	sum+=exp(-knownLevels[k].energy_/temperature/boltzConst)*(2.*knownLevels[k].J_+1.);
      }
      for(double j=0.;j<10.0;j+=1.) {	
	double jValue = spinOffset+j;
	LevelDensity* den = new RauscherLevelDensity(Z_,A_,jValue);
	gsl_integration_workspace * w 
	  = gsl_integration_workspace_alloc (1000);
	gsl_function F;
	F.function = &gsl_partfunc_integrand;  
	double error,value;
	struct gsl_partfunc_params params = {temperature,den};
	F.params=&params;
	gsl_integration_qagiu(&F,highestLevel,0.0,1e-3,1000,w,&value,&error);
	gsl_integration_workspace_free (w);
	sum+=2.*(2.*jValue+1.)*value;
	delete den; 
      }
      if(i==0) partFunc_.push_back(std::pair<double,double>(temperature,sum/(2.*groundStateJ_+1.)));
      else partFuncMACS_.push_back(std::pair<double,double>(temperature,sum/(2.*groundStateJ_+1.)));
    }
  }
}


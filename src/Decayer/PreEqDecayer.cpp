#include "Decayer/PreEqDecayer.h"
#include "Databases/NuclearMass.h"
#include "PreEqTransitionRateFunc.h"
#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#ifndef MPI_BUILD
#include <omp.h>
#endif

#ifndef MPI_BUILD
extern unsigned int randomSeed[12];
#endif

PreEqDecayer::PreEqDecayer(int initialNeutronNumber,
			   int initialNeutronHoleNumber,
			   int initialProtonNumber,
			   int initialProtonHoleNumber,
			   int Z, int A, double jInitial,
			   int piInitial, double energy) :
  initialNeutronNumber_(initialNeutronNumber),
  initialNeutronHoleNumber_(initialNeutronHoleNumber),
  initialProtonNumber_(initialProtonNumber),
  initialProtonHoleNumber_(initialProtonHoleNumber),
  Z_(Z), A_(A), piInitial_(piInitial), jInitial_(jInitial), 
  energy_(energy)
{
  totalIntegral_=0.;    

  initialExcitonNumber_ = initialNeutronNumber+
    initialNeutronHoleNumber+
    initialProtonNumber+
    initialProtonHoleNumber;

  double qValueProton,qValueNeutron;
  if(!NuclearMass::QValue(Z,A,Z,A-1,qValueNeutron)) {
    std::cout << "Unknown masses requested.  Aborting." << std::endl;
    exit(1);
  }
  if(!NuclearMass::QValue(Z,A,Z-1,A-1,qValueProton)) {
    std::cout << "Unknown masses requested.  Aborting." << std::endl;
    exit(1);
  }

  if(initialExcitonNumber_>1&&(qValueNeutron+energy>0||qValueProton+energy>0)) {
    for(double l=0;l<=maxL_;l+=1.) {
      int piFinal = (int(l)%2==0) ? piInitial : -1*piInitial;
      for(double s = fabs(l-0.5); s<=l+0.5; s+=1.) {
	for(double jFinal = fabs(s-jInitial); jFinal<=s+jInitial;jFinal+=1.) {
	  if(qValueNeutron+energy>0) {
	    bool exists = false;
	    for(int i =0;i<spinRatePairs_.size();i++) {
	      if(spinRatePairs_[i].Z_==Z&&
		 spinRatePairs_[i].A_==A-1&&
		 spinRatePairs_[i].spin_==jFinal&&
		 spinRatePairs_[i].parity_==piFinal) {
		exists=true;
		break;
	      }
	    }
	    if(!exists) {
	      PreEqTransitionRateFunc* newFunc = 
		new PreEqTransitionRateFunc(0,1,Z,A-1,
					    initialNeutronNumber,
					    initialNeutronHoleNumber,
					    initialProtonNumber,
					    initialProtonHoleNumber,
					    initialNeutronNumber-1,
					    initialNeutronHoleNumber,
					    initialProtonNumber,
					    initialProtonHoleNumber,
					    jInitial,piInitial,
					    jFinal,piFinal,0.5,1,maxL_,energy,
					    qValueNeutron);
	      spinRatePairs_.push_back(PreEqSpinRatePair(initialNeutronNumber-1,initialNeutronHoleNumber,
							 initialProtonNumber, initialProtonHoleNumber,
							 Z,A-1,jFinal,piFinal,qValueNeutron,
							 newFunc,newFunc->Integral()));
	      totalIntegral_+=newFunc->Integral();
	    }
	  }
	  if(qValueProton+energy>0) {
	    bool exists = false;
	    for(int i =0;i<spinRatePairs_.size();i++) {
	      if(spinRatePairs_[i].Z_==Z-1&&
		 spinRatePairs_[i].A_==A&&
		 spinRatePairs_[i].spin_==jFinal&&
		 spinRatePairs_[i].parity_==piFinal) {
		exists=true;
		break;
	      }
	    }
	    if(!exists) {
	      PreEqTransitionRateFunc* newFunc = 
		new PreEqTransitionRateFunc(1,1,Z-1,A-1,
					    initialNeutronNumber,
					    initialNeutronHoleNumber,
					    initialProtonNumber,
					    initialProtonHoleNumber,
					    initialNeutronNumber,
					    initialNeutronHoleNumber,
					    initialProtonNumber-1,
					    initialProtonHoleNumber,
					    jInitial,piInitial,
					    jFinal,piFinal,0.5,1,maxL_,energy,
					    qValueProton);
	      spinRatePairs_.push_back(PreEqSpinRatePair(initialNeutronNumber,initialNeutronHoleNumber,
							 initialProtonNumber-1,initialProtonHoleNumber,
							 Z-1,A-1,jFinal,piFinal,qValueProton,
							 newFunc,newFunc->Integral()));
	      totalIntegral_+=newFunc->Integral();
	    }
	  }
	}
      }
    }
  }
  //Neutron pair production
  PreEqTransitionRateFunc* newFunc = 
    new PreEqTransitionRateFunc(0,0,Z,A,
				initialNeutronNumber,
				initialNeutronHoleNumber,
				initialProtonNumber,
				initialProtonHoleNumber,
				initialNeutronNumber+1,
				initialNeutronHoleNumber+1,
				initialProtonNumber,
				initialProtonHoleNumber,
				jInitial,piInitial,
				jInitial,piInitial,0.0,1,maxL_,energy,
				0.);
  spinRatePairs_.push_back(PreEqSpinRatePair(initialNeutronNumber+1,initialNeutronHoleNumber+1,
					     initialProtonNumber,initialProtonHoleNumber,
					     Z,A,jInitial,piInitial,0.,
					     newFunc,newFunc->Integral()));
  totalIntegral_+=newFunc->Integral();
  //proton pair production
  newFunc = new PreEqTransitionRateFunc(0,0,Z,A,
					initialNeutronNumber,
					initialNeutronHoleNumber,
					initialProtonNumber,
					initialProtonHoleNumber,
					initialNeutronNumber,
					initialNeutronHoleNumber,
					initialProtonNumber+1,
					initialProtonHoleNumber+1,
					jInitial,piInitial,
					jInitial,piInitial,0.0,1,maxL_,energy,
					0.);
  spinRatePairs_.push_back(PreEqSpinRatePair(initialNeutronNumber,initialNeutronHoleNumber,
					     initialProtonNumber+1,initialProtonHoleNumber+1,
					     Z,A,jInitial,piInitial,0.,
					     newFunc,newFunc->Integral()));
  totalIntegral_+=newFunc->Integral();
  //proton pair conversion
  newFunc = new PreEqTransitionRateFunc(0,0,Z,A,
					initialNeutronNumber,
					initialNeutronHoleNumber,
					initialProtonNumber,
					initialProtonHoleNumber,
					initialNeutronNumber+1,
					initialNeutronHoleNumber+1,
					initialProtonNumber-1,
					initialProtonHoleNumber-1,
					jInitial,piInitial,
					jInitial,piInitial,0.0,1,maxL_,energy,
					0.);
  spinRatePairs_.push_back(PreEqSpinRatePair(initialNeutronNumber+1,initialNeutronHoleNumber+1,
					     initialProtonNumber-1,initialProtonHoleNumber-1,
					     Z,A,jInitial,piInitial,0.,
					     newFunc,newFunc->Integral()));
  totalIntegral_+=newFunc->Integral();
  //neutron pair conversion
  newFunc = new PreEqTransitionRateFunc(0,0,Z,A,
					initialNeutronNumber,
					initialNeutronHoleNumber,
					initialProtonNumber,
					initialProtonHoleNumber,
					initialNeutronNumber-1,
					initialNeutronHoleNumber-1,
					initialProtonNumber+1,
					initialProtonHoleNumber+1,
					jInitial,piInitial,
					jInitial,piInitial,0.0,1,maxL_,energy,
					0.);
  spinRatePairs_.push_back(PreEqSpinRatePair(initialNeutronNumber-1,initialNeutronHoleNumber-1,
					     initialProtonNumber+1,initialProtonHoleNumber+1,
					     Z,A,jInitial,piInitial,0.,
					     newFunc,newFunc->Integral()));
  totalIntegral_+=newFunc->Integral();

  if(!isCrossSection_) BuildCDF();
}

PreEqDecayer::~PreEqDecayer(){
  for(int i=0 ;i<spinRatePairs_.size();i++) {
    delete spinRatePairs_[i].rateFunc_;
  }
}

void PreEqDecayer::BuildCDF() {
  double offSet=0.;
  for(int pair = 0;pair<spinRatePairs_.size();pair++) {
    std::vector<XYPair> cumulativeSum = 
      spinRatePairs_[pair].rateFunc_->CumulativeSum();
    for(int e = 0;e<cumulativeSum.size();e++) {
      double cdfValue = (cumulativeSum[e].Y_+offSet)/totalIntegral_;
      cdf_.push_back(PreEqCDFEntry(pair,cumulativeSum[e].X_,cdfValue));
    }
    offSet+=spinRatePairs_[pair].rateFunc_->Integral();
  }
}

void PreEqDecayer::PrintCDF() {
  std::ofstream out("cdf.out");
  for(int i=0;i<cdf_.size();i++) {
    int pairIndex = cdf_[i].pairIndex_;
    out << std::fixed << spinRatePairs_[pairIndex].Z_ <<  ' ' << spinRatePairs_[pairIndex].A_ << ' '
	<< spinRatePairs_[pairIndex].neutronNumber_ << ' '
	<< spinRatePairs_[pairIndex].neutronHoleNumber_ << ' '
	<< spinRatePairs_[pairIndex].protonNumber_ << ' '
	<< spinRatePairs_[pairIndex].protonHoleNumber_ << ' '      
	<< spinRatePairs_[pairIndex].spin_ << ' ' << spinRatePairs_[pairIndex].parity_ << ' '
	<< std::scientific << cdf_[i].energy_ << ' ' <<  cdf_[i].value_ << std::endl;
  }
  out.flush();
  out.close();
}

bool PreEqDecayer::Decay(int& Z, int& A, double& jFinal, int& piFinal,
			 int& neutronNumber, int& neutronHoleNumber,
			 int& protonNumber, int& protonHoleNumber,
			 double& excitationEnergy, double& decayEnergy) {
  double randomNumber=0.;
#ifndef MPI_BUILD
  while(randomNumber==0.) randomNumber = double(rand_r(&randomSeed[omp_get_thread_num()]))/double(RAND_MAX);
#else
  while(randomNumber==0.) randomNumber = double(rand())/double(RAND_MAX);
#endif
  bool found = false;
  double previousValue = 0.;
  for(std::vector<PreEqCDFEntry>::const_iterator it = cdf_.begin();it<cdf_.end();++it) {
    if(previousValue<randomNumber&&randomNumber<=it->value_) {
      Z = spinRatePairs_[it->pairIndex_].Z_;
      A = spinRatePairs_[it->pairIndex_].A_;
      jFinal = spinRatePairs_[it->pairIndex_].spin_;
      piFinal = spinRatePairs_[it->pairIndex_].parity_;
      neutronNumber = spinRatePairs_[it->pairIndex_].neutronNumber_;
      neutronHoleNumber = spinRatePairs_[it->pairIndex_].neutronHoleNumber_;
      protonNumber = spinRatePairs_[it->pairIndex_].protonNumber_;
      protonHoleNumber = spinRatePairs_[it->pairIndex_].protonHoleNumber_;
      excitationEnergy = energy_+spinRatePairs_[it->pairIndex_].qValue_-it->energy_;
      decayEnergy = it->energy_;
      found = true;
      break;
    }
  }
  return found;
}

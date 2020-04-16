#include "DecayResults.h"
#include "DecayProduct.h"
#include "NuclearMass.h"
#include <fstream>
#include <iostream>
#include <chrono>
#include <thread>
#include "omp.h"

DecayResults::DecayResults(int Z, int A, double J, int Pi,
                           double initialEnergyLow, double initialEnergyHigh,
			   int suffixNo) :
 initialZ_(Z),initialA_(A),initialPi_(Pi),initialJ_(J),
 initialEnergyLow_(initialEnergyLow), initialEnergyHigh_(initialEnergyHigh) {
  char spin[10];
  if(initialPi_==1) sprintf(spin,"%.1f+",initialJ_);
  else sprintf(spin,"%.1f-",initialJ_);

  char filename[256];
  if(suffixNo!=0) {
    if(initialEnergyLow_==initialEnergyHigh_) 
      sprintf(filename,"Sapphire_%d%s_J=%s_E=%.3f_%d.root",
	      initialA_,NuclearMass::FindElement(initialZ_).c_str(),
	      spin,initialEnergyLow_,suffixNo);
    else sprintf(filename,"Sapphire_%d%s_J=%s_E=%.3f-%.3f_%d.root",
		 initialA_,NuclearMass::FindElement(initialZ_).c_str(),
		 spin,initialEnergyLow_,initialEnergyHigh_,suffixNo);
  } else {
    bool validFileName=false;
    int i  = 0;
    while(!validFileName) {
      if(i==0) {
	if(initialEnergyLow_==initialEnergyHigh_) 
	  sprintf(filename,"Sapphire_%d%s_J=%s_E=%.3f.root",
		  initialA_,NuclearMass::FindElement(initialZ_).c_str(),
		  spin,initialEnergyLow_);
	else sprintf(filename,"Sapphire_%d%s_J=%s_E=%.3f-%.3f.root",
		     initialA_,NuclearMass::FindElement(initialZ_).c_str(),
		     spin,initialEnergyLow_,initialEnergyHigh_);
      } else {
	if(initialEnergyLow_==initialEnergyHigh_) 
	  sprintf(filename,"Sapphire_%d%s_J=%s_E=%.3f_%d.root",
		  initialA_,NuclearMass::FindElement(initialZ_).c_str(),
		  spin,initialEnergyLow_,i+1);
	else sprintf(filename,"Sapphire_%d%s_J=%s_E=%.3f-%.3f_%d.root",
		     initialA_,NuclearMass::FindElement(initialZ_).c_str(),
		     spin,initialEnergyLow_,initialEnergyHigh_,i+1);
      }
      std::ifstream inTest(filename);
      if(!inTest) {
	validFileName=true;
      } else inTest.close();
      i++;
    }
  }
  std::cout << "Output filename is " << filename << std::endl;

  outputFile_ = new TFile(filename,"recreate");
  outputTree_ = new TTree("statDecay","Sapphire Results");
    
  outputTree_->Branch("initialEnergy",&initialEnergy_,"initialEnergy/D");
  outputTree_->Branch("neutronTotalWidth",&neutronTotalWidth_,"neutronTotalWidth/D");
  outputTree_->Branch("protonTotalWidth",&protonTotalWidth_,"protonTotalWidth/D");
  outputTree_->Branch("gammaTotalWidth",&gammaTotalWidth_,"gammaTotalWidth/D");
  outputTree_->Branch("alphaTotalWidth",&alphaTotalWidth_,"alphaTotalWidth/D");
  outputTree_->Branch("neutronEntranceWidth",&neutronEntranceWidth_,"neutronEntranceWidth/D");
  outputTree_->Branch("protonEntranceWidth",&protonEntranceWidth_,"protonEntranceWidth/D");
  outputTree_->Branch("gammaEntranceWidth",&gammaEntranceWidth_,"gammaEntranceWidth/D");
  outputTree_->Branch("alphaEntranceWidth",&alphaEntranceWidth_,"alphaEntranceWidth/D");
  outputTree_->Branch("numSteps",&numSteps_,"numSteps/I");
  outputTree_->Branch("fragmentEnergy",fragmentEnergy_,"fragmentEnergy[numSteps]/D");
  outputTree_->Branch("fragmentExcitation",fragmentExcitation_,"fragmentExcitation[numSteps]/D");
  outputTree_->Branch("fragmentMomentumX",fragmentMomentumX_,"fragmentMomentumX[numSteps]/D");
  outputTree_->Branch("fragmentMomentumY",fragmentMomentumY_,"fragmentMomentumY[numSteps]/D");
  outputTree_->Branch("fragmentMomentumZ",fragmentMomentumZ_,"fragmentMomentumZ[numSteps]/D");
  outputTree_->Branch("Z",Z_,"Z[numSteps]/I");
  outputTree_->Branch("A",A_,"A[numSteps]/I");
  outputTree_->Branch("numGammas",&numGammas_,"numGammas/I");
  outputTree_->Branch("gammaStepIndex",gammaStepIndex_,"gammaStepIndex[numGammas]/D");
  outputTree_->Branch("gammaEnergy",gammaEnergy_,"gammaEnergy[numGammas]/D");
  outputTree_->Branch("gammaMomentumX",gammaMomentumX_,"gammaMomentumX[numGammas]/D");
  outputTree_->Branch("gammaMomentumY",gammaMomentumY_,"gammaMomentumY[numGammas]/D");
  outputTree_->Branch("gammaMomentumZ",gammaMomentumZ_,"gammaMomentumZ[numGammas]/D");
  outputTree_->Branch("numNeutrons",&numNeutrons_,"numNeutrons/I");
  outputTree_->Branch("neutronStepIndex",neutronStepIndex_,"neutronStepIndex[numNeutrons]/D");
  outputTree_->Branch("neutronEnergy",neutronEnergy_,"neutronEnergy[numNeutrons]/D");
  outputTree_->Branch("neutronMomentumX",neutronMomentumX_,"neutronMomentumX[numNeutrons]/D");
  outputTree_->Branch("neutronMomentumY",neutronMomentumY_,"neutronMomentumY[numNeutrons]/D");
  outputTree_->Branch("neutronMomentumZ",neutronMomentumZ_,"neutronMomentumZ[numNeutrons]/D");
  outputTree_->Branch("numProtons",&numProtons_,"numProtons/I");
  outputTree_->Branch("protonStepIndex",protonStepIndex_,"protonStepIndex[numProtons]/D");
  outputTree_->Branch("protonEnergy",protonEnergy_,"protonEnergy[numProtons]/D");
  outputTree_->Branch("protonMomentumX",protonMomentumX_,"protonMomentumX[numProtons]/D");
  outputTree_->Branch("protonMomentumY",protonMomentumY_,"protonMomentumY[numProtons]/D");
  outputTree_->Branch("protonMomentumZ",protonMomentumZ_,"protonMomentumZ[numProtons]/D");
  outputTree_->Branch("numAlphas",&numAlphas_,"numAlphas/I");
  outputTree_->Branch("alphaStepIndex",alphaStepIndex_,"alphaStepIndex[numAlphas]/D");
  outputTree_->Branch("alphaEnergy",alphaEnergy_,"alphaEnergy[numAlphas]/D");
  outputTree_->Branch("alphaMomentumX",alphaMomentumX_,"alphaMomentumX[numAlphas]/D");
  outputTree_->Branch("alphaMomentumY",alphaMomentumY_,"alphaMomentumY[numAlphas]/D");
  outputTree_->Branch("alphaMomentumZ",alphaMomentumZ_,"alphaMomentumZ[numAlphas]/D");
}

DecayResults::~DecayResults() {
  std::cout << "Waiting on threads to be finished..." << std::endl;
  while(omp_get_num_threads()!=1){
    std::this_thread::sleep_for(std::chrono::seconds(1));
  }
  outputFile_->Write();
  outputFile_->Close();
  delete outputFile_;
}

void DecayResults::AddResults(std::vector<std::pair<DecayData, std::vector<DecayProduct> > >& results) {
  for(int i = 0;i<results.size();i++) {
    if(results[i].second.size()==0) continue;
    initialEnergy_ = results[i].first.energy();
    neutronTotalWidth_ = results[i].first.neutronTotalWidth();
    gammaTotalWidth_ = results[i].first.gammaTotalWidth();
    protonTotalWidth_ = results[i].first.protonTotalWidth();
    alphaTotalWidth_ = results[i].first.alphaTotalWidth();
    neutronEntranceWidth_ = results[i].first.neutronEntranceWidth();
    gammaEntranceWidth_ = results[i].first.gammaEntranceWidth();
    protonEntranceWidth_ = results[i].first.protonEntranceWidth();
    alphaEntranceWidth_ = results[i].first.alphaEntranceWidth();

    numSteps_=numNeutrons_=numGammas_=numProtons_=numAlphas_=0;
    for(std::vector<DecayProduct>::const_iterator it = results[i].second.begin();
	it<results[i].second.end();++it) {
      Z_[numSteps_] = it->Z_;
      A_[numSteps_] = it->A_;
      fragmentEnergy_[numSteps_] = it->fragmentEnergy_;
      fragmentExcitation_[numSteps_]= it->excitationEnergy_;
      fragmentMomentumX_[numSteps_] = it->fragmentMomentumX_;
      fragmentMomentumY_[numSteps_] = it->fragmentMomentumY_;
      fragmentMomentumZ_[numSteps_] = it->fragmentMomentumZ_; 
      if(it->particleType_==0) {
	gammaStepIndex_[numGammas_] = numSteps_;
	gammaEnergy_[numGammas_]=it->particleEnergy_;
	gammaMomentumX_[numGammas_]=it->particleMomentumX_;
	gammaMomentumY_[numGammas_]=it->particleMomentumY_;
	gammaMomentumZ_[numGammas_]=it->particleMomentumZ_;
	++numGammas_;
      } else if(it->particleType_==1) {
	neutronStepIndex_[numNeutrons_] = numSteps_;
	neutronEnergy_[numNeutrons_]=it->particleEnergy_;
	neutronMomentumX_[numNeutrons_]=it->particleMomentumX_;
	neutronMomentumY_[numNeutrons_]=it->particleMomentumY_;
	neutronMomentumZ_[numNeutrons_]=it->particleMomentumZ_;
	++numNeutrons_;
      } else if(it->particleType_==2) {
	protonStepIndex_[numProtons_] = numSteps_;
	protonEnergy_[numProtons_]=it->particleEnergy_;
	protonMomentumX_[numProtons_]=it->particleMomentumX_;
	protonMomentumY_[numProtons_]=it->particleMomentumY_;
	protonMomentumZ_[numProtons_]=it->particleMomentumZ_;
	++numProtons_;
      } else if(it->particleType_==3) {
	alphaStepIndex_[numAlphas_] = numSteps_;
	alphaEnergy_[numAlphas_]=it->particleEnergy_;
	alphaMomentumX_[numAlphas_]=it->particleMomentumX_;
	alphaMomentumY_[numAlphas_]=it->particleMomentumY_;
	alphaMomentumZ_[numAlphas_]=it->particleMomentumZ_;
	++numAlphas_;
      }
      ++numSteps_;
    }
    outputTree_->Fill();
  }
}

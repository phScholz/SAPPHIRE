#include "DecayResults.h"
#include "DecayProduct.h"
#include "NuclearMass.h"
#include <fstream>
#include <iostream>
#include <chrono>
#include <thread>
#include "omp.h"
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>



DecayResults::DecayResults(int Z, int A, double initialEnergyLow, double initialEnergyHigh, int suffixNo) :
                          initialZ_(Z),initialA_(A), initialJ_(-1.0), initialPi_(0),initialEnergyLow_(initialEnergyLow), initialEnergyHigh_(initialEnergyHigh) {
  
  /**
   * 1. Create File Name
   */
  std::string filename = CreateFileName(suffixNo);
  std::cout << "Output filename is " << filename << std::endl;

  /** 2. Open outFile as TFILE*/
  outputFile_ = new TFile(filename.c_str(),"recreate");

  /** 3. Create a new TTree "statDecay"*/
  outputTree_ = new TTree("statDecay","Sapphire Results");

  /** 4. Creating histograms*/
  CreateHists();

  /** 5. Creating branches for all decay results*/   
  CreateBranches();
}

DecayResults::DecayResults(int Z, int A, double J, int Pi, double initialEnergyLow, double initialEnergyHigh, int suffixNo) :
                          initialZ_(Z),initialA_(A),initialPi_(Pi),initialJ_(J), initialEnergyLow_(initialEnergyLow), initialEnergyHigh_(initialEnergyHigh) {
  
  std::string filename = CreateFileName(suffixNo);

  std::cout << "Output filename is " << filename << std::endl;

  /** 2. Open outFile as TFILE*/
  outputFile_ = new TFile(filename.c_str(),"recreate");

  /** 3. Create a new TTree "statDecay"*/
  outputTree_ = new TTree("statDecay","Sapphire Results");

  /** 4. Creating histograms*/
  CreateHists();

  /** 5. Creating branches for all decay results*/   
  CreateBranches();
}

void DecayResults::CreateHists(){
  ggMatrix_ = new TH2F("ggMatrix", "2D Histo; ggMatrix", 
                        4000, 0, 16.0,
                        4000, 0, 16.0);
                        
  ngMatrix_ = new TH2F("ngMatrix", "2D Histo; ngMatrix", 
                        4000, 0, 16.0,
                        4000, 0, 16.0);

  pgMatrix_ = new TH2F("pgMatrix", "2D Histo; pgMatrix", 
                        4000, 0, 16.0,
                        4000, 0, 16.0);

  agMatrix_ = new TH2F("agMatrix", "2D Histo; agMatrix", 
                        4000, 0, 16.0,
                        4000, 0, 16.0);
  
  TSCMatrix_ = new TH2F("TSCMatrix", "2D Histo; TSCMatrix", 
                        4000, 0, 16.0,
                        4000, 0, 16.0);
  
  gammaEnergyHist_ = new TH1F("GammaEnergyHist", "1D Gamma Energy Hist", 8000, 0, 16.0);
  protonEnergyHist_ = new TH1F("ProtonEnergyHist", "1D Proton Energy Hist", 8000, 0, 16.0);
  neutronEnergyHist_ = new TH1F("NeutronEnergyHist", "1D Neutron Energy Hist", 8000, 0, 16.0);
  alphaEnergyHist_ = new TH1F("AlphaEnergyHist", "1D Alpha Energy Hist", 8000, 0, 16.0);
}

void DecayResults::CreateBranches(){
  outputTree_->Branch("Events",&event,"Events");
  outputTree_->Branch("initialEnergy",&initialEnergy_,"initialEnergy/D");
  outputTree_->Branch("initialSpin",&initialSpin_,"initialSpin/D");
  outputTree_->Branch("initialParity",&initialParity_,"initialParity/D");
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

std::string DecayResults::CreateFileName(int suffixNo){ 

  char filename[256];

  if(initialJ_>=0){
    /** Construct file name either with SuffixNumber or not*/
    char spin[10];
    if(initialPi_==1) sprintf(spin,"%.1f+",initialJ_);
    else sprintf(spin,"%.1f-",initialJ_);
    
    if(suffixNo!=0) {
      if(initialEnergyLow_==initialEnergyHigh_){
        sprintf(filename,"output/Sapphire_%d%s_J=%s_E=%.3f_%d.root", initialA_,NuclearMass::FindElement(initialZ_).c_str(), spin,initialEnergyLow_,suffixNo);
      }
      else{
        sprintf(filename,"output/Sapphire_%d%s_J=%s_E=%.3f-%.3f_%d.root", initialA_,NuclearMass::FindElement(initialZ_).c_str(), spin,initialEnergyLow_,initialEnergyHigh_,suffixNo);
      }
    } 
    else
    {
      bool validFileName=false;
      int i  = 0;

      while(!validFileName) {
        if(i==0) 
        {
	        if(initialEnergyLow_==initialEnergyHigh_)
          {
	          sprintf(filename,"output/Sapphire_%d%s_J=%s_E=%.3f.root", initialA_,NuclearMass::FindElement(initialZ_).c_str(), spin,initialEnergyLow_);
          }
	        else
          {
            sprintf(filename,"output/Sapphire_%d%s_J=%s_E=%.3f-%.3f.root", initialA_,NuclearMass::FindElement(initialZ_).c_str(), spin,initialEnergyLow_,initialEnergyHigh_);
          }
        } 
        else
        {
	        if(initialEnergyLow_==initialEnergyHigh_)
          {
            sprintf(filename,"output/Sapphire_%d%s_J=%s_E=%.3f_%d.root", initialA_,NuclearMass::FindElement(initialZ_).c_str(), spin,initialEnergyLow_,i+1);
          }	
          else
          {
            sprintf(filename,"output/Sapphire_%d%s_J=%s_E=%.3f-%.3f_%d.root", initialA_,NuclearMass::FindElement(initialZ_).c_str(),spin,initialEnergyLow_,initialEnergyHigh_,i+1);
          }
        }

        std::ifstream inTest(filename);
        if(!inTest)
        {
	          validFileName=true;
        } 
        else
        {
          inTest.close();
        }
        i++;
      }
    }
  }
  else{    
    
    if(suffixNo!=0) {
      if(initialEnergyLow_==initialEnergyHigh_){
        sprintf(filename,"output/Sapphire_%d%s_E=%.3f_%d.root", initialA_,NuclearMass::FindElement(initialZ_).c_str(), initialEnergyLow_,suffixNo);
      }
      else{
        sprintf(filename,"output/Sapphire_%d%s_E=%.3f-%.3f_%d.root", initialA_,NuclearMass::FindElement(initialZ_).c_str(), initialEnergyLow_,initialEnergyHigh_,suffixNo);
      }
    } 
    else
    {
      bool validFileName=false;
      int i  = 0;
    
      while(!validFileName) {
        if(i==0) 
        {
	        if(initialEnergyLow_==initialEnergyHigh_)
          {
	          sprintf(filename,"output/Sapphire_%d%s_E=%.3f.root", initialA_,NuclearMass::FindElement(initialZ_).c_str(), initialEnergyLow_);
          }
	        else
          {
            sprintf(filename,"output/Sapphire_%d%s_E=%.3f-%.3f.root", initialA_,NuclearMass::FindElement(initialZ_).c_str(), initialEnergyLow_,initialEnergyHigh_);
          }
        } 
        else
        {
	        if(initialEnergyLow_==initialEnergyHigh_)
          {
            sprintf(filename,"output/Sapphire_%d%s_J=E=%.3f_%d.root", initialA_,NuclearMass::FindElement(initialZ_).c_str(), initialEnergyLow_,i+1);
          }	
          else
          {
            sprintf(filename,"output/Sapphire_%d%s_J=E=%.3f-%.3f_%d.root", initialA_,NuclearMass::FindElement(initialZ_).c_str(), initialEnergyLow_,initialEnergyHigh_,i+1);
          }
        }

        std::ifstream inTest(filename);
        if(!inTest)
        {
	          validFileName=true;
        } 
        else
        {
          inTest.close();
        }
        i++;
      }
    }
  }

  std::string fileN(filename);
  return fileN;
  
}

void DecayResults::WriteNCloseFile(){
  outputFile_->Write();
  outputFile_->Close();
}

DecayResults::~DecayResults() {}

void DecayResults::AddResults(std::vector<std::pair<DecayData, std::vector<DecayProduct> > >& results) {
  
  for(int i = 0;i<results.size();i++) {
    if(results[i].second.size()==0) continue;
    event.Data = results[i].first;
    event.Products = results[i].second;
    
    initialEnergy_ = results[i].first.energy();
    initialParity_ = results[i].first.parity();
    initialSpin_ = results[i].first.spin();
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
        gammaEnergyHist_->Fill(it->particleEnergy_);
        ++numGammas_;
      } else if(it->particleType_==1) {
	      neutronStepIndex_[numNeutrons_] = numSteps_;
	      neutronEnergy_[numNeutrons_]=it->particleEnergy_;
	      neutronMomentumX_[numNeutrons_]=it->particleMomentumX_;
	      neutronMomentumY_[numNeutrons_]=it->particleMomentumY_;
	      neutronMomentumZ_[numNeutrons_]=it->particleMomentumZ_;
        neutronEnergyHist_->Fill(it->particleEnergy_);
	      ++numNeutrons_;
      } else if(it->particleType_==2) {
	      protonStepIndex_[numProtons_] = numSteps_;
	      protonEnergy_[numProtons_]=it->particleEnergy_;
	      protonMomentumX_[numProtons_]=it->particleMomentumX_;
	      protonMomentumY_[numProtons_]=it->particleMomentumY_;
	      protonMomentumZ_[numProtons_]=it->particleMomentumZ_;
        protonEnergyHist_->Fill(it->particleEnergy_);
	      ++numProtons_;
      } else if(it->particleType_==3) {
	      alphaStepIndex_[numAlphas_] = numSteps_;
	      alphaEnergy_[numAlphas_]=it->particleEnergy_;
	      alphaMomentumX_[numAlphas_]=it->particleMomentumX_;
	      alphaMomentumY_[numAlphas_]=it->particleMomentumY_;
	      alphaMomentumZ_[numAlphas_]=it->particleMomentumZ_;
        alphaEnergyHist_->Fill(it->particleEnergy_);
	      ++numAlphas_;
      }
      ++numSteps_;
    }
    
    for(int i =0; i < numGammas_; i++){
      for(int j=i+1; j <numGammas_; j++){
        ggMatrix_->Fill(gammaEnergy_[i], gammaEnergy_[j], 1.);
        TSCMatrix_->Fill(gammaEnergy_[i], gammaEnergy_[i]+gammaEnergy_[j], 1.);
        TSCMatrix_->Fill(gammaEnergy_[j], gammaEnergy_[i]+gammaEnergy_[j], 1.);
      }
      
      for(int j=0; j < numNeutrons_; j++){
        ngMatrix_->Fill(gammaEnergy_[i], neutronEnergy_[j], 1.);
      }

      for(int j=0; j < numProtons_; j++){
        pgMatrix_->Fill(gammaEnergy_[i], protonEnergy_[j], 1.);
      }

      for(int j=0; j < numAlphas_; j++){
        agMatrix_->Fill(gammaEnergy_[i], alphaEnergy_[j], 1.);
      }
    }

    outputTree_->Fill();    
  }
}

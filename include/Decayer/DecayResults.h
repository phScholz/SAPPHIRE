/**
 * @file Decayer/DecayResults.h
 * @brief Declaration of DecayResults and DecayEvent
 */
#pragma once
#include <TTree.h>
#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TObject.h>
#include <vector>
#include "Decayer/DecayProduct.h"
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

/**
 * @brief Class to handle the results of the decayer module
 */
class DecayResults {
 public:
 /**
  * @brief Constructor for DecayResults for a single resonance
  * @param Z  charge number
  * @param A  mass number
  * @param J  initial spin
  * @param Pi initial parity
  * @param energyLow lower limit of initial energy
  * @param energyHigh upper limit of initial energy
  * @param suffixNo Suffix Number for outputfile
  */
  DecayResults(int Z, int A, double J, int Pi, double energyLow, double energyHigh, int suffixNo);

  /**
  * @brief Constructor for DecayResults for the decay of a distribution
  * @param Z  charge number
  * @param A  mass number
  * @param energyLow lower limit of initial energy
  * @param energyHigh upper limit of initial energy
  * @param suffixNo Suffix Number for outputfile
  */
  DecayResults(int Z, int A, double energyLow, double energyHigh, int suffixNo);



  /**
   * @brief Destructor for DecayResults
   */
  ~DecayResults();


  void AddResults(std::vector<std::pair<DecayData, std::vector<DecayProduct> > >&);
  
  void WriteNCloseFile();


private:
  /**
  * @brief Function to create the output file name.
  * @return A std::string of the file name
  */
  std::string CreateFileName(int suffixNo);

  /**
   * @brief Create branches in outputTree_
   * 
   */
  void CreateBranches();

  /**
   * @brief Create histograms in outputTree_
   * 
   */
  void CreateHists();

 private:
  int initialZ_;
  int initialA_;
  int initialPi_;
  double initialJ_;
  double initialEnergyLow_;
  double initialEnergyHigh_;
  DecayEvent event;
  TH2F* ngMatrix_;
  TH2F* pgMatrix_;
  TH2F* agMatrix_;
  TH2F* ggMatrix_;
  TH2F* TSCMatrix_;
  TH1F* gammaEnergyHist_;
  TH1F* neutronEnergyHist_;
  TH1F* protonEnergyHist_;
  TH1F* alphaEnergyHist_;
  TFile* outputFile_;
  TTree* outputTree_;
  Int_t Z_[100];
  Int_t A_[100];
  Int_t initialParity_;
  Int_t numNeutrons_;
  Int_t numGammas_;
  Int_t numProtons_;
  Int_t numAlphas_;
  Int_t numSteps_;
  Int_t protonStepIndex_[100];
  Int_t gammaStepIndex_[100];
  Int_t neutronStepIndex_[100];
  Int_t alphaStepIndex_[100];
  Double_t initialEnergy_;
  Double_t initialSpin_;
  Double_t fragmentEnergy_[100];
  Double_t fragmentExcitation_[100];
  Double_t fragmentMomentumX_[100];
  Double_t fragmentMomentumY_[100];
  Double_t fragmentMomentumZ_[100];
  Double_t protonTotalWidth_;
  Double_t neutronTotalWidth_;
  Double_t gammaTotalWidth_;
  Double_t alphaTotalWidth_;
  Double_t neutronEntranceWidth_;
  Double_t alphaEntranceWidth_;
  Double_t gammaEntranceWidth_;
  Double_t protonEntranceWidth_;
  Double_t protonEnergy_[100];
  Double_t protonMomentumX_[100];
  Double_t protonMomentumY_[100];
  Double_t protonMomentumZ_[100];
  Double_t gammaEnergy_[100];
  Double_t gammaMomentumX_[100];
  Double_t gammaMomentumY_[100];
  Double_t gammaMomentumZ_[100];
  Double_t neutronEnergy_[100];
  Double_t neutronMomentumX_[100];
  Double_t neutronMomentumY_[100];
  Double_t neutronMomentumZ_[100];
  Double_t alphaEnergy_[100];
  Double_t alphaMomentumX_[100];
  Double_t alphaMomentumY_[100];
  Double_t alphaMomentumZ_[100];
};

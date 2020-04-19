#ifndef DECAYRESULTS_H
#define DECAYRESULTS_H
#pragma once
#include <TTree.h>
#include <TFile.h>
#include <vector>

class DecayProduct;
class DecayData;

class DecayResults {
 public:
  DecayResults(int, int, double, int, double, double, int);
  ~DecayResults();
  void AddResults(std::vector<std::pair<DecayData, std::vector<DecayProduct> > >&);
 private:
  int initialZ_;
  int initialA_;
  int initialPi_;
  double initialJ_;
  double initialEnergyLow_;
  double initialEnergyHigh_;
  TFile* outputFile_;
  TTree* outputTree_;
  Int_t Z_[100];
  Int_t A_[100];
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
  Double_t gammaEnergy_[16000];
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

#endif

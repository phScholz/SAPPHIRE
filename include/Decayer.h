#ifndef DECAYER_H
#define DECAYER_H
#pragma once
#include <vector>
#include <cstdlib>
#include "SpinRatePair.h"
#include "TransitionRateFunc.h"

class CDFEntry {
 public:
 CDFEntry(int pairIndex,double energy,double value) :
  pairIndex_(pairIndex), energy_(energy), value_(value) {};
  int pairIndex_;
  double energy_;
  double value_;
};

class Level;

class Decayer{

  public:
    Decayer(int Z, int A, double jInitial, 
  	  int piInitial, double energy, double totalWidthForCorrection = 0.,
  	  double uncorrTotalWidthForCorrection = 0.,
  	  double uncorrTotalWidthSqrdForCorrection = 0.,
  	  Decayer* widthCorrectedDecayer=NULL);
    ~Decayer();
    bool Decay(int&,int&,double&,int&,double&,double&);
    void PrintFunctions();
    void PrintCDF();
    void CorrectWidthFluctuations();
    static void SetCrossSection(bool isCrossSection) {isCrossSection_=isCrossSection;};
    static void SetMaxL(double maxL) {maxL_=maxL;};
    static double GetMaxL() {return maxL_;};
    friend class CrossSection;
    double NeutronEntranceWidth() const {
      return neutronEntrance_;
    }
    double ProtonEntranceWidth() const {
      return protonEntrance_;
    }
    double AlphaEntranceWidth() const {
      return alphaEntrance_;
    }
    double GammaEntranceWidth() const {
      return gammaEntrance_;
    }
    double GammaTotalWidth() const {
      return gammaTotalWidth_;
    }
    double NeutronTotalWidth() const {
      return neutronTotalWidth_;
    }
    double AlphaTotalWidth() const {
      return alphaTotalWidth_;
    }
    double ProtonTotalWidth() const {
      return protonTotalWidth_;
    }

  private:
    void BuildCDF();
    bool BuildKnownCDF(int,std::vector<Level>&);
   private:
    static bool isCrossSection_;
    static double maxL_;
    int Z_;
    int A_;
    int piInitial_;
    double jInitial_;
    double energy_;
    double totalIntegral_;
    double totalIntegralSqrd_;
    double totalWidthForCorrection_;
    double uncorrTotalWidthForCorrection_;
    double uncorrTotalWidthSqrdForCorrection_;
    double neutronEntrance_;
    double protonEntrance_;
    double gammaEntrance_;
    double alphaEntrance_;
    double neutronTotalWidth_;
    double alphaTotalWidth_;
    double protonTotalWidth_;
    double gammaTotalWidth_;
    std::vector<SpinRatePair> spinRatePairs_;
    std::vector<CDFEntry> cdf_;
    Decayer* widthCorrectedDecayer_;
};


#endif

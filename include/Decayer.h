#ifndef DECAYER_H
#define DECAYER_H
#pragma once
#include <vector>
#include <cstdlib>
#include "SpinRatePair.h"
#include "TransitionRateFunc.h"
#include "SapphireInput.h"

/**
 * 
 */
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
    /**
     * @brief Constructor of the Decayer class
     * @param Z Charge/atomic number of the decaying nucleus
     * @param A Mass number of the decaying nucleus
     * @param jInitial The initial spin of the decaying nucleus
     * @param piInitial The initial parity of the decaying nucleus
     * @param energy The excitation energy of the decaying nucleus
     * @param totalWidthForCorrection By default 0.
     * @param uncorrTotalWidthForCorrection By default 0.
     * @param uncorrTotalWidthSqrdForCorrection By default 0.
     * @param widthCorrectedDecayer By default NULL.
     * @details 
     * 1. The following member variables will be initialized in the constructor:
     *  + Z_
     *  + A_
     *  + piInitial_
     *  + jInitial_
     *  + energy_
     *  + totalWidthForCorrection_
     *  + uncorrTotalWidthForCorrection_
     *  + uncorrTotalWidthSqrdForCorrection_
     *  + widthCorrectedDecayer_
     * 2. Calling of InitializeWidths() to initialize entrance and total widths for neutron, gamma, proton, and alpha are set to 0.
     * 3. 
     */
    Decayer(int Z, int A, double jInitial, 
  	  int piInitial, double energy, double totalWidthForCorrection = 0.,
  	  double uncorrTotalWidthForCorrection = 0.,
  	  double uncorrTotalWidthSqrdForCorrection = 0.,
  	  Decayer* widthCorrectedDecayer=NULL);
    
    /**
     * @brief Constructor of Decayer on the basis of a SapphireInput onject
     * @param input A SapphireInput object with input-parameters
     * @details See Decayer()
     */
    Decayer(SapphireInput & input);
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
    /**
     * @brief Initializes the member variables for partial and total widths.
     */
    void InitializeWidths();
    void BuildCDF();
    /**
     * @brief Building the cummulative distribution function (CDF) for known levels of the decaying nucleus.
     * @param levelIndex The level number of a bound level until which the levelscheme should be built.
     * @param knownLevels A std::vector of Level entries of the known levels of the decaying nucleus
     * @return True or false
     */
    bool BuildKnownCDF(int levelIndex,std::vector<Level>& knownLevels);
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

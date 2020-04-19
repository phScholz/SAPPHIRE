#pragma once

#include <vector>

class PreEqTransitionRateFunc;

class PreEqSpinRatePair {
 public:
 PreEqSpinRatePair(int neutronNumber, int neutronHoleNumber, 
		   int protonNumber, int protonHoleNumber,
		   int Z, int A, double spin, int parity, double qValue,
		   PreEqTransitionRateFunc* rateFunc, double integral) :
  neutronNumber_(neutronNumber), neutronHoleNumber_(neutronHoleNumber),
    protonNumber_(protonNumber), protonHoleNumber_(protonHoleNumber),
    Z_(Z), A_(A), rateFunc_(rateFunc),parity_(parity),spin_(spin),
    qValue_(qValue),integral_(integral) {};
  PreEqTransitionRateFunc* rateFunc_;
  int neutronNumber_;
  int neutronHoleNumber_;
  int protonNumber_;
  int protonHoleNumber_;
  int Z_;
  int A_;
  int parity_;
  double spin_;
  double qValue_;
  double integral_;
};

class PreEqCDFEntry {
 public:
 PreEqCDFEntry(int pairIndex,double energy,double value) :
  pairIndex_(pairIndex), energy_(energy), value_(value) {};
  int pairIndex_;
  double energy_;
  double value_;
};

class PreEqDecayer {
 public:
  PreEqDecayer(int initialNeutronNumber,
	       int initialNeutronHoleNumber,
	       int initialProtonNumber,
	       int initialProtonHoleNumber,
	       int Z, int A, double jInitial,
	       int piInitial, double energy);
  ~PreEqDecayer();
  bool Decay(int&, int&, double& , int& ,
	     int&, int&, int&, int&,
	     double&, double&);
  static void SetCrossSection(bool isCrossSection) {isCrossSection_=isCrossSection;};
  static void SetMaxL(double maxL) {maxL_=maxL;};
  static double GetMaxL() {return maxL_;};
  void PrintCDF();
 private: 
  void BuildCDF();
 private:
  static bool isCrossSection_;
  static double maxL_;
  int initialNeutronNumber_;
  int initialNeutronHoleNumber_;
  int initialProtonNumber_;
  int initialProtonHoleNumber_;
  int initialExcitonNumber_;
  int Z_;
  int A_;
  int piInitial_;
  double jInitial_;
  double energy_;
  double totalIntegral_;
  std::vector<PreEqSpinRatePair> spinRatePairs_;
  std::vector<PreEqCDFEntry> cdf_;
};


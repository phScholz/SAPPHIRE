#ifndef DECAYCONTROLLER_H
#define DECAYCONTROLLER_H

#include <vector>
#include "Decayer.h"
#include "DecayProduct.h"

class TVector3;

class DecayController {
 public:
 DecayController(int Z, int A, double jInitial, int piInitial, double energy,
		 int initialNeutronNumber = -1, int initialNeutronHoleNumber = -1,
		 int initialProtonNumber = -1, int initialProtonHoleNumber = -1) :
  Z_(Z), A_(A), piInitial_(piInitial),
    initialNeutronNumber_(initialNeutronNumber), initialNeutronHoleNumber_(initialNeutronHoleNumber),
    initialProtonNumber_(initialProtonNumber), initialProtonHoleNumber_(initialProtonHoleNumber),
    jInitial_(jInitial), energy_(energy) {};
  bool Decay(double&,double&,double&,double&,double&,double&,double&,double&);
 std::vector<DecayProduct> DecayProducts() const { return decayProducts_; };
 void PrintDecays();
 private:
 void CalcKinematics(int,int,int,int,double,double,int,double,TVector3&);
 private:
  static const double thresholdEnergy_=0.05;
  int Z_;
  int A_;
  int piInitial_;
  int initialNeutronNumber_;
  int initialNeutronHoleNumber_;
  int initialProtonNumber_;
  int initialProtonHoleNumber_;
  double jInitial_;
  double energy_;
  std::vector<DecayProduct> decayProducts_;
};

#endif

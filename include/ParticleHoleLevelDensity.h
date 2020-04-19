#pragma once

class ParticleHoleLevelDensity {
 public:
  ParticleHoleLevelDensity(int, int, double, int, int, int, int);
  double operator()(double energy,bool correct=true, bool spin=false);
  double PauliCorrection(int,int,int,int);
  double PairingCorrection(double energy);
  double gNu() const {return gNu_;};
  double gPi() const {return gPi_;};
  double FiniteDepth(double);
 private:
  bool invalid_;
  int totalExcitonNumber_;
  int totalHoleNumber_;
  int incident_;
  int A_;
  int Z_;
  double separationEnergy_;
  double preFactor_;
  double pauliCorrection_;
  double pairingEnergy_;
  double nCritical_;
  double spinFactor_;
  double gPi_;
  double gNu_;
};



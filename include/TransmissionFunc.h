#ifndef TRANSMISSIONFUNC_H
#define TRANSMISSIONFUNC_H

class TransmissionFunc {
 public:
  TransmissionFunc (int z2, int m2, double jInitial, int piInitial,
		    double jFinal, int piFinal, double maxL,
		    double totalWidthForCorrection,
		    double uncorrTotalWidthForCorrection,
		    double uncorrTotalWidthSqrdForCorrection,
		    TransmissionFunc* previous) :
  z2_(z2),m2_(m2),piInitial_(piInitial),piFinal_(piFinal),
    jInitial_(jInitial),jFinal_(jFinal), maxL_(maxL),
    totalWidthForCorrection_(totalWidthForCorrection),
    uncorrTotalWidthForCorrection_(uncorrTotalWidthForCorrection),
    uncorrTotalWidthSqrdForCorrection_(uncorrTotalWidthSqrdForCorrection),
    previous_(previous) {};  
  virtual ~TransmissionFunc() {};
  virtual double operator()(double) = 0;
  virtual bool IsValid() = 0;
 protected:
  int z2_;
  int m2_;
  int piInitial_;
  int piFinal_;
  double jInitial_;
  double jFinal_;
  double maxL_;
  double totalWidthForCorrection_;
  double uncorrTotalWidthForCorrection_;
  double uncorrTotalWidthSqrdForCorrection_;
  TransmissionFunc* previous_;
};

#endif

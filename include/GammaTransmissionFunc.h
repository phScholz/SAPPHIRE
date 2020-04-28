/**
 * @file GammaTransmissionFunc.h
 * @date 2020-04-27
 * @brief Declaration of the base class for Gamma Strength Functions
 */

#include "TransmissionFunc.h"
#include "Constants.h"
#include "NuclearMass.h"
#include "LevelDensity/LevelDensityFormula.h"
#include <math.h>

#pragma once
class GDRParameters {
 public:
  GDRParameters() : eta_(0.) {
    E_[0]=0.;
    E_[1]=0.;
    W_[0]=0.;
    W_[1]=0.;
  };
 GDRParameters(double eta, double E1, double W1, double kSigmaGamma1,
	       double E2, double W2, double kSigmaGamma2) : 
   eta_(eta) {
    E_[0]=E1;
    E_[1]=E2;
    W_[0]=W1;
    W_[1]=W2;
    kSigmaGamma_[0] = kSigmaGamma1;
    kSigmaGamma_[1] = kSigmaGamma2;
  };
  double eta_;
  double E_[2];
  double W_[2];
  double kSigmaGamma_[2];
};

typedef std::tr1::unordered_map<MassKey, GDRParameters > GDRTable;

/**
 * @brief Base class for gamma-ray strength functions. Child of TransmissionFunc
 */
class GammaTransmissionFunc : public TransmissionFunc {
 public:
  GammaTransmissionFunc(int,int,double,int,double,int,
			double,double,double,double,
			TransmissionFunc*); 
  virtual ~GammaTransmissionFunc() {};
  static GammaTransmissionFunc* 
    CreateGammaTransmissionFunc(int,int,double,int,double,int,
				double,LevelDensityFormula*,double,double,double,
				TransmissionFunc*,double); 
  bool IsValid() {
    return true;
  };
  double operator()(double);
  virtual double CalcStrengthFunction(double) = 0;
  static void InitializeGDRParameters(std::string);
  static void SetEGDRType(int);
  static void SetPorterThomas(bool);
 protected:
  static GDRTable gdrTable_;
  static int egdrType_;
  static const int mgdrType_=0;
  static const int egqrType_=0;
  static bool porterThomas_;
  GDRParameters gdrParameters_;
};


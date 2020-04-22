#pragma once
#include "TransitionRateFunc.h"

/**
 * @brief 
 */
class SpinRatePair {
 public:
 SpinRatePair(int Z, int A, double spin, int parity, double qValue,
	      TransitionRateFunc* rateFunc, double integral){
          Z_=Z;
          A_=A;
          parity_=parity;
          spin_=spin;
          qValue_=qValue;
          rateFunc_=rateFunc;
          integral_=integral;
        };
 TransitionRateFunc* rateFunc_;
 int Z_;
 int A_;
 int parity_;
 double spin_;
 double qValue_;
 double integral_;
};
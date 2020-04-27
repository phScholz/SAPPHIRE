/**
 * @file SpinRatePair.h
 * @date 2020-04-21
 * @brief Declaration of the class SpinRatePair
 */
#pragma once
#include "TransitionRateFunc.h"

/**
 * @brief  Class for a SpinRatePair
 * @details Combines the properties of the state with its transition rate
 */
class SpinRatePair {
 public:

 /**
  * @brief Contructor of a SpinRatePair which initializes its properties
  */
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

 TransitionRateFunc* rateFunc_; /**< Pointer to a TransitionRateFunc object*/
 int Z_; /**< Charge number of the nucleus*/
 int A_; /**< Mass number of the nucleus*/
 int parity_; /**< parity of the state*/
 double spin_; /**< Spin of the state*/
 double qValue_; /**< Qvalue for the decay*/
 double integral_; /**< Integral over the TransitionRateFunc*/
};
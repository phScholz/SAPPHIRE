/**
 * @file SLPair.h
 * @brief Declaration of the SLPair
 */
#pragma once
/**
 * @brief Class for a pair of spin and angular momentum
 */
class SLPair {
 public: 
  SLPair(double s, int l) : s_(s), l_(l) {};
  bool operator<(const SLPair& right) const {
    if(s_ == right.s_ ) return l_ < right.l_;
    else return s_ < right.s_;
  };
  double s_;
  int l_;
};
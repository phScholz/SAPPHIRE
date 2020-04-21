#pragma once
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
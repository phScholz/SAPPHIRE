/**
 * @file CDFEntry.h
 * @brief Declaration of class CDFEntry
 * @date 2020-04-22
 */
#pragma once
class CDFEntry {
 public:
  /**
   * @brief Construct a CDFEntry object
   * @param pairIndex
   * @param energy
   * @param value
   */
  CDFEntry(int pairIndex,double energy,double value) :
  pairIndex_(pairIndex), energy_(energy), value_(value) {};
  int pairIndex_; /**< The index of a SpinRatePair*/
  double energy_; /**< The energy of a transition*/
  double value_; /**< The CDF value for the transition*/
};
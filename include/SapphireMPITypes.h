/* #ifndef SAPPHIREMPITYPES_H
#define SAPPHIREMPITYPES_H

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/mpi/datatype.hpp>
#include <vector>
#include "Decayer/DecayProduct.h"

enum SapphireTags_t {SapphireTagProcess,SapphireTagDone,SapphireTagResults};

class InitialNucleusData {
 private:
  friend class boost::serialization::access;
  
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & preEq_;
    ar & Z_;
    ar & A_;
    ar & Pi_;
    ar & J_;
    ar & lowEnergy_;
    ar & highEnergy_;
  }
  bool preEq_;
  int Z_;
  int A_;
  int Pi_;
  double J_;
  double lowEnergy_;
  double highEnergy_;
 public:
  InitialNucleusData(){};
 InitialNucleusData(int Z, int A, double J, int Pi, double lowEnergy, double highEnergy, bool preEq) :
  preEq_(preEq), Z_(Z),A_(A),Pi_(Pi),J_(J),lowEnergy_(lowEnergy),highEnergy_(highEnergy) {};
  bool preEq() const { return preEq_; } ;
  int Z() const { return Z_; } ;
  int A() const { return A_; } ;
  int Pi() const { return Pi_; } ;
  double J() const { return J_; } ;
  double lowEnergy() const { return lowEnergy_; } ;
  double highEnergy() const { return highEnergy_; } ;
};
#endif */
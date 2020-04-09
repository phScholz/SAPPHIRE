#ifndef STATDECAYMPITYPES_H
#define STATDECAYMPITYPES_H

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/mpi/datatype.hpp>
#include <vector>
#include "DecayProduct.h"

enum SapphireTags_t {SapphireTagProcess,SapphireTagDone,SapphireTagResults};

class InitialNucleusData {
 private:
  friend class boost::serialization::access;
  
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
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

BOOST_IS_MPI_DATATYPE(InitialNucleusData);

namespace boost {
namespace serialization {

void serialize(Archive & ar, DecayData& g, const unsigned int version) {
    ar & g.energy_;
    ar & g.neutronEntranceWidth_;
    ar & g.protonEntranceWidth_;
    ar & g.alphaEntranceWidth_;
    ar & g.gammaEntranceWidth_;
    ar & g.neutronTotalWidth_;
    ar & g.protonTotalWidth_;
    ar & g.alphaTotalWidth_;
    ar & g.gammaTotalWidth_;
}

void serialize(Archive & ar, DecayProduct& g, const unsigned int version)
{
  ar &  g.Z_;
  ar &  g.A_;
  ar &  g.Pi_;
  ar &  g.particleType_;
  ar &  g.J_;
  ar &  g.excitationEnergy_;
  ar &  g.fragmentEnergyCM_;
  ar &  g.fragmentEnergy_;
  ar &  g.fragmentMomentumX_;
  ar &  g.fragmentMomentumY_;
  ar &  g.fragmentMomentumZ_;
  ar &  g.particleThetaCM_;
  ar &  g.particlePhiCM_;
  ar &  g.particleEnergyCM_;
  ar &  g.particleEnergy_;
  ar &  g.particleMomentumX_;
  ar &  g.particleMomentumY_;
  ar &  g.particleMomentumZ_;
}
} // namespace serialization
} // namespace boost

BOOST_IS_MPI_DATATYPE(DecayData);
BOOST_IS_MPI_DATATYPE(DecayProduct);

#endif

/* 
#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/mpi/datatype.hpp>
#include <vector>
#include "DecayProduct.h"
#include "SapphireMPITypes.h"

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

 */
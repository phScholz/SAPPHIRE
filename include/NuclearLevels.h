/**
 * @file NuclearLevels.h
 * @author Mary Beard, Philipp Scholz
 * @date 2020
 * @brief It contains classes and methods, handling the nuclear level scheme
 */

#ifndef NUCLEARLEVELS_H
#define NUCLEARLEVELS_H

#include "NuclearMass.h"
#include <vector>
#include <string>

class GammaTransition {
 public:
 GammaTransition(int levelIndex, double energy, double probability) :
  levelIndex_(levelIndex), energy_(energy), probability_(probability) {};
 int levelIndex_;
 double energy_;
 double probability_;
};

class Level {
 public:
 Level(double J, int Pi, double energy) :
  Pi_(Pi), J_(J), energy_(energy) {};
 int Pi_;
 double J_;
 double energy_;
 std::vector<GammaTransition> gammas_;
};

/**
 * @class LevelsContainer
 * @brief A container class for nuclear levels
 */
class LevelsContainer {
 public:
  LevelsContainer() {};                   /**< Simple constructor */
  LevelsContainer(std::istream&,int,int); /**< Constructructor from levels file */
  std::vector<Level> levels_;             /**< One std::vector object which for levels */
};

typedef std::tr1::unordered_map<MassKey, LevelsContainer> LevelsTable; /**< One map object which maps MassKey to LevelsContainer */

class NuclearLevels {
 public:
  static void InitializeLevels(std::string levelsDirectory,
			       std::string spinFile);
  static void PrintLevels(int,int);
  static std::vector<Level> FindLevels(int,int);
 private:
  static LevelsTable levelsTable_;
};
#endif

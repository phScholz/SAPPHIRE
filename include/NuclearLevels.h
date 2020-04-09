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

class LevelsContainer {
 public:
  LevelsContainer() {};
  LevelsContainer(std::istream&,int,int);
  std::vector<Level> levels_;
};

typedef std::tr1::unordered_map<MassKey, LevelsContainer> LevelsTable;

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

/**
 * @file Databases/NuclearLevels.h
 * @date 2020-04-21
 * @brief It contains classes and methods, handling the nuclear level scheme: GammaTransition, Level, LevelsContainer and NuclearLevels
 */

#pragma once
#include "Databases/NuclearMass.h"
#include <vector>
#include <string>

/**
 * @class GammaTransition
 * @brief Class which represents a gamma-ray transition, with probability, energy and a level index.
 */
class GammaTransition {
 public:
 GammaTransition(int levelIndex, double energy, double probability) :
  levelIndex_(levelIndex), energy_(energy), probability_(probability) {};
 int levelIndex_;
 double energy_;
 double probability_;
};

/**
 * @class Level
 * @brief Class which represents a nuclear level with spin, parity, energy, and gamma transitions.
 */
class Level {
 public:
 Level(double J, int Pi, double energy) :
  Pi_(Pi), J_(J), energy_(energy) {};
  
  /**
   * @brief Operator to contruct a less comparison between level energies
   */
    bool operator < (const Level& str) const
    {
        return (energy_ < str.energy_);
    }
  
  /**
   * @brief Operator to contruct a greater comparison between level energies
   */
    bool operator > (const Level& str) const
    {
        return (energy_ > str.energy_);
    }

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
  
  /**
   * @brief LevelsContainer constructer as called by NuclerLevels::InitializeLevels()
   * @param in  The reference to a std::ifstream of the isotopeFile at the beginning of the levels list.
   * @param numLevels The number of levels which follow.
   * @param numComplete The number of levels which have complete information in the level scheme.
   * @details
   * 
   * When constructed like that, LevelsContainer() will read numLevels level from the std::ifstream and
   * for each will extract information about the energy, spind and parity and number of gammas.
   * Then it will proceed to read the lines for the gammas and parse information about levelindex, energy and probability.
   * The levels, read in as class Level, belong then to the LevelsContainer which is used for the LevelsTable
   * NuclearLevels::levelsTable_.
   * 
   * For the structure of the levels files please look at:
   * - [RIPL-3 Website](https://www-nds.iaea.org/RIPL-3/)
   * - [README RIPL-3 2002](https://www-nds.iaea.org/RIPL-3/levels/levels.html)
   * - [README RIPL-3 2015](https://www-nds.iaea.org/RIPL-3/levels/levels-readme.html)
   * - R. Capote et al., *RIPL – Reference Input Parameter Library for Calculation of Nuclear Reactions and Nuclear Data Evaluations*, Nuclear Data Sheets **110** (2009) 3107, DOI: [10.1016/j.nds.2009.10.004](https://doi.org/10.1016/j.nds.2009.10.004).
   * 
   */
  LevelsContainer(std::istream& in,int numLevels,int numComplete); /**< Constructructor from levels file */
  std::vector<Level> levels_;             /**< One std::vector object which for levels */
};

typedef std::tr1::unordered_map<MassKey, LevelsContainer> LevelsTable; /**< One map object which maps MassKey to LevelsContainer */

/**
 * @brief The NuclearLevels class represents the levels information of all nuclei.
 */
class NuclearLevels {
 public:
 /**
  * @brief Read in levels from the RIPL files in /levels into levelsTable_.
  * @param levelsDirectory std::string with the path to the directory containing the RIPL-3 level files.
  * @param spinFile std::string with the path to the custon spinFile which is used to estimate the groundstate spin of unknown nuclei.
  * @details 
  * 
  * The known level scheme is read here from RIPL-3 files of known nuclear levels. 
  * For the structure of the levels files please look at:
  * - [RIPL-3 Website](https://www-nds.iaea.org/RIPL-3/)
  * - [README RIPL-3 2002](https://www-nds.iaea.org/RIPL-3/levels/levels.html)
  * - [README RIPL-3 2015](https://www-nds.iaea.org/RIPL-3/levels/levels-readme.html)
  * - R. Capote et al., *RIPL – Reference Input Parameter Library for Calculation of Nuclear Reactions and Nuclear Data Evaluations*, Nuclear Data Sheets **110** (2009) 3107, DOI: [10.1016/j.nds.2009.10.004](https://doi.org/10.1016/j.nds.2009.10.004).
  */
  static void InitializeLevels(std::string levelsDirectory, std::string spinFile);

  /**
   * @brief Print all known levels of a given nucleus to std::cout.
   * @param Z Charge/atomic number of the nucleus
   * @param A Mass number of the nucleus
   */
  static void PrintLevels(int,int);

  /**
   * @brief Find the levels of a given nucleus in levelsTable_ and returns a std::vector<Level>.
   * @param Z Charge/atomic number of the nucleus.
   * @param A Mass number of the nucleus.
   * @return A std::vector<Level> which contains the levels of a nucleus as provided by the LevelsContainer.
   */
  static std::vector<Level> FindLevels(int,int);

 private:
  static LevelsTable levelsTable_; /**< A map which contains of a LevelsContainer for each MassKey.*/
};


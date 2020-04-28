/**
 * @file RandomScheme.h
 * @author Philipp Scholz <pscholz@outlook.de>
 * @date 2020-04-25
 * @brief File for the RandomScheme class
 * 
 */

#include "GammaTransmissionFunc.h"
#include "NuclearLevels.h"
#include "LevelDensityFormula.h"
#include "TransmissionFunc.h"
#include <vector>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"


/**
 * @brief A class to build or fill levelschemes randomly via LevelDensity and GammaStrengthFunctions
 * @todo Everything
 * - GetRandomScheme
 * - FillRandomScheme
 * - ExtendRandomScheme
 */
class RandomScheme{
    public:
        /**
         * @brief Constructor
         */
        RandomScheme();

        /**
         * @brief Constructor where the user can set the maxJ and eStep attributes.
         */
        RandomScheme(int maxJ, double eStep);

        /**
         * @brief Destructor
         */
        ~RandomScheme();
    
    public:
        /**
         * @brief Method to build a random level scheme for one isotope
         * @param Z Charge number
         * @param A Mass number
         * @param eStart Starting energy from which on the scheme should be built.
         * @param energy Upper energy limit for the level scheme in MeV
         * @return A std::vector<Level> containing the randomly built levels
         * @todo Everything
         */
        void CreateRandomScheme(int Z, int A, double eStart, double energy);

        /**
         * @brief Method to fill a random level scheme for one isotope
         * @param Z Charge number
         * @param A Mass number
         * @param eMin Lower energy limit for the level scheme in MeV
         * @param eMax Upper energy limit for the level scheme in MeV
         * @return A std::vector<Level> containing the randomly built levels
         * @todo Everything
         */
        std::vector<Level> FillRandomScheme(int Z, int A, double eMin, double eMax);

        /**
         * @brief Method to extend an existing level scheme for one isotope with random levels
         * @param Z Charge number
         * @param A Mass number
         * @param eMax Upper energy limit for the level scheme in MeV
         * @return A std::vector<Level> containing the randomly built levels
         * @details
         * 1. knownLevels are found via NuclearLevels::FindLevels()
         * 2. CreateRandomScheme() is called to start the extension of the level
         * scheme at the energy of the last known level.
         */
        void ExtendRandomScheme(int Z, int A, double eMax);


        /**
         * @brief print the randomly generated level scheme.
         * @param maxE Upper energy limit of the printed level scheme.
         * @exception Throws out of range exception if RandomScheme is empty.
         */
        void PrintRandomScheme(double maxE);

        /** Getter for randomScheme pointer*/
        std::vector<Level> * GetScheme(){
            return randomScheme;
        }

    private:
        /** Getter leveldensity at excitation energy**/
        double CalcLevelDensity(double energy) {
            return levelDensity_->operator()(energy);
        };

          /** Getter transmission coefficient for specific energy**/
        double CalcTransmissionFunc(double energy) {
            return transmissionFunc_->operator()(energy);
        };

        /**
         * @brief Create the spin and parity of a groundstate from the spinOD file
         * @param Z_ charge number of the isotope
         * @param A_ mass number of the isotope
         * @return A Level object with the groundstate
         */
        Level GetGroundState(int Z_, int A_);

        /**
         * @brief A method to draw levels from a Poisson distribution
         * @param Z charge number of the isotope
         * @param A mass number of the isotope
         * @param eStart Energy from which levels should be created
         * @param eStop Energy until which levels should be created
         * @details
         * 1. Setup of the random number generator using boost::random
         * 2. If start energy is 0.0 MeV, a groundstate is created via GetGroundState()
         * 3. Loop over energies and spins (even mass: integer spins, odd mass: half-numbered spins)
         *  - Define a Level density for Z, A, and J
         *  - Create a Poisson distribution with expectation value from level density model
         *  - Draw number of levels from Poisson distribution and create as many levels in the energy interval
         *  - Energies are drawn from a uniform distribution iver the energy interval
         *  - Push back levels to randomScheme
         */
        void CreateLevels(int Z, int A, double eStart, double eStop);

        /**
         * @brief Create Gamma Transitions between the levels in randomScheme
         * @param eStart Start energy from which gamma transitions should be inserted
         * @details
         * 1. In a loop from the highest level to the lowest, gamma transitions are created to each lower lying level,
         * if the spin difference is not larger than 2. The partial strength are drawn from a GSF model via GammaTransmissionFunc.
         * 2. The probabilites are normalized by summing up all partial strengths to a totalstrength and consequently
         * deviding all partial strengths by the total strength for each level above eStart.
         * @todo 
         * - Porter-Thomas distribution for partial strengths
         */
        void CreateGammaTransitions(double eStart);


        int nldModel_;                   /**< Integer to select the level density model*/
        int gsfModel_;                     /**< Integer to select the gamma strength model*/
        LevelDensityFormula* levelDensity_;          /**< LevelDensity model*/
        TransmissionFunc* transmissionFunc_; /**< Particle or GammaTransmissionFunc*/
        TransmissionFunc* previous; /**< Particle or previous*/
        std::vector<Level> *randomScheme;

        int Z_;
        int A_;
        int maxL_ =2;

        double tWFC_=0;
        double uTWFC_=0;
        double uTWSFC_=0;
        
        double maxJ_; /**< Double to select the maximum number of spin*/
        double maxE_; /**< Double for the maximum excitation energy*/
        double eStep_; /**< Energy binning*/
        bool verbose_=false;

};
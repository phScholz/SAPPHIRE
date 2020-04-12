/**
 * @file SapphireInput.h
 * @brief Option class for Sapphire.
 * @date 2020
 * @author Philipp Scholz, <pscholz@outlook.de>
 * 
 */

#ifndef MODULE_SAPPHIREINPUT_H
#define MODULE_SAPPHIREINPUT_H
#endif

#include <string>

class SapphireInput{
    public:
    /**
     * @brief Simple constructor
     */
        SapphireInput();
    
    /**
     * @brief Constructor, directly passing the input file
     * @param fileName Path to the inputFile Name
     */
        SapphireInput(std::string fileName);

    /**
     * @brief Simple destructor
     */
        ~SapphireInput();

    /**
     * @brief Initialize default values
     */
        void Initialize();

    
    private:
        bool calcRates;             /**< Bool if rates should be calculated*/
        bool calcAverageWidth;      /**< Bool if average widths should be calculated*/
        bool energyFile;            /**< Bool if an energy file has been given*/
        bool reactionFile;          /**< Bool if an reactionFile has been given*/
        bool residualGamma_;        /**< Bool if residual cross section for gamma should be calculated*/
        bool residualNeutron_;      /**< Bool if residual cross section for neutron should be calculated*/
        bool residualProton_;       /**< Bool if residual cross section for proton should be calculated*/
        bool residualAlpha_;        /**< Bool if residual cross section for alpha should be calculated*/
        bool calculateGammaCutoff_; /**< */

        int entranceState;          /**< Int for the number of level, which should be the entrance State*/

        std::string energyFile;     /**< String with the path to the energyFile*/
        std::string reactionFile;   /**< String with the path to the reactionFile*/





}

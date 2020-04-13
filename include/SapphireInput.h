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
         * @brief Initialize default values
         */
        void Initialize();
        
        /**
         * @brief Read ini file and set the variables
         */
        void ReadInputFile(std::string filename);

        //Setter
        void CalcRates(bool x){calcRates=x;}                         /**< Setter for bool calcRates*/
        void CalcAverageWidth(bool x){calcAverageWidth=x;}           /**< Setter for bool calcAverageWidth*/
        void ResiudalGamma(bool x){residualGamma=x;}                 /**< Setter for bool residualGamma*/
        void ResiudalNeutron(bool x){residualNeutron=x;}             /**< Setter for bool residualNeutron*/
        void ResiudalProton(bool x){residualProton=x;}               /**< Setter for bool residualProton*/
        void ResiudalAlpha(bool x){residualAlpha=x;}                 /**< Setter for bool residualAlpha*/
        void CalculateGammaCutoff(bool x){calculateGammaCutoff=x;}   /**< Setter for bool calculateGammaCutoff*/

        void EntranceState(int x){entranceState=x;}                  /**< Setter for entranceState*/

        void EnergyFile(std::string x){energyFile=x;}                       /**< Setter for string energyFile*/
        void ReactionFile(std::string x){reactionFile=x;}                   /**< Setter for string reactionFile*/

        
        //Getter
        bool CalcRates(){return calcRates;}                     /**<Getter for calcrates*/    
        bool CalcAverageWidht(){return calcAverageWidth;}       /**<Getter for calcAverageWidth*/    
        bool ResidualGamma(){return residualGamma;}             /**<Getter for residualGamma*/    
        bool ResidualNeutron(){return residualNeutron;}         /**<Getter for residualNeutron*/    
        bool ResidualProton(){return residualProton;}           /**<Getter for residualProton*/
        bool ResidualAlpha(){return residualAlpha;}             /**<Getter for residualAlpha*/
        bool CalculateGammaCutoff(){return calculateGammaCutoff;} /**<Getter for calculateGammaCutoff*/

        int EntranceState(){return entranceState;}              /**<Getter for entranceState*/ 

        std::string EnergyFile(){return energyFile;}            /**<Getter for energyFile*/ 
        std::string ReactionFile(){return reactionFile;}        /**<Getter for reactionFile*/

        
    private:
        bool calcRates;             /**< Bool if rates should be calculated*/
        bool calcAverageWidth;      /**< Bool if average widths should be calculated*/
        bool residualGamma;        /**< Bool if residual cross section for gamma should be calculated*/
        bool residualNeutron;      /**< Bool if residual cross section for neutron should be calculated*/
        bool residualProton;       /**< Bool if residual cross section for proton should be calculated*/
        bool residualAlpha;        /**< Bool if residual cross section for alpha should be calculated*/
        bool calculateGammaCutoff; /**< */

        int entranceState;          /**< Int for the number of level, which should be the entrance State*/

        std::string energyFile;     /**< String with the path to the energyFile*/
        std::string reactionFile;   /**< String with the path to the reactionFile*/
};

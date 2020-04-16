#pragma once

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
#include <vector>


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
         * @param filename File name of the ini-File
         */
        void ReadInputFile(std::string filename);
        
        /**
         * @brief Print InputFile to std::cout
         * @param filename File name of the ini-File
         */
        void printIntputFile(std::string filename);

        /**
         * @brief Print Inputparameters to std::cout
         */
        void printIntputParameters();

        //Setter
        void CalcRates(bool x){calcRates=x;}                         /**< Setter for bool calcRates*/
        void CalcAverageWidth(bool x){calcAverageWidth=x;}           /**< Setter for bool calcAverageWidth*/
       
        void ResidualGamma(bool x){residualGamma=x;}                 /**< Setter for bool residualGamma*/
        void ResidualNeutron(bool x){residualNeutron=x;}             /**< Setter for bool residualNeutron*/
        void ResidualProton(bool x){residualProton=x;}               /**< Setter for bool residualProton*/
        void ResidualAlpha(bool x){residualAlpha=x;}                 /**< Setter for bool residualAlpha*/
        
        void CalculateGammaCutoff(bool x){calculateGammaCutoff=x;}   /**< Setter for bool calculateGammaCutoff*/
        void PorterThomas_g(bool x){porterThomas_g=x;}                  /**<Setter for bool porterThomas_g*/
        void PorterThomas_p(bool x){porterThomas_p=x;}                  /**<Setter for bool porterThomas_p*/

        void EntranceState(int x){entranceState=x;}                  /**< Setter for entranceState*/
        void g_ExitStates(int x){g_exitStates=x;}              /**<Setter for g_exitStates*/ 
        void n_ExitStates(int x){n_exitStates=x;}              /**<Setter for n_exitStates*/ 
        void p_ExitStates(int x){p_exitStates=x;}              /**<Setter for p_exitStates*/ 
        void a_ExitStates(int x){a_exitStates=x;}              /**<Setter for a_exitStates*/ 

        void g_Formalism(int x){g_formalism=x;}              /**<Setter for g_formalism*/
        void n_Formalism(int x){n_formalism=x;}              /**<Setter for n_formalism*/
        void p_Formalism(int x){p_formalism=x;}              /**<Setter for p_formalism*/
        void a_Formalism(int x){a_formalism=x;}              /**<Setter for a_formalism*/
        
        void DecayerMaxL(double x){decayerMaxL=x;}              /**<Setter for decayerMaxL*/ 
        void PreEqMaxL(double x){preEqMaxL=x;}              /**<Setter for preEqMaxL*/ 
        void g_CutoffEnergy(double x){g_cutoffEnergy=x;}    /**<Setter for g_CutoffEnergy*/ 

        void EnergyFile(std::string x){energyFile=x;}                       /**< Setter for string energyFile*/
        void Energies(std::string x){energies=x;}                       /**< Setter for string energies*/
        void ReactionFile(std::string x){reactionFile=x;}                   /**< Setter for string reactionFile*/
        void Suffix(std::string x){suffix=x;}                               /**< Setter for string module*/
        void Reaction(std::string x){reaction=x;}                    /**<Setter for reaction*/
        
        void MassTable(std::string x){massTable=x;}                               /**< Setter for string module*/
        void GdrParams(std::string x){gdrParams=x;}                               /**< Setter for string module*/
        void LevelDir(std::string x){levelDir=x;}                               /**< Setter for string module*/
        void SpinFile(std::string x){spinFile=x;}                               /**< Setter for string module*/

        
        //Getter
        bool CalcRates(){return calcRates;}                     /**<Getter for calcrates*/    
        bool CalcAverageWidth(){return calcAverageWidth;}       /**<Getter for calcAverageWidth*/    
        bool ResidualGamma(){return residualGamma;}             /**<Getter for residualGamma*/    
        bool ResidualNeutron(){return residualNeutron;}         /**<Getter for residualNeutron*/    
        bool ResidualProton(){return residualProton;}           /**<Getter for residualProton*/
        bool ResidualAlpha(){return residualAlpha;}             /**<Getter for residualAlpha*/
        bool CalculateGammaCutoff(){return calculateGammaCutoff;} /**<Getter for calculateGammaCutoff*/
        bool PorterThomas_g(){return porterThomas_g;}             /**<Getter for bool porterThomas_g*/
        bool PorterThomas_p(){return porterThomas_p;}             /**<Getter for bool porterThomas_p*/

        int EntranceState(){return entranceState;}              /**<Getter for entranceState*/ 
        int g_ExitStates(){return g_exitStates;}              /**<Getter for g_exitStates*/ 
        int n_ExitStates(){return n_exitStates;}              /**<Getter for n_exitStates*/ 
        int p_ExitStates(){return p_exitStates;}              /**<Getter for p_exitStates*/ 
        int a_ExitStates(){return a_exitStates;}              /**<Getter for a_exitStates*/ 

        int g_Formalism(){return g_formalism;}              /**<Getter for g_formalism*/ 
        int n_Formalism(){return n_formalism;}              /**<Getter for n_formalism*/ 
        int p_Formalism(){return p_formalism;}              /**<Getter for p_formalism*/ 
        int a_Formalism(){return a_formalism;}              /**<Getter for a_formalism*/ 

        double DecayerMaxL(){return decayerMaxL;}              /**<Getter for decayerMaxL*/
        double PreEqMaxL(){return preEqMaxL;}              /**<Getter for preEqMaxL*/ 
        double g_CutoffEnergy(){return g_cutoffEnergy;}    /**<Getter for g_CutoffEnergy*/ 

        std::string EnergyFile(){return energyFile;}            /**<Getter for energyFile*/ 
        std::string Energies(){return energies;}            /**<Getter for energies*/ 
        std::string ReactionFile(){return reactionFile;}        /**<Getter for reactionFile*/
        std::string Suffix(){return suffix;}                    /**<Getter for suffix*/
        std::string Reaction(){return reaction;}                    /**<Getter for reaction*/
        
        std::string MassTable(){return massTable;}                    /**<Getter for massTable*/
        std::string GdrParams(){return gdrParams;}                    /**<Getter for gdrParams*/
        std::string LevelDir(){return levelDir;}                    /**<Getter for levelDir*/
        std::string SpinFile(){return spinFile;}                    /**<Getter for spinFile*/

        
    private:
        bool calcRates;             /**< Bool if rates should be calculated*/
        bool calcAverageWidth;      /**< Bool if average widths should be calculated*/
        bool residualGamma;        /**< Bool if residual cross section for gamma should be calculated*/
        bool residualNeutron;      /**< Bool if residual cross section for neutron should be calculated*/
        bool residualProton;       /**< Bool if residual cross section for proton should be calculated*/
        bool residualAlpha;        /**< Bool if residual cross section for alpha should be calculated*/
        bool calculateGammaCutoff; /**< */
        bool porterThomas_g;             /**< Bool for PorterThomas usage gamma*/  
        bool porterThomas_p;           /**< Bool for PorterThomas usage particle*/
        bool printTrans;            /**< Bool if transmission should be printed*/
        bool preEq;                 /**< Bool for preEq.*/

        int entranceState;          /**< Int for the number of level, which should be the entrance State*/
        int g_exitStates;           /**< Int for the number of exitStates for gamma residual**/
        int n_exitStates;           /**< Int for the number of exitStates for neutron residual**/
        int p_exitStates;           /**< Int for the number of exitStates for proton residual**/
        int a_exitStates;           /**< Int for the number of exitStates for alpha residual**/

        int n_formalism;            /**< Choose neutron OMP*/
        int a_formalism;            /**< Choose alpha OMP*/
        int p_formalism;            /**< Choose proton OMP*/
        int g_formalism;            /**< Choose gamma strength function*/

        double decayerMaxL;            /**< Maximum l-value for the decayer*/
        double preEqMaxL;            /**< Maximum l-value for the preEq*/

        double g_cutoffEnergy;       /**< Double for the Gamma Cutoff energy*/

        
        std::string energyFile;     /**< String with the path to the energyFile*/
        std::string reactionFile;   /**< String with the path to the reactionFile*/
        std::string suffix;         /**< Suffix for output*/
        std::string reaction;         /**< String for reaction*/
        std::string energies;         /**< String for reaction*/

        std::string massTable;      /**< String for the path to the mass tables*/
        std::string gdrParams;      /**< String for the path to the GDR parameter file*/
        std::string levelDir;       /**< String for the path to the levels directory*/
        std::string spinFile;       /**< String for the path tot the spinFile*/
};

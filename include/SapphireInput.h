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
        *   @brief Function which can be called from Sapphire.cpp to print out an example InputFile.
        */
        void Go(int argc, char* argv[]);
        
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
        void PreEq(bool x){preEq=x;}                                    /**<Setter for bool preEq*/

        void EntranceState(int x){entranceState=x;}                  /**< Setter for entranceState*/
        void g_ExitStates(int x){g_exitStates=x;}              /**<Setter for g_exitStates*/ 
        void n_ExitStates(int x){n_exitStates=x;}              /**<Setter for n_exitStates*/ 
        void p_ExitStates(int x){p_exitStates=x;}              /**<Setter for p_exitStates*/ 
        void a_ExitStates(int x){a_exitStates=x;}              /**<Setter for a_exitStates*/ 
        void g_Formalism(int x){g_formalism=x;}              /**<Setter for g_formalism*/
        void n_Formalism(int x){n_formalism=x;}              /**<Setter for n_formalism*/
        void p_Formalism(int x){p_formalism=x;}              /**<Setter for p_formalism*/
        void a_Formalism(int x){a_formalism=x;}              /**<Setter for a_formalism*/
        void Events(int x){events=x;}              /**<Setter for events*/
        void ChunkSize(int x){chunkSize=x;}              /**<Setter for chunkSize*/
        void Parity(int x){parity=x;}              /**<Setter for parity*/
        void PType(int x){pType=x;}              /**<Setter for pType*/
        void MassNumber(int x){massNumber=x;}              /**<Setter for massNumber*/
        void ChargeNumber(int x){chargeNumber=x;}              /**<Setter for chargeNumber*/
        void Suffix(int x){suffix=x;}                               /**< Setter for string module*/

        void DecayerMaxL(double x){decayerMaxL=x;}              /**<Setter for decayerMaxL*/ 
        void PreEqMaxL(double x){preEqMaxL=x;}                  /**<Setter for preEqMaxL*/ 
        void g_CutoffEnergy(double x){g_cutoffEnergy=x;}        /**<Setter for g_CutoffEnergy*/ 
        void LowEnergy(double x){lowEnergy=x;}                 /**<Setter for lowEnergy*/ 
        void HighEnergy(double x){highEnergy=x;}                /**<Setter for highEnergy*/ 
        void Spin(double x){spin=x;}                            /**<Setter for spin*/        

        void EnergyFile(std::string x){energyFile=x;}                       /**< Setter for string energyFile*/
        void Energies(std::string x){energies=x;}                       /**< Setter for string energies*/
        void ReactionFile(std::string x){reactionFile=x;}                   /**< Setter for string reactionFile*/
        void Reaction(std::string x){reaction=x;}                    /**<Setter for reaction*/
        void Isotope(std::string x){isotope=x;}                    /**<Setter for isotope*/
        void PreEqConf(std::string x){preEqConf=x;}                    /**<Setter for preEqConf*/        
        void MassTable(std::string x){massTable=x;}                               /**< Setter for string module*/
        void GdrParams(std::string x){gdrParams=x;}                               /**< Setter for string module*/
        void LevelDir(std::string x){levelDir=x;}                               /**< Setter for string module*/
        void SpinFile(std::string x){spinFile=x;}                               /**< Setter for string module*/

        
        //Getter
        bool CalcRates() const {return calcRates;}                     /**<Getter for calcrates*/    
        bool CalcAverageWidth() const {return calcAverageWidth;}       /**<Getter for calcAverageWidth*/    
        bool ResidualGamma() const {return residualGamma;}             /**<Getter for residualGamma*/    
        bool ResidualNeutron() const {return residualNeutron;}         /**<Getter for residualNeutron*/    
        bool ResidualProton() const {return residualProton;}           /**<Getter for residualProton*/
        bool ResidualAlpha() const {return residualAlpha;}             /**<Getter for residualAlpha*/
        bool CalculateGammaCutoff() const {return calculateGammaCutoff;} /**<Getter for calculateGammaCutoff*/
        bool PorterThomas_g() const {return porterThomas_g;}             /**<Getter for bool porterThomas_g*/
        bool PorterThomas_p() const {return porterThomas_p;}             /**<Getter for bool porterThomas_p*/
        bool PreEq() const {return preEq;}                                /**<Getter for bool preEq*/

        int EntranceState() const {return entranceState;}              /**<Getter for entranceState*/ 
        int g_ExitStates() const {return g_exitStates;}              /**<Getter for g_exitStates*/ 
        int n_ExitStates() const {return n_exitStates;}              /**<Getter for n_exitStates*/ 
        int p_ExitStates() const {return p_exitStates;}              /**<Getter for p_exitStates*/ 
        int a_ExitStates() const {return a_exitStates;}              /**<Getter for a_exitStates*/ 
        int g_Formalism() const {return g_formalism;}              /**<Getter for g_formalism*/ 
        int n_Formalism() const {return n_formalism;}              /**<Getter for n_formalism*/ 
        int p_Formalism() const {return p_formalism;}              /**<Getter for p_formalism*/ 
        int a_Formalism() const {return a_formalism;}              /**<Getter for a_formalism*/ 
        int Events() const {return events;}              /**<Getter for events*/ 
        int ChunkSize() const {return chunkSize;}              /**<Getter for chunkSize*/ 
        int Parity() const {return parity;}                    /**<Getter for parity*/
        int PType() const {return pType;}                      /**<Getter for pType*/
        int ChargeNumber() const {return chargeNumber;}        /**<Getter for chargeNumber*/
        int MassNumber() const {return massNumber;}        /**<Getter for massNumber*/
        int Suffix() const {return suffix;}                    /**<Getter for suffix*/

        double DecayerMaxL() const {return decayerMaxL;}              /**<Getter for decayerMaxL*/
        double PreEqMaxL() const {return preEqMaxL;}              /**<Getter for preEqMaxL*/ 
        double g_CutoffEnergy() const {return g_cutoffEnergy;}    /**<Getter for g_CutoffEnergy*/ 
        double LowEnergy() const {return lowEnergy;}              /**<Getter for lowEnergy*/ 
        double HighEnergy() const {return highEnergy;}             /**<Getter for highEnergy*/ 
        double Spin() const {return spin;}                         /**<Getter for spin*/

        std::string EnergyFile() const {return energyFile;}            /**<Getter for energyFile*/ 
        std::string Energies() const {return energies;}            /**<Getter for energies*/ 
        std::string ReactionFile() const {return reactionFile;}        /**<Getter for reactionFile*/
        
        std::string Reaction() const {return reaction;}                    /**<Getter for reaction*/
        std::string Isotope() const {return isotope;}                    /**<Getter for isotope*/
        std::string PreEqConf() const {return preEqConf;}                    /**<Getter for preEqConf*/        
        std::string MassTable() const {return massTable;}                    /**<Getter for massTable*/
        std::string GdrParams() const {return gdrParams;}                    /**<Getter for gdrParams*/
        std::string LevelDir() const {return levelDir;}                    /**<Getter for levelDir*/
        std::string SpinFile() const {return spinFile;}                    /**<Getter for spinFile*/

        
    private:
        bool calcRates;             /**< Bool if rates should be calculated*/
        bool calcAverageWidth;      /**< Bool if average widths should be calculated*/
        bool residualGamma;        /**< Bool if residual cross section for gamma should be calculated*/
        bool residualNeutron;      /**< Bool if residual cross section for neutron should be calculated*/
        bool residualProton;       /**< Bool if residual cross section for proton should be calculated*/
        bool residualAlpha;        /**< Bool if residual cross section for alpha should be calculated*/
        bool calculateGammaCutoff; /**< Bool if the GammaCutoffenergy should be calculated or not.*/
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
        int events;                 /**< Number of decays*/
        int chunkSize;              /**< Portion of the total numbers of decays which is calculated at once.*/
        int parity;                     /**< Parity of the high energy resonance*/
        int pType;                  /**< Type of projectile*/
        int chargeNumber;                      /**< Charge number*/
        int massNumber;                      /**< Mass number*/
        int suffix;         /**< Suffix for output*/

        double decayerMaxL;            /**< Maximum l-value for the decayer*/
        double preEqMaxL;            /**< Maximum l-value for the preEq*/
        double g_cutoffEnergy;       /**< Double for the Gamma Cutoff energy*/
        double lowEnergy;           /**< lower limit for initial energy*/
        double highEnergy;          /**< higher limit for initial energy*/
        double spin;                /**< Spin of the decaying resonance*/
                
        std::string energyFile;     /**< String with the path to the energyFile*/
        std::string reactionFile;   /**< String with the path to the reactionFile*/
        std::string reaction;         /**< String for reaction*/
        std::string energies;         /**< String for reaction*/
        std::string isotope;        /**< String for the istope for the decay simulation*/
        std::string preEqConf;      /**< Preequillibrium exciton configuration of the initial state (pp,ph,np,nh)*/
        std::string massTable;      /**< String for the path to the mass tables*/
        std::string gdrParams;      /**< String for the path to the GDR parameter file*/
        std::string levelDir;       /**< String for the path to the levels directory*/
        std::string spinFile;       /**< String for the path tot the spinFile*/
};

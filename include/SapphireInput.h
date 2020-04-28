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
#include "NuclearMass.h"
#include "NuclearLevels.h"


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

        /**
         * @brief Get the massNumberString from the reactionString.
         * @param reactionString A string from the cmd line which represents the reaction, e.g., 60Fe+p
         * @returns The mass number part of the reactionString, e.g. "60"
         */
        std::string massNumberStringFromReactionString(std::string reactionString);

        /**
         * @brief Get the massNumberInt from the reactionString.
         * @param reactionString A string from the cmd line which represents the reaction, e.g., 60Fe+p
         * @returns The mass number part of the reactionString, e.g. 60, if successful. Otherwise returns 0.
         */
        int massNumberIntFromReactionString(std::string reactionString);

        /**
         * @brief Get the pTypeString from the reactionString.
         * @param reactionString A string from the cmd line which represents the reaction, e.g., 60Fe+p
         * @returns The pType part of the reactionString, e.g. "p"
         */
        std::string pTypeStringFromReactionString(std::string reactionString);

        /**
         * @brief Get a projectile integer from a projectile string
         * @param reactionString The part of the reactionString which defines the projectile.
         * @return Integer which represents the pType: 0 = gamma, 1 = neutron, 2 = proton, 3 = alpha.
         */
        int pTypeIntFromReactionString(std::string reactionString);

        /**
         * @brief Get the atomic number as string from reactionString.
         * @param reactionString A string from the cmd line which represents the reaction, e.g., 60Fe+p
         * @returns The atomic number part of the reactionString, e.g. Fe
        */
        std::string atomicNumberStringFromReactionString(std::string reactionString);

        /**
         * @brief Get the atomic number as int from reactionString.
         * @param reactionString A string from the cmd line which represents the reaction, e.g., 60Fe+p
         * @returns The atomic number part of the reactionString, e.g. 28 other wise returns 0.
         */
        int atomicNumberIntFromReactionString(std::string reactionString);

        /**
         * @brief Get the atomic number as string from isotopeString.
         * @param isotopeString A string from the cmd line which represents the reaction, e.g., 60Fe+p
         * @returns The atomic number part of the isotopeString, e.g. Fe
         */
        std::string atomicNumberStringFromIsotopeString(std::string isotopeString);

        /**
         * @brief Get the atomic number as int from isotopeString.
         * @param isotopeString A string from the cmd line which represents the reaction, e.g., 60Fe+p
         * @returns The atomic number part of the isotopeString, e.g. 28 other wise returns 0.
         */
        int atomicNumberIntFromIsotopeString(std::string isotopeString);

        /**
         * @brief Get the massNumberString from the isotopeString.
         * @param isotopeString A string from the cmd line which represents the reaction, e.g., 60Fe+p
         * @returns The mass number part of the isotopeString, e.g. "60"
         */
        std::string massNumberStringFromIsotopeString(std::string isotopeString);

        /**
         * @brief Get the massNumberInt from the isotopeString.
         * @param isotopeString A string from the cmd line which represents the reaction, e.g., 60Fe+p
         * @returns The mass number part of the isotopeString, e.g. 60, if successful. Otherwise returns 0.
         */
        int massNumberIntFromIsotopeString(std::string isotopeString);

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
        void PrintTrans(bool x){printTrans=x;}                      /**<Setter for bool printTrans*/
        void PreEq(bool x){preEq=x;}                                    /**<Setter for bool preEq*/
        void AlphaChannel(bool x){alphaChannel=x;}                /**<Setter for bool alphaChannel*/
        void ProtonChannel(bool x){protonChannel=x;}                /**<Setter for bool protonChannel*/
        void NeutronChannel(bool x){neutronChannel=x;}                /**<Setter for bool neutronChannel*/
        void GammaChannel(bool x){gammaChannel=x;}                /**<Setter for bool gammaChannel*/

        void EntranceState(int x){entranceState=x;}                  /**< Setter for entranceState*/
        void g_ExitStates(int x){               /**<Setter for g_exitStates*/ 
            g_exitStates=x;
            exitStates[0]=x;
        }              
        void n_ExitStates(int x){              /**<Setter for n_exitStates*/ 
            n_exitStates=x;
            exitStates[1]=x;
        }
        void p_ExitStates(int x){               /**<Setter for p_exitStates*/ 
            p_exitStates=x;
            exitStates[2]=x;
        }             
        void a_ExitStates(int x){              /**<Setter for a_exitStates*/ 
            a_exitStates=x;
            exitStates[3]=x;
        }
        void g_Formalism(int x){g_formalism=x;}              /**<Setter for g_formalism*/
        void n_Formalism(int x){n_formalism=x;}              /**<Setter for n_formalism*/
        void p_Formalism(int x){p_formalism=x;}              /**<Setter for p_formalism*/
        void a_Formalism(int x){a_formalism=x;}              /**<Setter for a_formalism*/
        void LevelDensity(int x){leveldensity=x;}              /**<Setter for a_formalism*/
        void Events(int x){events=x;}              /**<Setter for events*/
        void ChunkSize(int x){chunkSize=x;}              /**<Setter for chunkSize*/
        void Parity(int x){parity=x;}              /**<Setter for parity*/
        void PType(int x){pType=x;}              /**<Setter for pType*/
        void CompoundZ(int x){compoundZ=x;}              /**<Setter for compoundZ*/
        void CompoundA(int x){compoundA=x;}              /**<Setter for compoundA*/
        void GroundstatePi(int x){groundstatePi=x;}          /**<Setter for groundstatePi*/
        void Suffix(int x){suffix=x;}                               /**< Setter for string module*/
        void XsZ(int x) {xsZ=x;}                    /**<Setter for xsZ*/
        void XsA(int x) {xsA=x;}                    /**<Setter for xsA*/
        void DcZ(int x) {dcZ=x;}                    /**<Setter for dcZ*/
        void DcA(int x) {dcA=x;}                    /**<Setter for dcA*/

        void DecayerMaxL(double x){decayerMaxL=x;}              /**<Setter for decayerMaxL*/ 
        void PreEqMaxL(double x){preEqMaxL=x;}                  /**<Setter for preEqMaxL*/ 
        void g_CutoffEnergy(double x){g_cutoffEnergy=x;}        /**<Setter for g_CutoffEnergy*/ 
        void LowEnergy(double x){lowEnergy=x;}                 /**<Setter for lowEnergy*/ 
        void HighEnergy(double x){highEnergy=x;}                /**<Setter for highEnergy*/ 
        void Spin(double x){spin=x;}                            /**<Setter for spin*/        
        void GroundstateJ(double x){groundstateJ=x;}            /**<Setter for groundstateJ*/        
        void QValue(double x){qValue=x;}                        /**<Setter for qValue*/  


        void EnergyFile(std::string x){energyFile=x;}                       /**< Setter for string energyFile*/
        void Energies(std::string x){energies=x;}                       /**< Setter for string energies*/
        void ReactionFile(std::string x){reactionFile=x;}                   /**< Setter for string reactionFile*/
        void Reaction(std::string x){                                   /**<Setter for reaction*/
            reaction=x; 
            PType(pTypeIntFromReactionString(x));
            XsZ(atomicNumberIntFromReactionString(x));
            XsA(massNumberIntFromReactionString(x));
        }
        void Isotope(std::string x){                    /**<Setter for isotope*/
            isotope=x;
            DcZ(atomicNumberIntFromIsotopeString(x));
            DcA(massNumberIntFromIsotopeString(x));
        }
        void PreEqConf(std::string x){preEqConf=x;}                /**<Setter for preEqConf*/        
        void MassTable(std::string x){massTable=x;}                /**< Setter for string module*/
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
        bool PrintTrans() const {return printTrans;}                       /**Getter for bool printTrans*/
        bool PreEq() const {return preEq;}                                /**<Getter for bool preEq*/
        bool AlphaChannel() const {return alphaChannel;}              /**<Getter for bool alphaChannel*/
        bool ProtonChannel() const {return protonChannel;}              /**<Getter for bool protonChannel*/
        bool NeutronChannel() const {return neutronChannel;}              /**<Getter for bool neutronChannel*/
        bool GammaChannel() const {return gammaChannel;}              /**<Getter for bool gammaChannel*/

        int EntranceState() const {return entranceState;}              /**<Getter for entranceState*/ 
        int g_ExitStates() const {return g_exitStates;}              /**<Getter for g_exitStates*/ 
        int n_ExitStates() const {return n_exitStates;}              /**<Getter for n_exitStates*/ 
        int p_ExitStates() const {return p_exitStates;}              /**<Getter for p_exitStates*/ 
        int a_ExitStates() const {return a_exitStates;}              /**<Getter for a_exitStates*/ 
        int g_Formalism() const {return g_formalism;}              /**<Getter for g_formalism*/ 
        int n_Formalism() const {return n_formalism;}              /**<Getter for n_formalism*/ 
        int p_Formalism() const {return p_formalism;}              /**<Getter for p_formalism*/ 
        int a_Formalism() const {return a_formalism;}              /**<Getter for a_formalism*/ 
        int LevelDensity() const {return leveldensity;}              /**<Getter for a_formalism*/ 
        int Events() const {return events;}              /**<Getter for events*/ 
        int ChunkSize() const {return chunkSize;}              /**<Getter for chunkSize*/ 
        int Parity() const {return parity;}                    /**<Getter for parity*/
        int PType() const {return pType;}                      /**<Getter for pType*/
        int CompoundA() const {return compoundA;}                 /**<Getter for compoundA*/
        int CompoundZ() const {return compoundZ;}                /**<Getter for compoundZ*/
        int GroundstatePi() const {return groundstatePi;}       /**<Getter for groundstatePi*/
        int Suffix() const {return suffix;}                    /**<Getter for suffix*/
        int XsZ() const {return xsZ;}                    /**<Getter for xsZ*/
        int XsA() const {return xsA;}                    /**<Getter for xsA*/
        int DcZ() const {return dcZ;}                    /**<Getter for dcZ*/
        int DcA() const {return dcA;}                    /**<Getter for dcA*/

        double DecayerMaxL() const {return decayerMaxL;}              /**<Getter for decayerMaxL*/
        double PreEqMaxL() const {return preEqMaxL;}              /**<Getter for preEqMaxL*/ 
        double g_CutoffEnergy() const {return g_cutoffEnergy;}    /**<Getter for g_CutoffEnergy*/ 
        double LowEnergy() const {return lowEnergy;}              /**<Getter for lowEnergy*/ 
        double HighEnergy() const {return highEnergy;}             /**<Getter for highEnergy*/ 
        double Spin() const {return spin;}                         /**<Getter for spin*/
        double GroundstateJ() const {return groundstateJ;}              /**<Getter for groundstateJ*/
        double QValue() const {return qValue;}                         /**<Getter for qValue*/

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

        std::vector<int> exitStates;  /**< Vector of integer for the number of exitStates*/
    private:
        bool calcRates;             /**< Bool if rates should be calculated*/
        bool calcAverageWidth;      /**< Bool if average widths should be calculated*/
        bool residualGamma;        /**< Bool if residual cross section for gamma should be calculated*/
        bool residualNeutron;      /**< Bool if residual cross section for neutron should be calculated*/
        bool residualProton;       /**< Bool if residual cross section for proton should be calculated*/
        bool residualAlpha;        /**< Bool if residual cross section for alpha should be calculated*/
        bool calculateGammaCutoff; /**< Bool if the GammaCutoffenergy should be calculated or not.*/
        bool porterThomas_g;       /**< Bool for PorterThomas usage gamma*/  
        bool porterThomas_p;       /**< Bool for PorterThomas usage particle*/
        bool printTrans;            /**< Bool if transmission should be printed*/
        bool preEq;                 /**< Bool for preEq.*/
        bool alphaChannel;          /**< Toggle alpha channel*/
        bool protonChannel;          /**< Toggle proton channel*/
        bool neutronChannel;          /**< Toggle neutron channel*/
        bool gammaChannel;          /**< Toggle gamma channel*/

        int entranceState;          /**< Int for the number of level, which should be the entrance State*/
        int g_exitStates;           /**< Int for the number of exitStates for gamma residual**/
        int n_exitStates;           /**< Int for the number of exitStates for neutron residual**/
        int p_exitStates;           /**< Int for the number of exitStates for proton residual**/
        int a_exitStates;           /**< Int for the number of exitStates for alpha residual**/
        int n_formalism;            /**< Choose neutron OMP*/
        int a_formalism;            /**< Choose alpha OMP*/
        int p_formalism;            /**< Choose proton OMP*/
        int g_formalism;            /**< Choose gamma strength function*/
        int leveldensity;           /**< Choose a level density model*/
        int events;                 /**< Number of decays*/
        int chunkSize;              /**< Portion of the total numbers of decays which is calculated at once.*/
        int parity;                 /**< Parity of the high energy resonance*/
        int pType;                  /**< Type of projectile*/
        int compoundA;             /**< Int for the mass of the compound nuclei*/
        int compoundZ;             /**< Int for the charge of the compound nuclei*/
        int groundstatePi;          /**< Int for the parity of the groundstate*/
        int suffix;                 /**< Suffix for output*/
        int xsZ;                    /** Integer for the charge number of the isotope in reaction module*/
        int xsA;                    /** Integer for the mass number of the isotope in reaction module*/
        int dcZ;                    /** Integer for the charge number of the isotope in decayer module*/
        int dcA;                    /** Integer for the mass number of the isotope in decayer module*/

        double decayerMaxL;            /**< Maximum l-value for the decayer*/
        double preEqMaxL;            /**< Maximum l-value for the preEq*/
        double g_cutoffEnergy;       /**< Double for the Gamma Cutoff energy*/
        double lowEnergy;           /**< lower limit for initial energy*/
        double highEnergy;          /**< higher limit for initial energy*/
        double spin;                /**< Spin of the decaying resonance*/
        double groundstateJ;        /**< Spin of the groundstate*/
        double qValue;              /**< double for the qvalue of the reaction*/
                
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


    public:
        //NuclearMass masses;
        //NuclearLevels levels;
        
};

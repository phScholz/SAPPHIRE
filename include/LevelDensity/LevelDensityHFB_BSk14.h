/**
 * @file LevelDensityHFB_BSk14.h
 * @brief contains the LevelDensityHFB_BSk14 class
 * @date 2020-04-27
 */

#pragma once
#include "LevelDensityTable.h"
#include <tr1/unordered_map>


/**
 * @brief A Class which represents one line in the *.tab files of the HFB density
 * @details
 * For odd mass nuclei, the spins begin with 1/2 and end with 99/2.
 */
class HFBTabRow{

    public:
        HFBTabRow(){};
        void PrintRow();

    public:
        double U; /**< Energy in MeV*/
        double T; /**< Temperature in MeV*/
        double NCUMUL; /**< Cumulated number of levels*/
        double RHOOBS; /**< Observered Level Density*/
        double RHOTOT; /**< Total level density*/
        double J0;
        double J1;
        double J2;
        double J3;
        double J4;
        double J5;
        double J6;
        double J7;
        double J8;
        double J9;
        double J10;
        double J11;
        double J12;
        double J13;
        double J14;
        double J15;
        double J16;
        double J17;
        double J18;
        double J19;
        double J20;
        double J21;
        double J22;
        double J23;
        double J24;
        double J25;
        double J26;
        double J27;
        double J28;
        double J29;
        double J30;
        double J31;
        double J32;
        double J33;
        double J34;
        double J35;
        double J36;
        double J37;
        double J38;
        double J39;
        double J40;
        double J41;
        double J42;
        double J43;
        double J44;
        double J45;
        double J46;
        double J47;
        double J48;
        double J49;
};

typedef std::tr1::unordered_map<ChargeMassParityKey, std::vector<HFBTabRow> > HFBTable; /**< One map object which maps MassKey to LevelsContainer */

/**
 * @brief The HFB level density tables as class in sapphire
 */
class LevelDensityHFB_BSk14 : public LevelDensityTable {
    public:
        static HFBTable densityTable;
    private:
        HFBTabRow dummyRow;        
        void GetFileName();
        void ReadFile();
        void FillVector();
        std::vector<HFBTabRow> rows;
        bool verbose_=false;    

    public:
        LevelDensityHFB_BSk14(int Z, int A, double J, int parity) : LevelDensityTable(Z,A,J,parity){
            
            if(verbose_) std::cout << "LevelDensityHFB Constructor ... " << std::endl;
            SetTables(true);
            if(!FindDensities(Z_,A_,parity_))
            {
                if(verbose_) std::cout << "Getting Filename ... " << std::endl;
                GetFileName();
                if(verbose_) std::cout << "Reading File ... " << filename << std::endl;
                ReadFile();
            }
                if(verbose_) std::cout << "Filling Vector ... " << filename << std::endl;
                FillVector();
        }

        LevelDensityHFB_BSk14(int Z, int A, double J) : LevelDensityTable(Z,A,J,1){
            parity_=1;
            if(verbose_) std::cout << "LevelDensityHFB Constructor ... " << std::endl;
            SetTables(true);
            if(!FindDensities(Z_,A_,parity_))
            {
                if(verbose_) std::cout << "Getting Filename ... " << std::endl;
                GetFileName();
                if(verbose_) std::cout << "Reading File ... " << filename << std::endl;
                ReadFile();
            }
            if(verbose_) std::cout << "Filling Vector ... " << filename << std::endl;
            FillVector();
        }
    
    void PrintRows();

    /**
     * @brief Look if the densities are already stored in the HFBTable or not
     * @param Z charge number
     * @param A mass number
     * @param P parity
     * @return True or False depending on, if the specific std::vector<HFBTabRow> is already in HFBTable or not
     * @details
     * If the densities were found in the DensityTable then the std::vector<HFBTabRow> is copied to the member rows. 
     * Basically this is a skipping of the reading procedure, if the densities have been previously read.
     * With that, various calls for different spins should be accelerated.
     */
    bool FindDensities(int Z, int A, int P);
    void PrintDensityTable();
    double CalculateDensity(double E);
    double TotalLevelDensity(double E);


};


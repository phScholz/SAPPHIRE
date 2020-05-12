/**
 * @file LevelDensityHFB_BSk14.h
 * @brief contains the LevelDensityHFB_BSk14 class
 * @date 2020-04-27
 * @details
 * Quote from the RIPL-3 website:
 * 
 *       Directory: DENSITIES
 *       File:      level-densities-hfb.readme (28 November 2008)
 *       *******************************************************
 *       
 *                                level-densities-hfb/zxxx.dat  
 *       HFB plus combinatorial nuclear level densities at ground state deformations
 *                    (Provided by S. Goriely, November 28, 2008)
 *              *******************************************************
 *       
 *       Contents
 *       --------
 *       The files contains the HFB plus combinatorial nuclear level densities
 *       at ground state deformations [1]. The nuclear level density is coherently
 *       obtained on the basis of the single-particle level scheme and pairing energy
 *       derived at the ground state deformation based on the BSk14 Skyrme force [2].
 *       
 *       The nuclear level density obtained within the HFB plus combinatorial model [1]
 *       at the ground state deformation are given in table format (in a energy, spin 
 *       and parity grid). Each isotopic chain is included in a zxxx file (where xxx 
 *       corresponds to the atomic number Z) .
 *       
 *       The nuclear level densities are provided for about 9000 nuclei.
 *       
 *       Format
 *       ------
 *       Each record of the file, per isotope, contains 2 sets of data including the
 *       positive-parity and negative-parity nuclear level density in a spin (50 first spins)
 *       and energy (60 rows for U=0.25 to 200MeV) grid. Also included are an estimate
 *       of the nuclear temperature, the cumulative number of levels, the total level
 *       density and the total state density
 *       - title line including the saddle/isomer properties
 *            Z     : charge number
 *            A     : mass number
 *            B     : saddle or isomer energy above ground state (in MeV)
 *            hw    : saddle or isomer estimated width assuming a parabolic shape (in MeV)
 *            beta2 : quadrupole deformation parameter
 *            beta3 : octupole deformation parameter
 *            beta4 : hexadecapole deformation parameter
 *       
 *       The corresponding fortran format is (25x,i3,3x,i3,67x,f6.2,4x,f6.2,3(4x,f6.3))
 *       
 *       - 60 data lines including the energy-, spin-dependent level density for
 *       positive parities
 *            U     : excitation energy in MeV
 *            T     : nuclear temperature in MeV
 *            Ncumul: cumulative number of positive-parity levels
 *            Rhoobs: total level density in 1/MeV (positive parity)
 *            Rhotot: total state density in 1/MeV (positive parity)
 *            Rho   : spin-dependent level density (in 1/MeV) for the first 50 spins (positive parity)
 *       
 *       The corresponding fortran format is (f7.2,f7.3,1x,1p,53e9.2)
 *       
 *       - 60 data lines including the energy-, spin-dependent level density for
 *       negative parities
 *            U     : excitation energy in MeV
 *            T     : nuclear temperature in MeV
 *            Ncumul: cumulative number of negative-parity levels
 *            Rhoobs: total level density in 1/MeV (negative parity)
 *            Rhotot: total state density in 1/MeV (negative parity)
 *            Rho   : spin-dependent level density (in 1/MeV) for the first 50 spins (negative parity)
 *       
 *       The corresponding fortran format is (f7.2,f7.3,1x,1p,53e9.2)
 *       
 *       References
 *       ----------
 *       [1] S. Goriely, S. Hilaire, A.J. Koning, Improved microscopic nuclear level densities within
 *       the Hartree-Fock-Bogoliubov plus combinatorial method, Phys. Rev. C78 (2008) 064307. 
 *       [2] S. Goriely, M. Samyn, J.M. Pearson, Phys. Rev. C75 (2007) 064312

 */

#pragma once
#include "LevelDensityTable.h"
#include "NuclearMass.h"
#include <tr1/unordered_map>


/**
 * @brief A Class which represents one line in the *.tab files of the HFB density
 * @details
 * For odd mass nuclei, the spins begin with 1/2 and end with 99/2.
 */
class HFBTabRow{

    public:
        /**
         * @brief Constructor
         */
        HFBTabRow(){};

        /**
         * @brief Method to print a row to std::cout.
         */
        void PrintRow();

    public:
        double U; /**< Energy in MeV*/
        double T; /**< Temperature in MeV*/
        double NCUMUL; /**< Cumulated number of levels*/
        double RHOOBS; /**< Observered Level Density*/
        double RHOTOT; /**< Total level density*/
        double J0;  /**< Density for spin 0 or 1/2 depending on if odd or even mass nucleus.*/
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

/**
 * @brief A Class which represents one line in the *.ld files of the HFB density
 * @details
 * For odd mass nuclei, the spins begin with 1/2 and end with 99/2.
 */
class HFBLdRow{
    public:
        int Z;
        int A;
        double c;
        double d;

};

typedef std::tr1::unordered_map<MassKey, HFBLdRow > HFBCorrTable;
typedef std::tr1::unordered_map<ChargeMassParityKey, std::vector<HFBTabRow> > HFBTable; /**< One map object which maps MassKey to LevelsContainer */

/**
 * @brief The HFB level density tables as class in sapphire
 */
class LevelDensityHFB_BSk14 : public LevelDensityTable {
    public:
        /**
         * @brief A static HFBTable map.
         * @details
         * This was introduced to reduce the I/O for this NLD model. 
         * The densities for each nucleus are now read only once and then stored completely 
         * in this object.
         */
        static HFBTable densityTable;

        /**
         * @brief A static HFBCorrTable map.
         * @details
         * This was introduced to reduce the I/O for this NLD model. 
         * The correction parameters for each nucleus are now read only once and then stored completely 
         * in this object.
         */
        static HFBCorrTable corrTable;

    private:
        HFBTabRow dummyRow;    /**< Dummy row object which is used while reading the density files*/    

        /**
         * @brief Method to construct the path of the level density file and store it in LevelDensityTable::fileName.
         */
        void GetFileName();     

        /**
         * @brief Logic to read the density files and store it in the std::vector<HFBTabRow> rows.
         */
        void ReadFile();

        /**
         * @brief Logic to read the correction files.
         */
        void ReadCorr();

        /**
         * @brief Method to get the density of one spin and parity from std::vector<HFBTabRow> rows into DensityVector.
         */
        void FillVector();

        std::vector<HFBTabRow> rows; /**< std::vector object to store the rows from the density files*/
        bool verbose_=false;    /**< Boolean which can be set true for debugging purposes*/

    public:
        /**
         * @brief Constructor
         * @param Z Charge number
         * @param A Mass number
         * @param J Spin
         * @param parity Parity
         * @details
         * 1. SetTables() is set to true.
         * 2. It will be checked via FindDensities() if the densities of this nucleus has been read before.
         *  - If yes, the densities are taken from the DensityTable object
         *  - If not, the Filename is constructet via GetFileName() and ReadFile() will be called.
         * 3. FillVector() fills the DensityVector object with the densities for the respective spin and parity.
         */
        LevelDensityHFB_BSk14(int Z, int A, double J, int parity) : LevelDensityTable(Z,A,J,parity){
            
            if(verbose_) std::cout << "LevelDensityHFB Constructor ... " << std::endl;
            SetTables(true);
            if(J_>=0){
                if(!FindDensities())
                {
                    if(verbose_) std::cout << "Getting Filename ... " << std::endl;
                    GetFileName();
                    if(verbose_) std::cout << "Reading File ... " << filename << std::endl;
                    ReadFile();
                    ReadCorr();
                }
                    if(verbose_) std::cout << "Filling Vector ... " << filename << std::endl;
                    FillVector();
            }
        }

        /**
         * @brief Constructor w/o parity. This was needed for compatibility with older code fragments. See LevelDensityHFB_BSk14().
         */
        LevelDensityHFB_BSk14(int Z, int A, double J) : LevelDensityTable(Z,A,J,1){
            parity_=1;
            if(J_>=0){
                if(!FindDensities())
                {
                    if(verbose_) std::cout << "Getting Filename ... " << std::endl;
                    GetFileName();
                    if(verbose_) std::cout << "Reading File ... " << filename << std::endl;
                    ReadFile();
                    ReadCorr();
                }
                    if(verbose_) std::cout << "Filling Vector ... " << filename << std::endl;
                    FillVector();
            }
        }


        /**
         * @brief Printing the rows in the rows object. Calls HFBTabRow::PrintRow().
         */
        void PrintRows();

        /**
         * @brief Look if the densities are already stored in the HFBTable or not
         * @return True or False depending on, if the specific std::vector<HFBTabRow> is already in HFBTable or not
         * @details
         * If the densities were found in the DensityTable then the std::vector<HFBTabRow> is copied to the member rows. 
         * Basically this is a skipping of the reading procedure, if the densities have been previously read.
         * With that, various calls for different spins should be accelerated.
         */
        bool FindDensities();

        /**
         * @brief Print the content of DensityTable to std::out.
         */
        void PrintDensityTable();

        /**
         * @brief Main method which will be called by other objects.
         * @param E Energy.
         * @return The level density at a respective energy in [1/MeV].
         * @details
         * Returns the level density for one spin and parity.
         * Densities are interpolated by a linear function.
         * Adjustment flexibility is included via *c_* and *d_* which are 
         * by standard the corrections from 
         * [S. Goriely et al., Phys. Rev. C *78*, 064307 (2008)](https://doi.org/10.1103/PhysRevC.78.064307)
         * Adjustment via:
         * \f[\rho(U,J,P)_{adj.} = e^{c\sqrt{U-d}} \rho(U-d,J,P) \f]
         * @todo 
         * - Move from linear interpolation to something like exponential.
         */
        double CalculateDensity(double E);

        /**
         * @brief Get the total level density at a specific energy.
         * @param E Energy.
         * @returns Total level density in [1/MeV]
         * @details
         * Returns the entries of the column RHOTOT. Linear interpolation.
         * Adjustment flexibility is included via *c_* and *d_* which are 
         * by standard the corrections from 
         * [S. Goriely et al., Phys. Rev. C *78*, 064307 (2008)](https://doi.org/10.1103/PhysRevC.78.064307) 
         * Adjustment via:
         * \f[\rho(U,J,P)_{adj.} = e^{c\sqrt{U-d}} \rho(U-d,J,P) \f]
         * @todo
         * - Move from linear interpolation to something like exponential.
         */
        double TotalLevelDensity(double E);



};


/**
 * @file D1MQRPA.h
 * @brief Declaration of the D1M_QRPA GSF model.
 * @date 2020-04-30
 * @details
 * References:
 * 1. S. Goriely et al., The European Physical Journal A 55, 172 (2019)
 * 2. [S. Goriely, S. Hilaire, S. Péru, K. Sieja, Phys. Rev. C 98, 014327 (2018)](https://doi.org/10.1103/PhysRevC.98.014327)
 */

#include "GammaTransmissionFunc.h"
#include "NuclearMass.h"
#include <vector>


class QRPAE1row{
    public:
        double energy;
        std::vector<double> strength;
};

class QRPAM1row{
    public:
        double energy;
        std::vector<double> strength;
};

typedef std::tr1::unordered_map<MassKey, std::vector<QRPAE1row> > QRPA_E1_Table; /**< One map object which maps MassKey to QRPAE1row */
typedef std::tr1::unordered_map<MassKey, std::vector<QRPAM1row> > QRPA_M1_Table; /**< One map object which maps MassKey to QRPAM1row */


class D1MQRPA : public GammaTransmissionFunc {
    public:
        D1MQRPA(int z2, int m2, double jInitial, int piInitial,
			   double jFinal, int piFinal, double maxL, 
			   double totalWidthForCorrection,
			   double uncorrTotalWidthForCorrection,
			   double uncorrTotalWidthSqrdForCorrection,
			   TransmissionFunc* previous, double exEnergy):
               GammaTransmissionFunc(z2,m2,jInitial,piInitial,jFinal,piFinal,maxL,
			totalWidthForCorrection,uncorrTotalWidthForCorrection,
			uncorrTotalWidthSqrdForCorrection,previous), exEnergy_(exEnergy) {
                SetUpperLimit();
                ReadFile();
        };

        //virtual ~D1MQRPA();

        /**
         * @brief
         */
        double CalcStrengthFunction(double energy);

    private:
    void ReadFile();
    bool ReadE1();
    bool ReadM1();
    bool FindE1();
    bool FindM1();
    void PrintE1();
    void PrintM1();
    /**
     * @brief Set the parameters of the lower limit correction
     * @details
     * See [S. Goriely et al., Phys. Rev. C 98 014327 (2018)](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.98.014327) for details.
     */
    void SetLowerLimit();
    
    /**
     * @brief Set the parameters of the upper limit correction
     * @details
     * See [S. Goriely et al., Phys. Rev. C 98 014327 (2018)](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.98.014327) for details.
     */
    void SetUpperLimit();
    double CalcE1Strength(double energy);
    double CalcM1Strength(double energy);

    QRPAE1row dummyRowE1;
    QRPAM1row dummyRowM1;
    std::vector<QRPAE1row> e1Rows;
    std::vector<QRPAM1row> m1Rows;
    static QRPA_E1_Table e1Table;
    static QRPA_M1_Table m1Table;
    bool verbose_ =false;
    bool e1_;
    bool m1_;
    
    double exEnergy_; /**< Excitation energy.*/

    /* Parameters for the analytical correction*/
    double f0_;
    double e0_;
    double eta_;
    double C_;
};




#ifndef CROSSSECTION_H
#define CROSSSECTION_H
#include <vector>
#include <map>
#include <string>

class Decayer;
class SpinRatePair;

typedef std::vector<std::pair<Decayer*,std::vector<SpinRatePair*> > > DecayerVector;
typedef std::pair<int,double> int_double_pair;
struct int_double_pair_compare {
  bool operator()(const int_double_pair &lhs, int_double_pair const  &rhs) {
        return lhs.second < rhs.second;
  }
};

class CrossSectionValues {
 public:
 CrossSectionValues(double gamma, double neutron, double proton, double alpha,
		    double gammaStellar, double neutronStellar, double protonStellar, double alphaStellar) :
  gamma_(gamma),neutron_(neutron),proton_(proton),alpha_(alpha),
  gammaStellar_(gammaStellar),neutronStellar_(neutronStellar),protonStellar_(protonStellar),alphaStellar_(alphaStellar) {};
  double gamma_;
  double neutron_;
  double proton_;
  double alpha_;
  double gammaStellar_;
  double neutronStellar_;
  double protonStellar_;
  double alphaStellar_;
};

class CrossSection {
 public:
  CrossSection(int,int,int,std::string,bool,int entranceState = 0, std::vector<int> exitStates = std::vector<int>(4,-1));
  CrossSection_simple(int,int,int,std::string,bool,int entranceState = 0, std::vector<int> exitStates = std::vector<int>(4,-1));

  bool IsValid() const {
    return isValid_;
  }

  void Calculate();
  void PrintCrossSections();
  void PrintTransmissionTerms();
  std::pair<double,double> CalcAverageSWaveResWidth();
  std::pair<double,double> CalcAveragePWaveResWidth();
  std::pair<double,double> CalcAverageDWaveResWidth();
  void CalculateReactionRates(bool);
  void PrintReactionRates(bool);
  static void SetResidualGamma(bool residual) {residualGamma_=residual;};
  static void SetResidualNeutron(bool residual) {residualNeutron_=residual;};
  static void SetResidualProton(bool residual) {residualProton_=residual;};
  static void SetResidualAlpha(bool residual) {residualAlpha_=residual;};
  static void SetCalculateGammaCutoff(bool calc) {calculateGammaCutoff_=calc;};
  static void CreateTempVector();
  static void CreateMACSEnergiesVector();

 private:
  bool FillEnergies(std::string);
  bool CalcAllowedJPi(bool);
  bool CalcDecayerVector(double,DecayerVector&,bool forAverageWidth=false);
  void CalculateEnergyGrid();
  void CalcPartitionFunc();
 private:
  static bool residualGamma_;
  static bool residualNeutron_;
  static bool residualProton_;
  static bool residualAlpha_;
  static bool calculateGammaCutoff_;
  constexpr static double minEnergy_ = 0.001;
  constexpr static double maxEnergy_ = 15.0;
  constexpr static double dE_ = .1;
  bool gammaCutoffSet_;
  bool isValid_;
  bool energiesGiven_;
  int Z_;
  int A_;
  int pType_;
  int compoundA_;
  int compoundZ_;
  int groundStatePi_;
  int entranceState_;
  double preFactor_;
  double groundStateJ_;
  double seperationEnergy_;
  double skipEnergy_;
  std::vector<int> exitStates_;
  std::vector<double> specifiedExitSepE_;
  std::vector<double> specifiedExitJ_;
  std::vector<int> specifiedExitPi_;
  std::vector<bool> skipped_;
  std::vector<std::pair<double,int> > allowedJPi_;
  std::vector<std::pair<double,CrossSectionValues> > crossSections_;
  std::vector<std::pair<double,CrossSectionValues> > reactionRates_;
  std::map<int, std::vector<double> > entranceTrans_;
  std::map<int, std::vector<double> > gExitTrans_;
  std::map<int, std::vector<double> > nExitTrans_;
  std::map<int, std::vector<double> > pExitTrans_;
  std::map<int, std::vector<double> > aExitTrans_;
  std::vector<std::pair<double,double> > partFunc_;
  std::vector<std::pair<double,double> > partFuncMACS_;
  static std::vector<double> rateTemps_;
  static std::vector<double> macsEnergies_;
};


#endif

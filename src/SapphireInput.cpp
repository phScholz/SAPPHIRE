/**
 * @file SapphireInput.h
 * @brief Option class for Sapphire.
 * @date 2020
 * @author Philipp Scholz, <pscholz@outlook.de>
 * 
 */

#include "SapphireInput.h"
#include "NuclearMass.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/exception/all.hpp>
#include <boost/throw_exception.hpp>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>


    extern std::string sourceDirectory();
   
    SapphireInput::SapphireInput(){
        SapphireInput::Initialize();
    }

    void SapphireInput::Go(int argc, char* argv[]){
        SapphireInput::PrintIntputParameters("all");
    }

    void SapphireInput::Initialize(){
        /** Setting default configurations for the SapphireInput class. */
        std::cout << "Setting default values..." << std::endl;
        SapphireInput::exitStates = std::vector<int>(4,-1);
        SapphireInput::g_ExitStates(-1);
        SapphireInput::n_ExitStates(-1);
        SapphireInput::p_ExitStates(-1);
        SapphireInput::a_ExitStates(-1);   
        SapphireInput::CalcRates(false);           
        SapphireInput::CalcAverageWidth(false);
        SapphireInput::ResidualGamma(true);               
        SapphireInput::ResidualNeutron(false);           
        SapphireInput::ResidualProton(false);
        SapphireInput::ResidualAlpha(false);
        SapphireInput::CalculateGammaCutoff(false);   
        SapphireInput::PorterThomas_g(false);
        SapphireInput::PorterThomas_p(false);
        SapphireInput::EntranceState(0);
        SapphireInput::a_Formalism(0);
        SapphireInput::p_Formalism(0);
        SapphireInput::n_Formalism(0);
        SapphireInput::g_Formalism(1);
        SapphireInput::LevelDensity(1);
        SapphireInput::DecayerMaxL(8.);
        SapphireInput::PreEqMaxL(8.);
        SapphireInput::g_CutoffEnergy(10000.);
        SapphireInput::Reaction("25Mg+a");        
        SapphireInput::EnergyFile("/examples/energyFile");
        SapphireInput::ReactionFile("/examples/reactionFile");
        SapphireInput::MassTable(sourceDirectory()+"/tables/masses/masses.dat");
        SapphireInput::GdrParams(sourceDirectory()+"/tables/gamma/ripl3_gdr_parameters.dat");
        SapphireInput::LevelDir(sourceDirectory()+"/tables/levels/");
        SapphireInput::SpinFile(sourceDirectory()+"/tables/spinod.dat");        
        SapphireInput::Isotope("60Ni");
        SapphireInput::PorterThomas_p(false);
        SapphireInput::PorterThomas_g(false);        
        SapphireInput::PreEq(false);
        SapphireInput::PreEqConf("");
        SapphireInput::Spin(1.0);
        SapphireInput::Parity(-1);
        SapphireInput::LowEnergy(6.0);
        SapphireInput::HighEnergy(6.0);
        SapphireInput::Events(100000);
        SapphireInput::ChunkSize(10000);
        SapphireInput::AlphaChannel(true);
        SapphireInput::ProtonChannel(true);
        SapphireInput::NeutronChannel(true);
        SapphireInput::GammaChannel(true);
        SapphireInput::Suffix(0);
        SapphireInput::CTable(0);
        SapphireInput::RdZ(50);
        SapphireInput::RdA(120);
        SapphireInput::RdEmin(0);
        SapphireInput::RdEmax(0);
        SapphireInput::RdOutputFile("levelScheme.dat");
        SapphireInput::RdMode("create");
    }

    void SapphireInput::PrintIntputFile(std::string InputFile){
        /** 
        * Printing the inputfile as read by Sapphire. 
        * The InputFile is read in the same way here as for reading
        * the inputFile. So, if something is off with the input,
        * one can possibly see it from the output of this function.
        */
        std::cout << std::endl;
        std::cout << "INPUT File" << std::endl;
        std::cout << std::endl;
        std::string line;
        std::ifstream myfile (InputFile);
        if (myfile.is_open())
        {
            while ( getline (myfile,line) )
            {
                std::cout << line << std::endl;
            }
            
            myfile.close();
        }
    }

    void SapphireInput::PrintIntputParameters(std::string module){
        /**
        * This function prints out the current status of the input parameters,
        * in a way, which can be used as input file template for other calculations.
        * In comparison to a given input file, it can be double-checked if the 
        * parameters have been set as intended. 
        */
        std::cout << std::endl;
        std::cout << "INPUT PARAMETERS" << std::endl;
        std::cout << std::endl;

        std::cout << "[General]" << std::endl;
        std::cout << "\tMassTable            = "             << SapphireInput::MassTable() << std::endl;
        std::cout << "\tGDRParams            = "             << SapphireInput::GdrParams() << std::endl;
        std::cout << "\tLeveldir             = "             << SapphireInput::LevelDir() << std::endl;
        std::cout << "\tSpinFile             = "             << SapphireInput::SpinFile() << std::endl;
        std::cout << "\tProtonModel          = "             << SapphireInput::p_Formalism() << std::endl;
        std::cout << "\tNeutronModel         = "             << SapphireInput::n_Formalism() << std::endl;
        std::cout << "\tAlphaModel           = "             << SapphireInput::a_Formalism() << std::endl;
        std::cout << "\tGammaModel           = "             << SapphireInput::g_Formalism() << std::endl;
        std::cout << "\tPorterThomasParticle = "             << SapphireInput::PorterThomas_p() << std::endl;
        std::cout << "\tPorterThomasGamma    = "             << SapphireInput::PorterThomas_g() << std::endl;
        std::cout << "\tGammaChannel         = "             << SapphireInput::GammaChannel() << std::endl;
        std::cout << "\tNeutronChannel       = "             << SapphireInput::NeutronChannel() << std::endl;
        std::cout << "\tProtonChannel        = "             << SapphireInput::ProtonChannel() << std::endl;
        std::cout << "\tAlphaChannel         = "             << SapphireInput::AlphaChannel() << std::endl;
        std::cout << "\tLevelDensity         = "             << SapphireInput::LevelDensity() << std::endl;
        std::cout << "\tcTable               = "             << SapphireInput::CTable() << std::endl;
        std::cout << std::endl;
        
        if(module=="CrossSection" || module=="all"){
            
            std::cout << "[CrossSection]" << std::endl;
            std::cout << "\tReaction             = "             << SapphireInput::Reaction() << std::endl;
            std::cout << "\tEnergyFile           = "             << SapphireInput::EnergyFile() << std::endl;
            std::cout << "\tReactionFile         = "             << SapphireInput::ReactionFile() << std::endl;
            std::cout << "\tCalcRates            = "             << SapphireInput::CalcRates() << std::endl;
            std::cout << "\tCalcAverageWidth     = "             << SapphireInput::CalcAverageWidth() << std::endl;
            std::cout << "\tResidualGamma        = "             << SapphireInput::ResidualGamma() << std::endl;
            std::cout << "\tResidualNeutron      = "             << SapphireInput::ResidualNeutron() << std::endl;
            std::cout << "\tResidualProton       = "             << SapphireInput::ResidualProton() << std::endl;
            std::cout << "\tResidualAlpha        = "             << SapphireInput::ResidualAlpha() << std::endl; 
            std::cout << "\tCalculateGammaCutoff = "             << SapphireInput::CalculateGammaCutoff() << std::endl;
            std::cout << "\tEntranceState        = "             << SapphireInput::EntranceState() << std::endl;
            std::cout << "\tg_ExitStates         = "             << SapphireInput::g_ExitStates() << std::endl;
            std::cout << "\tn_ExitStates         = "             << SapphireInput::n_ExitStates() << std::endl;
            std::cout << "\tp_ExitStates         = "             << SapphireInput::p_ExitStates() << std::endl;
            std::cout << "\ta_ExitStates         = "             << SapphireInput::a_ExitStates() << std::endl;
            std::cout << std::endl;
        }
        if(module=="Decayer" || module=="all"){
            std::cout << "[Decayer]" << std::endl;    
            std::cout << "\tSuffix               = "             << SapphireInput::Suffix() << std::endl;
            std::cout << "\tIsotope              = "             << SapphireInput::Isotope() << std::endl;
            std::cout << "\tSpin                 = "             << SapphireInput::Spin() << std::endl;
            std::cout << "\tParity               = "             << SapphireInput::Parity() << std::endl;
            std::cout << "\tEnergyLow            = "             << SapphireInput::LowEnergy() << std::endl;
            std::cout << "\tEnergyHigh           = "             << SapphireInput::HighEnergy() << std::endl;
            std::cout << "\tEvents               = "             << SapphireInput::Events() << std::endl;
            std::cout << "\tChunkSize            = "             << SapphireInput::ChunkSize() << std::endl;        
            std::cout << "\tPreequillibrium      = "             << SapphireInput::PreEq() << std::endl;
            std::cout << "\tPreEqConfiguration   = "             << SapphireInput::PreEqConf() << std::endl;
            std::cout << std::endl;
        }
        if(module=="Random" || module=="all"){
            std::cout << "[Random]" << std::endl;
            std::cout << "\tZ                    = "             << SapphireInput::RdZ() << std::endl;
            std::cout << "\tA                    = "             << SapphireInput::RdA() << std::endl;
            std::cout << "\tMode                 = "             << SapphireInput::RdMode() << std::endl;
            std::cout << "\tOutputFile           = "             << SapphireInput::RdOutputFile() << std::endl;
            std::cout << "\tEmin                 = "             << SapphireInput::RdEmin() << std::endl;
            std::cout << "\tEmax                 = "             << SapphireInput::RdEmax() << std::endl;
            std::cout << std::endl;
        }
    }

    void SapphireInput::ReadInputFile(std::string InputFile){
        /**
        */
        std::cout << "Reading input file ..." << InputFile << std::endl;
        boost::property_tree::ptree pt;

        try{
            boost::property_tree::ini_parser::read_ini(InputFile, pt);
        }
        catch (const boost::exception& e)
        {
            std::string diag = diagnostic_information(e);                   
            std::cout << "Can't init settings." << diag << std::endl;
            exit(1);
        }

        //Reading General Input
        SapphireInput::MassTable(pt.get<std::string>("General.MassTable", SapphireInput::MassTable()));
        SapphireInput::GdrParams(pt.get<std::string>("General.GDRParams", SapphireInput::GdrParams()));
        SapphireInput::LevelDir(pt.get<std::string>("General.LevelDir", SapphireInput::LevelDir()));
        SapphireInput::SpinFile(pt.get<std::string>("General.SpinFile", SapphireInput::SpinFile()));
        SapphireInput::p_Formalism(pt.get<int>("General.ProtonModel", SapphireInput::p_Formalism()));
        SapphireInput::n_Formalism(pt.get<int>("General.NeutronModel", SapphireInput::n_Formalism()));
        SapphireInput::a_Formalism(pt.get<int>("General.AlphaModel", SapphireInput::a_Formalism()));
        SapphireInput::g_Formalism(pt.get<int>("General.GammaModel", SapphireInput::g_Formalism()));
        SapphireInput::LevelDensity(pt.get<int>("General.LevelDensity", SapphireInput::LevelDensity()));
        SapphireInput::CTable(pt.get<double>("General.cTable", SapphireInput::CTable()));
        SapphireInput::PorterThomas_p(pt.get<bool>("General.PorterThomasParticle", SapphireInput::PorterThomas_p()));
        SapphireInput::PorterThomas_g(pt.get<bool>("General.PorterThomasGamma", SapphireInput::PorterThomas_g()));
        SapphireInput::GammaChannel(pt.get<bool>("General.GammaChannel", SapphireInput::GammaChannel()));
        SapphireInput::AlphaChannel(pt.get<bool>("General.AlphaChannel", SapphireInput::AlphaChannel()));
        SapphireInput::ProtonChannel(pt.get<bool>("General.ProtonChannel", SapphireInput::ProtonChannel()));
        SapphireInput::NeutronChannel(pt.get<bool>("General.NeutronChannel", SapphireInput::NeutronChannel()));

        //Reading CrossSection Input
        SapphireInput::Reaction(pt.get<std::string>("CrossSection.Reaction", SapphireInput::Reaction()));
        SapphireInput::Energies(pt.get<std::string>("CrossSection.Energies", SapphireInput::Energies()));
        SapphireInput::EnergyFile(pt.get<std::string>("CrossSection.EnergyFile", SapphireInput::EnergyFile()));
        SapphireInput::ReactionFile(pt.get<std::string>("CrossSection.ReactionFile", SapphireInput::ReactionFile()));
        SapphireInput::CalcRates(pt.get<bool>("CrossSection.CalcRates", SapphireInput::CalcRates()));
        SapphireInput::CalcAverageWidth(pt.get<bool>("CrossSection.CalcAverageWidth", SapphireInput::CalcAverageWidth()));
        SapphireInput::ResidualGamma(pt.get<bool>("CrossSection.ResidualGamma", SapphireInput::ResidualGamma()));               
        SapphireInput::ResidualNeutron(pt.get<bool>("CrossSection.ResidualNeutron", SapphireInput::ResidualNeutron()));           
        SapphireInput::ResidualProton(pt.get<bool>("CrossSection.ResidualProton", SapphireInput::ResidualProton()));
        SapphireInput::ResidualAlpha(pt.get<bool>("CrossSection.ResidualAlpha", SapphireInput::ResidualAlpha()));
        SapphireInput::CalculateGammaCutoff(pt.get<bool>("CrossSection.CalculateGammaCutoff", SapphireInput::CalculateGammaCutoff()));

        //Reading Decayer Input
        SapphireInput::Suffix(pt.get<int>("Decayer.Suffix", SapphireInput::Suffix()));
        SapphireInput::Isotope(pt.get<std::string>("Decayer.Isotope", SapphireInput::Isotope()));
        SapphireInput::Spin(pt.get<double>("Decayer.Spin", SapphireInput::Spin()));
        SapphireInput::Parity(pt.get<double>("Decayer.Parity", SapphireInput::Parity()));
        SapphireInput::LowEnergy(pt.get<double>("Decayer.EnergyLow", SapphireInput::LowEnergy()));
        SapphireInput::HighEnergy(pt.get<double>("Decayer.EnergyHigh", SapphireInput::HighEnergy()));
        SapphireInput::Events(pt.get<int>("Decayer.Events", SapphireInput::Events()));
        SapphireInput::ChunkSize(pt.get<int>("Decayer.ChunkSize", SapphireInput::ChunkSize()));
        SapphireInput::PreEq(pt.get<bool>("Decayer.Preequillibrium", SapphireInput::PreEq()));
        SapphireInput::PreEqConf(pt.get<std::string>("Decayer.PreEqConfiguration", SapphireInput::PreEqConf()));

        //Reading Random Input
        SapphireInput::RdZ(pt.get<int>("Random.Z", SapphireInput::RdZ()));
        SapphireInput::RdA(pt.get<int>("Random.A", SapphireInput::RdA()));
        SapphireInput::RdMode(pt.get<std::string>("Random.Mode", SapphireInput::RdMode()));
        SapphireInput::RdOutputFile(pt.get<std::string>("Random.OutputFile", SapphireInput::RdOutputFile()));
        SapphireInput::RdEmax(pt.get<double>("Random.Emax", SapphireInput::RdEmax()));
        SapphireInput::RdEmin(pt.get<double>("Random.Emin", SapphireInput::RdEmin()));

        //Initialize Datatables
        //SapphireInput::masses.InitializeElements();
        //SapphireInput::masses.InitializeMasses(SapphireInput::MassTable());
        //SapphireInput::levels.InitializeLevels(SapphireInput::LevelDir(), SapphireInput::SpinFile());
    }

    std::string SapphireInput::pTypeStringFromReactionString(std::string reactionString){
        std::string massNumberString=SapphireInput::massNumberStringFromReactionString(reactionString);
        std::string atomicNumberString=SapphireInput::atomicNumberStringFromReactionString(reactionString);

        reactionString.erase(0,massNumberString.length());
        reactionString.erase(0,atomicNumberString.length());
      
        std::string projectileString;
        for(unsigned int i = 0; i<reactionString.length(); i++) {
            std::string nextChar(reactionString,i,1);
            if(nextChar=="+") continue;
            else projectileString+=nextChar;
        }

        return projectileString;
    }

    int SapphireInput::pTypeIntFromReactionString(std::string reactionString){
        std::string projectileString = SapphireInput::pTypeStringFromReactionString(reactionString);

        if(projectileString=="g")
            return 0;
        else if(projectileString=="n") 
            return 1;      
        else if(projectileString=="p")
            return 2;
        else if(projectileString=="a")
            return 3;
        
        /** Return "neutron" by default*/
        return 1;
    }

    std::string SapphireInput::massNumberStringFromReactionString(std::string reactionString){
        std::string massNumberS;
        for(unsigned int i = 0; i<reactionString.length(); i++) {
            std::string nextChar(reactionString,i,1);
            std::istringstream stm(nextChar);
            int nextDigit;
            if(!(stm>>nextDigit)) break;
                else massNumberS+=nextChar;
        }

        return massNumberS;
    }

    int SapphireInput::massNumberIntFromReactionString(std::string reactionString){
        std::string massNumberS = SapphireInput::massNumberStringFromReactionString(reactionString);
        if(massNumberS.length()>0)
            return atoi(massNumberS.c_str());
        else
            return 0;
    }

    std::string SapphireInput::atomicNumberStringFromReactionString(std::string reactionString){
        std::string massNumberS = SapphireInput::massNumberStringFromReactionString(reactionString);
        reactionString.erase(0,massNumberS.length());
        std::string atomicNumberString;
        for(unsigned int i = 0; i<reactionString.length(); i++) {
            std::string nextChar(reactionString,i,1);
            if(nextChar=="+") break;
                else atomicNumberString+=nextChar;
            }        
        return atomicNumberString;
    }

    int SapphireInput::atomicNumberIntFromReactionString(std::string reactionString){
        std::string atomicNumberString = SapphireInput::atomicNumberStringFromReactionString(reactionString);
        if(NuclearMass::FindZ(atomicNumberString) != -1) 
            return NuclearMass::FindZ(atomicNumberString);
        else
            return 0;
    }

    std::string SapphireInput::massNumberStringFromIsotopeString(std::string isotopeString){
        std::string massNumberString;

        for(unsigned int i = 0; i<isotopeString.length(); i++) {
      	    std::string nextChar(isotopeString,i,1);
            std::istringstream stm(nextChar);
            int nextDigit;
            if(!(stm>>nextDigit)) break;
            else massNumberString+=nextChar;
      }
      return massNumberString;
    }

    int SapphireInput::massNumberIntFromIsotopeString(std::string isotopeString){
        std::string massNumberString = SapphireInput::massNumberStringFromIsotopeString(isotopeString);
        if(massNumberString.length()>0)
            return atoi(massNumberString.c_str());
        else
            return 0;
    }

    std::string SapphireInput::atomicNumberStringFromIsotopeString(std::string isotopeString){
        std::string massNumberString = SapphireInput::massNumberStringFromIsotopeString(isotopeString);
        std::string atomicNumberString = isotopeString.substr(massNumberString.length());
        return atomicNumberString;
    }

    int SapphireInput::atomicNumberIntFromIsotopeString(std::string isotopeString){
        std::string atomicNumberString = SapphireInput::atomicNumberStringFromIsotopeString(isotopeString);
        int Z = NuclearMass::FindZ(atomicNumberString);
        return Z;
    }

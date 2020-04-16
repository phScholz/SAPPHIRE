/**
 * @file SapphireInput.h
 * @brief Option class for Sapphire.
 * @date 2020
 * @author Philipp Scholz, <pscholz@outlook.de>
 * 
 */

#include "SapphireInput.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/exception/all.hpp>
#include <boost/throw_exception.hpp>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <vector>


extern std::string sourceDirectory();
   
    SapphireInput::SapphireInput(){
        SapphireInput::Initialize();
    }

    void SapphireInput::Initialize(){
        std::cout << "Setting default values..." << std::endl;
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

        SapphireInput::DecayerMaxL(8.);
        SapphireInput::PreEqMaxL(8.);
        SapphireInput::g_CutoffEnergy(10000.);

        SapphireInput::Reaction("25Mg+a");
        SapphireInput::Energies("25Mg+a");
        SapphireInput::EnergyFile("");
        SapphireInput::ReactionFile("");
        SapphireInput::MassTable(sourceDirectory()+"/tables/masses.dat");
        SapphireInput::GdrParams(sourceDirectory()+"/tables/ripl3_gdr_parameters.dat");
        SapphireInput::LevelDir(sourceDirectory()+"/levels/");
        SapphireInput::SpinFile(sourceDirectory()+"/tables/spinod.dat");
        
    }

    void SapphireInput::printIntputFile(std::string InputFile){
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

    void SapphireInput::printIntputParameters(){
        std::cout << std::endl;
        std::cout << "INPUT PARAMETERS" << std::endl;
        std::cout << std::endl;

        std::cout << "\t[General]" << std::endl;
        std::cout << "\tMassTable            = "             << SapphireInput::MassTable() << std::endl;
        std::cout << "\tGDRParams            = "             << SapphireInput::GdrParams() << std::endl;
        std::cout << "\tLeveldir             = "             << SapphireInput::LevelDir() << std::endl;
        std::cout << "\tSpinFile             = "             << SapphireInput::SpinFile() << std::endl;
        std::cout << "\tprotonOMP            = "             << SapphireInput::p_Formalism() << std::endl;
        std::cout << "\tneutronOMP           = "             << SapphireInput::n_Formalism() << std::endl;
        std::cout << "\talphaOMP             = "             << SapphireInput::a_Formalism() << std::endl;
        std::cout << "\tstrength             = "             << SapphireInput::g_Formalism() << std::endl;
        
        std::cout << std::endl;
        std::cout << "\t[CrossSection]" << std::endl;
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
        std::cout << "\t[Decayer]" << std::endl;
        std::cout << "\tSuffix               = "             << SapphireInput::Suffix() << std::endl;
    }

    void SapphireInput::ReadInputFile(std::string InputFile){
        std::cout << "Reading input file ..." << InputFile << std::endl;
        boost::property_tree::ptree pt;

        try{
            boost::property_tree::ini_parser::read_ini(InputFile, pt);
        }
        catch (const boost::exception& e)
        {
            std::string diag = diagnostic_information(e);
            // display your error message here, then do whatever you need to, e.g.        
            std::cout << "Can't init settings." << diag << std::endl;
            exit(1);
        }

        //Reading General Input
        SapphireInput::MassTable(pt.get<std::string>("General.MassTable", SapphireInput::MassTable()));
        SapphireInput::GdrParams(pt.get<std::string>("General.GDRParams", SapphireInput::GdrParams()));
        SapphireInput::LevelDir(pt.get<std::string>("General.LevelDir", SapphireInput::LevelDir()));
        SapphireInput::SpinFile(pt.get<std::string>("General.SpinFile", SapphireInput::SpinFile()));
        SapphireInput::SpinFile(pt.get<std::string>("General.SpinFile", SapphireInput::SpinFile()));
        
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
        SapphireInput::Suffix(pt.get<std::string>("Decayer.Suffix", SapphireInput::Suffix()));
    }


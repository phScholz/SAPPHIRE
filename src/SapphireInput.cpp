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
        std::cout << "\t\tMassTable            = "             << SapphireInput::MassTable() << std::endl;
        std::cout << "\t\tGDRParams            = "             << SapphireInput::GdrParams() << std::endl;
        std::cout << "\t\tLeveldir             = "             << SapphireInput::LevelDir() << std::endl;
        std::cout << "\t\tSpinFile             = "             << SapphireInput::SpinFile() << std::endl;
        std::cout << "\t\tSuffix               = "             << SapphireInput::Suffix() << std::endl;
        std::cout << std::endl;
        std::cout << "\t[CrossSection]" << std::endl;
        std::cout << "\t\tReaction             = "             << SapphireInput::Reaction() << std::endl;
        std::cout << "\t\tEnergyFile           = "             << SapphireInput::EnergyFile() << std::endl;
        std::cout << "\t\tReactionFile         = "             << SapphireInput::ReactionFile() << std::endl;
        std::cout << "\t\tCalcRates            = "             << SapphireInput::CalcRates() << std::endl;
        std::cout << "\t\tCalcAverageWidth     = "             << SapphireInput::CalcAverageWidth() << std::endl;
        std::cout << "\t\tResidualGamma        = "             << SapphireInput::ResidualGamma() << std::endl;
        std::cout << "\t\tResidualNeutron      = "             << SapphireInput::ResidualNeutron() << std::endl;
        std::cout << "\t\tResidualProton       = "             << SapphireInput::ResidualProton() << std::endl;
        std::cout << "\t\tResidualAlpha        = "             << SapphireInput::ResidualAlpha() << std::endl;
        std::cout << "\t\tCalculateGammaCutoff = "             << SapphireInput::CalculateGammaCutoff() << std::endl;
        std::cout << "\t\tEntranceState        = "             << SapphireInput::EntranceState() << std::endl;
        std::cout << "\t\tg_ExitStates         = "             << SapphireInput::g_ExitStates() << std::endl;
        std::cout << "\t\tn_ExitStates         = "             << SapphireInput::n_ExitStates() << std::endl;
        std::cout << "\t\tp_ExitStates         = "             << SapphireInput::p_ExitStates() << std::endl;
        std::cout << "\t\ta_ExitStates         = "             << SapphireInput::a_ExitStates() << std::endl;
    }

    void SapphireInput::ReadInputFile(std::string InputFile){
        SapphireInput::printIntputFile(InputFile);

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
        SapphireInput::Suffix(pt.get<std::string>("General.Suffix", SapphireInput::Suffix()));
        //Reading CrossSection Input
        SapphireInput::Reaction(pt.get<std::string>("CrossSection.Reaction", SapphireInput::Reaction()));
        SapphireInput::EnergyFile(pt.get<std::string>("CrossSection.EnergyFile", SapphireInput::EnergyFile()));
        SapphireInput::ReactionFile(pt.get<std::string>("CrossSection.ReactionFile", SapphireInput::ReactionFile()));
        SapphireInput::CalcRates(pt.get<bool>("CrossSection.CalcRates", SapphireInput::CalcRates()));
        SapphireInput::CalcAverageWidth(pt.get<bool>("CrossSection.CalcAverageWidth", SapphireInput::CalcAverageWidth()));
        SapphireInput::ResidualGamma(pt.get<bool>("CrossSection.ResidualGamma", SapphireInput::ResidualGamma()));               
        SapphireInput::ResidualNeutron(pt.get<bool>("CrossSection.ResidualNeutron", SapphireInput::ResidualNeutron()));           
        SapphireInput::ResidualProton(pt.get<bool>("CrossSection.ResidualProton", SapphireInput::ResidualProton()));
        SapphireInput::ResidualAlpha(pt.get<bool>("CrossSection.ResidualAlpha", SapphireInput::ResidualAlpha()));
        SapphireInput::CalculateGammaCutoff(pt.get<bool>("CrossSection.CalculateGammaCutoff", SapphireInput::CalculateGammaCutoff()));

    }


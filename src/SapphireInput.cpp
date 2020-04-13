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
#include <iostream>
#include <fstream>


   
    SapphireInput::SapphireInput(){
        SapphireInput::Initialize();
    }

    void SapphireInput::Initialize(){

        SapphireInput::CalcRates(false);           
        SapphireInput::CalcAverageWidth(false);
        SapphireInput::ResiudalGamma(false);               
        SapphireInput::ResiudalNeutron(false);           
        SapphireInput::ResiudalProton(false);
        SapphireInput::ResiudalAlpha(false);
        SapphireInput::CalculateGammaCutoff(false);

        SapphireInput::EntranceState(0);

        SapphireInput::EnergyFile("");
        SapphireInput::ReactionFile("");
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

        std::cout << "\tCrossSection.CalcRates            = "             << SapphireInput::CalcRates() << std::endl;
        std::cout << "\tCrossSection.CalcAverageWidth     = "      << SapphireInput::CalcAverageWidth() << std::endl;
        std::cout << "\tCrossSection.ResidualGamma        = "         << SapphireInput::ResidualGamma() << std::endl;
        std::cout << "\tCrossSection.ResidualNeutron      = "       << SapphireInput::ResidualNeutron() << std::endl;
        std::cout << "\tCrossSection.ResidualProton       = "        << SapphireInput::ResidualProton() << std::endl;
        std::cout << "\tCrossSection.ResidualAlpha        = "         << SapphireInput::ResidualAlpha() << std::endl;
        std::cout << "\tCrossSection.CalculateGammaCutoff = "  << SapphireInput::CalculateGammaCutoff() << std::endl;
    }

    void SapphireInput::ReadInputFile(std::string InputFile){
        SapphireInput::printIntputFile(InputFile);

        boost::property_tree::ptree pt;
        boost::property_tree::ini_parser::read_ini(InputFile, pt);

        SapphireInput::CalcRates(pt.get<bool>("CrossSection.CalcRates", SapphireInput::CalcRates()));
        SapphireInput::CalcAverageWidth(pt.get<bool>("CrossSection.CalcAverageWidth", SapphireInput::CalcAverageWidth()));
        SapphireInput::ResiudalGamma(pt.get<bool>("CrossSection.ResidualGamma", SapphireInput::ResidualGamma()));               
        SapphireInput::ResiudalNeutron(pt.get<bool>("CrossSection.ResidualNeutron", SapphireInput::ResidualNeutron()));           
        SapphireInput::ResiudalProton(pt.get<bool>("CrossSection.ResidualProton", SapphireInput::ResidualProton()));
        SapphireInput::ResiudalAlpha(pt.get<bool>("CrossSection.ResidualAlpha", SapphireInput::ResidualAlpha()));
        SapphireInput::CalculateGammaCutoff(pt.get<bool>("CrossSection.CalculateGammaCutoff", SapphireInput::CalculateGammaCutoff()));

    }


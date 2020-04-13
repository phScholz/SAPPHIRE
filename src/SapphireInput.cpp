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

    void SapphireInput::ReadInputFile(std::string InputFile){
        SapphireInput::printIntputFile(InputFile);

        std::cout << std::endl;
        std::cout << "INPUT PARAMETERS" << std::endl;
        std::cout << std::endl;
        boost::property_tree::ptree pt;
        boost::property_tree::ini_parser::read_ini(InputFile, pt);

        SapphireInput::CalcRates(pt.get<bool>("CrossSection.CalcRates", SapphireInput::CalcRates()));
        std::cout << "CrossSection.CalcRates = " << SapphireInput::CalcRates() << std::endl;
        SapphireInput::CalcAverageWidth(pt.get<bool>("CrossSection.CalcRates", SapphireInput::CalcAverageWidth()));
        std::cout << "CrossSection.CalcAverageWidth = " << pt.get<bool>("CrossSection.CalcAverageWidth", SapphireInput::CalcAverageWidth()) << std::endl;
        std::cout << "CrossSection.ResidualNeutron = " << pt.get<bool>("CrossSection.ResidualNeutron", SapphireInput::ResidualNeutron()) << std::endl;                                
    }


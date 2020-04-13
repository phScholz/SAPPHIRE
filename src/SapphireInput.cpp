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

    void SapphireInput::ReadInputFile(std::string InputFile){
        std::cout << std::endl;
        std::cout << "INPUT PARAMETERS" << std::endl;
        std::cout << std::endl;
        boost::property_tree::ptree pt;
        boost::property_tree::ini_parser::read_ini(InputFile, pt);
        std::cout << "CrossSection.CalcRates" << pt.get<bool>("CrossSection.CalcRates", SapphireInput::CalcRates()) << std::endl;
        std::cout << "CrossSection.CalcAverageWidth" << pt.get<bool>("CrossSection.CalcAverageWidth", SapphireInput::CalcAverageWidth()) << std::endl;
        std::cout << "CrossSection.ResidualNeutron" << pt.get<bool>("CrossSection.ResidualNeutron", SapphireInput::ResidualNeutron()) << std::endl;                                
    }


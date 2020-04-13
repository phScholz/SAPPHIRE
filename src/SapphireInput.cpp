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
        boost::property_tree::ptree pt;
        boost::property_tree::ini_parser::read_ini(InputFile, pt);
        std::cout << pt.get<std::string>("CrossSection.CalcRates","False") << std::endl;
        std::cout << pt.get<std::string>("CrossSection.CalcAverageWidth", "False") << std::endl;
        std::cout << pt.get<std::string>("CrossSection.ResidualNeutron", "False") << std::endl;                                
    }


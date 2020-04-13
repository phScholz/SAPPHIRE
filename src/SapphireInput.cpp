/**
 * @file SapphireInput.h
 * @brief Option class for Sapphire.
 * @date 2020
 * @author Philipp Scholz, <pscholz@outlook.de>
 * 
 */

#include "SapphireInput.h"

namespace SapphireInput{
    
    void SapphireInput(){
        Initialize();
    }

    void Initialize(){

        CalcRates(false);           
        CalcAverageWidth(false);
        ResiudalGamma(false);               
        ResiudalNeutron(false);           
        ResiudalProton(false);
        ResiudalAlpha(false);
        CalculateGammaCutoff(false);

        EntranceState(0);
        
        EnergyFile("");
        ReactionFile("");
    }


}

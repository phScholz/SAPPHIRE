#include "gtest/gtest.h"
#include "SapphireInput.h"

#include "Databases/NuclearMass.h"
#include "Decayer/Decayer.h"
#include "RandomScheme.h"
#include "CrossSection.h"
#include "GammaStrength/GammaTransmissionFunc.h"
#include "ParticleTransmissionFunc.h"
#include "TransitionRateFunc.h"
#include "LevelDensity/LevelDensityTable.h"
#include "LevelDensity/LevelDensityHFB_BSk14.h"

#include <iostream>
#include <iomanip>

extern std::string sourceDirectory();
extern void Initialize();

TEST(NuclearMass, FindMass){
    double m=0;
    for(int charge=1; charge <= 90; charge++){
        for(int mass=2.5*charge; mass<=3*charge; mass++){
            EXPECT_EQ(true, NuclearMass::FindMass(charge, mass, m));
            //std::cout << std::setw(10) << charge << std::setw(10) <<mass << std::setw(10) << m << std::endl;
        }
    }
}

TEST(NuclearMass, MassDifference){
    double diff1=0;
    double diff2=0;
    for(int charge1=1; charge1 <= 90; charge1++){
        for(int charge2=1; charge2 <= 90; charge2++){
            for(int mass1=2.5*charge1; mass1<=3*charge1; mass1++){
                for(int mass2=2.5*charge2; mass2<=3*charge2; mass2++){
                    EXPECT_EQ(true, NuclearMass::MassDifference(charge1,mass1,charge2,mass2,diff1));
                    EXPECT_EQ(true, NuclearMass::MassDifference(charge2,mass2,charge1,mass1,diff2));
                    EXPECT_FLOAT_EQ(diff1,-diff2);
                }
            }
        }
    }
}


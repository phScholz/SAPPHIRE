#include "LevelDensity/LevelDensityHFB_BSk14.h"
#include "LevelDensity/RauscherLevelDensity.h"
#include "LevelDensity/LevelDensityTable.h"
#include "gtest/gtest.h"

extern void Initialize();

TEST(LevelDensity, HFB_BSk14_Reading){
    Initialize();

    //AC
    auto * NLD = new LevelDensityHFB_BSk14(89,184,0,1);    
    EXPECT_EQ(round(NLD->CalculateDensity(0.5)*100)/100.0, 1.35);
    
    NLD = new LevelDensityHFB_BSk14(89,288,3,-1);    
    EXPECT_FLOAT_EQ(NLD->CalculateDensity(9.0), 4.27e8);
    
    NLD = new LevelDensityHFB_BSk14(89,184,1,1);
    EXPECT_EQ(round(NLD->CalculateDensity(0.5)*1000)/1000.0, 7.27e-1);
    
    //AG
    NLD = new LevelDensityHFB_BSk14(47,84,1,1);
    EXPECT_EQ(round(NLD->CalculateDensity(0.5)*100)/100.0, 0);

    NLD = new LevelDensityHFB_BSk14(47,85,0.5,-1);
    EXPECT_EQ(round(NLD->CalculateDensity(1.5)*100)/100.0, 1.2);    
    EXPECT_EQ(NLD->operator()(1.5), 1.2);

    //AL
    NLD = new LevelDensityHFB_BSk14(13,21,0.5,1);    
    EXPECT_FLOAT_EQ(NLD->CalculateDensity(9.0), 4.44);
    
    NLD = new LevelDensityHFB_BSk14(13,48,3,-1);    
    EXPECT_FLOAT_EQ(NLD->CalculateDensity(9.0), 4.57e2);
        
}

TEST(LevelDensity, RauscherLevelDensity){
    for(int charge=8; charge<=100; charge++){
        for(int mass=charge; mass <= charge*3; mass++){\
            if(mass%2){
                for(double spin=0.5; spin<=7.5; spin = spin + 0.5){
                    auto * NLD1 = new RauscherLevelDensity(charge,mass,spin,1);
                    auto * NLD2 = new RauscherLevelDensity(charge,mass,spin,-1);
                    
                    for(double energy = 0.5; energy <=20.0; energy=energy+1.0){
                        EXPECT_GT(NLD1->operator()(energy), 0.0);
                        EXPECT_EQ(NLD1->operator()(energy),NLD2->operator()(energy));
                    }
                    
                }
            }
            if(!mass%2){
                for(double spin=0; spin<=8.0; spin = spin + 1){
                    auto * NLD1 = new RauscherLevelDensity(charge,mass,spin,1);
                    auto * NLD2 = new RauscherLevelDensity(charge,mass,spin,-1);
                    for(double energy = 0.5; energy <=20.0; energy=energy+1.0){
                        EXPECT_GT(NLD1->operator()(energy), 0.0);
                        EXPECT_EQ(NLD1->operator()(energy),NLD2->operator()(energy));
                    }
                }
            }
        }
    }
}
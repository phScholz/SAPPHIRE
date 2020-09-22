#include "LevelDensity/LevelDensityHFB_BSk14.h"
#include "LevelDensity/LevelDensityTable.h"
#include "gtest/gtest.h"

TEST(LevelDensity, HFB_BSk14_Reading){

    auto * NLD = new LevelDensityHFB_BSk14(89,184,0,1);
    
    //XPECT_EQ(NLD->CalculateDensity(0.25), 9.03e-03);
    //
    //LD = new LevelDensityHFB_BSk14(89,184,1,1);
   
    //XPECT_EQ(NLD->CalculateDensity(0.25), 4.42e-02);
    //
    //LD = new LevelDensityHFB_BSk14(47,84,1,1);
   
    //XPECT_EQ(NLD->CalculateDensity(0.25), 0);

    NLD = new LevelDensityHFB_BSk14(47,85,0.5,-1);

    EXPECT_EQ(NLD->CalculateDensity(1.5), 1.2);    

}
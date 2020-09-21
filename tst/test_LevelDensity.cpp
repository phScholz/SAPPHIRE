#include "LevelDensity/LevelDensityHFB_BSk14.h"
#include "LevelDensity/LevelDensityTable.h"
#include "gtest/gtest.h"

TEST(LevelDensity, HFB_BSk14_Reading){

    auto * NLD1 = new LevelDensityHFB_BSk14(89,184,0,1);
    
    EXPECT_EQ(NLD1->CalculateDensity(0.25), 9.03e-03);

    delete NLD1;

    auto * NLD2 = new LevelDensityHFB_BSk14(89,184,1,1);
   
    EXPECT_EQ(NLD2->CalculateDensity(0.25), 4.42e-02);

    delete NLD2;

    auto * NLD3 = new LevelDensityHFB_BSk14(47,84,1,1);
   
    EXPECT_EQ(NLD3->CalculateDensity(0.25), 0);

    NLD3 = new LevelDensityHFB_BSk14(47,85,0.5,-1);

    EXPECT_EQ(NLD3->CalculateDensity(1.5), 1.2);

    delete NLD3;

}
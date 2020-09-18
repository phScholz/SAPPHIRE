#include "gtest/gtest.h"
#include "SapphireInput.h"


TEST(SapphireInput, ExitStates){
    SapphireInput * input = new SapphireInput();
    input->a_ExitStates(2);
    EXPECT_EQ(2, input->a_ExitStates());

    input->g_ExitStates(2);
    EXPECT_EQ(2, input->a_ExitStates());

    input->n_ExitStates(2);
    EXPECT_EQ(2, input->a_ExitStates());

    input->p_ExitStates(2);
    EXPECT_EQ(2, input->a_ExitStates());
}

TEST(SapphireInput, PrintBool){
    SapphireInput * input = new SapphireInput();
    
    input->PrintTrans(true);
    EXPECT_EQ(true, input->PrintTrans());

    input->PrintXs(true);
    EXPECT_EQ(true, input->PrintXs());

    input->PrintRate(true);
    EXPECT_EQ(true, input->PrintRate());

    input->PrintMACS(true);
    EXPECT_EQ(true, input->PrintMACS());
}

TEST(SapphireInput, CalcBool){
    SapphireInput * input = new SapphireInput();
    
    input->CalcRates(true);
    EXPECT_EQ(true, input->CalcRates());

    input->CalcMACS(true);
    EXPECT_EQ(true, input->CalcMACS());

    input->CalcXS(true);
    EXPECT_EQ(true, input->CalcXS());

    input->CalcAverageWidth(true);
    EXPECT_EQ(true, input->CalcAverageWidth());
}

TEST(SapphireInput, ResidualBool){
    SapphireInput * input = new SapphireInput();
    
    input->ResidualGamma(true);
    EXPECT_EQ(true, input->ResidualGamma());

    input->ResidualNeutron(true);
    EXPECT_EQ(true, input->ResidualNeutron());

    input->ResidualProton(true);
    EXPECT_EQ(true, input->ResidualProton());

    input->ResidualAlpha(true);
    EXPECT_EQ(true, input->ResidualAlpha());
}
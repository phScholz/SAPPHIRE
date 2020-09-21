
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

extern std::string sourceDirectory();
extern void Initialize();

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

    delete input;
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

    delete input;
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

    input->CalculateGammaCutoff(true);
    EXPECT_EQ(true, input->CalculateGammaCutoff());

    delete input;
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

    delete input;
}

TEST(SapphireInput, PorterThomasBool){
    SapphireInput * input = new SapphireInput();
    
    input->PorterThomas_g(true);
    EXPECT_EQ(true, input->PorterThomas_g());

    input->PorterThomas_p(true);
    EXPECT_EQ(true, input->PorterThomas_p());

    delete input;
}

TEST(SapphireInput, Formalism){
    auto * input = new SapphireInput();

    input->a_Formalism(0);
    input->p_Formalism(0);
    input->n_Formalism(0);
    input->E1_Formalism(0);
    input->E2_Formalism(0);
    input->M1_Formalism(0);
    input->LevelDensity(0);

    EXPECT_EQ(0, input->a_Formalism());
    EXPECT_EQ(0, input->p_Formalism());
    EXPECT_EQ(0, input->n_Formalism());
    EXPECT_EQ(0, input->E1_Formalism());
    EXPECT_EQ(0, input->E2_Formalism());
    EXPECT_EQ(0, input->M1_Formalism());
    EXPECT_EQ(0, input->LevelDensity());

    input->a_Formalism(1);
    input->p_Formalism(2);
    input->n_Formalism(3);
    input->E1_Formalism(4);
    input->E2_Formalism(5);
    input->M1_Formalism(6);
    input->LevelDensity(7);

    EXPECT_EQ(1, input->a_Formalism());
    EXPECT_EQ(2, input->p_Formalism());
    EXPECT_EQ(3, input->n_Formalism());
    EXPECT_EQ(4, input->E1_Formalism());
    EXPECT_EQ(5, input->E2_Formalism());
    EXPECT_EQ(6, input->M1_Formalism());
    EXPECT_EQ(7, input->LevelDensity());

    delete input;
}

TEST(SapphireInput, Files){
    auto * input = new SapphireInput();


    delete input;
}

TEST(SapphireInput, Initialize){
    auto * input = new SapphireInput();
    input->Initialize();

    EXPECT_EQ(input->ResidualGamma(), false);
    EXPECT_EQ(input->ResidualProton(), false);
    EXPECT_EQ(input->ResidualAlpha(), false);
    EXPECT_EQ(input->ResidualNeutron(), false);
    
    EXPECT_EQ(input->g_ExitStates(), -1);
    EXPECT_EQ(input->n_ExitStates(), -1);
    EXPECT_EQ(input->p_ExitStates(), -1);
    EXPECT_EQ(input->a_ExitStates(), -1);

    EXPECT_EQ(input->PrintTrans(),  false);
    EXPECT_EQ(input->PrintXs(),  true);
    EXPECT_EQ(input->PrintRate(),  false);
    EXPECT_EQ(input->PrintMACS(),  false);

    EXPECT_EQ(input->CalcRates(),  false);
    EXPECT_EQ(input->CalcMACS(),  false);
    EXPECT_EQ(input->CalcXS(),  true);
    EXPECT_EQ(input->CalcAverageWidth(),  false);
    EXPECT_EQ(input->CalculateGammaCutoff(),  false);

    EXPECT_EQ(input->ResidualGamma(), false);
    EXPECT_EQ(input->ResidualNeutron(), false);
    EXPECT_EQ(input->ResidualAlpha(), false);
    EXPECT_EQ(input->ResidualProton(), false);

    EXPECT_EQ(input->PorterThomas_g(), false);
    EXPECT_EQ(input->PorterThomas_p(), false);
    
    EXPECT_EQ(input->a_Formalism() , 0);
    EXPECT_EQ(input->p_Formalism() , 0);
    EXPECT_EQ(input->n_Formalism() , 0);
    EXPECT_EQ(input->E1_Formalism(), 3);
    EXPECT_EQ(input->M1_Formalism(), 3);
    EXPECT_EQ(input->E2_Formalism(), 0);
    EXPECT_EQ(input->LevelDensity(), 1);
    
    EXPECT_EQ(input->EntranceState(), 0);
    EXPECT_EQ(input->DecayerMaxL(), 8.);
    EXPECT_EQ(input->PreEqMaxL(), 8.);
    EXPECT_EQ(input->g_CutoffEnergy(), 10000);

    EXPECT_EQ(input->Reaction(), "25Mg+a");
    EXPECT_EQ(input->ReactionFile(), "");
    EXPECT_EQ(input->EnergyFile(), "");

    EXPECT_EQ(input->MassTable(), sourceDirectory()+"/tables/masses/masses.dat");
    EXPECT_EQ(input->GdrParams(), sourceDirectory()+"/tables/gamma/ripl3_gdr_parameters.dat");
    EXPECT_EQ(input->LevelDir(), sourceDirectory()+"/tables/levels/");
    EXPECT_EQ(input->SpinFile(), sourceDirectory()+"/tables/spinod.dat");
    EXPECT_EQ(input->Isotope(), "60Ni");
        
    EXPECT_EQ(input->PreEq(),false);
    EXPECT_EQ(input->PreEqConf(),"");
    EXPECT_EQ(input->Spin(),1.0);
    EXPECT_EQ(input->Parity(), -1);
    EXPECT_EQ(input->LowEnergy(), 6.0);
    EXPECT_EQ(input->HighEnergy(), 6.0);
    EXPECT_EQ(input->Events(), 100000);
    EXPECT_EQ(input->ChunkSize(), 10000);
    EXPECT_EQ(input->AlphaChannel(),true);
    EXPECT_EQ(input->ProtonChannel(), true);
    EXPECT_EQ(input->NeutronChannel(), true);
    EXPECT_EQ(input->GammaChannel(), true);
    EXPECT_EQ(input->Suffix(), 0);
    EXPECT_EQ(input->CTable(), 0);
    EXPECT_EQ(input->RdZ(), 50);
    EXPECT_EQ(input->RdA(), 120);
    EXPECT_EQ(input->RdEmin(), 0);
    EXPECT_EQ(input->RdEmax(), 0);
    EXPECT_EQ(input->RdOutputFile(), "levelScheme.dat");
    EXPECT_EQ(input->RdMode(), "create");
    EXPECT_EQ(input->Gnorm(), 1.0);
    EXPECT_EQ(input->Pnorm(), 1.0);
    EXPECT_EQ(input->Nnorm(), 1.0);
    EXPECT_EQ(input->Anorm(), 1.0);
    EXPECT_EQ(input->DisFile(), "");
    EXPECT_EQ(input->EBinning(), 100);
    
    delete input;
}

TEST(SapphireInput, SetInputTransitionRate){
    auto * input = new SapphireInput;
    
    int nldmodel = 1;
    int eBinning = 50;
    double gCutOff = 100000;

    input->LevelDensity(nldmodel);
    input->EBinning(eBinning);
    input->g_CutoffEnergy(gCutOff);

    input->SetInputTransitionRate();

    auto * TRF = new TransitionRateFunc();
    EXPECT_EQ(TRF->NLDmodel(), nldmodel);
    EXPECT_EQ(TRF->Binning(), eBinning);
    EXPECT_EQ(TRF->GetGammaCutoffEnergy(), gCutOff);
}

TEST(SapphireInput, SetInputLevelDensity){
    auto * input = new SapphireInput();
    
    double cTable = 1;    
    double dTable = 1;

    input->CTable(cTable);
    input->DTable(dTable);    

    input->SetInputLevelDensity();

    auto * NLD = new LevelDensityHFB_BSk14(50,120,1,1);
    EXPECT_EQ(NLD->GetDtable(), dTable);
    EXPECT_EQ(NLD->GetCtable(), cTable);
}

TEST(SapphireInput, SetInputParticleTransmission){
    auto * input = new SapphireInput();

    int aF=0;
    int pF=0;
    int nF=0;

    bool pTp = true;
    
    double pNorm=2.0;
    double aNorm=2.0;
    double nNorm=2.0;

    input->a_Formalism(aF);
    input->p_Formalism(pF);
    input->n_Formalism(nF);

    input->Anorm(aNorm);
    input->Pnorm(pNorm);
    input->Nnorm(nNorm);

    input->PorterThomas_p(pTp);

    input->SetInputParticleTransmission();

    EXPECT_EQ(ParticleTransmissionFunc::GetAlphaFormalism(),aF);
    EXPECT_EQ(ParticleTransmissionFunc::GetProtonFormalism(),pF);
    EXPECT_EQ(ParticleTransmissionFunc::GetNeutronFormalism(),nF);

    EXPECT_EQ(ParticleTransmissionFunc::GetPorterThomas(),pTp);

    EXPECT_EQ(ParticleTransmissionFunc::GetAnorm(),aNorm);
    EXPECT_EQ(ParticleTransmissionFunc::GetPnorm(),pNorm);
    EXPECT_EQ(ParticleTransmissionFunc::GetNnorm(),nNorm);   
}

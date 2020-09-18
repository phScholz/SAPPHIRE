set(BINARY test_sapphire)


set(TEST_SOURCE ${TEST_SOURCE} ./main.cpp)
set(TEST_SOURCE ${TEST_SOURCE} ./test_NuclearMass.cpp)
set(TEST_SOURCE ${TEST_SOURCE} ./test_SapphireInput.cpp)

set(TEST_SOURCE ${TEST_SOURCE} ../src/LevelDensity/LevelDensityFormula.cpp ../src/LevelDensity/LevelDensityTable.cpp ../src/LevelDensity/RauscherLevelDensity.cpp ../src/LevelDensity/LevelDensityHFB_BSk14.cpp)
set(TEST_SOURCE ${TEST_SOURCE} ../src/Databases/NuclearLevels.cpp ../src/Databases/NuclearMass.cpp)
set(TEST_SOURCE ${TEST_SOURCE} ../src/Potentials/ParticleTransmissionFunc.cpp ../src/Potentials/JLMPotential.cpp ../src/Potentials/EquivSquareWell.cpp ../src/Potentials/McFaddenSatchlerPotential.cpp ../src/Potentials/Potential.cpp)
set(TEST_SOURCE ${TEST_SOURCE} ../src/GammaStrength/GammaTransmissionFunc.cpp ../src/GammaStrength/BrinkAxelGSF.cpp ../src/GammaStrength/D1MQRPA.cpp ../src/GammaStrength/McCullaghGSF.cpp ../src/GammaStrength/KopeckyUhlGSF.cpp)
set(TEST_SOURCE ${TEST_SOURCE} ../src/Modules/Module_OldSapphire.cpp ../src/Modules/Module_BreitWigner.cpp ../src/Modules/Module_Compound.cpp ../src/Modules/Module_NLD.cpp ../src/Modules/Module_GSF.cpp ../src/Modules/Module_RandomScheme.cpp ../src/Modules/Module_Decayer.cpp ../src/Modules/Module_CrossSection.cpp)
set(TEST_SOURCE ${TEST_SOURCE} ../src/Decayer/Decayer.cpp ../src/Decayer/DecayResults.cpp ../src/Decayer/DecayController.cpp ../src/Decayer/PreEqDecayer.cpp)
set(TEST_SOURCE ${TEST_SOURCE} ../src/CrossSection.cpp ../src/SapphireInput.cpp ../src/CoulFunc.cpp ../src/TransitionRateFunc.cpp ../src/Setup.cpp ../src/ParticleHoleLevelDensity.cpp ../src/PreEqTransitionRateFunc.cpp ../src/SapphireMPITypes.cpp ../src/RandomScheme.cpp)


#message(\n"TEST_SOURCE ${TEST_SOURCE}\n")

add_executable(${BINARY} ${TEST_SOURCE})
include_directories(${Sapphire_SOURCE_DIR}/include)


add_test(NAME ${BINARY} COMMAND ${BINARY})

target_link_libraries(${BINARY} gtest ${GSL_LIBRARIES} ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${EXTRA_LIBS})

install(TARGETS test_sapphire DESTINATION ${Sapphire_BINARY_DIR})
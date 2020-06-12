set(Sapphire_SOURCES ${Sapphire_SOURCES} LevelDensity/LevelDensityFormula.cpp LevelDensity/LevelDensityTable.cpp LevelDensity/RauscherLevelDensity.cpp LevelDensity/LevelDensityHFB_BSk14.cpp)
set(Sapphire_SOURCES ${Sapphire_SOURCES} Databases/NuclearLevels.cpp Databases/NuclearMass.cpp)
set(Sapphire_SOURCES ${Sapphire_SOURCES} Potentials/ParticleTransmissionFunc.cpp Potentials/JLMPotential.cpp Potentials/EquivSquareWell.cpp Potentials/McFaddenSatchlerPotential.cpp Potentials/Potential.cpp)
set(Sapphire_SOURCES ${Sapphire_SOURCES} GammaStrength/GammaTransmissionFunc.cpp GammaStrength/BrinkAxelGSF.cpp GammaStrength/D1MQRPA.cpp GammaStrength/McCullaghGSF.cpp GammaStrength/KopeckyUhlGSF.cpp)
set(Sapphire_SOURCES ${Sapphire_SOURCES} Modules/Module_NLD.cpp Modules/Module_GSF.cpp Modules/Module_RandomScheme.cpp Modules/Module_OldSapphire.cpp Modules/Module_Decayer.cpp Modules/Module_CrossSection.cpp)
set(Sapphire_SOURCES ${Sapphire_SOURCES} Decayer/Decayer.cpp Decayer/DecayResults.cpp Decayer/DecayController.cpp Decayer/PreEqDecayer.cpp)
set(Sapphire_SOURCES ${Sapphire_SOURCES} CrossSection.cpp SapphireInput.cpp CoulFunc.cpp TransitionRateFunc.cpp Sapphire.cpp Setup.cpp ParticleHoleLevelDensity.cpp PreEqTransitionRateFunc.cpp SapphireMPITypes.cpp RandomScheme.cpp)

include_directories(${Sapphire_SOURCE_DIR}/include)

add_executable (sapphire ${Sapphire_SOURCES})

target_link_libraries (sapphire ${GSL_LIBRARIES} ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${EXTRA_LIBS})

install(TARGETS sapphire DESTINATION ${Sapphire_BINARY_DIR})
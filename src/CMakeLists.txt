set(Sapphire_SOURCES Module_OldSapphire.cpp Module_Decayer.cpp SapphireInput.cpp CoulFunc.cpp GammaTransmissionFunc.cpp ParticleTransmissionFunc.cpp TransitionRateFunc.cpp DecayController.cpp LevelDensity.cpp RauscherLevelDensity.cpp Decayer.cpp NuclearMass.cpp Sapphire.cpp NuclearLevels.cpp DecayResults.cpp Setup.cpp ParticleHoleLevelDensity.cpp PreEqTransitionRateFunc.cpp PreEqDecayer.cpp SapphireMPITypes.cpp EquivSquareWell.cpp Potential.cpp BrinkAxelGSF.cpp KopeckyUhlGSF.cpp McCullaghGSF.cpp McFaddenSatchlerPotential.cpp JLMPotential.cpp)

if (NOT MPI_BUILD) 
set(Sapphire_SOURCES ${Sapphire_SOURCES} CrossSection.cpp Module_CrossSection.cpp)
endif (NOT MPI_BUILD)

include_directories(${Sapphire_SOURCE_DIR}/include)

add_executable (sapphire ${Sapphire_SOURCES})

#if (MPI_BUILD) 
target_link_libraries (sapphire ${GSL_LIBRARIES} ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${EXTRA_LIBS})
#else (MPI_BUILD)
#target_link_libraries (sapphire ${GSL_LIBRARIES} ${ROOT_LIBRARIES} ${EXTRA_LIBS})
#endif (MPI_BUILD)

install(TARGETS sapphire DESTINATION ${Sapphire_BINARY_DIR})

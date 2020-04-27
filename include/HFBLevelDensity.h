/**
 * @file HFBLevelDensity.h
 * @brief contains the HFBLevelDensity class
 * @date 2020-04-27
 */

#pragma once
#include "LevelDensityTable.h"

/**
 * @brief The HFB level density tables as class in sapphire
 */
class HFBLevelDensity : public LevelDensityTable {
    HFBLevelDensity(int Z, int A, double J, int parity) : LevelDensityTable(Z,A,J){

    }
};
/**
 * @file CompoundStates.h
 * @author Philipp Scholz
 * @brief Contains the class CompoundStates.h
 * @date 2020-06-19
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#pragma once
#include "CompoundState.h"
#include <vector>

class CompoundStates{
    public:
        CompoundStates(){};

        ~CompoundStates(){};

        std::vector<CompoundState> states;
};
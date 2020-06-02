/**
 * @file SPDistribution.h
 * @brief Home of SPDistribution
 * @date 2020-05
 */

#pragma once
#include <vector>
#include "SPPopulation.h"

/**
 * @brief Class which includes differen SPPopulation members
 */
class SPDistribution{
    public:
        /**
         * @brief Simple Constructor
         */
        SPDistribution();

        std::vector<SPPopulation> distribution; /**< Vector which stores the population data*/
};


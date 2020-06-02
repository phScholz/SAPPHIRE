/**
 * @file SPPopulation.h
 * @brief Home of SPPopulation
 * @date 2020-05
 */

#pragma once

/**
 * @brief Simple class to store a population value for a spin parity pair
 */
class SPPopulation{
    public:
        /**
         * @brief Simple Constructor
         */
        SPPopulation():spin(0.0), parity(1), pop(0.0){
        };

        /**
         * @brief Constructor which initializes the member variables
         * @param spin Double value for the spin
         * @param parity integer value for parity, 1 == +; -1 == -
         */
        SPPopulation(double spin, int parity, double pop):spin(spin),parity(parity), pop(pop){};
        
        double spin; /**< Spin*/
        int parity; /**< parity*/
        double pop; /**< population*/
};
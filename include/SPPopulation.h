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
        SPPopulation():spin(0.0), parity(1), pop(0.0), cdf(0.0){
        };

        ~SPPopulation(){};

        /**
         * @brief Constructor which initializes the member variables
         * @param spin Double value for the spin
         * @param parity integer value for parity, 1 == +; -1 == -
         */
        SPPopulation(double spin, int parity, double pop):spin(spin),parity(parity), pop(pop), cdf(0.0){};
        
        double Spin() const {return spin;};        /**< Getter for spin*/
        double Spin(double x){spin = x;};   /**< Setter for spin*/

        int Parity() const{return parity;};        /**< Getter for parity*/
        int Parity(int x){parity = x;};   /**< Setter for parity*/

        double Pop() const {return pop;};        /**< Getter for pop*/
        double Pop(double x){pop = x;};   /**< Setter for pop*/

        double Cdf() const {return cdf;};        /**< Getter for cdf*/
        double Cdf(double x){cdf = x;};   /**< Setter for cdf*/

    private:
        double spin; /**< Spin*/
        int parity; /**< parity*/
        double pop; /**< population*/
        double cdf; /**< cdf*/
};
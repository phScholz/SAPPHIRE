/**
 * @file SPDistribution.h
 * @brief Home of SPDistribution
 * @date 2020-05
 */

#pragma once
#include <vector>
#include <iostream>
#include "SPPopulation.h"

/**
 * @brief Class which includes differen SPPopulation members
 */
class SPDistribution{
    public:
        /**
         * @brief Simple Constructor
         */
        SPDistribution(): total(0.0){};

        std::vector<SPPopulation> *distribution; /**< Vector which stores the population data*/
        double Total() const {return total;}; /**<Getter for total*/
        double Total(double x){total = x;}; /**< Setter for total*/
        
        /**
         * @brief Normalize the distribution.
         */
        void Normalize(){
            if(verbose_) std::cout << std::endl << "Normalizing spin distribution ... " << std::endl;

            if(this->distribution->size()>0){
                /**
                 * 1. Initialize total
                 */
                this->Total(0);

                /**
                 * 2. Calculate the total probability
                 */
                for(std::vector<SPPopulation>::iterator it = distribution->begin(); it !=distribution->end(); ++it){
                    double pop = it->Pop();
                    this->Total(this->Total()+pop);
                }

                /**
                 * 3. Normalize the single probabilities on the total probability
                 */
                for(std::vector<SPPopulation>::iterator it = distribution->begin(); it !=distribution->end(); ++it){
                    double pop = it->Pop();
                    double total = this->Total();
                    it->Pop(pop/total);
                }

                /**
                 * 4. Calculate cdf values
                 */

                double offset = 0;
                for(std::vector<SPPopulation>::iterator it = distribution->begin(); it !=distribution->end(); ++it){
                    it->Cdf(offset+it->Pop());
                    offset = it->Cdf();
                }
            }
            else{
                std::cout << std::endl << "Cannot normalize empty object!" << std::endl;
            }
        };

        void PrintPopulation(){
            for(std::vector<SPPopulation>::iterator it = distribution->begin(); it !=distribution->end(); ++it){
                std::cout << "\t" << it->Spin() << "\t" << it->Parity() << "\t" << it->Pop() << std::endl;
            }
        }

    private:    
        double total;
        bool verbose_ = true;

};


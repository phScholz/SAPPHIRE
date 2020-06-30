/**
 * @file CompoundState.h
 * @author Philipp Scholz 
 * @brief Contains the class CompoundState.h
 * @date 2020-06-19
 * 
 * 
 */
#pragma once
#include <string>

class CompoundState{
    public:
        CompoundState(int A, int Z, double J, int parity, double energy): A_(A), Z_(Z), J_(J), P_(parity), E_(energy){};
        ~CompoundState(){};

        /**
         * @brief Function to obtain the total width for a specific channel.
         * 
         * @param channel "gamma", "neutron", "proton", "alpha"
         * @return double Total decay channel width
         */
        double TotalWidth(std::string channel) const{
            if(channel == "gamma") return gTotalWidth_;
            if(channel == "neutron") return nTotalWidth_;
            if(channel == "proton") return pTotalWidth_;
            if(channel == "alpha") return aTotalWidth_;
            if(channel == "all") return gTotalWidth_+nTotalWidth_+pTotalWidth_+aTotalWidth_;
            return -1;
        };

        /**
         * @brief Function to get the groundstate partiel width for a specfic channel.
         * 
         * @param channel "gamma", "neutron", "proton", "alpha"
         * @return double Groundstate decay width.
         */
        double GroundStateWidth(std::string channel) const{
            if(channel == "gamma")   {return   gGSwidth_;}
            if(channel == "neutron") {return nGSwidth_;}
            if(channel == "proton")  {return  pGSwidth_;}
            if(channel == "alpha")   {return   aGSwidth_;}
            return -1;
        };

        /**
         * @brief Method to set the total width of a respective channel
         * 
         * @param channel "gamma", "neutron", "proton", "alpha"
         * @param width Double for the width
         */
        void TotalWidth(std::string channel, double width){
            if(channel == "gamma")   {gTotalWidth_=width;}
            if(channel == "neutron") {nTotalWidth_=width;}
            if(channel == "proton")  {pTotalWidth_=width;}
            if(channel == "alpha")   {aTotalWidth_=width;}
                        
        };

        /**
         * @brief Method to set the groundstate width of a respective channel
         * 
         * @param channel "gamma", "neutron", "proton", "alpha"
         * @param width Double for the width
         */
        void GroundStateWidth(std::string channel, double width){
            if(channel == "gamma")   gGSwidth_=width;
            if(channel == "neutron") nGSwidth_=width;
            if(channel == "proton")  pGSwidth_=width;
            if(channel == "alpha")   aGSwidth_=width;
        };

        int    A() const {return A_;}
        int    Z() const {return Z_;}
        int    P() const {return P_;}
        double E() const {return E_;}
        double J() const {return J_;}
        double Density() const {return D_;}

        void Density(double x) {D_ = x;}



    private:
        int A_; /**< Mass number*/
        int Z_; /**< Charge number*/
        int P_; /**< Parity number*/
        
        double J_; /**< Spin number*/
        double E_; /**< Excitation energy*/
        double D_; /**< LevelDensity*/
        
        double gGSwidth_; /**< Groundstate gamma width*/
        double nGSwidth_; /**< Groundstate neutron width*/
        double pGSwidth_; /**< Groundstate proton width*/
        double aGSwidth_; /**< Groundstate alpha width*/

        double gTotalWidth_;    	/**< Total gamma width*/ 
        double nTotalWidth_;    	/**< Total neutron width*/
        double pTotalWidth_;    	/**< Total proton width*/
        double aTotalWidth_;    	/**< Total alpha width*/

};
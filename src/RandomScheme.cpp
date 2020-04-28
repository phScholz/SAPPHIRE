/**
 * @file RandomScheme.cpp
 * @author Philipp Scholz <pscholz@outlook.de>
 * @date 2020-04-25
 * @brief File for the RandomScheme class
 * 
 */

#include "RandomScheme.h"
#include "LevelDensity/LevelDensityHFB_BSk14.h"
#include "LevelDensity/RauscherLevelDensity.h"
#include "NuclearLevels.h"
#include <stdexcept>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
#include <ctime>
#include <vector>
#include <string>
#include "omp.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>

extern std::string sourceDirectory();

template <class NLD>
RandomScheme<NLD>::RandomScheme():maxJ_(5), maxE_(5.0), eStep_(0.01), previous(NULL){
    randomScheme = new std::vector<Level>();
    nldModel_=0;
    gsfModel_=0;
}

template <class NLD>
RandomScheme<NLD>::RandomScheme(int maxJ, double eStep) : maxJ_(maxJ), maxE_(5.0), eStep_(eStep), previous(NULL){
    randomScheme = new std::vector<Level>();    
}

template <class NLD>
RandomScheme<NLD>::~RandomScheme(){
    /*try{
        delete levelDensity_;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }

    try
    {
        delete randomScheme;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }*/
}

template <class NLD>
void RandomScheme<NLD>::CreateLevels(int Z, int A, double eStart, double eStop){
    if(verbose_) std::cout << "Creating Levels..." << std::endl;
    maxE_=eStop;

    //Select random number generator
    typedef boost::mt19937 base_generator_type;
    
    //Use the system time as seed
    base_generator_type generator(time(0));
    
    // Define a uniform random number distribution which produces "double"
    // values between 0 and 1 (0 inclusive, 1 exclusive).
    boost::uniform_real<> uni_dist(0,1);
    
    boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);

    if(eStart == 0.0 && randomScheme->size() == 0) randomScheme->push_back(GetGroundState(Z,A));

    if(verbose_) std::cout << "Begin Loop over energies." << std::endl;
    for(double energy=eStart+eStep_/2.; energy<=maxE_; energy+=eStep_){
        if(verbose_) std::cout << "Begin Loop over spins." << std::endl;
        for(double j=0; j<=maxJ_; j+=0.5){
            //If even A and the spin float is an integer
            if(A%2==0 && abs(j-floor(j))==0){
                //Create levelDensity object
                levelDensity_ = new NLD(Z,A,j);

                double levelsPerEnergyStep = CalcLevelDensity(energy)*eStep_;
                boost::poisson_distribution<> poisson_dist(levelsPerEnergyStep);
                boost::variate_generator<base_generator_type&, boost::poisson_distribution<>> poisson(generator, poisson_dist);
                double i=0;    
                double levelnum=poisson();       
                if(verbose_) std::cout << energy <<"\t"<<j<<"\t" << levelsPerEnergyStep << "\t" << levelnum << std::endl; 
                while(i < levelnum){                
                    int pi= (uni() >= 0.5) ? 1 : -1;
                    double randomEnergy = uni()*eStep_ + (energy-eStep_/2.);
                    randomScheme->push_back(Level(j,pi, randomEnergy));
                    i=i+1;
                }
            }
            
            //If odd A and non integer spin
            if(verbose_) std::cout << A%2 << " " << abs(j-floor(j)) << std::endl;
            if(A%2!=0 && abs(j-floor(j))!=0){
                //Create levelDensity object
                levelDensity_ = new NLD(Z,A,j);
                double levelsPerEnergyStep = CalcLevelDensity(energy)*eStep_;
                boost::poisson_distribution<> poisson_dist(levelsPerEnergyStep);
                boost::variate_generator<base_generator_type&, boost::poisson_distribution<>> poisson(generator, poisson_dist);
                double i=0;    
                double levelnum=poisson();       
                if(verbose_) std::cout << energy <<"\t"<<j<<"\t" << levelsPerEnergyStep << "\t" << levelnum << std::endl; 
                while(i < levelnum){                
                    int pi= (uni() >= 0.5) ? 1 : -1;
                    double randomEnergy = uni()*eStep_ + (energy-eStep_/2.);
                    randomScheme->push_back(Level(j,pi, randomEnergy));
                    i=i+1;
                }
            }
            
        }
    }

}

template <class NLD>
void RandomScheme<NLD>::CreateGammaTransitions(double eStart){
    if(verbose_) std::cout << "Creating Transitions..." << std::endl;
    size_t index = randomScheme->size() - 1;
    //Loop from highest level to lowest level
    if(verbose_) std::cout << "Begin Loop for gamma rays" << std::endl;
    for(std::vector<Level>::reverse_iterator rit = randomScheme->rbegin(); rit != randomScheme->rend(); ++rit, --index){
        if (rit->energy_ > eStart){
            int i = 1;
            double partialStrength = 0;
            double totalStrength = 0;
            //loop from lowest level to highes level
            for(std::vector<Level>::iterator it = randomScheme->begin(); it != randomScheme->end(); ++it) {
                if(abs(rit->J_-it->J_) <= 2 && index > i){
                    transmissionFunc_= GammaTransmissionFunc::CreateGammaTransmissionFunc(Z_, A_, 
                                rit->J_, rit->Pi_, it->J_, it->Pi_,maxL_,
		    					 levelDensity_, tWFC_,uTWFC_,uTWSFC_, previous, rit->energy_);

                    if(!transmissionFunc_->IsValid()) {
                        throw std::invalid_argument("Input for transmissionfunc is invalid.");
                        delete transmissionFunc_;
                    }

                    if(abs(rit->energy_-it->energy_)!=0){
                        partialStrength=CalcTransmissionFunc(abs(rit->energy_-it->energy_));                    
                        totalStrength+=partialStrength;
                        rit->gammas_.push_back(GammaTransition(i,abs(rit->energy_-it->energy_), partialStrength));
                    }

                    if(transmissionFunc_->IsValid()) delete transmissionFunc_;
                }

                if(index==1 && abs(rit->J_-it->J_)>2 && rit->gammas_.size()==0){
                    rit->gammas_.push_back(GammaTransition(0,abs(rit->energy_-it->energy_),1.0));
                    totalStrength=1;
                }

                i++;
            }

            if(verbose_) PrintRandomScheme(maxE_);

            if(verbose_) std::cout << "Renormalization" << std::endl;
            if(rit->gammas_.size()>0 && rit->energy_>eStart){
                totalStrength=0;
                //calculate total strength
                for(std::vector<GammaTransition>::iterator it = rit->gammas_.begin(); it !=rit->gammas_.end(); ++it){

                    totalStrength += it->probability_;
                }
                //normalize
                for(std::vector<GammaTransition>::iterator it = rit->gammas_.begin(); it !=rit->gammas_.end(); ++it){
                        it->probability_ = it->probability_/totalStrength;
                }
            }
        }
    }
}

template <class NLD>
void RandomScheme<NLD>::CreateRandomScheme(int Z, int A, double eStart, double eStop){
    
    CreateLevels(Z,A,eStart,eStop);
    CreateGammaTransitions(eStart);
}

template <class NLD>
void RandomScheme<NLD>::PrintRandomScheme(double maxE){
    std::cout.precision(3);
    if(verbose_) std::cout << "Random level scheme generated: " << std::endl;

    std::cout << std::endl;

    int index = 1;

    if(randomScheme->size()>0){
        for(std::vector<Level>::iterator it = randomScheme->begin(); it != randomScheme->end(); ++it) {
            if(it->energy_<=maxE){
                std::cout << index << "\t" << std::fixed << it->energy_ << "\t" << std::fixed << it->J_ << "\t" << it->Pi_ << std::endl;
                if(it->gammas_.size()>0){
                    for(std::vector<GammaTransition>::iterator g = it->gammas_.begin(); g !=it->gammas_.end(); ++g){
                        std::cout << "\t" << g->levelIndex_ << "\t" << std::fixed << g->energy_ << "\t" << std::fixed << g->probability_ << std::endl;
                    }
                }
                index++;
            }
        }
    }else throw std::out_of_range("RandomScheme is empty.");
    std::cout << std::endl;
}

template <class NLD>
Level RandomScheme<NLD>::GetGroundState(int Z_, int A_){
    std::string spinFile(sourceDirectory()+"/tables/spinod.dat");
    std::ifstream in(spinFile.c_str());
    int N_=A_-Z_;
    
    if(!in){/**std::cout << "!!! Cannot read spinFile !!!" << std::endl;*/ throw std::invalid_argument("Cannot open/read spinFile.") ;}
  
    //Is this only to find the ground state spins, if no groundstate is known?
    while(!in.eof()) {
        std::string line;
        std::getline(in,line);
    
        if(!in.eof()) {
            int Z,N,A,parityZ,parityN;
            int spinZ,spinN;
            std::istringstream lineStream(line);
            lineStream >> Z;
            char spinZString[8];
            lineStream.read(spinZString,7);
            spinZString[7]='\0';
      
            if(Z%2==0) {
	            spinZ=0;
	            parityZ=1;
            } else {
	            char* pch = strtok(spinZString,"/");
	            std::istringstream spinZStringStream(pch);
	            spinZStringStream>>spinZ;
	            pch = strtok(NULL,"/");
	            if(pch[1]=='+') parityZ=1;
	            else parityZ=-1;
            }

            lineStream >> N;
            char spinNString[8];
            lineStream.read(spinNString,7);
            spinNString[7]='\0';
      
            if(N%2==0) {
	            spinN=0;
	            parityN=1;
            } else {
	            char* pch = strtok(spinNString,"/");
	            std::istringstream spinNStringStream(pch);
	            spinNStringStream>>spinN;
	            pch = strtok(NULL,"/");
	            if(pch[1]=='+') parityN=1;
	            else parityN=-1;
            }

            lineStream >> A;
            double spin;
      
            if(spinZ==0||spinN==0) spin = double(spinZ+spinN)/2.; 
            else {
	            bool parallelZ = false;
	            bool parallelN = false;
	
                if(parityZ>0) {
	                if(((spinZ+1)/2)%2!=0) parallelZ = true;
	                } else if(((spinZ+1)/2)%2==0){
                        parallelZ = true;
    	            }
    
                if(parityN>0) {
	                if(((spinN+1)/2)%2!=0) parallelN = true;
	            } else  {
	                if(((spinN+1)/2)%2==0) parallelN = true;
	            }
    
                spin = (parallelZ&&parallelN) ? double(spinZ+spinN)/2. :
	            fabs(double(spinZ-spinN))/2.;
            }

            int parity = parityZ*parityN;

            if(Z==Z_ && N==N_){
                return Level(spin,parity,0.);
            }
        }
    }

}

template <class NLD>
void RandomScheme<NLD>::ExtendRandomScheme(int Z, int A, double eMax){
    std::vector<Level> knownLevels = NuclearLevels::FindLevels(Z,A);
    //std::cout << knownLevels.size() << std::endl;
    if(knownLevels.size()>0){
        for(std::vector<Level>::iterator it = knownLevels.begin(); it != knownLevels.end(); ++it){
            //std::cout << "Pushing back" << std::endl;
            randomScheme->push_back(*it);
        }
        CreateRandomScheme(Z,A,knownLevels.at(knownLevels.size()-1).energy_, eMax);
    }else{
        CreateRandomScheme(Z,A,0.0,eMax);
    }
}

template class RandomScheme<RauscherLevelDensity>;
//template class RandomScheme<LevelDensityHFB_BSk14>;

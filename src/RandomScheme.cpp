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
#include "Databases/NuclearLevels.h"
#include <stdexcept>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
#include <ctime>
#include <vector>
#include <string>
#include "omp.h"
#include "Progressbar.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>

extern std::string sourceDirectory();

template <class NLD>
RandomScheme<NLD>::RandomScheme():maxJ_(4), maxE_(5.0), eStep_(0.01), previous(NULL){
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
    std::cout << "Creating Levels..." << std::endl;
    maxE_=eStop;
    Z_=Z;
    A_=A;


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
    
    int numSteps = (int) (eStop-eStart)/eStep_;

    ProgressBar pg;
    pg.start(numSteps);
    for(int i=0; i<=numSteps; i++){
        pg.update(i);
        double energy=eStart+eStep_/2.0+i*eStep_;
        if(verbose_) std::cout << "Begin Loop over spins." << std::endl;
        for(double j=randomScheme->at(0).J_-1; j<=randomScheme->at(0).J_+1; j+=0.5){
            //If even A and the spin float is an integer
            if(A%2==0 && abs(j-floor(j))==0){
                //Create levelDensity object
                levelDensity_ = new NLD(Z,A,j);

                //Mean of poisson distribution isn't allowed to be <=0 ... thats why we are setting the level density in these cases to 1e-10
                double nld = CalcLevelDensity(energy);
                double levelsPerEnergyStep = (nld >0) ? nld*eStep_ : 1e-10;
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
                double nld = CalcLevelDensity(energy);
                double levelsPerEnergyStep = (nld >0) ? nld*eStep_ : 1e-10;
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
    //sorting level scheme by energy_
    //std::sort(randomScheme->begin(), randomScheme->end());

}

template <class NLD>
void RandomScheme<NLD>::AverageGroundStateBranching(){
    std::cout << std::endl << "\t" << "Energy" << "\t" << "spin" << "\t" << "groundstate/total" << std::endl;
    for(std::vector<Level>::iterator it=randomScheme->begin(); it != randomScheme->end(); ++it){
        double total = 0;
        double groundstate=0;
        for(std::vector<GammaTransition>::iterator git = it->gammas_.begin(); git !=it->gammas_.end(); ++git){
            total += git->probability_;
            
            if(git->levelIndex_==1){
                groundstate=git->probability_;
            }
        }

        std::cout.precision(3);
        if(abs(it->J_-randomScheme->at(0).J_) == 1){
                std::cout << std::fixed << "\t" << it->energy_ << "\t" << it->J_ << "\t" << groundstate/total << std::endl;
        }
    }
}

template <class NLD>
void RandomScheme<NLD>::CreateGammaTransitions(double eStart){
    //Select random number generator
    typedef boost::mt19937 base_generator_type;
    
    //Use the system time as seed
    base_generator_type generator(time(0));
    
    //Define a normal random number distribution which produces "double" for Porter-Thomas
    boost::normal_distribution<> normal_dist(0.0,1.0);
    
    boost::variate_generator<base_generator_type&, boost::normal_distribution<> > gauss(generator, normal_dist);

    std::cout << std::endl<< "Creating Transitions..." << std::endl;
    size_t index = randomScheme->size() - 1;
    //Loop from highest level to lowest level
    if(verbose_) std::cout << "Begin Loop for gamma rays" << std::endl;
    ProgressBar pg;
    int size=randomScheme->size();
    pg.start(size-1);
    int j=0;
    if(size>1){
        for(unsigned int k = size-1; k>=0; k--){
        //for(std::vector<Level>::reverse_iterator rit = randomScheme->rbegin(); rit != randomScheme->rend(); ++rit, --index){
            pg.update(j);
            if (randomScheme->at(k).energy_ > eStart){
                double partialStrength = 0;
                double totalStrength = 0;

                //loop from lowest level to highes level
                //for(std::vector<Level>::iterator it = randomScheme->begin(); it != randomScheme->end(); ++it) {
                for(unsigned int l=0; l<k; l++){
                    if(abs(randomScheme->at(k).J_-randomScheme->at(l).J_) <= 1 && k > l){
                        transmissionFunc_= GammaTransmissionFunc::CreateGammaTransmissionFunc(Z_, A_, 
                                    randomScheme->at(k).J_, randomScheme->at(k).Pi_, randomScheme->at(l).J_, randomScheme->at(l).Pi_,maxL_,
	    	    					tWFC_,uTWFC_,uTWSFC_, previous, randomScheme->at(k).energy_);

                        if(!transmissionFunc_->IsValid()) {
                            delete transmissionFunc_;
                            throw std::invalid_argument("Input for transmissionfunc is invalid.");

                        }

                        if(abs(randomScheme->at(k).energy_-randomScheme->at(l).energy_)>=0.4){
                            levelDensity_ = new NLD(Z_,A_,randomScheme->at(k).J_);
                            double nld = CalcLevelDensity(randomScheme->at(k).energy_);
                            double transmission=CalcTransmissionFunc(abs(randomScheme->at(k).energy_-randomScheme->at(l).energy_));
                            partialStrength=pow(gauss(),2)*transmission/nld;                    
                            totalStrength+=partialStrength;
                            randomScheme->at(k).gammas_.push_back(GammaTransition(l+1,abs(randomScheme->at(k).energy_-randomScheme->at(l).energy_), partialStrength));
                        }

                        if(transmissionFunc_->IsValid()) delete transmissionFunc_;
                    }

                    if(k==1 && randomScheme->at(k).gammas_.size()==0){
                        randomScheme->at(k).gammas_.push_back(GammaTransition(1,abs(randomScheme->at(k).energy_-randomScheme->at(l).energy_),1.0));
                        totalStrength=1;
                    }
                }

                if(verbose_) PrintRandomScheme(maxE_);

                if(verbose_) std::cout << "Renormalization" << std::endl;
                if(randomScheme->at(k).gammas_.size()>0 && randomScheme->at(k).energy_>eStart){
                    totalStrength=0;
                    //calculate total strength
                    for(std::vector<GammaTransition>::iterator it = randomScheme->at(k).gammas_.begin(); it !=randomScheme->at(k).gammas_.end(); ++it){

                        totalStrength += it->probability_;
                    }

                    //normalize
                    for(std::vector<GammaTransition>::iterator it = randomScheme->at(k).gammas_.begin(); it !=randomScheme->at(k).gammas_.end(); ++it){
                            it->probability_ = it->probability_/totalStrength;
                    }
                }
            }
            j++;
        }
        pg.update(size-1);
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

template <class NLD>
void RandomScheme<NLD>::WriteRandomScheme(std::string file){
    std::ofstream out(file, std::ios::out);
    out.precision(3);

    if(verbose_) std::cout << "Random level scheme generated: " << std::endl;
    
    ProgressBar pg;
    pg.start(randomScheme->size());
    if(randomScheme->size()>0){
        int index = 1;
        for(std::vector<Level>::iterator it = randomScheme->begin(); it != randomScheme->end(); ++it) {
            out << index << "\t" << std::fixed << it->energy_ << "\t" << std::fixed << it->J_ << "\t" << it->Pi_ << std::endl;
            if(it->gammas_.size()>0){
                for(std::vector<GammaTransition>::iterator g = it->gammas_.begin(); g !=it->gammas_.end(); ++g){
                    if(g->probability_ >= 0.001)
                    out << "\t" << g->levelIndex_ << "\t" << std::fixed << g->energy_ << "\t" << std::fixed << g->probability_ << std::endl;
                }
            }
        index++;
        pg.update(index);
        }
    }
    else throw std::out_of_range("RandomScheme is empty.");
}

template class RandomScheme<RauscherLevelDensity>;
template class RandomScheme<LevelDensityHFB_BSk14>;

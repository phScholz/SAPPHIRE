/**
 * @file Sapphire.cpp
 * @brief Entry point for the Module_OldSapphire 
 * @author Mary Beard, Philipp Scholz <pscholz@outlook.de>
 * @date 2020
 * 
 * The content of this file is basically the original version of Sapphire.cpp.
 * 
 */

#include "Module_OldSapphire.h"
#include "Progressbar.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include "DecayController.h"
#include "NuclearMass.h"
#include "DecayResults.h"

#include "CrossSection.h"
#include "omp.h" /** Currently only used for the Decayer*/

#include "SapphireMPITypes.h"
//#include <boost/mpi.hpp>

#include "TransitionRateFunc.h"
#include "ParticleTransmissionFunc.h"
#include "GammaTransmissionFunc.h"
/* #include "Setup.cpp"

extern void Initialize(); */

unsigned int randomSeed[12];


namespace Module_OldSapphire{
    
    typedef struct EntrancePairs {
        EntrancePairs(int Z,int A,int pType) {
            Z_=Z;
            A_=A;
            pType_=pType;
        };
        int Z_;
        int A_;
        int pType_;
    } EntrancePairs;


    void Go(int argc,char *argv[]){
       //std::cout << "Not yet working."  << std::endl;
        oldSapphire(argc, argv);
    }
    
    

    
    /**
     * @brief CMD line parameters are parsed for the Decay module
     * @param args cmd line string
     * @param Z nuclear charge number
     * @param A nuclear mass number
     * @param J Reference to a spin double
     * @param Pi Reference to a parity int
     * @param lowEnergy Reference to a lowEnergy double
     * @param highEnergy  Reference to a highEnergy double
     * @param events  Reference to the events int.
     * 
     * @returns True or False
     * 
     * @note This has to move to another class in future. This does not belong in a main file but in the class file of the respective module.
     * 
     */
    bool parseCommandLineForDecay(std::vector<std::string>& args, 
    			      int& Z, int& A, double& J, int& Pi, 
    			      double& lowEnergy, double& highEnergy,
    			      int& events) {

      if(args.size()<4) return false;
    
      bool goodA=false;
      bool goodZ=false;
      bool goodPi=false;
      bool goodJ=false;
      bool goodEnergy=false;
    
      std::string isotopeString(args[1]);
      std::string massNumberString;
      for(int i = 0; i<isotopeString.length(); i++) {
      	std::string nextChar(isotopeString,i,1);
        std::istringstream stm(nextChar);
        int nextDigit;
        if(!(stm>>nextDigit)) break;
        else massNumberString+=nextChar;
      }
      if(massNumberString.length()>0) {
        A = atoi(massNumberString.c_str());
        goodA = true;
      }
      if(NuclearMass::FindZ(isotopeString.substr(massNumberString.length()))!=-1) {
        Z = NuclearMass::FindZ(isotopeString.substr(massNumberString.length()));
        goodZ = true;
      }
    
      std::string jPiString(args[2]);
      std::string firstString,secondString,parityString,delimiterString;
      bool foundDelimiter = false;
      for(int i = 0;i<jPiString.length();i++) {
      	std::string nextChar(jPiString,i,1);
      	if(nextChar=='/'||nextChar=='.') {
      	  foundDelimiter=true;
      	  delimiterString=nextChar;
      	  continue;
      	} else if(nextChar=='+'||nextChar=='-') {
      	  parityString = nextChar;
          break;
        }
        std::istringstream stm(nextChar);
        int digit;
        if(stm>>digit) {
          if(!foundDelimiter) firstString+=nextChar;
          else secondString+=nextChar;
        }
      }
      if(parityString.length()>0) {
        if(parityString=="-") Pi=-1;
        else if(parityString=="+") Pi=1;
        goodPi=true;
      }
      if(firstString.length()>0) {
      	if(foundDelimiter&&delimiterString=="/") {
      	  if(secondString.length()>0) J = atof(firstString.c_str())/atof(secondString.c_str());
      	  else J = atof(firstString.c_str());
      	} else if(foundDelimiter&&delimiterString==".") {
      	  firstString += '.'+secondString;
      	  J = atof(firstString.c_str());
      	} else J = atof(firstString.c_str());
      	double intPart;
        if(modf(J*2.,&intPart)==0.) goodJ=true; 
      }
    
      std::string energyString(args[3]); 
      if(energyString.length()>0) {
        std::string lowEnergyString;
        std::string highEnergyString;
        int dashPos = energyString.find_first_of('-');
        if(dashPos!=std::string::npos) {
          lowEnergyString = energyString.substr(0,dashPos+1);
          highEnergyString = energyString.substr(dashPos+1);
        } else {
          lowEnergyString = highEnergyString = energyString;
        }
        std::istringstream lowEnergyStream(lowEnergyString);
        std::istringstream highEnergyStream(highEnergyString);
        if((lowEnergyStream>>lowEnergy)&&
           (highEnergyStream>>highEnergy)) goodEnergy=true;
      }
    
      if(args.size()>4) {
        std::string eventsString(args[4]);
        if(eventsString.length()>0) {
          std::istringstream stm(eventsString);
          if(!(stm>>events)) events = 1;
        }
      } else events = 1;

      return (goodA&&goodZ&&goodPi&&goodJ&&goodEnergy);       
    }

    
    /**
    * @brief CMD line parameters are parsed for the Cross section module
    * @param args cmd line string
    * @param Z nuclear charge number
    * @param A nuclear mass number
    * @param pType Projectile
    * @param energyFile Reference to the string for the path to the EnergyFile
    * @param asciiIn Boolean which shows if there is a asciiFile or not
    * @param highEnergy  Reference to a highEnergy double
    * @param events  Reference to the events int.
    * 
    * @returns True or False
    * 
    * @note This has to move to another class in future. This does not belong in a main file but in the class file of the respective module.
    * 
    */
    bool parseCommandLineForXS(std::vector<std::string>& args,int& Z, int&A, 
    			   int& pType, std::string& energyFile, bool asciiIn) {

      if(asciiIn) {
        if(args.size()==2) energyFile = std::string(args[1]);
        return true;
      }
    
      if(args.size()<2) return false;

      bool goodA=false;
      bool goodZ=false;
      bool goodPType=false;

      std::string reactionString(args[1]);
      std::string massNumberString;
      for(int i = 0; i<reactionString.length(); i++) {
        std::string nextChar(reactionString,i,1);
        std::istringstream stm(nextChar);
        int nextDigit;
        if(!(stm>>nextDigit)) break;
        else massNumberString+=nextChar;
      }
      if(massNumberString.length()>0) {
        A = atoi(massNumberString.c_str());
        goodA = true;
      }
      reactionString.erase(0,massNumberString.length());
      std::string atomicNumberString;
      for(int i = 0; i<reactionString.length(); i++) {
        std::string nextChar(reactionString,i,1);
        if(nextChar=="+") break;
        else atomicNumberString+=nextChar;
      }
      if(NuclearMass::FindZ(atomicNumberString) != -1) {
        Z = NuclearMass::FindZ(atomicNumberString);
        goodZ = true;
      }
      reactionString.erase(0,atomicNumberString.length());
      std::string projectileString;
      for(int i = 0; i<reactionString.length(); i++) {
        std::string nextChar(reactionString,i,1);
        if(nextChar=="+") continue;
        else projectileString+=nextChar;
      }
      if(projectileString=="g") {
        pType = 0;
        goodPType=true;
      } else if(projectileString=="n") {
        pType = 1;
        goodPType=true;
      } else if(projectileString=="p") {
        pType = 2;
        goodPType=true;
      } else if(projectileString=="a") {
        pType = 3;
        goodPType=true;
      }

      if(args.size()==3) energyFile = std::string(args[2]);

      return (goodA&&goodZ&&goodPType);
    }

    

    
    void parseCommandLineForOptions(std::vector<std::string>& args, int& suffixNo, bool &preEq,
				int& numPiParticles, int& numPiHoles, int& numNuParticles, int& numNuHoles,
                                bool& calcAverageWidth, bool& calcRates, bool& asciiIn,
				std::string& inFile, int& entranceState, std::vector<int>& exitStates,
				bool& printTrans) {
    
    /* void parseCommandLineForOptions(std::vector<std::string>& args, int& suffixNo, bool &preEq,
				int& numPiParticles, int& numPiHoles, int& numNuParticles, int& numNuHoles) {
     */
        std::vector<std::string>::iterator it = args.begin();
        while(it!=args.end()) {
            if(it->compare(0,8,"--l-max=")==0) {
                std::istringstream stm(it->substr(8));
            
            double maxL;
      
            if(stm>>maxL) {
	            Decayer::SetMaxL(maxL);
	            std::cout << "Maximum l-value set to " << int(Decayer::GetMaxL()) 
		        << "." << std::endl;
            }
      
            it = args.erase(it);
            

            } else if(it->compare(0,15,"--gamma-cutoff=")==0) {
                std::istringstream stm(it->substr(15));
                double cutoffEnergy;
                if(stm>>cutoffEnergy) {
	                CrossSection::SetCalculateGammaCutoff(false);
	                TransitionRateFunc::SetGammaCutoffEnergy(cutoffEnergy);
	                std::cout << "Gamma channel integration cutoff energy set to " 
		                << TransitionRateFunc::GetGammaCutoffEnergy() << "." 
		                << std::endl;
                }
      
                it = args.erase(it);
            } else if(it->compare(0,17,"--residual-gamma=")==0) {
                std::string residual = it->substr(17);
      
                if(residual=="Y") {
                    CrossSection::SetResidualGamma(true);
                    std::cout << "Setting calculation of residual gamma cross section." << std::endl;
                } else if(residual=="N") {
                    CrossSection::SetResidualGamma(false);
                    std::cout << "Unsetting calculation of residual gamma cross section." << std::endl;
                }
      
                it = args.erase(it);
            
            } else if(it->compare(0,19,"--residual-neutron=")==0) {
                std::string residual = it->substr(19);
                if(residual=="Y") {
                    CrossSection::SetResidualNeutron(true);
                    std::cout << "Setting calculation of residual neutron cross section." << std::endl;
                } else if(residual=="N") {
                    CrossSection::SetResidualNeutron(false);
                    std::cout << "Unsetting calculation of residual neutron cross section." << std::endl;
                }
                
                it = args.erase(it);
            
            } else if(it->compare(0,18,"--residual-proton=")==0) {
                std::string residual = it->substr(18);
      
                if(residual=="Y") {
                    CrossSection::SetResidualProton(true);
                    std::cout << "Setting calculation of residual proton cross section." << std::endl;
                } else if(residual=="N") {
                    CrossSection::SetResidualProton(false);
                    std::cout << "Unsetting calculation of residual proton cross section." << std::endl;
                }
      
                it = args.erase(it);
            
            } else if(it->compare(0,17,"--residual-alpha=")==0) {
                std::string residual = it->substr(17);
      
                if(residual=="Y") {
                    CrossSection::SetResidualAlpha(true);
                    std::cout << "Setting calculation of residual alpha cross section." << std::endl;
                } else if(residual=="N") {
                    CrossSection::SetResidualAlpha(false);
                    std::cout << "Unsetting calculation of residual alpha cross section." << std::endl;
                }
      
                it = args.erase(it);
    
            } else if(it->compare(0,19,"--average-rad-width")==0) {
                calcAverageWidth = true;
                it = args.erase(it);
    
            } else if(it->compare(0,7,"--rates")==0) {
                calcRates = true;
                it = args.erase(it);
    
            } else if(it->compare(0,5,"--in=")==0)  {
                asciiIn = true;
                inFile = it->substr(5);
                it = args.erase(it);
    
            } else if(it->compare(0,11,"--entrance=")==0)  {
                std::istringstream stm(it->substr(11));
       
                if(!(stm>>entranceState)) entranceState = 0;
                it = args.erase(it);
    
            } else if(it->compare(0,13,"--gamma-exit=")==0)  {
                std::istringstream stm(it->substr(13));
                int exit;
      
                if(stm>>exit) {
	                exitStates[0]=exit;
                }
      
                it = args.erase(it);
            
            } else if(it->compare(0,15,"--neutron-exit=")==0)  {
                std::istringstream stm(it->substr(15));
                int exit;
                
                if(stm>>exit) {
	                exitStates[1]=exit;
                }
                
                it = args.erase(it);
    
            } else if(it->compare(0,14,"--proton-exit=")==0)  {
                std::istringstream stm(it->substr(14));
                int exit;
      
                if(stm>>exit) {
	                exitStates[2]=exit;
                }
      
                it = args.erase(it);
    
            } else if(it->compare(0,13,"--alpha-exit=")==0)  {
                std::istringstream stm(it->substr(13));
                int exit;
      
                if(stm>>exit) {
	                exitStates[3]=exit;
                }

                it = args.erase(it);
    
            } else if(it->compare(0,14,"--transmission")==0)  {
                printTrans=true;
                it = args.erase(it);

            } else if(it->compare(0,9,"--suffix=")==0) {
                std::istringstream stm(it->substr(9));
                if(!(stm>>suffixNo)) suffixNo = 0;
                it = args.erase(it);
    
            } else if(it->compare(0,9,"--pre-eq=")==0) {
                preEq=true;
                std::string particleHoleString = it->substr(9);
                std::string temp;
                int argNum = 0;
                
                for(unsigned char i = 0;i<=particleHoleString.length();++i) {
                    if(particleHoleString[i]==','||i==particleHoleString.length()) {
                        std::istringstream stm(temp);
                        if(argNum == 0) {
                            stm>>numPiParticles;
                        } else if(argNum == 1) {
                            stm>>numPiHoles; 
                        } else if(argNum == 2) {
                            stm>>numNuParticles;
                        } else if(argNum == 3) {
                            stm>>numNuHoles;
                        }
                    
                        temp.clear();
                        argNum++; 
                    } else {
	                    temp.push_back(particleHoleString[i]);
                    }
                }

                if(numPiParticles<0||numPiHoles<0||numNuParticles<0||numNuHoles<0) 
	                std::cout << "Error parsing pre-eq options." << std::endl;
                it = args.erase(it);
            
            } else if(it->compare(0,11,"--opt-alpha")==0) {
                std::cout << "Using optical model for alpha transmission term." << std::endl;
                ParticleTransmissionFunc::SetAlphaFormalism(1);
                it = args.erase(it);
            
            }else if(it->compare(0,13,"--opt-neutron")==0) {
                std::cout << "Using optical model for neutron transmission term." << std::endl;
                ParticleTransmissionFunc::SetNeutronFormalism(1);
                it = args.erase(it);
            
            } else if(it->compare(0,12,"--opt-proton")==0) {
                std::cout << "Using optical model for proton transmission term." << std::endl;
                ParticleTransmissionFunc::SetProtonFormalism(1);
                it = args.erase(it);
            } else if(it->compare(0,7,"--EGDR=")==0) {
                std::string type = it->substr(7);
                std::transform(type.begin(), type.end(),
		        type.begin(), ::tolower);
      
                if(type=="edslo") {
	                std::cout << "Using ED-SLO EGDR shape." << std::endl;
	                GammaTransmissionFunc::SetEGDRType(2);
                } else if(type=="slo") {
	                std::cout << "Using SLO EGDR shape." << std::endl;
	                GammaTransmissionFunc::SetEGDRType(0);
                } else if(type=="glo") {
	                std::cout << "Using GLO EGDR shape." << std::endl;
	                GammaTransmissionFunc::SetEGDRType(1);
                }
                it = args.erase(it);
            } else if(it->compare(0,15,"--porter-thomas")==0) {
                std::cout << "Enabling Porter-Thomas distrbutions." << std::endl;
                GammaTransmissionFunc::SetPorterThomas(true); 
                ParticleTransmissionFunc::SetPorterThomas(true); 
                it = args.erase(it);
            } else ++it;
        }
    }

/* #ifdef MPI_BUILD
int oldSapphire_MPI(int argc, char *argv[]) {
  for(int i=1;i<argc;i++) 
    if(strcmp(argv[i],"--help")==0) {
      printHelp();
      return 0;
    }

  boost::mpi::environment env(argc,argv);
  boost::mpi::communicator world;
  
  if(world.size()<2) {
    std::cout << "Must run with at least 2 cores." << std::endl;
    env.abort(-1);
  }
  
  if(world.rank()==0) std::cout << "Initializing Databases..." << std::endl;
  Initialize();
  
  int events = 0, suffixNo=0;
  InitialNucleusData initialNucleus;

  std::vector<std::string> args;	
  for(int 1 = 0;i<argc;i++) args.push_back(std::string(argv[i]));

  bool preEq = false;
  int numPiParticles=-1,numPiHoles=-1,numNuParticles=-1,numNuHoles=-1;
  parseCommandLineForOptions(args,suffixNo,preEq,numPiParticles,numPiHoles,numNuParticles,numNuHoles);

  if(world.rank()==0) {
    if(argc<2) {
      std::cout << "Too few arguments given."  << std::endl;
      env.abort(-1);
    }
           
    int A,Z,Pi;
    double J,lowEnergy,highEnergy;
    if(!parseCommandLineForDecay(args,Z,A,J,Pi,lowEnergy,highEnergy,events)) {
      std::cout << "Isotope, J-Pi, and Energy must be specified for decay (i.e. 60Fe 1- 15)." << std::endl;
      env.abort(-1);
    }
    std::cout << "Input Values For Parent Nucleus:" << std::endl
	      << std::setw(15) << "Z:"             << std::setw(12) << Z          
	      << std::setw(0) << std::endl
	      << std::setw(15) << "A:"             << std::setw(12) << A          
	      << std::setw(0) << std::endl
	      << std::setw(15) << "J:"             << std::setw(12) << J          
	      << std::setw(0) << std::endl
	      << std::setw(15) << "Pi:"            << std::setw(12) << Pi         
	      << std::setw(0) << std::endl
	      << std::setw(15) << "energy (low):"  << std::setw(12) << lowEnergy  
	      << std::setw(0) << std::endl
	      << std::setw(15) << "energy (high):" << std::setw(12) << highEnergy 
	      << std::setw(0) << std::endl
	      << std::setw(15) << "events:"        << std::setw(12) << events     
	      << std::setw(0) << std::endl;
    initialNucleus = InitialNucleusData(Z,A,J,Pi,lowEnergy,highEnergy,preEq);
  }
  boost::mpi::broadcast(world,initialNucleus,0);
  
  if(world.rank()==0) {
    std::cout << "Starting Decay Simulation With " << world.size()-1 
	      << " Processes..." << std::endl;
    masterProcess(world,initialNucleus,suffixNo,events);
  } else {
    slaveProcess(world,initialNucleus);
  }

  return 0;
}
#endif */

    
/* #ifdef MPI_BUILD
void masterProcess(boost::mpi::communicator& world,InitialNucleusData initalNucleus,
                   int suffixNo,int events) {
  
  int maxChunkSize = 5000;
  int chunkSize = int(double(events)/double(world.size()-1));
  if(chunkSize>maxChunkSize) chunkSize = maxChunkSize;
  std::cout << "Using Event Chunk Size: " << chunkSize << std::endl;
  
  DecayResults* results = (events==1) ? NULL:
    new DecayResults(initalNucleus.Z(),
		     initalNucleus.A(),
		     initalNucleus.J(),
		     initalNucleus.Pi(),
		     initalNucleus.lowEnergy(),
		     initalNucleus.highEnergy(),
		     suffixNo);
  
  int eventsDispatched=0;
  int runningProcesses=0;
  for(int i = 1;i<world.size();i++) {
    if(eventsDispatched+chunkSize<=events) {
      world.send(i,SapphireTagProcess,chunkSize);
      ++runningProcesses;
      eventsDispatched+=chunkSize;
    } else if(eventsDispatched<events) {
      int eventsRemaining = events-eventsDispatched;
      world.send(i,SapphireTagProcess,eventsRemaining);
      ++runningProcesses;
      eventsDispatched=events;
    } else {
      world.send(i,SapphireTagDone,0);
    }    
  }
  
  bool done = (eventsDispatched >= events) ;
  int completedEvents = 0;
  while(!done||runningProcesses) {
    boost::optional<boost::mpi::status> stat = world.iprobe(boost::mpi::any_source,boost::mpi::any_tag);
    if(stat) {
      std::vector<std::pair<DecayData,std::vector<DecayProduct> > > r;
      world.recv(stat->source(),boost::mpi::any_tag,r);
      completedEvents+=r.size();
      if(completedEvents>0&&completedEvents%1000==0) std::cout << "Decayed " << completedEvents
							       << " of " << events << " nuclei..."
							       << std::endl;
      if(events>1) {
	results->AddResults(r);
      }
      if(!done) {
	if(eventsDispatched+chunkSize<=events) {
	  world.send(stat->source(),SapphireTagProcess,chunkSize);
	  eventsDispatched+=chunkSize;
	} else if(eventsDispatched<events) {
	  int eventsRemaining = events-eventsDispatched;
	  world.send(stat->source(),SapphireTagProcess,eventsRemaining);
	  eventsDispatched=events;
	}
	done = (eventsDispatched >= events);
      } else {
	world.send(stat->source(),SapphireTagDone,0);
	--runningProcesses;
      }
    }
  }
  if(results) delete results;
}
#endif */

/* #ifdef MPI_BUILD
void slaveProcess(boost::mpi::communicator& world,InitialNucleusData initialNucleus) {
  srand ( time(NULL) + world.rank());
  
  bool done=false;
  while (!done) {
    int numToDecay;
    boost::mpi::status stat = world.recv(0,boost::mpi::any_tag,numToDecay);
    if(stat.tag()==SapphireTagProcess) {
      std::vector<std::pair<DecayData,std::vector<DecayProduct> > > chunkResults;
      for(int j = 0;j<numToDecay;j++) {
	double energy = (initialNucleus.lowEnergy()==initialNucleus.highEnergy()) ? initialNucleus.lowEnergy() :
	  initialNucleus.lowEnergy()+(initialNucleus.highEnergy()-initialNucleus.lowEnergy())*
	  double(rand())/double(RAND_MAX);
	DecayController* controller;
	if(initialNucleus.preEq()) {
	    controller= new DecayController(initialNucleus.Z(),
					    initialNucleus.A(),
					    initialNucleus.J(),
					    initialNucleus.Pi(),
					    energy,numNuParticles,numNuHoles,numPiParticles,numPiHoles);
	} else controller = new DecayController(initialNucleus.Z(),
						initialNucleus.A(),
						initialNucleus.J(),
						initialNucleus.Pi(),
						energy);
	double neutronEntranceWidth = 0.;
	double protonEntranceWidth = 0.;
	double gammaEntranceWidth = 0.;
	double alphaEntranceWidth = 0.;
	double neutronTotalWidth = 0.;
	double protonTotalWidth = 0.;
	double gammaTotalWidth = 0.;
	double alphaTotalWidth = 0.;
       	controller->Decay(neutronEntranceWidth,protonEntranceWidth,alphaEntranceWidth,gammaEntranceWidth,
			  neutronTotalWidth,protonTotalWidth,alphaTotalWidth,gammaTotalWidth);
	//chunkResults.push_back(std::pair<DecayData,std::vector<DecayProduct> >(DecayData(),std::vector<DecayProduct>()));
	chunkResults.push_back(std::pair<DecayData,std::vector<DecayProduct> >(DecayData(energy,neutronEntranceWidth,protonEntranceWidth,
											 alphaEntranceWidth,gammaEntranceWidth,
											 neutronTotalWidth,protonTotalWidth,
											 alphaTotalWidth,gammaTotalWidth),controller->DecayProducts()));
	delete controller;
      }
      world.send(0,SapphireTagResults,chunkResults);
    } else {
      done = true;
    }
  }
}
#endif */

void printHelp() {
    std::cout << "Module: old" << std::endl;
    std::cout << std::endl;
    std::cout  << "\tSyntax (Cross-Section): sapphire old <options> entrance-pair [energies-file]" << std::endl;
    std::cout  << "\tSyntax  (Monte-Carlo) : sapphire old <options> nucleus spin-parity energy-range events" << std::endl << std::endl
   << "Options:" << std::endl
   << std::setw(25) << std::left << "\t--opt-alpha" << std::setw(0) << "Sets the alpha transmission coefficients to be" << std::endl
   << std::setw(25)  << '\t' << std::setw(0) << "calculated with the optical model from" << std::endl
   << std::setw(25)  << '\t' << std::setw(0) << "McFadden-Satchler parameters.  Default calculates" << std::endl
   << std::setw(25)  << '\t' << std::setw(0) << "using equivalent square well parametrization." << std::endl
   << std::setw(25) << std::left << "\t--l-max=VALUE:" << std::setw(0) << "Sets the maximum value of orbital angular momentum to" << std::endl
   <<std::setw(25)  << '\t' << std::setw(0) << "VALUE." << std::endl
   << std::setw(25) << std::left << "\t--suffix=VALUE:" << std::setw(0) << "Monte-Carlo mode only.  Sets the suffix for the output" <<std::endl
   << std::setw(25)  << '\t' << std::setw(0) << "file name (i.e. Sapphire_*_VALUE.root)." << std::endl
   << std::setw(25) << std::left << "\t--pre-eq=pp,ph,np,nh" << std::setw(0) << "Monte-Carlo mode only.  Sets the initial exciton" << std::endl
   << std::setw(25)  << '\t' << std::setw(0) << "configuration and simulates the equilibrium" << std::endl
   << std::setw(25)  << '\t' << std::setw(0) << "process allowing for pre-equilibrium particle decay." << std::endl
   << std::setw(25) << std::left << "\t--gamma-cutoff=VALUE" << std::setw(0) << "Cross section mode only.  Sets the maximum excitation" << std::endl
   << std::setw(25)  << '\t' << std::setw(0) << "energy in the compound nucleus where gamma emission is" << std::endl
   << std::setw(25)  << '\t' << std::setw(0) << "exclusive.  Transitions to final states above this energy" << std::endl
   << std::setw(25)  << '\t' << std::setw(0) << "will not contribute to residual capture cross sections." << std::endl
   << std::setw(25) << std::left << "\t--residual-gamma=Y/N" << std::setw(0) << "Cross section mode only.  Toggles if residual or total" << std::endl
   << std::setw(25)  << '\t' << std::setw(0) << "capture cross sections are calculated." << std::endl
   << std::setw(25) << std::left << "\t--residual-neutron=Y/N" << std::setw(0) << "Cross section mode only.  Toggles if residual or total" << std::endl
   << std::setw(25)  << '\t' << std::setw(0) << "neutron cross sections are calculated." << std::endl
   << std::setw(25) << std::left << "\t--residual-proton=Y/N" << std::setw(0) << "Cross section mode only.  Toggles if residual or total" << std::endl
   << std::setw(25)  << '\t' << std::setw(0) << "proton cross sections are calculated." << std::endl
   << std::setw(25) << std::left << "\t--residual-alpha=Y/N" << std::setw(0) << "Cross section mode only.  Toggles if residual or total" << std::endl
   << std::setw(25)  << '\t' << std::setw(0) << "alpha cross sections are calculated." << std::endl
   << std::setw(25) << std::left << "\t--average-rad-width" << std::setw(0) << "Cross section mode only.  Calculates the average s-wave" << std::endl
   << std::setw(25)  << '\t' << std::setw(0) << "radiative width at threshold and exits." << std::endl
   << std::setw(25) << std::left << "\t--rates" << std::setw(0) << "Cross section mode only.  Calculates the astrophysical" << std::endl 
   << std::setw(25)  << '\t' << std::setw(0) << "reaction rates in addition to cross sections.  If" << std::endl
   << std::setw(25)  << '\t' << std::setw(0) << "projectile is a neutron, Maxwellian averaged capture" << std::endl
   << std::setw(25)  << '\t' << std::setw(0) << "cross sections are also calculated." << std::endl 
   << std::setw(25) << std::left << "\t--EGDR=type" << std::setw(0) << "Sets the EGDR shape.  Variable 'type' can either be" << std::endl 
   << std::setw(25)  << '\t' << std::setw(0) << "edslo (SMOKER), slo (BRINK-AXEL), or glo (KOPECKY-UHL)." << std::endl;
}

    int oldSapphire(int argc, char *argv[]){
        

        for(int i=1;i<argc;i++){ 
            if(strcmp(argv[i],"--help")==0) {
                printHelp();
                exit(0);
            }
        }

        if(argc<3) {
            std::cout << "Too few arguments given."  << std::endl;
            printHelp();
            return -1;
        }

        /** Calling the Initialize() function from Setup.cpp*/
        /* Initialize(); */
  
        /** Create a vector which containts the cmd parameter. 
         *  Now int i starts at 1, the first parameter is "old" now.
        */
        std::vector<std::string> args; 
        for(int i = 1;i<argc;i++) args.push_back(std::string(argv[i])); 

        /** Define sum variables which will be needed for the cmd line parsing*/
        int suffixNo=0;
        int entranceState=0;
        bool preEq = false;
        bool calcAverageWidth = false;
        bool calcRates = false;
        bool asciiIn = false;
        bool printTrans = false;
        int numPiParticles=-1,numPiHoles=-1,numNuParticles=-1,numNuHoles=-1;
        std::string inFile;
        std::vector<int> exitStates(4,-1);
        
        /** Parse optional cmd line arguments*/
        parseCommandLineForOptions( args,
                                    suffixNo,
                                    preEq,
                                    numPiParticles,
                                    numPiHoles,
                                    numNuParticles,
                                    numNuHoles,
                                    calcAverageWidth,
                                    calcRates,
                                    asciiIn,
                                    inFile,
                                    entranceState,
                                    exitStates,
                                    printTrans);

        /** Random number generator seeding?*/
        time_t seed = time(NULL);
        seed+=473879*suffixNo;
        srand ( seed );

        /** Fill the randomSeed array with 12 random numbers*/
        for(unsigned char i = 0;i<12;i++) randomSeed[i]=rand();

        /*
        * This here was in the old Sapphire code the test which mode should be run.
        * Calculations of cross sections or decay?
        * In the new version of this code, this will be handled by the choice of the module.
        * So, now ... this should start a cross section calculation.
        * CROSS_SECTION CALCULATION STARTS HERE
        */
        
        if(args.size()<4) {
            int Z,A,pType;
            std::string energyFile;
            std::vector<EntrancePairs> entrancePairs;
            if(!parseCommandLineForXS(args,Z,A,pType,energyFile,asciiIn)) {
                std::cout << "Initial pair must be given for cross section (i.e. 60Fe+n)." << std::endl;
                return -1;
            }
        /**
         * Apparently here is an undocumented option to make calculations for more than one reaction. 
         * This is probably enabled by the "--in=" keyword. See parseCommandLineForOptions().
         */
            if(!asciiIn) {
                /*
                * If no ascii file was given, cross section is only calculated for one reaction.
                * Target nucleus and projectile are written into std::cout
                */
                std::cout << "Input Values For Cross Section:"   << std::endl
		            << std::setw(14) << "Z:"               << std::setw(12) 
		            << Z      << std::setw(0) << std::endl
		            << std::setw(14) << "A:"               << std::setw(12) 
		            << A      << std::setw(0) << std::endl;
                
                
                if(pType==0) 
	                std::cout << std::setw(14) << "projectile:"  << std::setw(12) << "g" << std::setw(0) << std::endl;
                else if(pType==1) 
	                std::cout << std::setw(14) << "projectile:"  << std::setw(12) << "n" << std::setw(0) << std::endl;
                else if(pType==2) 
	                std::cout << std::setw(14) << "projectile:"  << std::setw(12) << "p" << std::setw(0) << std::endl;
                else if(pType==3) 
	                std::cout << std::setw(14) << "projectile:"  << std::setw(12) << "a" << std::setw(0) << std::endl;
                
                std::cout << "Starting Cross Section Calculation..." << std::endl;
            } else {
                /** 
                * Apparently, a file was given with for than one reaction.
                * So, now it will be checked if the file can be read. 
                * If not, the program exits.
                * If, then EntrancePairs will be generated and pushed back in entrancePairs.
                */
                std::ifstream in(inFile.c_str());
                if(!in) {
	                std::cout << "Could not open " << inFile << " for reading." << std::endl;
	                return -1;
                } else {
	                std::cout << "Reading nuclei from " << inFile << "." << std::endl;

	                while(!in.eof()) {
	                    std::string line;
	                    std::getline(in,line);
	                    if(!in.eof()) {
	                        std::istringstream stm(line);
	                        if(stm >> Z >> A >> pType)
	                        entrancePairs.push_back(EntrancePairs(Z,A,pType));
	                    }
	                }
	                
                    in.close();
                }
            }
    
            /**
             * The Decayer is informed, that a CrossSection calculation will happen.
             */
            Decayer::SetCrossSection(true);  
            /** 
             * If no ascii-file was given, then the length of the following loop is 1.
             */
            int loopMax = (asciiIn) ? entrancePairs.size() : 1;
    
            for(int i = 0;i<loopMax;++i) {
                if(asciiIn) std::cout << "Calculating for Z: "
			                    << entrancePairs[i].Z_ << " A: " << entrancePairs[i].A_ 
			                    << " and Particle Type: " << entrancePairs[i].pType_ << std::endl;
                /**
                 * If an ascii-file was given
                 */
                CrossSection* crossSection = (asciiIn) ? new CrossSection(entrancePairs[i].Z_,entrancePairs[i].A_, entrancePairs[i].pType_,energyFile,calcRates) : new CrossSection(Z,A,pType,energyFile,calcRates,entranceState,exitStates);
                if(crossSection->IsValid()) {
	                if(calcAverageWidth) {
	                    std::pair<double,double> sWave = crossSection->CalcAverageSWaveResWidth();
	                    std::pair<double,double> pWave = crossSection->CalcAveragePWaveResWidth();
	                    std::pair<double,double> dWave = crossSection->CalcAverageDWaveResWidth();
	                    std::cout << "  Average S-Wave Radiative Width [meV]: " 
	                    	    << 1.e9*sWave.first
	                    	    << std::endl;
	                    std::cout << "Average S-Wave Resonance Spacing [keV]: " 
	                    	    << 1.e3*sWave.second
	                    	    << std::endl;	
	                    std::cout << "  Average P-Wave Radiative Width [meV]: " 
	                    	    << 1.e9*pWave.first 
	                    	    << std::endl;
	                    std::cout << "Average P-Wave Resonance Spacing [keV]: " 
	                    	    << 1.e3*pWave.second
	                    	    << std::endl;	
	                    std::cout << "  Average D-Wave Radiative Width [meV]: " 
	                    	    << 1.e9*dWave.first 
	                    	    << std::endl;
	                    std::cout << "Average D-Wave Resonance Spacing [keV]: " 
		                << 1.e3*dWave.second
		                << std::endl;	
	                } else {
	                    crossSection->Calculate();
	                    crossSection->PrintCrossSections();
	                    if(printTrans) crossSection->PrintTransmissionTerms();
	                    
                        if(calcRates) {
	                        if(pType==1) {
	                            crossSection->CalculateReactionRates(true);
	                            crossSection->PrintReactionRates(true);
	                        }
	                        crossSection->CalculateReactionRates(false);
	                        crossSection->PrintReactionRates(false);	  
	                    }
	                }
                } else {
	            std::cout << "Could not calculate cross section." << std::endl;
                }
                delete crossSection;
            }
            return 0;
        }

        /** 
         * If there are more arguments, then we are apparently in
         * the statistical decay calculations, right?
         * This needs to change. This can be made much more readable and maintainable.
         */
        auto start = std::chrono::steady_clock::now();
        int chunkSize = 10000; /**< Portion of events which are written to the root tree.*/

        int A,Z,Pi,events; /** Apparently placeholders for cmd line parameters*/
        double J,lowEnergy,highEnergy; /** Apparently placeholders for cmd line parameters*/
        
        /**
         * Cmd line parameters for the decay mode will be read.
         */
        if(!parseCommandLineForDecay(args,Z,A,J,Pi,lowEnergy,highEnergy,events)) {
            std::cout << "Isotope, J-Pi, and Energy must be specified for decay (i.e. 60Fe 1- 15)." << std::endl;
            return -1;
        } else {
  	        std::cout << "Input Values For Parent Nucleus:" << std::endl
  	          << std::setw(15) << "Z:"             << std::setw(12) << Z          
  	          << std::setw(0) << std::endl
  	          << std::setw(15) << "A:"             << std::setw(12) << A          
  	          << std::setw(0) << std::endl
  	          << std::setw(15) << "J:"             << std::setw(12) << J          
  	          << std::setw(0) << std::endl
  	          << std::setw(15) << "Pi:"            << std::setw(12) << Pi         
  	          << std::setw(0) << std::endl
  	          << std::setw(15) << "energy (low):"  << std::setw(12) << lowEnergy  
  	          << std::setw(0) << std::endl
  	          << std::setw(15) << "energy (high):" << std::setw(12) << highEnergy 
  	          << std::setw(0) << std::endl
  	          << std::setw(15) << "events:"        << std::setw(12) << events     
  	          << std::setw(0) << std::endl
              << std::setw(15) << "chunk size:"    << std::setw(12) << chunkSize     
  	          << std::setw(0) << std::endl;
        }

        /**
         * Ah, I see. The total events will be divided in chunks.
         */
        int remainder = events%chunkSize;
        int chunks = (events-remainder)/chunkSize;

        std::cout << "Starting Decay Simulation..." << std::endl;

        int numDecayed = 0; /** Initialize the counting of decays.*/
        DecayResults* results = NULL;
        if(events>1) results = new DecayResults(Z,A,J,Pi,lowEnergy,highEnergy,suffixNo);
        
        
        omp_lock_t writelock;
        omp_init_lock(&writelock);
  

        for(int i = 0;i<=chunks;i++) {

          
            int numInChunk = (i==chunks) ? remainder : chunkSize;
            if(numInChunk==0) continue;
            if(events>=chunkSize)
              std::cout << "Decay chunk " << i+1 << " of " << chunks << " started ..." << std::endl;
            else
              std::cout << "Decay of " << numInChunk << " nuclei started ..." << std::endl;
            std::vector<std::pair<DecayData,std::vector<DecayProduct> > > chunkResults;
            chunkResults.resize(numInChunk);
            /**
             * Using Open Multiprocessing (OMP) for the parallel execution of the following for-loop.
             * Nice.
             */
            



            ProgressBar pg;
            pg.start(numInChunk); 

            
            #pragma omp parallel for
            for(int j = 0;j<numInChunk;j++) {
                
                //std::this_thread::sleep_for(std::chrono::seconds(1));
                int localNumDecayed = numDecayed++;

                //if(localNumDecayed%(events/20)==0&&localNumDecayed>0){
                  //if(j%(numInChunk/20)==0)
                    pg.update(j);
                //std::cout << "Decayed " << localNumDecayed 
	    	        //  << " of " << events << " nuclei..." << std::endl;
                //}
                double energy = (lowEnergy==highEnergy) ? lowEnergy :
    
                lowEnergy+(highEnergy-lowEnergy)*double(rand_r(&randomSeed[omp_get_thread_num()]))/double(RAND_MAX);

                DecayController* controller;

                if(preEq) {
	               controller= new DecayController(Z,A,J,Pi,energy,numNuParticles,numNuHoles,numPiParticles,numPiHoles);
                } else controller = new DecayController(Z,A,J,Pi,energy);

                double neutronEntranceWidth = 0.;
                double protonEntranceWidth = 0.;
                double gammaEntranceWidth = 0.;
                double alphaEntranceWidth = 0.;
                double neutronTotalWidth = 0.;
                double protonTotalWidth = 0.;
                double gammaTotalWidth = 0.;
                double alphaTotalWidth = 0.;

                controller->Decay(neutronEntranceWidth,protonEntranceWidth,alphaEntranceWidth,gammaEntranceWidth,
	    		    neutronTotalWidth,protonTotalWidth,alphaTotalWidth,gammaTotalWidth); 

                chunkResults[j] = 
	            std::pair<DecayData,std::vector<DecayProduct> >(DecayData(energy,neutronEntranceWidth,protonEntranceWidth,
	    							       alphaEntranceWidth,gammaEntranceWidth,
	    							       neutronTotalWidth,protonTotalWidth,
	    							       alphaTotalWidth,gammaTotalWidth),controller->DecayProducts());

                if(events==1) controller->PrintDecays();
                delete controller;
            }
            
            omp_set_lock(&writelock);            
            if(events>1){
                std::cout << std::endl << "Writing ROOT Tree..." << std::endl;
                
                results->AddResults(chunkResults);
            }
            omp_unset_lock(&writelock);
        }   
        omp_destroy_lock(&writelock);
        

        if(results) delete results;

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << std::endl << "Total calculation time: " << elapsed_seconds.count() << "s\n";

        
        while(omp_get_num_threads()>1){
          std::this_thread::sleep_for(std::chrono::seconds(1));
          std::cout << "Waiting ..." << std::endl;
        }
        return 0;
    }
}
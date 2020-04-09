#include "NuclearLevels.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <string.h>

LevelsContainer::LevelsContainer(std::istream& in, int numLevels, int numComplete) {
  for(int i=0;i<numLevels;i++) {
    std::string level;
    std::getline(in,level);

    std::istringstream levelStream;
    levelStream.str(level);

    std::string dummy;
    double energy;
    double J;
    int Pi;
    int numGammas;
    levelStream >> dummy >> energy >> J >> Pi;
    levelStream.ignore(11);
    levelStream >> numGammas; 
    Level tempLevel(J,Pi,energy);
    for(int j = 0;j<numGammas;j++) {
      std::string gammaLine;
      std::getline(in,gammaLine);
      std::istringstream gammaStream;
      gammaStream.str(gammaLine);
      int levelIndex;
      double energy,probability;
      gammaStream >> levelIndex >> energy >> probability;
      tempLevel.gammas_.push_back(GammaTransition(levelIndex,energy,probability));
    }
    if(i<numComplete&&(Pi!=0||energy==0.)&&J>=0.) levels_.push_back(tempLevel);
    //if((Pi!=0||energy==0.)&&J>=0.) levels_.push_back(tempLevel);
  }
}

void NuclearLevels::InitializeLevels(std::string levelsDirectory,
				     std::string spinFile) {
#ifndef MPI_BUILD
  std::cout << "Reading known nuclear levels..." << std::endl;
#endif
  for(int Z = 0; Z<=118; Z++) {
    char isotopeFile[256];
    sprintf(isotopeFile,"%sz%03d.dat",levelsDirectory.c_str(),Z);
    std::ifstream in(isotopeFile);
    if(!in) continue;
    while(!in.eof()) {
      std::string line;
      std::getline(in,line);
      if(!in.eof()) {
	int A, numLevels, numComplete;
	std::string dummy;
	std::istringstream stm;
	stm.str(line);
	stm >> dummy >> A >> dummy >> numLevels >> dummy >> numComplete;
	if(numLevels!=0) levelsTable_[MassKey(Z,A)]=LevelsContainer(in,numLevels,numComplete);
      }
    }
    in.close();
  }
  std::ifstream in(spinFile.c_str());
  if(!in) return;
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
	} else  {
	  if(((spinZ+1)/2)%2==0) parallelZ = true;
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
      LevelsTable::iterator it = levelsTable_.find(MassKey(Z,A));
      if(it==levelsTable_.end()) {
	LevelsContainer levelsContainer;
	levelsContainer.levels_.push_back(Level(spin,parity,0.));
	levelsTable_[MassKey(Z,A)] = levelsContainer;
      } else if(it->second.levels_.size() == 0) {
	it->second.levels_.push_back(Level(spin,parity,0.));
      } else if(it->second.levels_.size() > 0 &&
		it->second.levels_[0].energy_ != 0.) {
	it->second.levels_.insert(it->second.levels_.begin(),
				  Level(spin,parity,0.));
      }
      it = levelsTable_.find(MassKey(Z,A));
      if(it->second.levels_[0].Pi_==0) it->second.levels_[0].Pi_=parity;
    }
  }
}

void NuclearLevels::PrintLevels(int Z, int A) {
  LevelsTable::const_iterator it = levelsTable_.find(MassKey(Z,A));
  if(it!=levelsTable_.end()) {
    for(std::vector<Level>::const_iterator level = it->second.levels_.begin();
	level<it->second.levels_.end();level++) {
      std::cout << std::setw(12) << level->energy_ 
		<< std::setw(5)  << level->J_  
		<< std::setw(5)  << level->Pi_
		<< std::setw(0)  << std::endl;
      for(std::vector<GammaTransition>::const_iterator gamma = level->gammas_.begin();
	  gamma<level->gammas_.end();gamma++)
	std::cout << '\t' << std::setw(5) << gamma->levelIndex_
		  << std::setw(12) << gamma->energy_
		  << std::setw(12) << gamma->probability_
		  << std::setw(0)  << std::endl;
    }
  }
}

std::vector<Level> NuclearLevels::FindLevels(int Z, int A) {
  std::vector<Level> levelsVector;
  
  LevelsTable::const_iterator it = levelsTable_.find(MassKey(Z,A));
  if(it!=levelsTable_.end()) levelsVector=it->second.levels_;
    
  return levelsVector;
}

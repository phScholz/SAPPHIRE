/**
 * @file WidthEntry.h
 * @author Philipp Scholz
 * @brief Class to store entries in the partial widths file for Module_BreitWigner
 * @date 2020-06-29
 * 
 * 
 */

#pragma once
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

class WidthEntry{
    private:
        void parseJPiString(double & J, int & Pi, std::string & jPiString){
            bool foundDelimiter = false;
            bool goodPi = false;
            bool goodJ = false;

            std::string firstString,secondString,parityString,delimiterString;

            for(int i = 0;i<jPiString.length();i++) {

            	std::string nextChar(jPiString,i,1);


            	if(nextChar == "/" || nextChar == ".") {
            	  foundDelimiter=true;
            	  delimiterString=nextChar;
            	  continue;

            	} else if(nextChar=="+" || nextChar=="-") {
            	    parityString = nextChar;
                    break;
                }

                std::istringstream stm(nextChar);
                int digit;

                if(stm>>digit) {
                    if(!foundDelimiter){
                        firstString+=nextChar;
                    }
                    else{
                        secondString+=nextChar;
                    }
                }
            }

            if(parityString.length()>0) {
              if(parityString=="-") Pi=-1;
              else if(parityString=="+") Pi=1;
              goodPi=true;
            }

            if(firstString.length()>0) {
              	if(foundDelimiter&&delimiterString=="/") {
                  
              	    if(secondString.length()>0){
                        J = atof(firstString.c_str())/atof(secondString.c_str());

                    }
              	    else{
                        J = atof(firstString.c_str());
                    }

              	} else if(foundDelimiter&&delimiterString==".") {
              	  firstString += '.'+secondString;
              	  J = atof(firstString.c_str());

              	} else {
                    J = atof(firstString.c_str());
                }

              	double intPart;

                if(modf(J*2.,&intPart)==0.) goodJ=true; 

              }
        };

        double spinDoubleFromString(std::string &jPiString){
      
          double J;
          int Pi;

          parseJPiString(J,Pi, jPiString);

          return J;
        }

        int parityIntFromString(std::string &jPiString){
      
            double J;
            int Pi;

            parseJPiString(J,Pi, jPiString);

            return Pi;
        };

    public:
        /**
         * @brief Contructor initializing member variables
         */
        WidthEntry(double energy, double spin, int parity; double width; double widtherror) : energy_(energy), spin_(spin), parity_(parity), width_(width), widthError_(widtherror){
        };

        /**
         * @brief Construct a new Width Entry object from a line in a partial width file
         * 
         * @param line 
         */
        WidthEntry(std::string line){
            std::istringstream lineStream(line);
            lineStream >> energy_ >> jPiString_ >> width_ >> widthError_;
            spin_ = spinDoubleFromString(jPiString_);
            parity_ = parityIntFromString(jPiString_);
        };

        double energy_; /*<< Energy of the particle*/
        double spin_;   /*<< Spin of the compound states*/
        int parity_;    /*<< Parity of the compound states*/

        std::string jPiString_; /*<< String of the JPI input in the partial width file*/

        double width_;      /*<< Partial width*/
        double widthError_; /*<< Uncertainty of the Partial Width*/

        

};
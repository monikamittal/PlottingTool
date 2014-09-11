#ifndef GenerateRandomNumber_H
#define GenerateRandomNumber_H

#include "TRandom.h"

class GenerateRandomNumber {
 public :
   GenerateRandomNumber(){};
   std::vector<double>  mygauss(int n_randomnumbers, double randomseed, double mean, double sigma);
   ~GenerateRandomNumber(){};

};

#endif

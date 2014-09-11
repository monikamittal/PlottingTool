#define ReducibleBackgroundEstimate_cxx
#include "GenerateRandomNumber.h"
#include <iostream>

// Generate Gaussian distributed Random numbers
std::vector<double>  GenerateRandomNumber::mygauss(int n_randomnumbers, double randomseed, double mean, double sigma){
  TRandom rndm;
  std::vector<double>  vec;
  rndm.SetSeed(randomseed);
  for (int ii =0; ii<n_randomnumbers; ii++ ){
    //std::cout<<rndm.Gaus(mean,sigma)<<std::endl;
    vec.push_back(rndm.Gaus(mean,sigma));
  }
  return(vec);
}



// Generate Exponential distributed Random numbers
std::vector<double>  GenerateRandomNumber::mygauss(int n_randomnumbers, double randomseed, double mean, double sigma){
  TRandom rndm;
  std::vector<double>  vec;
  rndm.SetSeed(randomseed);
  for (int ii =0; ii<n_randomnumbers; ii++ ){
    //std::cout<<rndm.Gaus(mean,sigma)<<std::endl;
    vec.push_back(rndm.Gaus(mean,sigma));
  }
  return(vec);
}

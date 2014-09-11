#include <cstdlib>
#include <iostream>
#include <string>
#include <stdio.h>
#include <fstream>
#include "TMath.h"
#include <vector>
#include "TH1F.h"
#include "TRandom.h"

#include "../PlottingTools/interface/CreateCanvas.h"

double findLogLikelihood(double bi  , double bi_bs , double data ){
  double LogLi  ;
  LogLi = (bi - bi_bs) + (data * log(bi_bs/bi)) ;
  return LogLi;
}


int main(){
  CreateCanvas canvas;
  
  TRandom t;
  const int ndataentries = 100000;
  const int nbins = 1000;
  const int nbin_likelihood = 1000;
  const int nExperiment = 10000;
  
  double mean[ndataentries];
  double xlow, xhigh;
  xlow = -10.0;
  xhigh = -8.0; 
  TH1F *h1 =new TH1F("exp", "exponential",nbins,0,25);
  TH1F *likelihood =new TH1F("likelihood", "likelihood", nbin_likelihood , xlow , xhigh);
  
  //----------finding bi for hypothesis Ho = exp(-3t)
  //Hypothesis H0 corresponding to background only
  double integral = 0.0;
  std::vector<double>  bi,x;
  bi.clear();
  x.clear();
  double binwidth = h1->GetBinWidth(2);
  
  //for(int i=1;i<nbins;i++){
  for(int i=1;i<=nbins;i++){
    x.push_back( double ((i- 0.5)*binwidth + 0.0)) ;
    bi.push_back( 3*( exp ( -3.0 * (x[i-1]) ) )*binwidth );
    //std::cout<<" bi = "<<bi[i-1]<<std::endl;
    integral += 3*( exp ( -3.0 * (x[i-1]) ) )*binwidth ;
  }
  std::cout<<" integral of H0 function = "<<integral<<std::endl;
  
  
  //----------finding bi + si for hypothesis H` = exp(-5t)
  double integral_bs = 0.0;
  std::vector<double> bi_bs, x_bs;
  x_bs.clear();
  bi_bs.clear();
  
  for(int i=1;i<=nbins;i++){
    x_bs.push_back( double ((i- 0.5)*binwidth + 0.0)) ;
    bi_bs.push_back( 5.0 * (exp ( -5.0 * x_bs[i-1] )) * binwidth );
    integral_bs += bi_bs[i-1] ;
  }
  std::cout<<" integral of H1 function = "<<integral_bs<<std::endl;
  
  //---------Data exponentially distributed
  double LogLikelihood,data_for_p_value;
  double data[nbins];
  double alpha=5.0;
  double LogLi ;
  for(int iExp = 0; iExp < nExperiment+1 ; iExp++){
    t.SetSeed(iExp+6);
    for(int i=0;i<ndataentries;i++){
      double q=t.Exp(alpha); 
      h1->Fill(q);
    }
    h1->Scale(1/h1->Integral());
    
    LogLikelihood = 0.0 ;
    for(int i=1;i<nbins;i++){
      data[i] = h1->GetBinContent(i);
      LogLi = findLogLikelihood ( bi[i-1] , bi_bs[i-1] , data[i]);
      LogLikelihood = LogLi + LogLikelihood ;
      //      std::cout<<" LogLikelihood = "<<LogLikelihood<<std::endl;
    }
    if(iExp<nExperiment){
      likelihood->Fill(LogLikelihood);
    }
    if(iExp==nExperiment) {
      data_for_p_value = LogLikelihood ;
    }
  }
  if(likelihood->Integral()!=0)
    likelihood->Scale(1.0/likelihood->Integral());
  TCanvas* c = canvas.declareCanvas();
  likelihood->Draw();
  c->SaveAs("loglikelyhood.png");

  std::cout<<" observed data point = "<<data_for_p_value<<std::endl;
  int obs_bin = int(abs((xlow-data_for_p_value)/likelihood->GetBinWidth(2))); //tobs is at this bin
  std::cout<<" obs_bin = "<<obs_bin<<std::endl;
  double pValue = likelihood->Integral(obs_bin,nbin_likelihood);
  //if( LogLikelihood < likelihood->GetMean())
  // pValue = 1 - pValue ;
  std::cout<<" pValue =  "<<pValue<<" integral   "<<likelihood->Integral()<<std::endl;
 
  return 0;
}

#include"TRandom2.h"
#include "TH1F.h"
void exp(){
  TRandom2 t ;
  const int n = 1000;
  double mean[n];
  const int nbins = 1000;
  double xlow, xhigh;
  const int nbin_exp = 100;
  xlow = -10.0;
  xhigh = -0.0; 
  TH1F *h1 =new TH1F("exp", "exponential",nbins,0,25);
  TH1F *likelihood =new TH1F("likelihood", "likelihood", nbin_exp , xlow , xhigh);
  //----------finding bi for hypothesis Ho = exp(-3t)
  double integral = 0.0;
  double bi[nbins]; 
  double binwidth = h1->GetBinWidth(2);
  //std::cout<<"binwidth ="<<binwidth<<std::endl;
  double x[nbins];
  for(int i=1;i<=nbins;i++){
    x[i] = double ((i- 0.5)*binwidth + 0.0) ;
    bi[i] = 3*( exp ( -3.0 * (x[i]) ) )*binwidth ;
    //cout<<"bin center["<<i<<"] (calculated by hand)= "<<x[i]<<endl;
    //x[i] = double (h1->GetXaxis()->GetBinCenter(i)) ;
    //std::cout<<" x = "<<x[i] << std::endl;
    //bi[i] += 3*( exp ( -3.0 * (x[i]+0.5*binwidth) )) ;
    //bi[i] *= 0.5;
    integral += bi[i] ;
    //cout<<"bin center["<<i<<"] = "<<x[i]<<" and bi "<<bi[i]<<endl;
  }
  std::cout<<" integral of H0 function = "<<integral<<std::endl;
  //----------finding bi + si for hypothesis H` = exp(-5t)
  double integral_bs = 0.0;
  double bi_bs[nbins];
  double x_bs[nbins];
  for(int i=1;i<=nbins;i++){
    x_bs[i] = double ((i- 0.5)*binwidth + 0.0) ;
    //x_bs[i] = double (h1->GetXaxis()->GetBinCenter(i)) ;
    bi_bs[i] = 5.0 * (exp ( -5.0 * x_bs[i] )) * binwidth ;
    integral_bs += bi_bs[i] ;
    //cout<<"bin center["<<i<<"] = "<<x_bs[i]<<" and bi_bs "<<bi_bs[i]<<endl;
  }
  std::cout<<" integral of H1 function = "<<integral_bs<<std::endl;
  //---------Data exponentially distributed
  double data[nbins];
  double alpha=5.0;
  double LogLi ;
  const int nExperiment = 1000;
  for(int iExp = 0; iExp < nExperiment+1 ; iExp++){
    t.SetSeed(iExp+9);
    for(int i=0;i<n;i++){
      double q=t.Exp(alpha); 
      h1->Fill(q);
    }
    h1->Scale(1/h1->Integral());
    double LogLikelihood = 0.0 ;
    for(int i=1;i<nbins;i++){
      data[i] = h1->GetBinContent(i);
      //std::cout<<" bi[i] = "<<bi[i]<<"   bi_bs[i]   = "<<bi_bs[i]<<"   data[i]  =  "<<data[i]<<std::endl;
      LogLi = findLogLikelihood ( bi[i] , bi_bs[i] , data[i]);
      LogLikelihood = LogLi + LogLikelihood ;
      //std::cout<<"LogLikelihood = "<<LogLikelihood<<std::endl;
    }
    if(iExp<nExperiment){
      likelihood->Fill(LogLikelihood);
    }
    if(iExp==nExperiment) {
      //std::cout<<" observed data point = "<<LogLikelihood<<std::endl;
    }
    //std::cout<<" LogLikleihood =  "<<LogLikelihood<<std::endl;
  }
  if(likelihood->Integral()!=0)
  likelihood->Scale(1.0/likelihood->Integral());
  likelihood->Draw();
  std::cout<<" observed data point = "<<LogLikelihood<<std::endl;
  int obs_bin = int(abs((xlow-LogLikelihood)/likelihood->GetBinWidth(2))); //tobs is at this bin
  double pValue = likelihood->Integral(obs_bin,nbin_exp);
  if( LogLikelihood < likelihood->GetMean())
    pValue = 1 - pValue ;
  std::cout<<" pValue =  "<<pValue<<" integral   "<<likelihood->Integral()<<std::endl;
  
}

double findLogLikelihood(double bi  , double bi_bs , double data ){
  double LogLi  ;
  LogLi = (bi - bi_bs) + (data * log(bi_bs/bi)) ;
  //  std::cout<<" LogLi = "<<LogLi<<std::endl;
  return LogLi;
}

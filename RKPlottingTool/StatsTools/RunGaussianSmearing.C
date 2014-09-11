#include <cstdlib>
#include <iostream>
#include <string>
#include <stdio.h>
#include <fstream>
#include "TMath.h"
#include <vector>
#include "TH1F.h"
#include "GenerateRandomNumber.C"
#include "CreateCanvas.C"

int main(){
  GenerateRandomNumber generate_rndm;
  CreateCanvas canvas;

  
  TH1F* h_wide = new TH1F("h_wide","h_wide",100,-10.,15.);
  TH1F* h_narrow = new TH1F("h_narrow","h_narrow",100,-10.,15.);
  TH1F* h_smeared = new TH1F("h_smeared","h_smeared",100,-10.,15.);
  
    
  std::vector<double> gausdistributed_wide =  generate_rndm.mygauss(200000,1.0,2.,4);//int n_randomnumbers, double randomseed, double mean, double sigma
  std::vector<double> gausdistributed_narrow =  generate_rndm.mygauss(200000,2.0,2.,3.5);//int n_randomnumbers, double randomseed, double mean, double sigma
  double sigma_smear = TMath::Sqrt( (4.*4.) - (3.5*3.5) );
  std::vector<double> gausdistributed_smear =  generate_rndm.mygauss(200000,3.0,0.,sigma_smear);//int n_randomnumbers, double randomseed, double mean, double sigma

  
  for(int i=0; i<int(gausdistributed_wide.size());i++) {
    //std::cout<<gausdistributed_wide[i]<<" : "<<gausdistributed_narrow[i]<<" : "<<gausdistributed_smear[i]<<std::endl;
    h_wide->Fill(gausdistributed_wide[i]);
    h_narrow->Fill(gausdistributed_narrow[i]);
    h_smeared->Fill(gausdistributed_narrow[i] +  gausdistributed_smear[i]  );
  }    
  h_wide->SetLineColor(2);
  h_narrow->SetLineColor(4);
  h_smeared->SetLineColor(1);
  
  h_wide->SetLineWidth(2);
  h_narrow->SetLineWidth(2);
  h_smeared->SetLineWidth(2);
  
  h_wide->SetMaximum(10000.);
  //h_narrow->SetMaximum(10000.);
  //h_smeared->SetLineWidth(10000.0);
  TCanvas* c = canvas.declareCanvas();
  h_wide->Draw();
  h_narrow->Draw("sames");
  h_smeared->Draw("sames");
  c->SaveAs("plot2.png");
  return (0.0);
}

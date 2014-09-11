#include "TROOT.h"
#include <cstdlib>
#include <iostream>
#include <string>
#include <stdio.h>
#include <fstream>
#include "TMath.h"
#include <vector>
#include "TH1F.h"
#include "TRandom.h"
#include "TFile.h"
#include "TString.h"
#include <sstream>
#include "TLegend.h"
#include "TStyle.h"
#include "THStack.h"
#include <algorithm>
#include "../interface/PlotStackedHisto.h"
#include"../interface/CreateCanvas.h"


void PlotStackedHisto:: DrawHistograms(std::vector<FileInfo> fobjVec, std::vector<HistoInfo> hobjVec){
  CreateCanvas CC;
  TLegend* leg;
  double setmaximum_ = 0.;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetFrameLineWidth(3);

  THStack* hstack = new THStack("generalised_stack","");
  
  TCanvas* c = CC.declareCanvasFor1Plot();
  std::cout<<" Plot Stack  ---- finfiVec size ---"<<fobjVec.size()<<std::endl;
  
  // fonjVec.size gives the number of block in the text file
  for(int ifilevec = 0; ifilevec < (int) fobjVec.size(); ifilevec++){
    const int NFILES = fobjVec[ifilevec].filename.size();
    TFile* f[NFILES];
    
    // how many plots to be made (how many blocks)
    for(int ihisto =0; ihisto < (int) hobjVec.size(); ihisto++){
      const int NHISTO = hobjVec[ihisto].histname.size();
      double x1 = hobjVec[ihisto].legpos[0];
      double y1 = hobjVec[ihisto].legpos[1];
      double x2 = hobjVec[ihisto].legpos[2];
      double y2 = hobjVec[ihisto].legpos[3];
      
      std::cout<< " legend position = "<<x1 
	       <<"  "<<y1
	       <<"  "<<x2
	       <<"  "<<y2
	       <<std::endl;
      
      leg = new TLegend(x1,y1,x2,y2,NULL,"brNDC");
      leg->SetBorderSize(0);
      leg->SetLineColor(1);
      leg->SetLineStyle(1);
      leg->SetLineWidth(1);
      leg->SetFillColor(0);
      leg->SetFillStyle(0);
      leg->SetTextFont(22);
      
      const int nhistos = NHISTO;
      TH1F* h[nhistos] ;
      
      for (int ifile=0; ifile < NFILES; ifile++) {
	f[ifile]  = (TFile*)gROOT->GetListOfFiles()->FindObject(fobjVec[ifilevec].filename[ifile].c_str());
	if (!f[ifile])
	  f[ifile]  = new TFile(fobjVec[ifilevec].filename[ifile].c_str());
	f[ifile]->cd();
	std::cout<<" filename = "<<fobjVec[ifilevec].filename[ifile].c_str()<<std::endl;
	std::cout<<" dirname = "<<fobjVec[ifilevec].dirname.c_str()<<std::endl;
	
	// this is helpful when all histograms are in same root file
	for (int IHISTO=0; IHISTO < NHISTO; IHISTO++) {
	  std::cout<< "--------------- Inside NHISTO Loop "<<IHISTO<<std::endl;
	  TString dir = fobjVec[ifilevec].dirname.c_str();
	  std::cout<<" size = "<<dir.Sizeof()<<std::endl;
	  TString histo;
	  // check if input histo is in a directory
	  if(dir.Sizeof()==1)
	    histo = hobjVec[ihisto].histname[IHISTO].c_str();
	  if(dir.Sizeof()>1)
	    histo = dir+"/"+hobjVec[ihisto].histname[IHISTO].c_str();
	  
	  h[ifile] = (TH1F*) f[ifile]->Get(histo);
	  h[ifile]->GetXaxis()->SetRangeUser(hobjVec[ihisto].XRange[0], hobjVec[ihisto].XRange[1]);
	  double integral = h[ifile]->Integral();
	  //double scaleFactor = xsection * lumi / totaleventhist ;
	  double scaleFactor = 1.0;
	  h[ifile]->Scale(scaleFactor);
	  	  
	  //std::cout<<" added to stack "<<std::endl;
	  //h[ifile]->GetXaxis()->SetLabelFont(22);
	  //h[ifile]->GetYaxis()->SetLabelFont(22);
	  //TString xaxis = hobjVec[ihisto].XAXISTitle.c_str();
	  //h[ifile]->GetXaxis()->SetTitle(hobjVec[ihisto].XAXISTitle.c_str());
	  //h[ifile]->GetYaxis()->SetTitle(hobjVec[ihisto].YAXISTitle.c_str());
	  //h[ifile]->GetXaxis()->SetTitleFont(22);
	  //h[ifile]->GetYaxis()->SetTitleFont(22);
	  //h[ifile]->GetXaxis()->SetTitleSize(0.045);
	  //h[ifile]->GetYaxis()->SetTitleSize(0.045);
	  //h[ifile]->GetXaxis()->SetLabelSize(0.045);
	  //h[ifile]->GetYaxis()->SetLabelSize(0.045);
	  //h[ifile]->GetYaxis()->SetTitleOffset(1.45);
	  //std::cout<<" xmin, xmax "<<hobjVec[ihisto].XRange[0] <<"    "<<hobjVec[ihisto].XRange[1]<<std::endl;
	  //h[ifile]->GetXaxis()->SetRangeUser(hobjVec[ihisto].XRange[0], hobjVec[ihisto].XRange[1]);
	  
	  
	  if(ifile!=NFILES-1) {
	    h[ifile]->SetLineColor(kBlack);
	    h[ifile]->SetFillColor(fobjVec[ifilevec].color[ifile]); // if more than 1 input file
	    hstack->Add(h[ifile]);
	    if(ifile == (NFILES-2) ) hstack->Draw();
	    std::cout<<" added histogram to stack "<<std::endl;
	  }
	  else{
	    std::cout<<" inside data "<<std::endl;
	    h[ifile]->Draw("SAME,E1");
	    h[ifile]->SetMarkerColor(kBlack);
	    h[ifile]->SetMarkerStyle(kFullCircle);
	    h[ifile]->SetMarkerSize(1.2);
	  }
	  
	  //std::cout<<" int = "<<h[ifile]->Integral()<<std::endl;
	  TString legendtext ;
	  if(NFILES>1) legendtext = fobjVec[ifilevec].legend[ifile].c_str() ;
	  if(NHISTO>1) legendtext = hobjVec[ihisto].legend[IHISTO].c_str() ;
	  
	  if(ifile<NFILES-1) leg->AddEntry(h[ifile],legendtext,"F");
	  else leg->AddEntry(h[ifile],legendtext,"P");
	}//for (int IHISTO=0; IHISTO < NHISTO; IHISTO++) {
      }//for (int ifile=0; ifile < NFILES; ifile++) {
      leg->Draw();
      
      TPaveText *pt = new TPaveText(0.1577181,0.9,0.9580537,0.96,"brNDC");
      pt->SetBorderSize(0);
      pt->SetTextAlign(12);
      pt->SetFillStyle(0);
      pt->SetTextFont(22);
      pt->SetTextSize(0.045);
      TText *text = pt->AddText(0.01,0.5,"CMS Simulation");
      //text = pt->AddText(0.3,0.6,"#sqrt{s} = 8 TeV, L = 19.8 fb^{-1}");
      text = pt->AddText(0.4,0.5,"#sqrt{s} = 8 TeV");
      pt->Draw();
      
      TString plotname = hobjVec[ihisto].plotname.c_str();
      TString prefix = fobjVec[ifilevec].prefixname.c_str();
      c->SaveAs("plots/"+prefix+plotname+".pdf");
      c->SaveAs("plots/"+prefix+plotname+".root");
      c->SaveAs("plots/"+prefix+plotname+".C");
      c->SaveAs("plots/"+prefix+plotname+".png");
    }//for(int ihisto =0; ihisto < (int) hobjVec.size(); ihisto++){
  }//for(int ifilevec = 0; ifilevec < (int) fobjVec.size(); ifilevec++){
}

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
#include <algorithm>
#include "../interface/PlotPOverlappingHisto.h"
#include"../interface/CreateCanvas.h"


void PlotPOverlappingHisto:: DrawHistograms(std::vector<FileInfo> fobjVec, std::vector<HistoInfo> hobjVec){
  CreateCanvas CC;
  TLegend* leg;
  double setmaximum_ = 0.;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetFrameLineWidth(3);
  
  TCanvas* c = CC.declareCanvasFor1Plot();
  std::cout<<" PlotOverlap  ---- finfiVec size ---"<<fobjVec.size()<<std::endl;
  
  // fonjVec.size gives the number of block in the text file
  for(int ifilevec = 0; ifilevec < (int) fobjVec.size(); ifilevec++){
    const int NFILES = fobjVec[ifilevec].filename.size();
    TFile* f[NFILES];
    TH1F* h ;

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
      //leg->SetTextFont(60);
      //    SetLabelFont(22);
      
      // file information of One single block
      std::vector<double> maxi;
      maxi.clear();

      for (int ifile=0; ifile < NFILES; ifile++) {
	f[ifile]  = (TFile*)gROOT->GetListOfFiles()->FindObject(fobjVec[ifilevec].filename[ifile].c_str());
	if (!f[ifile])
	  f[ifile]  = new TFile(fobjVec[ifilevec].filename[ifile].c_str());
	f[ifile]->cd();
	std::cout<<" filename = "<<fobjVec[ifilevec].filename[ifile].c_str()<<std::endl;
	std::cout<<" dirname = "<<fobjVec[ifilevec].dirname.c_str()<<std::endl;
	
	// find maximum of all the histograms
	if(fobjVec[ifilevec].option == "overlap"){
	  for (int IHISTO=0; IHISTO < NHISTO; IHISTO++) {
	    TString dir = fobjVec[ifilevec].dirname.c_str();
	    std::cout<<" size = "<<dir.Sizeof()<<std::endl;
	    TString histo;
	    if(dir.Sizeof()==1) 	      histo = hobjVec[ihisto].histname[IHISTO].c_str();
	    if(dir.Sizeof()>1)    	      histo = dir+"/"+hobjVec[ihisto].histname[IHISTO].c_str();
	    h = (TH1F*) f[ifile]->Get(histo);
	    h->Rebin(hobjVec[ihisto].ReBin);
	    double integral = h->Integral();
	    if(fobjVec[ifilevec].normalizeUnity==1)   h->Scale(1.0/integral);
	    maxi.push_back(h->GetMaximum());
	    std::cout<<"----maxi-------"<<h->GetMaximum()<<std::endl;
	  }
	}
	setmaximum_ = *(std::max_element(maxi.begin(), maxi.end()));
	std::cout<<" ----------------------------------------------------------------- maximum of all is "<<1.2*setmaximum_<<std::endl;
      }
      

      for (int ifile=0; ifile < NFILES; ifile++) {
	f[ifile]  = (TFile*)gROOT->GetListOfFiles()->FindObject(fobjVec[ifilevec].filename[ifile].c_str());
	if (!f[ifile])
	  f[ifile]  = new TFile(fobjVec[ifilevec].filename[ifile].c_str());
	f[ifile]->cd();
	std::cout<<" filename = "<<fobjVec[ifilevec].filename[ifile].c_str()<<std::endl;
	std::cout<<" dirname = "<<fobjVec[ifilevec].dirname.c_str()<<std::endl;
	
	// number of histograms to be overlapped
	
	
	for (int IHISTO=0; IHISTO < NHISTO; IHISTO++) {
	  TString dir = fobjVec[ifilevec].dirname.c_str();
	  //std::cout<<" size = "<<dir.Sizeof()<<std::endl;
	  //f[ifile]->cd(fobjVec[ifilevec].dirname.c_str());
	  TString histo;
	  // check if input histo is in a directory
	  if(dir.Sizeof()==1)
	    histo = hobjVec[ihisto].histname[IHISTO].c_str();
	  if(dir.Sizeof()>1)
	    histo = dir+"/"+hobjVec[ihisto].histname[IHISTO].c_str();
	  
	  //std::cout<< " histo = "<<histo<<std::endl;
	  
	  h = (TH1F*) f[ifile]->Get(histo);
	  h->GetXaxis()->SetRangeUser(hobjVec[ihisto].XRange[0], hobjVec[ihisto].XRange[1]);
	  // Added for suvankar
	  double integral = h->Integral();
	  
	  if(fobjVec[ifilevec].normalizeUnity==1)   h->Scale(1.0/integral);
	  h->SetMaximum(1.2*setmaximum_);
	  // Added ended
	  //h->Rebin(hobjVec[ihisto].ReBin);
	  
	  if(ifile==0 && IHISTO==0) {
	    h->Draw();
	    h->SetMaximum(1.2*setmaximum_);
	  }
	  
	  else {
	    //std::cout<<" inside else "<<std::endl;
	    h->Draw("same");
	    h->SetMaximum(1.2*setmaximum_);
	  }
	  h->SetLineWidth(2);
	  //std::cout<<" linewidth"<<std::endl;
	  if(NFILES>1) h->SetLineColor(fobjVec[ifilevec].color[ifile]); // if more than 1 input file
	  if(NHISTO>1) h->SetLineColor(hobjVec[ihisto].color[IHISTO]);  // if 1 file but multiple histograms
	  h->GetXaxis()->SetLabelFont(22);
	  h->GetYaxis()->SetLabelFont(22);
	  TString xaxis = hobjVec[ihisto].XAXISTitle.c_str();
	  //std::cout<<" xaxis = "<<xaxis<<std::endl;
	  h->GetXaxis()->SetTitle(hobjVec[ihisto].XAXISTitle.c_str());
	  h->GetYaxis()->SetTitle(hobjVec[ihisto].YAXISTitle.c_str());
	  //std::cout<<" title"<<std::endl;
	  h->GetXaxis()->SetTitleFont(22);
	  h->GetYaxis()->SetTitleFont(22);
	  
	  h->GetXaxis()->SetTitleSize(0.045);
          h->GetYaxis()->SetTitleSize(0.045);
	  h->GetXaxis()->SetLabelSize(0.045);
          h->GetYaxis()->SetLabelSize(0.045);
	  
	  h->GetYaxis()->SetTitleOffset(1.45);
	  std::cout<<" xmin, xmax "<<hobjVec[ihisto].XRange[0] <<"    "<<hobjVec[ihisto].XRange[1]<<std::endl;
	  h->GetXaxis()->SetRangeUser(hobjVec[ihisto].XRange[0], hobjVec[ihisto].XRange[1]);

	  //std::cout<<" int = "<<h->Integral()<<std::endl;
	  TString legendtext ;
	  if(NFILES>1) legendtext = fobjVec[ifilevec].legend[ifile].c_str() ;
	  if(NHISTO>1) { 
	    legendtext = hobjVec[ihisto].legend[IHISTO].c_str() ;
	    //std::cout<<" legend = "<<hobjVec[ihisto].legend[IHISTO].c_str() <<std::endl;
	  }
	  leg->AddEntry(h,legendtext,"L");
	  
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
    }
  }
}

#define DrawGraph_cxx
#include "../interface/DrawGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"

TCanvas* DrawGraph::DrawGraph_9pad_x_6graph(TCanvas* c, std::vector<TObject*> num, std::vector<TObject*> den , std::vector<TString> legend_text, std::vector<double>& FR_value){
  int ipad=0;
  TH1F* NUM[54];
  TH1F* DEN[54];
  std::vector<int> icolor;
  icolor.clear();
  icolor.push_back(1); 
  icolor.push_back(2); 
  icolor.push_back(3); 
  icolor.push_back(4); 
  icolor.push_back(5); 
  icolor.push_back(6); 
  TLegend *leg;
  double maximum;
  for(int i=0;i<int(num.size());i++){
    
    TGraphAsymmErrors *gr = new TGraphAsymmErrors();
    
    NUM[i] = (TH1F*)num[i]->Clone();
    DEN[i] = (TH1F*)den[i]->Clone();
    std::cout<<" NUM->Integral() "<<NUM[i]->Integral()
	     <<" DEN->Integral() "<<DEN[i]->Integral()
	     <<" FR = "<<NUM[i]->Integral()/DEN[i]->Integral()
	     <<std::endl;
    FR_value.push_back(NUM[i]->Integral()/DEN[i]->Integral());
    
    int j = (i) % 6;
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.7);
    gr->SetMarkerColor(icolor[j]);
    
    std::cout<<" j = "<<j<<std::endl;
    if(i==0 || j == 0 ){
      maximum = 0.0;
      ipad++;
      leg = new TLegend(0.55,0.45,0.95,0.90,NULL,"brNDC");
      leg->SetBorderSize(0);
      leg->SetLineColor(1);
      leg->SetLineStyle(1);
      leg->SetLineWidth(1);
      leg->SetFillColor(0);
      leg->SetFillStyle(0);
      leg->SetTextFont(42);
      leg->SetTextFont(60);
      leg->AddEntry(gr,legend_text[i],"P");
      
      TPad *p1;
      p1 = (TPad*)c->cd(ipad);
      gr->BayesDivide(NUM[i], DEN[i] );
      gr->Draw("AP");
      maximum = TMath::Max(gr->GetHistogram()->GetMaximum(),maximum);
      std::cout<<" gr->GetMaximum() = "<<gr->GetHistogram()->GetMaximum()<<std::endl;
      std::cout<<" maximum = "<<maximum<<std::endl;
    }//f(i==0 || j<0.00001 ){

    else {
      //NUM->Draw("same");
      gr->BayesDivide(NUM[i], DEN[i] );
      //num[i]->Draw("P,same");
      gr->Draw("P,same");
      maximum = TMath::Max(gr->GetMaximum(),maximum);
      leg->AddEntry(gr,legend_text[i],"P");
    }//else {
    if(i%5==0) leg->Draw();
    gr->GetXaxis()->SetLabelSize(0.06);
    gr->GetYaxis()->SetLabelSize(0.06);
    gr->GetXaxis()->SetTitle("p_{T}^{jet}");
    gr->GetYaxis()->SetTitle("fake rate");
    gr->GetXaxis()->SetTitleSize(0.08);
    gr->GetYaxis()->SetTitleSize(0.08);
    gr->GetXaxis()->SetTitleOffset(0.86);
    gr->GetYaxis()->SetTitleOffset(0.86);

    gr->SetMaximum(0.06);

  }
  return c  ;
}





TCanvas* DrawGraph::DrawGraph_4(TCanvas* c, std::vector<TObject*> num, std::vector<TObject*> den , std::vector<TString> legend_text, 
				double xmin, double xmax, double ymin, double ymax,  TString xtitle, bool plotratio){
  
  
  std::cout<<" legend size = "<<legend_text.size()<<std::endl;
  int ipad=0;
  TH1F* NUM[4];
  TH1F* DEN[4];
  std::vector<int> icolor;
  icolor.clear();
  icolor.push_back(3); 
  icolor.push_back(1); 
  icolor.push_back(2); 
  icolor.push_back(4); 
  
  TLegend *leg;
  leg = new TLegend(0.6,0.65,0.9,0.90,NULL,"brNDC");
  std::cout<<" num = "<<num.size()<<std::endl;
  TPad *p1;
  p1 = (TPad*)c->cd(1);
  p1->SetLogy(1);
  for(int i=0;i<int(num.size());i++){
    TGraphAsymmErrors *gr = new TGraphAsymmErrors();
    NUM[i] = (TH1F*)num[i]->Clone();
    DEN[i] = (TH1F*)den[i]->Clone();
    std::cout<<" NUM->Integral() "<<NUM[i]->Integral()
	     <<" DEN->Integral() "<<DEN[i]->Integral()
	     <<" FR = "<<NUM[i]->Integral()/DEN[i]->Integral()
	     <<std::endl;
    
    
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.8);
    gr->SetMarkerColor(icolor[i]);
    //std::cout<<" j = "<<j<<std::endl;
    
    leg->SetBorderSize(0);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    //leg->SetTextFont(5);
    leg->AddEntry(gr,legend_text[i],"P");
    
    gr->BayesDivide(NUM[i], DEN[i] );
    if(i==0){
      gr->Draw("AP");
    }
    else{
      gr->Draw("P,same");
    }
    
    gr->GetXaxis()->SetLabelSize(0.04);
    gr->GetYaxis()->SetLabelSize(0.04);
    gr->GetXaxis()->SetTitle(xtitle);
    gr->GetYaxis()->SetTitle("fake rate");
    gr->GetXaxis()->SetTitleSize(0.06);
    gr->GetYaxis()->SetTitleSize(0.06);
    gr->GetXaxis()->SetTitleOffset(0.86);
    gr->GetYaxis()->SetTitleOffset(0.86);

    
    
    //gr->SetMaximum(0.06);
    //gr->SetMinimum(0.035);

    gr->SetMaximum(ymax);
    gr->SetMinimum(ymin);
    //gr->GetXaxis()->SetRangeUser(0.0,3000.);
    gr->GetXaxis()->SetLimits(xmin,xmax);
    //    gr->GetYaxis()->SetLimits(0.0004,0.1);
  }
  
  leg->Draw();
  
  
  //---------------------------------------------------
  //------Ratio canvas---------------------------------
  //---------------------------------------------------
  if(plotratio){
    p1 = (TPad*)c->cd(2);
    
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    TH1F* NUM0;
    TH1F* NUM1;
    TH1F* NUM2;
    TH1F* NUM3;
    
    TH1F* DEN0;
    TH1F* DEN1;
    TH1F* DEN2;
    TH1F* DEN3;
    
    TH1F* Ratio1;
    TH1F* Ratio2;
    TH1F* Ratio3;
    
    NUM0 = (TH1F*)num[0]->Clone();
    NUM1 = (TH1F*)num[1]->Clone();
    NUM2 = (TH1F*)num[2]->Clone();
    NUM3 = (TH1F*)num[3]->Clone();
    
    DEN0 = (TH1F*)den[0]->Clone();
    DEN1 = (TH1F*)den[1]->Clone();
    DEN2 = (TH1F*)den[2]->Clone();
    DEN3 = (TH1F*)den[3]->Clone();
    
    Ratio1 = (TH1F*)num[1]->Clone();
    Ratio2 = (TH1F*)num[2]->Clone();
    Ratio3 = (TH1F*)num[3]->Clone();
    
    Ratio1->Clear();
    Ratio2->Clear();
    Ratio3->Clear();
    
    
    std::cout<<" num of bins = "<<NUM1->GetNbinsX()<<std::endl;
    for(int ibin=1; ibin<=NUM1->GetNbinsX();ibin++){
      double den,num1,num2,num3;
      if (DEN0->GetBinContent(ibin)>0)   den = NUM0->GetBinContent(ibin)/DEN0->GetBinContent(ibin);
      if (DEN1->GetBinContent(ibin)>1)   num1 = NUM1->GetBinContent(ibin)/DEN1->GetBinContent(ibin);
      if (DEN2->GetBinContent(ibin)>2)   num2 = NUM2->GetBinContent(ibin)/DEN2->GetBinContent(ibin);
      if (DEN3->GetBinContent(ibin)>3)   num3 = NUM3->GetBinContent(ibin)/DEN3->GetBinContent(ibin);
      
      if(den>0){
	Ratio1->SetBinContent(ibin,num1/den);
	Ratio2->SetBinContent(ibin,num2/den);
	Ratio3->SetBinContent(ibin,num3/den);
      }
    }
    Ratio1->SetMinimum(1.);
    Ratio1->SetMaximum(1.5);
    
    Ratio1->GetXaxis()->SetRangeUser(xmin,xmax);
    Ratio1->Draw("P");
    Ratio1->SetMarkerStyle(20);
    Ratio1->SetMarkerSize(0.8);
    Ratio1->SetMarkerColor(icolor[1]);
    Ratio1->GetXaxis()->SetTitle(xtitle);
    Ratio1->GetYaxis()->SetTitle("XYZ/default");
    Ratio1->GetXaxis()->SetLabelSize(0.04);
    Ratio1->GetYaxis()->SetLabelSize(0.04);
    Ratio1->GetXaxis()->SetTitleSize(0.06);
    Ratio1->GetYaxis()->SetTitleSize(0.06);
    Ratio1->GetXaxis()->SetTitleOffset(0.86);
    Ratio1->GetYaxis()->SetTitleOffset(1.2);
    
    
    Ratio2->GetXaxis()->SetRangeUser(xmin,xmax);
    Ratio2->SetMarkerStyle(20);
    Ratio2->Draw("Psame");
    Ratio2->SetMarkerStyle(20);
    Ratio2->SetMarkerSize(0.8);
    Ratio2->SetMarkerColor(icolor[2]);
    Ratio2->GetXaxis()->SetRangeUser(xmin,xmax);
    
    Ratio3->GetXaxis()->SetRangeUser(xmin,xmax);
    Ratio3->SetMarkerStyle(20);
    Ratio3->Draw("Psame");
    Ratio3->SetMarkerStyle(20);
    Ratio3->SetMarkerSize(0.8);
    Ratio3->SetMarkerColor(icolor[3]);
    Ratio3->GetXaxis()->SetRangeUser(xmin,xmax);
    std::cout<<" xmax = "<<xmax<<std::endl;
    return c  ;
  }
}
  




TCanvas* DrawGraph::DrawGraph_6(TCanvas* c, std::vector<double> FR_value){
  std::vector<int> icolor;
  icolor.clear();
  icolor.push_back(1); 
  icolor.push_back(2); 
  icolor.push_back(3); 
  icolor.push_back(4); 
  icolor.push_back(6); 
  icolor.push_back(7); 
  
  int ipad=0;
  TH1F* FR[6] ;
  for(int i=0;i<6;i++){
    char histpostname[100];
    sprintf(histpostname,"_%d",i);
    std::string histname(histpostname);
    FR[i] = new TH1F(("FR_"+histname).c_str(),"FR",10,1,11);
  }
  
  double maximum = 0.0;
  double minimum = 1.0;
  for(int ihist = 0; ihist < 6; ihist++){
    for (int i = 0; i<9; i++){
      FR[ihist]->SetBinContent(i+1,FR_value[(i*6)+ihist]);
      //FR[ihist]->SetMarkerColor(icolor[i]);
    }
    
    if(ihist==0)     FR[ihist]->Draw("P");
    else FR[ihist]->Draw("P,same");
    FR[ihist]->GetXaxis()->SetLabelSize(0.04);
    FR[ihist]->GetYaxis()->SetLabelSize(0.04);
    FR[ihist]->GetXaxis()->SetTitle("dz");
    FR[ihist]->GetYaxis()->SetTitle("fake rate");
    FR[ihist]->GetXaxis()->SetTitleSize(0.06);
    FR[ihist]->GetYaxis()->SetTitleSize(0.06);
    FR[ihist]->GetXaxis()->SetTitleOffset(0.86);
    FR[ihist]->GetYaxis()->SetTitleOffset(1.2);
    FR[ihist]->SetMarkerStyle(19+icolor[ihist]);
    FR[ihist]->SetMarkerSize(0.8);
    //FR[ihist]->SetMinimum(0.0);
    //FR[ihist]->SetMaximum(0.03);
    FR[ihist]->SetMarkerColor(icolor[ihist]);
    FR[ihist]->SetStats(0);
    maximum = TMath::Max(FR[ihist]->GetMaximum(),maximum);
    //minimum = TMath::Min(FR[ihist]->GetMinimum(),minimum);
  }
  
  for(int ihist=0; ihist<6;ihist++){
    FR[ihist]->SetMaximum(maximum+(0.5*maximum));
    FR[ihist]->SetMinimum(maximum-(maximum*0.8));
    FR[ihist]->GetXaxis()->SetBinLabel(1,"dz < 0.2 ");
    FR[ihist]->GetXaxis()->SetBinLabel(2,"dz < 0.3 ");
    FR[ihist]->GetXaxis()->SetBinLabel(3,"dz < 0.4 ");
    FR[ihist]->GetXaxis()->SetBinLabel(4,"dz < 0.5 ");
    FR[ihist]->GetXaxis()->SetBinLabel(5,"dz < 1.0 ");
    FR[ihist]->GetXaxis()->SetBinLabel(6,"dz < 2.0 ");
    FR[ihist]->GetXaxis()->SetBinLabel(7,"dz < 3.0 ");
    FR[ihist]->GetXaxis()->SetBinLabel(8,"dz < 4.0 ");
    FR[ihist]->GetXaxis()->SetBinLabel(9,"dz < 5.0 ");

 
  }

  TLegend *leg;
  leg = new TLegend(0.6,0.65,0.95,0.9,NULL,"brNDC");
    
  leg->SetBorderSize(0);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->AddEntry(FR[0],"ip < 0.3","P");
  leg->AddEntry(FR[1],"ip < 0.5","P");
  leg->AddEntry(FR[2],"ip < 0.10","P");
  leg->AddEntry(FR[3],"ip < 0.25","P");
  leg->AddEntry(FR[4],"ip < 0.50","P");
  leg->AddEntry(FR[5],"ip < 1.00","P");
  leg->Draw();
  return c  ;
}
  






TCanvas* DrawGraph::DrawGraph_6(TCanvas* c, std::vector<TObject*> histvector,
				std::vector<TString> legend,
				TString XAxis, TString YAxis,
				double rebinx, double rebiny,
				double xlow, double xhigh,
				double ylow, double yhigh,
				double logx, double logy){
  
  int ipad=0;
  TH2F* histograms[6];
  TLegend *leg;
  double maximum;

  
  TPad *p1;
  for(int i=0;i<int(histvector.size());i++){
    
    
    ipad++;
    histograms[i] = (TH2F*)histvector[i]->Clone();
    p1 = (TPad*)c->cd(ipad);
    p1->SetLogy(logy);
    p1->SetLogx(logx);
    p1->SetLeftMargin(0.2);
    p1->SetBottomMargin(0.2);
    p1->SetBottomMargin(0.2);

    histograms[i]->Draw("COLZ");
    histograms[i]->Rebin2D(rebinx,rebiny);

    histograms[i]->GetXaxis()->SetRangeUser(xlow,xhigh);
    histograms[i]->GetYaxis()->SetRange(ylow,yhigh);
    
    
    gStyle->SetOptStat(0);

    histograms[i]->GetXaxis()->SetLabelSize(0.06);
    histograms[i]->GetYaxis()->SetLabelSize(0.06);
    histograms[i]->GetXaxis()->SetTitle(XAxis);
    histograms[i]->GetYaxis()->SetTitle(YAxis);
    histograms[i]->GetXaxis()->SetTitleSize(0.08);
    histograms[i]->GetYaxis()->SetTitleSize(0.08);
    histograms[i]->GetXaxis()->SetTitleOffset(0.86);
    histograms[i]->GetYaxis()->SetTitleOffset(0.86);
    
    
    

    TPaveText *pt = new TPaveText(0.2577181,0.9,0.9580537,0.9947552,"brNDC");
    pt->SetBorderSize(0);
    pt->SetTextAlign(12);
    pt->SetFillStyle(0);
    pt->SetTextFont(42);
    pt->SetTextSize(0.07);
    TText *text = pt->AddText(0.01,0.5,legend[i]);
    pt->Draw();


    //TLatex latex;
    //latex.SetTextSize(0.05);
    //latex.DrawLatex(0.05,110.,legend[i]);
  }
return c  ;
}
  





TCanvas* DrawGraph::DrawHistograms_4(TCanvas* c, std::vector<TObject*> histvector,std::vector<TString> legend_text, TString legendheader){
  Double_t Ymin[4]={0.80,0.60,0.40,0.20};
  Double_t Ymax[4]={0.9999,0.80,0.60,0.40};
  
  std::vector<int> icolor;
  icolor.clear();
  icolor.push_back(1); 
  icolor.push_back(2); 
  icolor.push_back(4); 
  icolor.push_back(6); 
  
  
  TLegend *leg;
  leg = new TLegend(0.55,0.6,0.83,0.90,NULL,"brNDC");
  leg->SetHeader(legendheader);
  TH1F* HIST[4];

  TF1* myfit  = new TF1("myfit"," ( ([0]) * TMath::Exp (-0.5*(x-[1])*(x-[1]) / ([2]*[2]))) ",0.8,1.2);
  myfit->SetParName(0,"c0");
  myfit->SetParName(1,"c1");
  myfit->SetParName(1,"c3");
  
  myfit->SetParameter(0,  0.1);
  myfit->SetParameter(1, 1. );
  myfit->SetParameter(2, 0.1 );
  myfit->SetLineColor(4);
  
  
  
  for(int i=0;i<int(histvector.size());i++){
  //for(int i=0;i<int(4);i++){
    HIST[i] = (TH1F*)histvector[i]->Clone();
    HIST[i]->Rebin(5);
    if(i==0){
      HIST[i]->Draw();
    }
    else {
      HIST[i]->Draw("sames");
    }
    
    
    if(legend_text[i]=="dr < 0.12") { 
      HIST[i]->Fit("myfit","R");
      myfit->SetRange(0.0,3.0);
      myfit->Draw("same");
      HIST[i-1]->Draw("sames");
      HIST[i-2]->Draw("sames");
      
    }
    gPad->Update();
    TPaveStats *st = (TPaveStats*)HIST[i]->FindObject("stats");
    st->SetY1NDC(Ymin[i]);                                                                                                                  
    st->SetY2NDC(Ymax[i]);                                                                                                                  
    //st->SetX1NDC(Ymin[i]);                                                                                                                  
    //st->SetX2NDC(Ymax[i]);                                                                                                                  
    st->SetLineColor(icolor[i]);
    st->SetLineWidth(2);
    st->SetTextColor(icolor[i]);
    gStyle->SetOptStat(1110);
    
    //HIST[i]->SetMarkerStyle(20);
    //HIST[i]->SetMarkerSize(0.8);
    HIST[i]->SetLineColor(icolor[i]);
    HIST[i]->SetLineWidth(2);
    
    leg->SetBorderSize(0);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->AddEntry(HIST[i],legend_text[i],"l");
    
    HIST[i]->GetXaxis()->SetLabelSize(0.04);
    HIST[i]->GetYaxis()->SetLabelSize(0.04);
    HIST[i]->GetXaxis()->SetTitle("p_{T}^{reco}_{#tau}/p_{T}^{gen}_{#tau}");
    HIST[i]->GetYaxis()->SetTitle("au");
    HIST[i]->GetXaxis()->SetTitleSize(0.06);
    HIST[i]->GetYaxis()->SetTitleSize(0.06);
    HIST[i]->GetXaxis()->SetTitleOffset(0.86);
    HIST[i]->GetYaxis()->SetTitleOffset(0.86);
    
  }
  leg->Draw();
  
  
  
  return c  ;
}
  





TCanvas* DrawGraph::DrawTProfile(TCanvas* c, std::vector<TObject*> histvector,std::vector<TString> legend_text,TString legendheader){
  
  std::vector<int> icolor;
  icolor.clear();
  icolor.push_back(1); 
  icolor.push_back(2); 
  icolor.push_back(4); 
  //icolor.push_back(3); 
  icolor.push_back(6); 
  
  
  TLegend *leg;
  //leg = new TLegend(0.55,0.6,0.83,0.90,NULL,"brNDC");
  leg = new TLegend(0.6,0.7,0.9,0.90,NULL,"brNDC");
  leg->SetHeader(legendheader);
  TProfile* HIST[4];
  for(int i=0;i<int(histvector.size());i++){
    HIST[i] = (TProfile*)histvector[i]->Clone();
    HIST[i]->Rebin(4);
    if(i==0){
      HIST[i]->Draw();
    }
    else {
      HIST[i]->Draw("same");
    }
    std::cout<<" integral = "<<HIST[i]->Integral()<<std::endl;
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    HIST[i]->SetLineColor(icolor[i]);
    HIST[i]->SetLineWidth(2);
    
    leg->SetBorderSize(0);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->AddEntry(HIST[i],legend_text[i],"l");
    
    HIST[i]->GetXaxis()->SetLabelSize(0.04);
    HIST[i]->GetYaxis()->SetLabelSize(0.04);
    HIST[i]->GetXaxis()->SetTitle("p_{T}^{gen}_{#tau}");
    HIST[i]->GetYaxis()->SetTitle("<p_{T}^{XYZ}/p_{T}^{gen}_{#tau}>");
    HIST[i]->GetXaxis()->SetTitleSize(0.06);
    HIST[i]->GetYaxis()->SetTitleSize(0.06);
    HIST[i]->GetXaxis()->SetTitleOffset(0.86);
    HIST[i]->GetYaxis()->SetTitleOffset(0.86);
    HIST[i]->GetXaxis()->SetRangeUser(0.0,1600.);
    HIST[i]->SetMinimum(0.2);
    HIST[i]->SetMaximum(1.5);
  }
  leg->Draw();
  return c  ;
}
  



TCanvas* DrawGraph::DrawGraph_NN_pad_x_MM_histo(TCanvas* c, std::vector<TObject*> histovector, std::vector<TString> legend_text , int NN, int MM){
  const int array_size = NN * MM ; 
  
  int ipad=0;
  TH1F* HISTOGRAMS[array_size];
  std::cout<< "array_size = "<<array_size
	   <<" histovector = "<<histovector.size()
	   <<std::endl;
  std::vector<int> icolor;
  icolor.clear();
  icolor.push_back(1); 
  icolor.push_back(2); 
  icolor.push_back(3); 
  icolor.push_back(4); 
  icolor.push_back(5); 
  icolor.push_back(6); 
  
  TLegend *leg;
  double maximum;
  for(int i=0;i<int(histovector.size());i++){
    
    HISTOGRAMS[i] = (TH1F*)histovector[i]->Clone();
    std::cout<<" i = "<<i<<std::endl;
    int j = (i) % MM;
    HISTOGRAMS[i]->SetMarkerStyle(20);
    HISTOGRAMS[i]->SetMarkerSize(0.7);
    HISTOGRAMS[i]->SetMarkerColor(icolor[j]);
    
    std::cout<<" j = "<<j<<std::endl;
    if(i==0 || j == 0 ){
      maximum = 0.0;
      ipad++;
      leg = new TLegend(0.55,0.45,0.95,0.90,NULL,"brNDC");
      leg->SetBorderSize(0);
      leg->SetLineColor(1);
      leg->SetLineStyle(1);
      leg->SetLineWidth(1);
      leg->SetFillColor(0);
      leg->SetFillStyle(0);
      leg->SetTextFont(42);
      leg->SetTextFont(60);
      leg->AddEntry(HISTOGRAMS[i],legend_text[i],"P");
      
      TPad *p1;
      p1 = (TPad*)c->cd(ipad);
      HISTOGRAMS[i]->Draw();
      std::cout<<" integral = "<<      HISTOGRAMS[i]->Integral()<<std::endl;
      //maximum = TMath::Max(HISTOGRAMS[i]->GetMaximum(),maximum);
      
      
    }//f(i==0 || j<0.00001 ){

    else {
      HISTOGRAMS[i]->Draw("same");
      //maximum = TMath::Max(HISTOGRAMS[i]->GetMaximum(),maximum);
      leg->AddEntry(HISTOGRAMS[i],legend_text[i],"l");
    }//else {
    
    //if(i%MM==0) leg->Draw();
    HISTOGRAMS[i]->GetXaxis()->SetLabelSize(0.06);
    HISTOGRAMS[i]->GetYaxis()->SetLabelSize(0.06);
    HISTOGRAMS[i]->GetXaxis()->SetTitle("p_{T}^{jet}");
    HISTOGRAMS[i]->GetYaxis()->SetTitle("fake rate");
    HISTOGRAMS[i]->GetXaxis()->SetTitleSize(0.08);
    HISTOGRAMS[i]->GetYaxis()->SetTitleSize(0.08);
    HISTOGRAMS[i]->GetXaxis()->SetTitleOffset(0.86);
    HISTOGRAMS[i]->GetYaxis()->SetTitleOffset(0.86);
    HISTOGRAMS[i]->SetMaximum(0.06);
    std::cout<<"working"<<std::endl;
  }
  std::cout<< "working till  now"<<std::endl;
  return c  ;
}




TCanvas* DrawGraph::DrawHistograms_n(TCanvas* c, std::vector<TObject*> histvector,std::vector<TString> legend_text, TString legendheader){
  Double_t Ymin[5]={0.80,0.65,0.50,0.35,0.20};
  Double_t Ymax[5]={0.95,0.80,0.65,0.50,0.35};
  
  const int arraysize = histvector.size();
  
  std::vector<int> icolor;
  icolor.clear();
  icolor.push_back(1); 
  icolor.push_back(2); 
  icolor.push_back(4); 
  icolor.push_back(6); 
  icolor.push_back(8); 
  
  
  TLegend *leg;
  leg = new TLegend(0.15,0.85,0.30,0.999,NULL,"brNDC");
  leg->SetTextSize(0.023); 
  //leg->SetHeader(legendheader);
  TH1F* HIST[arraysize];
  
  for(int i=0;i<int(histvector.size());i++){
    HIST[i] = (TH1F*)histvector[i]->Clone();
    HIST[i]->Rebin(5);
    if(i==0){
      HIST[i]->Draw();
    }
    else {
      HIST[i]->Draw("sames");
    }
    
    
    gPad->Update();
    TPaveStats *st = (TPaveStats*)HIST[i]->FindObject("stats");
    std::cout<<" Ymin = "<<Ymin[i]
	     <<" Ymax = "<<Ymax[i]
	     <<std::endl;
    st->SetY1NDC(Ymin[i]);                                                                                                                  
    st->SetY2NDC(Ymax[i]);                                                                                                                  
    st->SetLineColor(icolor[i]);
    st->SetLineWidth(2);
    st->SetTextColor(icolor[i]);
    gStyle->SetOptStat(1110);
    
    HIST[i]->SetLineColor(icolor[i]);
    HIST[i]->SetLineWidth(2);
    
    leg->SetBorderSize(0);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->AddEntry(HIST[i],legend_text[i],"l");
    
    HIST[i]->GetXaxis()->SetLabelSize(0.04);
    HIST[i]->GetYaxis()->SetLabelSize(0.04);
    HIST[i]->GetXaxis()->SetTitle("p_{T}^{reco}_{#tau}/p_{T}^{gen}_{#tau}");
    HIST[i]->GetYaxis()->SetTitle("au");
    HIST[i]->GetXaxis()->SetTitleSize(0.06);
    HIST[i]->GetYaxis()->SetTitleSize(0.06);
    HIST[i]->GetXaxis()->SetTitleOffset(0.86);
    HIST[i]->GetYaxis()->SetTitleOffset(0.86);
    
  }
  leg->Draw();
  
  
  
  return c  ;
}
  



TCanvas* DrawGraph::DrawGraph_NN_pad_x_MM_histo(TCanvas* c, RKSaveObject ObjToSave){
  const int array_size = ObjToSave.xx * ObjToSave.yy ;  //xx for pad,  yy for histograms
  std::cout<<" array_size = "<<array_size<<std::endl;
  
 int ipad=0;
 TH1F* numerator[array_size];
 TH1F* denominator[array_size];
 
 std::cout<< "array_size = "<<array_size
	  <<" histovector = "<<ObjToSave.num.size()
	  <<std::endl;
  std::vector<int> icolor;
  icolor.clear();
  icolor.push_back(1); 
  icolor.push_back(2); 
  icolor.push_back(3); 
  icolor.push_back(4); 
  icolor.push_back(5); 
  icolor.push_back(6); 
  
  TLegend *leg;
  double maximum;
  for(int i=0;i<int(ObjToSave.num.size());i++){
    TGraphAsymmErrors *gr = new TGraphAsymmErrors();

    numerator[i] = (TH1F*)ObjToSave.num[i]->Clone();
    denominator[i] = (TH1F*)ObjToSave.den[i]->Clone();
    
    std::cout<<" num int = "<<numerator[i]->Integral()
	     <<", den int = "<<denominator[i]->Integral()
	     <<", efficiency = "<<(double)numerator[i]->Integral()/(double)denominator[i]->Integral()
	     <<std::endl;
    int j = (i) % (int)ObjToSave.yy ;
        
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.7);
    gr->SetMarkerColor(icolor[j]);
    

    if(i==0 || j == 0 ){
      maximum = 0.0;
      ipad++;
      //leg = new TLegend(0.45,0.15,0.9,0.50,NULL,"brNDC");
      leg = new TLegend(0.45,0.65,0.9,0.90,NULL,"brNDC");
      leg->SetBorderSize(0);
      leg->SetLineColor(1);
      leg->SetLineStyle(1);
      leg->SetLineWidth(1);
      leg->SetFillColor(0);
      leg->SetFillStyle(0);
      leg->SetTextFont(42);
      leg->SetTextFont(60);
      leg->AddEntry(gr,ObjToSave.legend_text[i],"P");
      //c->cd(1);
      TPad *p1;
      p1 = (TPad*)c->cd(ipad);
      p1->SetLeftMargin(0.2);
      p1->SetTopMargin(0.1);
      p1->SetBottomMargin(0.15);
      p1->SetLogy(1);
      
      gr->Divide(numerator[i], denominator[i] );
      gr->Draw("AP");
      
    }//f(i==0 || j==0 ){

    else {
      gr->Divide(numerator[i], denominator[i] );
      gr->Draw("same,P");
      leg->AddEntry(gr,ObjToSave.legend_text[i],"P");
      }//else {
    
    //if(i%MM==0) leg->Draw();
    gr->GetXaxis()->SetLabelSize(0.06);
    gr->GetYaxis()->SetLabelSize(0.06);
    gr->GetXaxis()->SetTitle("p_{T}^{jet}");
    gr->GetYaxis()->SetTitle("fake rate");
    gr->GetXaxis()->SetTitleSize(0.08);
    gr->GetYaxis()->SetTitleSize(0.08);
    gr->GetXaxis()->SetTitleOffset(0.8);
    gr->GetYaxis()->SetTitleOffset(0.89);
    //gr->SetMaximum(0.06);
    leg->Draw();
    std::cout<<" done till end"<<std::endl;
    //    c->SaveAs("efficiency.png");
  }
  std::cout<< "working till  now"<<std::endl;
  return c  ;
}








TCanvas* DrawGraph::DrawGraph_1_pad_x_MM_histo(TCanvas* c, RKSaveObject ObjToSave){
  const int array_size = ObjToSave.num.size() ;  //xx for pad,  yy for histograms
  std::cout<<" array_size = "<<array_size<<std::endl;
  
 int ipad=0;
 TH1F* numerator[array_size];
 TH1F* denominator[array_size];
 
 std::cout<< "array_size = "<<array_size
	  <<" histovector = "<<ObjToSave.num.size()
	  <<" den size  = "<<ObjToSave.den.size()
	  <<" leg size = "<<ObjToSave.legendplace.size()
	  <<std::endl;
  std::vector<int> icolor;
  icolor.clear();
  icolor.push_back(1); 
  icolor.push_back(2); 
  icolor.push_back(3); 
  icolor.push_back(4); 
  icolor.push_back(5); 
  icolor.push_back(6); 
  
  TLegend *leg;
  double maximum;
  TPad *p1;
  p1 = (TPad*)c->cd(1);
  
  for(int i=0;i<int(ObjToSave.num.size());i++){
    
    TGraphAsymmErrors *gr = new TGraphAsymmErrors();

    numerator[i] = (TH1F*)ObjToSave.num[i]->Clone();
    denominator[i] = (TH1F*)ObjToSave.den[i]->Clone();
    
    std::cout<<" num int = "<<numerator[i]->Integral()
	     <<", den int = "<<denominator[i]->Integral()
	     <<", efficiency = "<<(double)numerator[i]->Integral()/(double)denominator[i]->Integral()
	     <<std::endl;
    int j = (i) % (int)ObjToSave.yy ;
        
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.7);
    gr->SetMarkerColor(icolor[j]);
    gr->GetXaxis()->SetRangeUser(0.,1000.);

    if(i==0 ){
      maximum = 0.0;
      ipad++;
      //leg = new TLegend(0.45,0.15,0.9,0.50,NULL,"brNDC");
      //leg = new TLegend( 0.45,0.65,0.9,0.90 ,NULL,"brNDC");
      leg = new TLegend( ObjToSave.legendplace[0], ObjToSave.legendplace[1], ObjToSave.legendplace[2], ObjToSave.legendplace[3] ,NULL,"brNDC");
      std::cout<<" working 1"<<std::endl;
      leg->SetBorderSize(0);
      leg->SetLineColor(1);
      leg->SetLineStyle(1);
      leg->SetLineWidth(1);
      leg->SetFillColor(0);
      leg->SetFillStyle(0);
      leg->SetTextFont(42);
      leg->SetTextFont(60);
      leg->AddEntry(gr,ObjToSave.legend_text[i],"P");

      gr->Divide(numerator[i], denominator[i] );
      gr->Draw("AP");
      gr->GetXaxis()->SetRangeUser(0.,1000.);
    }//f(i==0 ){

    else {
      gr->Divide(numerator[i], denominator[i] );
      gr->Draw("same,P");
      gr->GetXaxis()->SetRangeUser(0.,1000.);
      leg->AddEntry(gr,ObjToSave.legend_text[i],"P");
      }//else {
    
    //if(i%MM==0) leg->Draw();
    gr->GetXaxis()->SetLabelSize(0.055);
    gr->GetYaxis()->SetLabelSize(0.055);
    gr->GetXaxis()->SetTitle(ObjToSave.xaxistitle);
    gr->GetYaxis()->SetTitle(ObjToSave.yaxistitle);
    gr->GetXaxis()->SetTitleSize(0.07);
    gr->GetYaxis()->SetTitleSize(0.07);
    gr->GetXaxis()->SetTitleOffset(0.85);
    gr->GetYaxis()->SetTitleOffset(0.79);
    //gr->SetMaximum(0.06);
    leg->Draw();
    std::cout<<" done till end"<<std::endl;
    //    c->SaveAs("efficiency.png");
  }
  std::cout<< "working till  now"<<std::endl;



  //---------------------------------------------------
  //------Ratio canvas---------------------------------
  //---------------------------------------------------
  p1 = (TPad*)c->cd(2);
  p1->SetBottomMargin(0.20);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  TH1F* NUM0;
  TH1F* NUM1;
  TH1F* NUM2;
  
  TH1F* DEN0;
  TH1F* DEN1;
  TH1F* DEN2;
    
  TH1F* Ratio1;
  TH1F* Ratio2;
  
  NUM0 = (TH1F*)ObjToSave.num[0]->Clone();
  NUM1 = (TH1F*)ObjToSave.num[1]->Clone();
  NUM2 = (TH1F*)ObjToSave.num[2]->Clone();
  		
  DEN0 = (TH1F*)ObjToSave.den[0]->Clone();
  DEN1 = (TH1F*)ObjToSave.den[1]->Clone();
  DEN2 = (TH1F*)ObjToSave.den[2]->Clone();
  
  Ratio1 = (TH1F*)ObjToSave.num[1]->Clone();
  Ratio2 = (TH1F*)ObjToSave.num[2]->Clone();
  
  //Ratio1->Clear();
  //Ratio2->Clear();
  
  
  std::cout<<" num of bins = "<<NUM1->GetNbinsX()<<std::endl;
  for(int ibin=1; ibin<=NUM1->GetNbinsX();ibin++){
    double den,num1,num2,num3;
    if (DEN0->GetBinContent(ibin)>0)   den = NUM0->GetBinContent(ibin)/DEN0->GetBinContent(ibin);
    if (DEN1->GetBinContent(ibin)>0)   num1 = NUM1->GetBinContent(ibin)/DEN1->GetBinContent(ibin);
    if (DEN2->GetBinContent(ibin)>0)   num2 = NUM2->GetBinContent(ibin)/DEN2->GetBinContent(ibin);
    
    if(den>0){
      Ratio1->SetBinContent(ibin,num1/den);
      Ratio2->SetBinContent(ibin,num2/den);
    }
  }
  Ratio1->SetMinimum(0.9);
  Ratio1->SetMaximum(1.1);

  Ratio1->GetXaxis()->SetRangeUser(0.,1010.);
  Ratio1->Draw("P");
  Ratio1->SetMarkerStyle(20);
  Ratio1->SetMarkerSize(0.8);
  Ratio1->SetMarkerColor(icolor[1]);
  
  Ratio1->GetXaxis()->SetLabelSize(0.055);
  Ratio1->GetYaxis()->SetLabelSize(0.055);
  Ratio1->GetXaxis()->SetTitle(ObjToSave.xaxistitle);
  Ratio1->GetYaxis()->SetTitle("Ratio");
  Ratio1->GetXaxis()->SetTitleSize(0.07);
  Ratio1->GetYaxis()->SetTitleSize(0.07);
  Ratio1->GetXaxis()->SetTitleOffset(0.85);
  Ratio1->GetYaxis()->SetTitleOffset(0.79);
  
    
  
  //  Ratio2->GetXaxis()->SetRangeUser(0.9,1.1);
  Ratio2->SetMarkerStyle(20);
  Ratio2->Draw("Psame");
  Ratio2->SetMarkerStyle(20);
  Ratio2->SetMarkerSize(0.8);
  Ratio2->SetMarkerColor(icolor[2]);
  //Ratio2->GetXaxis()->SetRangeUser(0,1100.);
  Ratio2->GetXaxis()->SetRangeUser(0.,1010.);
  return c  ;
}


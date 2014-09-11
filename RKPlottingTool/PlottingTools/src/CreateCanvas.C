
#define CreateCanvas_cxx
#include "../interface/CreateCanvas.h"
#include "TCanvas.h"

TCanvas* CreateCanvas::declareCanvas(){
  //TCanvas* c = new TCanvas("canvas", "canvas",286,86,600,600);
  TCanvas* c = new TCanvas("canvas", "canvas",0,0,1200,750);
  
  c->Range(-68.75,-7.5,856.25,42.5);
  c->SetLogy(1);
  
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  c->SetTickx(1);
  c->SetTicky(1);
  c->SetLeftMargin(0.10);
  c->SetRightMargin(0.10);
  c->SetTopMargin(0.30);
  c->SetBottomMargin(0.10);
  c->SetFrameFillStyle(0);
  c->SetFrameBorderMode(0);
  c->SetFrameFillStyle(0);
  c->SetFrameBorderMode(0);
  return c;
}



TCanvas* CreateCanvas::declareCanvasFor1Plot(){
  TCanvas* c = new TCanvas("canvas", "canvas",286,86,600,600);
  //TCanvas* c = new TCanvas("canvas", "canvas",0,0,1300,700);
  c->Range(-68.75,-7.5,856.25,42.5);
  c->SetLogy(1);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  c->SetTickx(1);
  c->SetTicky(1);
  c->SetLeftMargin(0.15);
  c->SetRightMargin(0.10);
  c->SetTopMargin(0.10);
  c->SetBottomMargin(0.20);
  c->SetFrameFillStyle(0);
  c->SetFrameBorderMode(0);
  c->SetFrameFillStyle(0);
  c->SetFrameBorderMode(0);
  return c;
}



TCanvas* CreateCanvas::declareCanvasFor1Plot(int i){
  TCanvas* c = new TCanvas("canvas", "canvas",286,86,700,800);
  //TCanvas* c = new TCanvas("canvas", "canvas",0,0,1300,700);
  
  c->Range(-68.75,-7.5,856.25,42.5);
  c->SetLogy(1);
  
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  c->SetTickx(1);
  c->SetTicky(1);
  c->SetLeftMargin(0.15);
  c->SetRightMargin(0.10);
  c->SetTopMargin(0.15);
  c->SetBottomMargin(0.20);
  c->SetFrameFillStyle(0);
  c->SetFrameBorderMode(0);
  c->SetFrameFillStyle(0);
  c->SetFrameBorderMode(0);
  return c;
}

TCanvas* CreateCanvas::declareCanvasFor2Plot(){
  TCanvas* c = new TCanvas("canvas", "canvas",286,86,600,600);
  //TCanvas* c = new TCanvas("canvas", "canvas",0,0,1300,700);
  
  c->Range(-68.75,-7.5,856.25,42.5);
  c->SetLogy(1);
  
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  c->SetTickx(1);
  c->SetTicky(1);
  c->SetLeftMargin(0.15);
  c->SetRightMargin(0.10);
  c->SetTopMargin(0.10);
  c->SetBottomMargin(0.20);
  c->SetFrameFillStyle(0);
  c->SetFrameBorderMode(0);
  c->SetFrameFillStyle(0);
  c->SetFrameBorderMode(0);
  return c;
}


TCanvas* CreateCanvas::DivideManually(TCanvas* c, int row, int column){
  c->Divide(row,column);
  return c;
}


void CreateCanvas::SaveCanvas(TCanvas* c, TString canvasname, TString directoryname){
  c->SaveAs("Plots/"+directoryname+"/png/"+canvasname+".png");
  c->SaveAs("Plots/"+directoryname+"/pdf/"+canvasname+".pdf");
  c->SaveAs("Plots/"+directoryname+"/root/"+canvasname+".root");
  c->SaveAs("Plots/"+directoryname+"/macro/"+canvasname+".C");
  
}

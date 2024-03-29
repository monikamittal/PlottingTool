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
#include "../interface/ParseHistoInfo.h"


std::vector<HistoInfo> ParseHistoInfo::ExtractHistoInformation(char* histcardname){
  histInfoVector.clear();  
  std::ifstream histfile(histcardname);
  std::string Line;
  
  while( std::getline( histfile, Line ) )    {  // output of readline is stored in "line"
    std::istringstream issHist( Line );
    std::string resultHist;
    if( std::getline( issHist, resultHist , ' ') ){
      //Extract Histname
      if(resultHist=="Histname"){
	std::string histname;
	while( std::getline( issHist, histname, ' ') ){
	  Histname.push_back(histname);
	  if(debug) std::cout<<" Histname = "<<"  "<<Histname.size()<<std::endl;
      	}
	if(Histname.size() == 0) {
      	  std::cout<<" Code can not run without input files .. fix the plotting card"<<std::endl;
      	  return histInfoVector;
      	}
      }//if(resultHist=="Histname"){
  
      

      //Extract XAxis Title
      if(resultHist=="Xaxis"){
	std::string xaxis;
	while( std::getline( issHist, xaxis) ){
	  Xaxis = xaxis;
	}
	if(true) std::cout<<" xaxis = "<<Xaxis<<"  "<<Xaxis.size()<<std::endl;
      	if(Xaxis.size() == 0) {
      	  std::cout<<" You havenot provided X Axis title, if you want this then fix the plotting card"<<std::endl;
      	}
      }//if(resultHist=="Histname"){
  
      
      //Extract YAxis Title
      if(resultHist=="Yaxis"){
	std::string yaxis;
	while( std::getline( issHist, yaxis) ){
	  Yaxis = yaxis;
	}
	if(debug) std::cout<<" yaxis = "<<Yaxis<<"  "<<Yaxis.size()<<std::endl;
      	if(Yaxis.size() == 0) {
      	  std::cout<<" You havenot provided Y Axis title, if you want this then fix the plotting card"<<std::endl;
      	}
      }//if(resultHist=="Histname"){
  

      
      //Extract Rebin value
      if(resultHist=="Rebin"){
	std::string Rebin;
	while( std::getline( issHist, Rebin) ){
	  rebin = atoi(Rebin.c_str());
	}
	if(debug)	std::cout<<" rebin = "<<rebin<<std::endl;
	if(Rebin.size() == 0) {
      	  std::cout<<" You havenot provided Y Axis title, if you want this then fix the plotting card"<<std::endl;
      	}
      }//if(resultHist=="Histname"){
  
      
      //Extract XAxis Range
      if(resultHist=="XRange"){
	std::string XRange;
	double xrange;
	while( std::getline( issHist, XRange, ' ') ){
	  xrange = atof(XRange.c_str());
	  xRange.push_back(xrange);
      	  if(debug) std::cout<<" xrange = "<<xrange<<std::endl;
	  std::cout<<" xrange = "<<xrange<<std::endl;
	  if(xRange.size() == 0) {
	    std::cout<<" You have not provided te XRange, if you want to use this then fix the plotting card"<<std::endl;
	  }
	}
      }//if(resultHist=="XRange"){



      //Extract YAxis Range
      if(resultHist=="YRange"){
	std::string YRange;
	int yrange;
	while( std::getline( issHist, YRange, ' ') ){
	  yrange = atoi(YRange.c_str());
	  yRange.push_back(yrange);
      	  if(debug) std::cout<<" yrange = "<<yrange<<std::endl;
	  if(yRange.size() == 0) {
	    std::cout<<" You have not provided te YRange, if you want to use this then fix the plotting card"<<std::endl;
	  }
	}
      }//if(resultHist=="YRange"){


      
      //Extract Legend Pos
      if(resultHist=="legendplace"){
	std::string legendplace;
	float  position;
	while( std::getline( issHist, legendplace, ' ') ){
	  std::cout<<" legendplace "<<legendplace<<std::endl;
	  position = atof(legendplace.c_str());
	  LegendPlace.push_back(position);
      	  if(true) std::cout<<" leg pos = "<<position<<std::endl;
	  if(yRange.size() == 0) {
	    std::cout<<" You have not provided the legend place, fix the plotting card"<<std::endl;
	  }
	}
      }//if(resultHist=="legendplace"){


      //Extract filename
      if(resultHist=="PlotName"){
	std::string PLOTNAME;
	while( std::getline( issHist, PLOTNAME, ' ') ){
      	  if(debug) std::cout<<" PLOTNAME = "<<PLOTNAME<<std::endl;
	  PlotName = PLOTNAME;
	  if(PLOTNAME.size() == 0) {
	    std::cout<<" You have not provided the filename to save fix the plotting card"<<std::endl;
	  }
	}
      }//if(resultHist=="YRange"){
      
      
      //Extract Legend                                                                                                                                               
      if(resultHist=="legend"){
	std::string legend;
	while( std::getline( issHist, legend, ' ') ){
	  if(debug) std::cout<<" legend = "<<legend<<std::endl;
	  legends.push_back(legend);
	}
	if(legends.size() == 0 && Histname.size() > 1 )  {  // use legend from Rootfile info if there are >1 files
	  std::cout<<" Code can not run without legends, in case you dont need legend then put some dummy legend .. fix the plotting card"<<std::endl;
	  return histInfoVector;
	}
      }//if(result=="legend"){                                                                                                                                       

      
      //Extract color                                                                                                                                                
      if(resultHist=="color"){
	std::string color;
	int Color;
	while( std::getline( issHist, color, ' ') ){
	  Color = atoi(color.c_str());
	  if(debug) std::cout<<" color = "<<Color<<std::endl;
	  colors.push_back(Color);
	}
	if(colors.size() == 0 && Histname.size() > 1){  // use colors from Rootfile info if there are >1 files
	  std::cout<<" Please provide color Informations .. fix the plotting card"<<std::endl;
	  return histInfoVector;
	}
      }//if(result=="color"){                                                                                                                                        


      
      //Check  End of One Histo
      if(resultHist=="endofhisto"){
	histoinfo.histname = Histname;
	histoinfo.XRange = xRange;
	histoinfo.YRange = yRange;
	histoinfo.XAXISTitle = Xaxis;
	histoinfo.YAXISTitle = Yaxis;
	histoinfo.ReBin = rebin;
	histoinfo.plotname = PlotName ;
	histoinfo.legpos = LegendPlace;
	histoinfo.legend = legends;
        histoinfo.color = colors;
	//std::cout<<" histoinfo.YAXISTitle  = "<<histoinfo.YAXISTitle<<std::endl;
	//Push each histo info in the vector
	histInfoVector.push_back(histoinfo);
	Initialize();
	
	// call function to fill or plot and save histograms
      }//if(resultHist=="###"){
  
      
    }//if( std::getline( issHist, resultHist , ' ') ){
  }//while( std::getline( histfile, Line ) )    {

  return histInfoVector;
}


void ParseHistoInfo::Initialize(){
  debug = false;
  Histname.clear();
  legends.clear();
  colors.clear();
  PlotName.clear();
  LegendPlace.clear();
  xRange.clear();
  yRange.clear();
  rebin = 1;
}

void ParseHistoInfo::Print(std::vector<HistoInfo> hinfovector){
  for(int i=0; i<(int) hinfovector.size(); i++){
    std::cout<<" xaxis = "<<hinfovector[i].XAXISTitle
	     <<" yaxis = "<<hinfovector[i].YAXISTitle
	     <<std::endl;
  }
}

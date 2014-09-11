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

//--- My Files
#include "PlottingTools/interface/CreateCanvas.h"
#include "PlottingTools/interface/DrawGraph.h"
#include "PlottingTools/interface/FileInfo.h"
#include "PlottingTools/interface/HistoInfo.h"
#include "PlottingTools/interface/ParseHistoInfo.h"
#include "PlottingTools/interface/ParseRootFileInfo.h"
#include "PlottingTools/interface/PlotPOverlappingHisto.h"
#include "PlottingTools/interface/PlotStackedHisto.h"
int main(int argc, char *argv[]){
  bool debug = false;

  if(argc<3) {
    std::cout<<" please provide sufficient inputs"<<std::endl;
    return 0;
  }
  
  char* fileinfo  = argv[1];
  char* histoinfo = argv[2];
  
  // Fill file information here 
  std::vector<FileInfo> finfoVector ;
  finfoVector.clear();
  ParseRootFileInfo pfile;
  pfile.Initialize();
  char* plottingcardname = fileinfo;
  finfoVector = pfile.ExtractFileInformation(plottingcardname);
  //finfoVector = pfile.ExtractFileInformation(fileinfo);
  if(debug) pfile.Print(finfoVector);

  // Fill histogrm info here
  std::vector<HistoInfo> hinfovector;
  hinfovector.clear();
  ParseHistoInfo phisto;
  char* histocardname = histoinfo;
  phisto.Initialize();
  hinfovector = phisto.ExtractHistoInformation(histocardname);
  if(debug) phisto.Print(hinfovector);
  
  
 
  
  //
  //PlotPOverlappingHisto POverlap;
  //POverlap.DrawHistograms(finfoVector, hinfovector);
  
  PlotStackedHisto PStack;
  PStack.DrawHistograms(finfoVector, hinfovector);
  
  return 0;
}

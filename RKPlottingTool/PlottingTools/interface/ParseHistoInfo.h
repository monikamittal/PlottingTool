#ifndef ParseHistoInfo_H
#define ParseHistoInfo_H

#include <iostream>
#include "CommonFilesToInclude.h"
#include "HistoInfo.h"

class ParseHistoInfo {
 public :
  
  ParseHistoInfo(){};
  HistoInfo histoinfo;
  std::vector<HistoInfo> histInfoVector;
  bool debug ;
  std::vector<std::string> Histname;
  std::vector<double> xRange;
  std::vector<double> yRange;
  std::vector<double> LegendPlace;
  std::vector<std::string> legends;
  std::vector<int> colors;
  std::string Xaxis, Yaxis;
  std::string PlotName ;
  int rebin ;


  void Initialize();
  void Print(std::vector<HistoInfo> hinfovector);
  std::vector<HistoInfo> ExtractHistoInformation(char* plottingcardname);
  ~ParseHistoInfo(){};
  
};
#endif 

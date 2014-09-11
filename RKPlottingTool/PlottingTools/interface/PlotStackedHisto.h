#ifndef PlotStackedHisto_H
#define PlotStackedHisto_H

#include <iostream>
#include "CommonFilesToInclude.h"
#include "FileInfo.h"
#include "FileInfo.h"
#include "HistoInfo.h"

class PlotStackedHisto {
 public :
  
  PlotStackedHisto(){};
  void DrawHistograms(std::vector<FileInfo> fobj, std::vector<HistoInfo> hobjVec);
  ~PlotStackedHisto(){};
  
};
#endif 

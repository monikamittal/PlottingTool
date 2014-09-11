#ifndef PlotPOverlappingHisto_H
#define PlotPOverlappingHisto_H

#include <iostream>
#include "CommonFilesToInclude.h"
#include "FileInfo.h"
#include "FileInfo.h"
#include "HistoInfo.h"

class PlotPOverlappingHisto {
 public :
  
  PlotPOverlappingHisto(){};
  void DrawHistograms(std::vector<FileInfo> fobj, std::vector<HistoInfo> hobjVec);
  ~PlotPOverlappingHisto(){};
  
};
#endif 

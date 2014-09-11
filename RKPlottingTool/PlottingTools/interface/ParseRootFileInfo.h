#ifndef ParseRootFileInfo_H
#define ParseRootFileInfo_H

#include <iostream>
#include "CommonFilesToInclude.h"
#include "FileInfo.h"

class ParseRootFileInfo {
 public :
  
  ParseRootFileInfo(){};
  bool debug ;
  std::vector<std::string> options;
  int nFile;
  int Normalize;
  int nHist;
  double Lumi;
  
  std::vector<int> normalizeunity;
  std::vector<std::string> files;
  std::vector<std::string> legends;
  std::vector<int> colors;
  std::vector<double> xsections;
  std::string DirnameInRootfile;
  std::string PrefixToSave;
  TH1F* TotalEvents;
  FileInfo fileinfo ;
  void Initialize();
  void Print(std::vector<FileInfo> finfo);
  std::vector<FileInfo> fileinfoVector;
  std::vector<FileInfo> ExtractFileInformation(char* plottingcardname);
  ~ParseRootFileInfo(){};
  
};
#endif 

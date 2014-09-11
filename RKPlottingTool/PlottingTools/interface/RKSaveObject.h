#ifndef RKSaveObject_H
#define RKSaveObject_H

#include <iostream>
#include "CommonFilesToInclude.h"

class RKSaveObject {
 public :
  RKSaveObject(){
    num.clear();
    den.clear();
  }
  //modeofplot = 0; 
  /*    modeofplot 
	0 : efficiency 
	1 : fake rate
	2 : other
  */
  
  
  std::vector<TObject*> num;
  std::vector<TObject*> den;
  std::vector<TString> legend_text;
  TString Xaxistitle;
  TString Yaxistitle;
  double xmax;
  double ymax;
  TCanvas* c;
  double xx;
  double yy;
  TString pngfilename;
  TString  dirname;
  std::vector<double> legendplace;
  int modeofplot;
  TString xaxistitle;
  TString yaxistitle;
  void Initialize(int mode);
 

  ~RKSaveObject(){};
  
};
#endif 

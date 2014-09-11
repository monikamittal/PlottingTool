#ifndef DrawGraph_H
#define DrawGraph_H

#include "CommonFilesToInclude.h"
#include "RKSaveObject.h"
#include "CreateCanvas.h"
class DrawGraph {
 public :
  DrawGraph(){};
  TCanvas* DrawGraph_9pad_x_6graph(TCanvas* c, std::vector<TObject*> num, std::vector<TObject*> den, std::vector<TString> legend_text, std::vector<double>& FR_values);
  TCanvas* DrawGraph_4(TCanvas* c, std::vector<TObject*> num, std::vector<TObject*> den , std::vector<TString> legend_text,double xmin,double xmax,double ymin, 
		       double ymax, TString xtitle, bool plotratio);
  TCanvas* DrawGraph_6(TCanvas* c, std::vector<double> FR_value);
  TCanvas* DrawGraph_6(TCanvas* c, std::vector<TObject*> histvector,
		       std::vector<TString> legend,
		       TString XAxis, TString YAxis,
		       double rebinx, double rebiny,
		       double xlow, double xhigh,
		       double ylow, double yhigh,
		       double logx=0, double logy=1);

  TCanvas* DrawHistograms_4(TCanvas* c, std::vector<TObject*> histvector,std::vector<TString> legend_text, TString legendheader);
  TCanvas* DrawTProfile(TCanvas* c, std::vector<TObject*> histvector,std::vector<TString> legend_text, TString legendheader);
  TCanvas* DrawGraph_NN_pad_x_MM_histo(TCanvas* c, std::vector<TObject*> histovector, std::vector<TString> legend_text , int NN, int MM);
  TCanvas* DrawHistograms_n(TCanvas* c, std::vector<TObject*> histvector,std::vector<TString> legend_text, TString legendheader);
  TCanvas* DrawGraph_NN_pad_x_MM_histo(TCanvas* c ,RKSaveObject ObjToSave);
  TCanvas* DrawGraph_1_pad_x_MM_histo(TCanvas* c ,RKSaveObject ObjToSave);

  ~DrawGraph(){};
};
#endif

#ifndef CREATECANVAS_H
#define CREATECANVAS_H

#include "TCanvas.h"
class CreateCanvas {
 public :
  CreateCanvas(){};
  TCanvas* declareCanvas();
  TCanvas* DivideManually(TCanvas* c, int row, int column);
  TCanvas* declareCanvasFor1Plot();
  TCanvas* declareCanvasFor1Plot(int i);
  TCanvas* declareCanvasFor2Plot();
  void SaveCanvas(TCanvas* c, TString canvasname, TString directoryname);
  ~CreateCanvas(){};
};
#endif

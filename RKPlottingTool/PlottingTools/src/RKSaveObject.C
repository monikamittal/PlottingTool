#define RKSaveObject_cxx
#include "../interface/RKSaveObject.h"
#include "TCanvas.h"

void  RKSaveObject::Initialize(int modeofplot){
  if(modeofplot==0){
    xaxistitle = "p_{T}^{#tau}";
    yaxistitle = "efficiency";
    legendplace.push_back(0.45);
    legendplace.push_back(0.19);
    legendplace.push_back(0.9);
    legendplace.push_back(0.50);
  }
  if(modeofplot==1){
    xaxistitle = "p_{T}^{jet}";
    yaxistitle = "fake rate";
    legendplace.push_back(0.45);
    legendplace.push_back(0.65);
    legendplace.push_back(0.90);
    legendplace.push_back(0.90);
  }
  if(modeofplot>1){
  legendplace.push_back(0.45);
  legendplace.push_back(0.15);
  legendplace.push_back(0.9);
  legendplace.push_back(0.50);
  }
}


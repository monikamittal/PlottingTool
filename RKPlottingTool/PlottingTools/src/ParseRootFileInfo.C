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
#include "../interface/ParseRootFileInfo.h"


std::vector<FileInfo> ParseRootFileInfo::ExtractFileInformation(char* plottingcardname){
  fileinfoVector.clear();
  std::ifstream infile(plottingcardname);
  //std::ifstream infile("plottingcard.txt");
  std::string line;
  while( std::getline( infile, line ) )    {  // output of readline is stored in "line"
    std::istringstream iss( line );
    std::string result;
    if( std::getline( iss, result , ' ') ){
            
      // Extract Option
      if(result== "option"){
	std::string option;
	while( std::getline( iss, option, ' ') ){
	  if(debug) std::cout<<" option = "<<option<<std::endl;
	  std::cout<<" You are running code with "<<option<<" "<<result<<std::endl;
	  options.push_back(option);
	}
	
	if(options.size()>1) { // you can run code with one option at one time
	  std::cout<<" Code can not run for "<<options.size()<<" options .. Fix the plotting card"<<std::endl;
	  return fileinfoVector;
	}
      }//if(result== "option"){



      
      // Extract DirName
      if(result== "dirnameinrootfile"){
	std::string dirnameinrootfile;
	while( std::getline( iss, dirnameinrootfile, ' ') ){
	  if(debug) std::cout<<" dirnameinrootfile = "<<dirnameinrootfile<<std::endl;
	  DirnameInRootfile = dirnameinrootfile;
	}
	
	if(dirnameinrootfile.size()==0) {
	  std::cout<<" You have not provided any dirname "<<std::endl;
	}
      }//if(result== "dirnameinrootfile"){


      
      // Extract PrefixToSave
      if(result== "rootfileprefix"){
	std::string prefixforrootfile;
	while( std::getline( iss, prefixforrootfile, ' ') ){
	  if(debug) std::cout<<" prefixforrootfile = "<<prefixforrootfile<<std::endl;
	  PrefixToSave = prefixforrootfile;
	}
	
	if(PrefixToSave.size()==0) {
	  std::cout<<" You have not provided any PrefixToSave "<<std::endl;
	}
      }//if(result== "dirnameinrootfile"){

      

      //Extract nfiles
      if(result=="nfile"){
	std::string nfile;
	while( std::getline( iss, nfile, ' ') ){
	  nFile = atoi(nfile.c_str());
	}
	if(debug) std::cout<<" nfile = "<<nFile<<std::endl;
	if(nFile == 0) {
	  std::cout<<" Code can not run without input files .. fix the plotting card"<<std::endl;
	  return fileinfoVector;
	}
      }//if(result=="nfile"){



      //Extract Normalize Unity option
      if(result=="normalizedUnity"){
	std::string normalize;
	while( std::getline( iss, normalize, ' ') ){
	  Normalize = atoi(normalize.c_str());
	}
	if(debug) std::cout<<" Normalize = "<<Normalize<<std::endl;
	if(nFile == 0) {
	  std::cout<<" Code will not normalize plots to unity."<<std::endl;
	}
      }//if(result=="normalizedUnity"){

      
      
      //Extract files
      if(result=="files"){
	std::string file;
      	while( std::getline( iss, file, ' ') ){
      	  if(debug) std::cout<<" file = "<<file<<std::endl;
      	  files.push_back(file);
      	}
      	if(files.size() == 0) {
      	  std::cout<<" Code can not run without input files .. fix the plotting card"<<std::endl;
      	  return fileinfoVector;
      	}
      }//if(result=="nfile"){
      
      // skip reading following information if there is only file
      if(nFile>1){
	//Extract nhist
	if(result=="nhist"){
	  std::string nhist;
	  while( std::getline( iss, nhist, ' ') ){
	    nHist = atoi(nhist.c_str());
	    if(debug) std::cout<<" nHist = "<<nHist<<std::endl;
	  }
	  if(nHist == 0) {
	    std::cout<<" Please provide number hist information ... fix the plotting card"<<std::endl;
	    return fileinfoVector;
	  }
	}//if(result=="nfile"){
	
	
	
	//Extract Legend
	if(result=="legend"){
	  std::string legend;
	  while( std::getline( iss, legend, ';') ){
	    if(debug) std::cout<<" legend = "<<legend<<std::endl;
	    legends.push_back(legend);
	  }
	  if(legends.size() == 0) {
	    std::cout<<" Code can not run without legends, in case you dont need legend then put some dummy legend .. fix the plotting card"<<std::endl;
	    return fileinfoVector;
	  }
	}//if(result=="legend"){
	
	
	//Extract color
	if(result=="color"){
	  std::string color;
	  int Color;
	  while( std::getline( iss, color, ' ') ){
	    Color = atoi(color.c_str());
	    if(debug) std::cout<<" color = "<<Color<<std::endl;
	    colors.push_back(Color);
	  }
	  if(colors.size() == 0) {
	    std::cout<<" Please provide color Informations .. fix the plotting card"<<std::endl;
	    return fileinfoVector;
	  }
	}//if(result=="color"){
	
      }//if(nFile>1){
      
      if(result=="endoffile"){
	fileinfo.filename = files;
	fileinfo.legend = legends;
	fileinfo.color = colors;
	fileinfo.xsection = xsections;
	fileinfo.option = options[0];
	fileinfo.nhist = nHist;
	fileinfo.nfile = nFile;
	fileinfo.lumi = Lumi;
	fileinfo.dirname = DirnameInRootfile;
	fileinfo.normalizeUnity = Normalize;
	//fileinfo.totalevent = (TH1F*) TotalEvents->Clone("totalevents");
	fileinfo.prefixname = PrefixToSave;
	fileinfoVector.push_back(fileinfo);
	Initialize();
      }//if(result=="endofhisto"){
      
    }//if( std::getline( iss, result , ' ') ){
    
  }//while( std::getline( file, line ) )    {

  
  if(debug) std::cout<<" Extracted Information of Rootfiles "<<std::endl;
  if(debug) std::cout<<" size = "<<options.size()
		     <<" option = "<<options[0]
		     <<std::endl;
  return fileinfoVector;
}


void ParseRootFileInfo::Initialize(){
  debug = true;
  Lumi = 0.0;
  options.clear();
  files.clear();
  normalizeunity.clear();
  legends.clear();
  colors.clear();
  DirnameInRootfile="";
}


void ParseRootFileInfo::Print(std::vector<FileInfo> finfovector){
  for(int i=0; i<(int) finfovector.size(); i++){
    std::cout<<" option = "<<finfovector[i].option
	     <<" legend= "<<finfovector[i].legend[0]
	     <<" color= "<<finfovector[i].color[0]
	     <<" nhist= "<<finfovector[i].nhist
	     <<" nfile= "<<finfovector[i].nfile
	     <<std::endl;
  }
  
}

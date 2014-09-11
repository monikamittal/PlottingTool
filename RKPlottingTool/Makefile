OBJ=CreateCanvas.o DrawGraph.o RKSaveObject.o ParseRootFileInfo.o ParseHistoInfo.o PlotPOverlappingHisto.o PlotStackedHisto.o
CFLAGS = -Wall -Wno-deprecated -I
all: saveplots.exe
saveplots.exe : $(OBJ)
	g++ `root-config --cflags` Main.C $(OBJ) -o saveplots.exe `root-config --libs`
CreateCanvas.o : PlottingTools/src/CreateCanvas.C 
	g++ `root-config --cflags` PlottingTools/src/CreateCanvas.C -c 

DrawGraph.o : PlottingTools/src/DrawGraph.C
	g++ `root-config --cflags` PlottingTools/src/DrawGraph.C -c 

RKSaveObject.o : PlottingTools/src/RKSaveObject.C
	g++ `root-config --cflags` PlottingTools/src/RKSaveObject.C -c 

ParseRootFileInfo.o : PlottingTools/src/ParseRootFileInfo.C
	g++ `root-config --cflags` PlottingTools/src/ParseRootFileInfo.C -c 


ParseHistoInfo.o : PlottingTools/src/ParseHistoInfo.C
	g++ `root-config --cflags` PlottingTools/src/ParseHistoInfo.C -c 

PlotPOverlappingHisto.o : PlottingTools/src/PlotPOverlappingHisto.C
	g++ `root-config --cflags` PlottingTools/src/PlotPOverlappingHisto.C -c 

PlotStackedHisto.o : PlottingTools/src/PlotStackedHisto.C
	g++ `root-config --cflags` PlottingTools/src/PlotStackedHisto.C -c 

clean: 
	rm -rf CreateCanvas.o DrawGraph.o RKSaveObject.o saveplots.exe	ParseRootFileInfo.o ParseHistoInfo.o PlotPOverlappingHisto.o PlotStackedHisto.o
	echo Cleaning done	

c: 
	rm -rf *.o
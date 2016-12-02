# Makefile depends on ROOT.
# Define the variables ROOTSYS in your environment.


ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --glibs)

CXXFLAGS      = -Wall -fPIC  -g -W
CXXFLAGS     += $(ROOTCFLAGS)

LIBS          = $(ROOTLIBS) -lGenVector
CC            = g++

DIR 	      := $(PWD)/../
KERNELDIR     := $(DIR)../TbKernel/
ALGODIR       := $(DIR)../TbAlgorithms/

OBJS = $(KERNELDIR)bin/TbGeometrySvc.o $(ALGODIR)bin/TbDecoder.o $(ALGODIR)bin/TbClustering.o $(ALGODIR)bin/TbToyData.o


all:			$(DIR)bin/bach

$(DIR)bin/bach: 	$(OBJS)
			${CC} $(CXXFLAGS) $(LIBS) -o $(DIR)bin/bach $(DIR)src/bach.cpp $(OBJS)

$(KERNELDIR)bin/TbGeometrySvc.o:	$(KERNELDIR)src/TbGeometrySvc.cpp $(KERNELDIR)src/TbGeometrySvc.h
			${CC} $(CXXFLAGS) -c $(KERNELDIR)src/TbGeometrySvc.cpp -o $(KERNELDIR)bin/TbGeometrySvc.o

$(ALGODIR)bin/TbDecoder.o:	$(ALGODIR)src/TbDecoder.cpp $(ALGODIR)src/TbDecoder.h
			${CC} $(CXXFLAGS) -c $(ALGODIR)src/TbDecoder.cpp -o $(ALGODIR)bin/TbDecoder.o
$(ALGODIR)bin/TbClustering.o:	$(ALGODIR)src/TbClustering.cpp $(ALGODIR)src/TbClustering.h
			${CC} $(CXXFLAGS) -c $(ALGODIR)src/TbClustering.cpp -o $(ALGODIR)bin/TbClustering.o
$(ALGODIR)bin/TbToyData.o:	$(ALGODIR)src/TbToyData.cpp $(ALGODIR)src/TbToyData.h
			${CC} $(CXXFLAGS) -c $(ALGODIR)src/TbToyData.cpp -o $(ALGODIR)bin/TbToyData.o



clean:
	rm -fr $(DIR)bin/*
	rm -fr $(KERNELDIR)bin/*
	rm -fr $(ALGODIR)bin/*

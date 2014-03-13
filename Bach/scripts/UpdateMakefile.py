import os, sys

def UM():
    classlist = []
    for l in os.listdir("../../TbAlgorithms/src/"):
        if ".h" in l or "~" in l:
            continue
        classlist.append( l.strip(".cpp") )



    s = """
# Makefile depends on ROOT. 
# Define the variables ROOTSYS in your environment.


ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --glibs)

CXXFLAGS      = -Wall -fPIC  -g -W
CXXFLAGS     += $(ROOTCFLAGS) 

LIBS          = $(ROOTLIBS) -lGenVector
CC            = g++

DIR 	      := $(PWD)/
KERNELDIR     := $(DIR)../TbKernel/
ALGODIR       := $(DIR)../TbAlgorithms/



OBJS = $(KERNELDIR)bin/TbGeometrySvc.o $(KERNELDIR)bin/Millepede.o """

    for c in classlist:
        s += "$(ALGODIR)bin/"+c+".o "

    s += """

all:			$(DIR)bin/bach

$(DIR)bin/bach: 	$(OBJS)
			${CC} $(CXXFLAGS) $(LIBS) -o $(DIR)bin/bach $(DIR)src/bach.cpp $(OBJS) 

$(KERNELDIR)bin/TbGeometrySvc.o:	$(KERNELDIR)src/TbGeometrySvc.cpp $(KERNELDIR)src/TbGeometrySvc.h
			${CC} $(CXXFLAGS) -c $(KERNELDIR)src/TbGeometrySvc.cpp -o $(KERNELDIR)bin/TbGeometrySvc.o

$(KERNELDIR)bin/Millepede.o:	$(KERNELDIR)src/Alignment/Millepede.cpp $(KERNELDIR)src/Alignment/Millepede.h
			${CC} $(CXXFLAGS) -c $(KERNELDIR)src/Alignment/Millepede.cpp -o $(KERNELDIR)bin/Millepede.o
"""

    for c in classlist:
        s += """
$(ALGODIR)bin/%s.o:	$(ALGODIR)src/%s.cpp $(ALGODIR)src/%s.h
			${CC} $(CXXFLAGS) -c $(ALGODIR)src/%s.cpp -o $(ALGODIR)bin/%s.o

""" % (c,c,c,c,c)

    s += """



clean: 
	rm -fr $(DIR)bin/*
	rm -fr $(KERNELDIR)bin/*
	rm -fr $(ALGODIR)bin/*
"""
    f = open("../Makefile",'w')
    f.write(s)
    f.close()
    print "Makefile generated"

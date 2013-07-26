ROOTCFLAGS  = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS    = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS   = $(shell $(ROOTSYS)/bin/root-config --glibs)

CXX         = g++ -m32
#CXX         = g++
#CXXFLAGS    = -g -fPIC -Wno-deprecated -0 -ansi -D_GNU_SOURCE -g -02
CXXFLAGS    = -Wno-deprecated

CXXFLAGS   += $(ROOTCFLAGS)
LIBS        = $(ROOTLIBS)

NGLIBS      = $(ROOTGLIBS)
NGLIBS     += -lMinuit
GLIBS     = $(filter-out -lNew -lz , $(NGLIBS))

#INCLUDEDIR  = ./include/
#SRCDIR      = ./src/
#CXX        += -I$(INCLUDEDIR) -I.
#OUTLIB      = ./lib/


.SUFFIXES: .cc,.C,.cxx,.hh,.h
#.PREFIXES: ./lib/

all: runabcd

runabcd: runABCD.C \
	abcdMethod.o 
	$(CXX) $(CXXFLAGS) -o runabcd *.o $(GLIBS) $ $<
	touch runabcd

abcdMethod.o: abcdMethod.cc
	$(CXX) $(CXXFLAGS) -c -o abcdMethod.o $<


# clean
clean: 
	rm -f *.o 
	rm -f runabcd

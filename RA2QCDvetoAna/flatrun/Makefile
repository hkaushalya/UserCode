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

all: runsmear

runsmear: FactorizationBySmearing.cc \
	NtupleSelector.o \
        mydict.o \
	SmearFunction.o \
	Utils.o 
	$(CXX) $(CXXFLAGS) -o runsmear *.o $(GLIBS) $ $<
	touch runsmear

NtupleSelector.o: NtupleSelector.cc
	$(CXX) $(CXXFLAGS) -c -o NtupleSelector.o $<

SmearFunction.o: SmearFunction.cc
	$(CXX) $(CXXFLAGS) -c -o SmearFunction.o $<

Utils.o: Utils.cc
	$(CXX) $(CXXFLAGS) -c -o Utils.o $<

mydict.cxx: myLinkDef.h
	@echo "Generating dictionary ..."
	@rootcint mydict.cxx -c myLinkDef.h

mydict.o: mydict.cxx
	$(CXX) $(CXXFLAGS) -c -o mydict.o $<

# clean
clean: 
	rm -f *.o 
	rm -f runsmear

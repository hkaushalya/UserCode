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
NGLIBS     += -lMinuit2
GLIBS     = $(filter-out -lNew -lz , $(NGLIBS))

#INCLUDEDIR  = ./include/
#SRCDIR      = ./src/
#CXX        += -I$(INCLUDEDIR) -I.
#OUTLIB      = ./lib/


.SUFFIXES: .cc,.C,.cpp,.cxx,.hh,.h
#.PREFIXES: ./lib/

all: runsmear

runsmear: FactorizationBySmearing.cc \
	NtupleSelector.o \
        mydict.o \
	SmearFunction.o \
	Utils.o \
	ChengHanBisect_Mt2_332_Calculator.o \
	pluginUserCodeTopTagger.so	
	$(CXX) $(CXXFLAGS) -o runsmear *.o $(GLIBS) $ $<
	touch runsmear

NtupleSelector.o: NtupleSelector.cc
	$(CXX) $(CXXFLAGS) -c -o NtupleSelector.o $<

ChengHanBisect_Mt2_332_Calculator.o: ChengHanBisect_Mt2_332_Calculator.cpp \
	SolutionType.h \
	Mt2Units.h \
	Mt2Util.h \
	Mt2_332_Calculator.h \
	Mt2Calculator.o \
	mt2_bisect.o \
	Mt2Vectors.o \
	Mt2TwoVector.o \
	Mt2LorentzTransverseVector.o \
	Basic_Mt2_332_Calculator.o
	$(CXX) $(CXXFLAGS) -c -o ChengHanBisect_Mt2_332_Calculator.o $<

mt2_bisect.o: mt2_bisect.cpp
	$(CXX) $(CXXFLAGS) -c -o mt2_bisect.o $<

Mt2Vectors.o: Mt2Vectors.cpp
	$(CXX) $(CXXFLAGS) -c -o Mt2Vectors.o $<

Mt2TwoVectors.o: Mt2TwoVectors.cpp
	$(CXX) $(CXXFLAGS) -c -o Mt2TwoVectors.o $<

Mt2LorentzTransverseVector.o: Mt2LorentzTransverseVector.cpp
	$(CXX) $(CXXFLAGS) -c -o Mt2LorentzTransverseVector.o $<

Basic_Mt2_332_Calculator.o: Basic_Mt2_332_Calculator.cpp
	$(CXX) $(CXXFLAGS) -c -o Basic_Mt2_332_Calculator.o $<

Mt2Calculator.o: Mt2Calculator.cpp
	$(CXX) $(CXXFLAGS) -c -o Mt2Calculator.o $<

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

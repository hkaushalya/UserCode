#ifndef UTILS_HH
#define UTILS_HH
#include "TLorentzVector.h"
#include "TH1.h"
#include <iostream>

#define PRINTER(name) printer(#name, (name))

template <class T>
inline void printer(const char *name, const T& val)
{
	        std::cout << name << " = " << val << std::endl;
}
bool PtAComparator(const TLorentzVector& a,const TLorentzVector& b);
void DumpHist(const TH1* h);
void printProgBar( int percent );


//to get all vector like inputs
//vector<double> NumStringToVec(const string s);
#endif

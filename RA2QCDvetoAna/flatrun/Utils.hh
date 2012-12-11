#ifndef UTILS_HH
#define UTILS_HH
#include "TLorentzVector.h"
#include "TH1.h"
using namespace std;

bool PtAComparator(const TLorentzVector& a,const TLorentzVector& b);
void DumpHist(const TH1* h);

//to get all vector like inputs
//vector<double> NumStringToVec(const string s);
#endif

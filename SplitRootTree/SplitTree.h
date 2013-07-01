#ifndef SplitTree_H
#define SplitTree_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "NtupleSelector.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <sstream>
#include "IOColors.hh"
#include "TDirectory.h"
#include <utility>

/***************************************************************
 * See cc file for class description.
 * Sam Hewamanage, Florida International University
 **************************************************************/

class SplitTree : public NtupleSelector {

	public:
		SplitTree(const TString &inputFileList="foo.txt");
		~SplitTree();
		Bool_t   FillChain(TChain *chain, const TString &inputFileList);
		Long64_t LoadTree(Long64_t entry);
		void     EventLoop(const string outFileName, const int evt2Process, const unsigned evtsPerFile);


	private:
		TFile *outRootFile;
		double smearedJetPt_;
		unsigned uNTRIES;

		void DumpJets(const vector<TLorentzVector>& jets, const double& minPt=10.0) const ;
		void DumpJet(const TLorentzVector& jet) const ;
		void DumpJetsAndCSV(const vector<TLorentzVector>& jets, const vector<double>& csv,
							const TLorentzVector& met) const;
		void PrintEventNumber() const;
		void NewOutputFile(const string outFileName, TFile*& fnew, TTree*& newTree);	
		void CloseFile(TFile* f);
		unsigned fileCounter;  //suffix fore new output files
};
#endif

#ifdef SplitTree_cxx


Bool_t SplitTree::FillChain(TChain *chain, const TString &inputFileList) {

	ifstream infile(inputFileList, ifstream::in);
	std::string buffer;

	if(!infile.is_open()) {
		std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
		return kFALSE;
	}

	while(1) {
		infile >> buffer;
		if(!infile.good()) break;
		//std::cout << "Adding tree from " << buffer.c_str() << std::endl;                                                              
		chain->Add(buffer.c_str());
	}
	std::cout << "No. of Entries in this tree = " << chain->GetEntries() << std::endl;
	return kTRUE;
}

Long64_t SplitTree::LoadTree(Long64_t entry) {
	// Set the environment to read one entry                                                                                          
	if (!fChain) return -5;
	Long64_t centry = fChain->LoadTree(entry);
	if (centry < 0) return centry;
	if (!fChain->InheritsFrom(TChain::Class()))  return centry;
	TChain *chain = (TChain*)fChain;
	if (chain->GetTreeNumber() != fCurrent) {
		fCurrent = chain->GetTreeNumber();
		//    Notify();
	}
	return centry;
}

SplitTree::~SplitTree() { 

	if (!fChain) return;
	delete fChain->GetCurrentFile();
}

#endif

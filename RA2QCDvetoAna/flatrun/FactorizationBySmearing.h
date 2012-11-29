#ifndef FactorizationBySmearing_H
#define FactorizationBySmearing_H

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
#include "SmearFunction.h"
#include <sstream>
#include "IOColors.hh"
#include "TDirectory.h"

/***************************************************************
 * See cc file for class description.
 * Sam Hewamanage, Florida International University
 **************************************************************/

class FactorizationBySmearing : public NtupleSelector{

	public:
		FactorizationBySmearing(const TString &inputFileList="foo.txt", 
								const char *outFileName="histo.root",
								const int njet50min=3, const int njet50max=1000); 
		~FactorizationBySmearing();
		Bool_t   FillChain(TChain *chain, const TString &inputFileList);
		Long64_t LoadTree(Long64_t entry);
		void     EventLoop(const char * datasetname, const int evt2Process);
		void     BookHistogram(const char *);
		double   DeltaPhi(double, double);
		double   DeltaR(double eta1, double phi1, double eta2, double phi2);
		double   HT (const std::vector<TLorentzVector>&);
		TLorentzVector MHT(const std::vector<TLorentzVector>&);
		TVector3 GetMHT(vector<TLorentzVector>);
		void CreateRecoJetVec(std::vector<TLorentzVector>& vjets);
		void CreateGenJetVec(std::vector<TLorentzVector>& vjets);
		double JetResolutionHist_Pt_Smear(const double& pt, const double& eta, const int& i);
		void SmearingGenJets(const vector<TLorentzVector>& jets_gen, std::vector<TLorentzVector> &genJets_smeared);
		int GetIndex(const double& x, const std::vector<double>* vec);
		double GetHFProb(const double& pt, const double& eta, const int& i_jet);
		void JetDeltaMin(const string s);
		void PtBinEdges_scaling(const string s);
		void EtaBinEdges_scaling(const string s);
		void AdditionalSmearing(const string s);
		void LowerTailScaling(const string s);
		void UpperTailScaling(const string s);
		void PtBinEdges(const string s);
		void EtaBinEdges(const string s);
		bool PassCuts(const vector<TLorentzVector>& jets);
		vector<TLorentzVector> GetPt50Eta2p5Jets(const vector<TLorentzVector>& jets);

		TFile *oFile;
			
		struct JetHist_t {
			TH1D *h_Jet_pt;
			TH1D *h_Jet_eta;
			TH1D *h_Jet_phi;
			TH1D *h_Jet_dphi;
		};
		struct CommonHist_t {
			TH1D *h_Njet50eta2p5;
			TH1D *h_Njet30eta5p0;
			TH1D *h_Mht;      //this Mht/Ht plots are for cut confirmations. Total HT/MHT plots are separately made
 			TH1D *h_Ht;
			TH1D *h_DphiMin;
			TH2D *h_DphiMinVsMht;
			TH1D* pass[6];
			TH1D* passFineBin[6];
			TH1D* fail[6];
			TH1D* failFineBin[6];
			TH1D* signal;
			TH1D* signalFineBin;
			TH1D* sidebandSyst[2];
			TH1D* sidebandSystFineBin[2];
		};


		struct Hist_t {
			CommonHist_t hv_RecoEvt;
			CommonHist_t hv_GenEvt;
			CommonHist_t hv_SmearedEvt;
			vector<JetHist_t> hv_RecoJets;
			vector<JetHist_t> hv_GenJets;
			vector<JetHist_t> hv_SmearedJets;
		};


	private:
		bool bDEBUG;
		SmearFunction *smearFunc_;
		double smearedJetPt_;
		std::vector<double> PtBinEdges_scaling_;
		std::vector<double> EtaBinEdges_scaling_;
		std::vector<double> AdditionalSmearing_;
		std::vector<double> LowerTailScaling_;
		std::vector<double> UpperTailScaling_;
		double AdditionalSmearing_variation_;
		double LowerTailScaling_variation_;
		double UpperTailScaling_variation_;
		bool absoluteTailScaling_;
		bool applyDphiCut_; 

		double MHTcut_low_;
		double MHTcut_medium_;
		double MHTcut_high_;
		double HTcut_low_;
		double HTcut_medium_;
		double HTcut_high_;
		double HTcut_veryhigh_;
		double HTcut_extremehigh_;
		std::vector<double> PtBinEdges_, EtaBinEdges_, HtBins_, MhtBins_, NjetBins_;
		vector<vector<Hist_t> > Hist; //for each HT/MHT bins
		unsigned NJet50_min_, NJet50_max_;
		std::vector<float> vDphiVariations;
		double nRecoJetEvts, nGenJetEvts, nSmearedJetEvts;

		vector<vector<TH1*> > jerHist; //for each pt/eta bins
		void BookJerDebugHists();

		void DumpJets(const vector<TLorentzVector>& jets);
		void GetHist(TDirectory *dir, Hist_t& hist, 
				const float htmin, const float htmax,
				const float mhtmin, const float mhtmax);
		void GetJetHist(vector<JetHist_t>& Hist, const string jetcoll
							, const string htmhtrangelabel	);
		unsigned GetVectorIndex(const vector<double>& bins, const double& val);
		void FillHistogram(const vector<TLorentzVector>& jets, const int jetcoll, const double& wgt=1.0);
		int CountJets(const std::vector<TLorentzVector>& vjets, 
			const double minPt, const double maxEta);		
		void FillJetHistogram(const vector<TLorentzVector>& jets, 
				vector<JetHist_t> hist, const int jetcoll, const TLorentzVector& mhtvec, const double& wgt=1.0); 
		float DelPhiMin(const vector<TLorentzVector>& jets, 
					const TLorentzVector& mhtVec);
		bool PassRA2dphiCut(const vector<TLorentzVector>& jets, const TLorentzVector& mht);
		double GetLumiWgt(const string& datasetname, const double& dataLumi);
		void DivideByBinWidth(TH1* h);
		bool PassCleaning();
};
#endif

#ifdef FactorizationBySmearing_cxx

FactorizationBySmearing::FactorizationBySmearing(
				const TString &inputFileList, 
				const char *outFileName,
				const int njet50min,
				const int njet50max
				) {

	//isItMC = true;
	//isItZ  = false;

	TChain *tree = new TChain("treeMaker/tree");  

	if( ! FillChain(tree, inputFileList) ) {
		std::cerr << "Cannot get the tree " << std::endl;
		assert(false);
	}

	Init(tree);

	HtBins_.push_back(0);
	HtBins_.push_back(500);
	HtBins_.push_back(750);
	HtBins_.push_back(1000);
	HtBins_.push_back(1250);
	HtBins_.push_back(1500);
	HtBins_.push_back(8000);

	MhtBins_.push_back(0);
	//MhtBins_.push_back(200);
	//MhtBins_.push_back(350);
	//MhtBins_.push_back(500);
	MhtBins_.push_back(8000);

	//jet bins: 2, 3-5, 6-7,>=8
	NJet50_min_  = njet50min;
	NJet50_max_  = njet50max;
	nRecoJetEvts = 0;
	nGenJetEvts  = 0;
	nSmearedJetEvts = 0;

	BookHistogram(outFileName);
}

Bool_t FactorizationBySmearing::FillChain(TChain *chain, const TString &inputFileList) {

	ifstream infile(inputFileList, ifstream::in);
	std::string buffer;

	if(!infile.is_open()) {
		std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
		return kFALSE;
	}

	std::cout << "TreeUtilities : FillChain " << std::endl;
	while(1) {
		infile >> buffer;
		if(!infile.good()) break;
		//std::cout << "Adding tree from " << buffer.c_str() << std::endl;                                                              
		chain->Add(buffer.c_str());
	}
	std::cout << "No. of Entries in this tree : " << chain->GetEntries() << std::endl;
	return kTRUE;
}

Long64_t FactorizationBySmearing::LoadTree(Long64_t entry) {
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

FactorizationBySmearing::~FactorizationBySmearing() { 

	if (!fChain) return;
	delete fChain->GetCurrentFile();
	oFile->cd();
	oFile->Write();
	oFile->Close();

}

double FactorizationBySmearing::DeltaPhi(double phi1, double phi2) {
	double result = phi1 - phi2;
	while (result > M_PI)    result -= 2 * M_PI;
	while (result <= -M_PI)  result += 2 * M_PI;
	return result;
}

double FactorizationBySmearing::DeltaR(double eta1, double phi1, double eta2, double phi2) {
	double deta = eta1 - eta2;
	double dphi = DeltaPhi(phi1, phi2);
	return std::sqrt(deta*deta + dphi*dphi);
}

#endif

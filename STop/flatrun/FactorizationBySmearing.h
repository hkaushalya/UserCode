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
#include <utility>
#include "Type3TopTagger.h"

/***************************************************************
 * See cc file for class description.
 * Sam Hewamanage, Florida International University
 **************************************************************/

class FactorizationBySmearing : public NtupleSelector {

	public:
		FactorizationBySmearing(
								const std::vector<unsigned> vCutMasks, 
								const TString &inputFileList="foo.txt", 
								const char *outFileName="histo.root"
								);
		~FactorizationBySmearing();
		Bool_t   FillChain(TChain *chain, const TString &inputFileList);
		Long64_t LoadTree(Long64_t entry);
		void     EventLoop(const char * datasetname, const int evt2Process, const int systVariations);
		//void     BookHistogram(const char *);
		void     BookHistogram(TFile *oFile, const bool mcFlag);
		double   DeltaPhi(double, double);
		double   DeltaR(double eta1, double phi1, double eta2, double phi2);
		double   HT (const std::vector<TLorentzVector>&);
		TLorentzVector MHT(const std::vector<TLorentzVector>& jets,
							const double& minPt=30.0, const double& maxEta=5.0); 
		TLorentzVector GetMET(const vector<TLorentzVector>& jets,
							const double& minPt=10.0, const double& maxEta=5.0); 
		void CreateRecoJetVec(std::vector<TLorentzVector>& vjets, std::vector<double>& bDisc);
		void CreateGenJetVec(std::vector<TLorentzVector>& vjets);
		double JetResolutionHist_Pt_Smear(const double& pt, const double& eta, const int& i);
		double BJetResolutionHist_Pt_Smear(const double& pt, const double& eta, const int& i);
		void SmearingGenJets(const vector<TLorentzVector>& jets_gen, const vector<double>& bDisc_gen, 
						vector<unsigned>& bJetInds,
						std::vector<TLorentzVector> &genJets_smeared, std::vector<double>& bDisc_smeared);
		void SetGenJetBdiscriminators(const vector<TLorentzVector>& jets_reco, const std::vector<double>& bDisc_reco,
									const std::vector<TLorentzVector> &genJets, std::vector<double>& bDisc_gen);
		int GetIndex(const double& x, const std::vector<double>* vec);
		double GetHFProb(const double& pt, const double& eta, const int& i_jet);
		bool PassCuts(const vector<TLorentzVector>& jets);
		vector<TLorentzVector> GetPt50Eta2p5Jets(const vector<TLorentzVector>& jets);
		void Print4vec(const TLorentzVector& tl1, const TLorentzVector& tl2) const;
		void Print4vec(const TLorentzVector& tl1) const;
		TLorentzVector GetSmearUnclMet(const TLorentzVector& met);

		struct JetHist_t {
			TH1D *h_Jet_pt;
			TH1D *h_Jet_eta;
			TH1D *h_Jet_phi;
			TH1D *h_Jet_dphimet;
			TH1D *h_Jet_dphimht;
		};
		struct CommonHist_t {
			TH1D *h_evtWeight;
			TH1D *h_nVtx;
			TH1D *h_Njet50eta2p5;
			TH1D *h_Njet30eta5p0;
			TH1D *h_Njet70eta2p4;
			TH1D *h_Njet50eta2p4;
			TH1D *h_Njet30eta2p4;
			TH1D *h_Mht;
			TH1D *h_Met;
			TH1D *h_Met1;
			TH1D *h_Met2;
			TH1D *h_Met3;
			TH1D *h_Met4;
			TH1D *h_UnclMet;
 			TH1D *h_Ht;
 			TH1D *h_dPhiMin_met;
 			TH1D *h_dPhiMin_mht;
			TH1D *h_nbjets;
			TH1D *h_bjetMass;
			TH1D *h_bjetPt;
			TH1D *h_M123;
			TH1D *h_MT2;
			TH1D *h_MTb;
			TH1D *h_MTt;
			TH1D *h_MTb_p_MTt;
			TH1D* h_passDphi;
			TH1D* h_passDphiFineBin;
			TH1D* h_failDphi;
			TH1D* h_failDphiFineBin;
		};


		struct Hist_t {
			CommonHist_t hv_RecoEvt;
			CommonHist_t hv_GenEvt;
			CommonHist_t hv_SmearedEvt;
			vector<JetHist_t> hv_RecoJets;
			vector<JetHist_t> hv_GenJets;
			vector<JetHist_t> hv_SmearedJets;
		};

		enum CutMasks { NoCutMask_=0, NjetCutMask_= 1, metCutMask_=2, tripletCutMask_=4, topmassCutMask_=8, topplusbjetCutMask_=16, mt2CutMask_=32, minbjetCutMask_=64, dphiCutMask_=128, invdphiCutMask_=256};

	private:
		TFile *outRootFile;
		vector<string> vBadHcalLaserEvts, vBadEcalLaserEvts;
		bool bDEBUG, bNON_STD_MODE;
		string sNON_STD_MODE_EXPLAIN;
		bool bRUNNING_ON_MC;
		bool bUSE_BJET_SMEAR_FUNC; //control using b-jet reponse functions to smear b-jets
		bool bDO_TRIG_PRESCALING, bDO_TRIG_SELECTION, bDO_PU_WEIGHING, bDO_GENJET_SMEARING;
		bool bDO_LUMI_WEIGHING;
		SmearFunction *smearFunc_, *bjetsmearFunc_;
		double smearedJetPt_;
		unsigned uNTRIES;
		std::vector<double> PtBinEdges_scaling_;
		std::vector<double> EtaBinEdges_scaling_;
		std::vector<double> AdditionalSmearing_;
		std::vector<double> LowerTailScaling_;
		std::vector<double> UpperTailScaling_;
		double AdditionalSmearing_variation_;
		double LowerTailScaling_variation_;
		double UpperTailScaling_variation_;
		bool absoluteTailScaling_;
		unsigned uMinNjet70Eta2p4_, uMinNjet50Eta2p4_, uMinNjet30Eta2p4_, uMinTriplets_, uMinBjets_;
		double dMinMet_, dMaxMet_, dMinTopMass_, dMaxTopMass_, dMinTopPlusBjetMass_, dMinMt2_;
		unsigned nBadEcalLaserEvts;

		vector<double> HtBins_, MhtBins_;
		vector <unsigned> BitMaskBins_;
		vector<vector<vector<Hist_t> > > Hist; //for each njet/HT/MHT bins
		vector<double> vDphiVariations;
		double nRecoJetEvts, nGenJetEvts, nSmearedJetEvts, nVectorInexWarnings;
		vector<string> vTriggersToUse; 

		vector<vector<TH1*> > jerHist; //for each pt/eta bins
		void BookJerDebugHists();

		void DumpJets(const vector<TLorentzVector>& jets, const double& minPt=10.0) const ;
		void DumpJet(const TLorentzVector& jet) const ;
		void DumpJetsAndCSV(const vector<TLorentzVector>& jets, const vector<double>& csv,
							const TLorentzVector& met) const;
		vector<unsigned> FindBjets(const vector<TLorentzVector>& jets, const vector<double>& bDisc);
		void GetHist(TDirectory *dir, Hist_t& hist, 
					const string bitMaskString,
					const pair<unsigned, unsigned> htrange,
					const pair<unsigned, unsigned> mhtrange,
					const bool mcFlag
					);
		void GetJetHist(vector<JetHist_t>& Hist, const string jetcoll
							, const string htmhtrangelabel	);
		unsigned GetVectorIndex(const vector<double>& bins, const double& val);
		unsigned GetVectorIndex(const vector< pair<unsigned, unsigned> >& binEdges, const unsigned& val);
		bool FillHistogram(const vector<TLorentzVector>& jets, const vector<unsigned>& bJetInds, 
						const vector<double>& bDisc, 
						const TLorentzVector& reco_uncl_met, const int& jetcoll, const double& wgt=1.0);
		int CountJets(const std::vector<TLorentzVector>& vjets, 
			const double minPt, const double maxEta);		
		void FillJetHistogram(const vector<TLorentzVector>& jets, 
				vector<JetHist_t> hist, const int jetcoll, 
				const TLorentzVector& metvec, const TLorentzVector& mhtvec, const double& wgt=1.0); 
		double DelPhiMin(const vector<TLorentzVector>& jets, 
					const TLorentzVector& mhtVec, const unsigned& njet50min);
		bool PassDphiCut(const vector<TLorentzVector>& jets, 
					const TLorentzVector& mht, const unsigned njet50min, 
					const double& cut1, const double& cut2, const double& cut3);
		double GetLumiWgt(const string& datasetname, const double& dataLumi);
		void DivideByBinWidth(TH1* h);
		bool PassCleaning();
		void TrigPrescaleWeight(bool &failTrig, double &weight) const;
		void LoadBadHcalLaserEvents();
		void LoadBadEcalLaserEvents();
		void PrintEventNumber() const;
		bool PassHOfilter();
   	void TripletSelector(const std::vector<TLorentzVector> & , const std::vector<double> bDisc, std::vector<TLorentzVector> & , std::vector<TLorentzVector> & , std::vector<TLorentzVector> & , double & , double & );
		topTagger::type3TopTagger* topt;
		unsigned nreco, ngen, nsmear;
		string GetMaskString(const unsigned& bitMask);
		unsigned GetMaskBin(const unsigned& passBitMask);
		unsigned tnjetonly, tmetonly, tnjetmetonly, tnjetmetdphionly;
};
#endif

#ifdef FactorizationBySmearing_cxx


Bool_t FactorizationBySmearing::FillChain(TChain *chain, const TString &inputFileList) {

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
	outRootFile->cd();
	outRootFile->Write();
	outRootFile->Close();
   if (bRUNNING_ON_MC && bDO_GENJET_SMEARING && smearFunc_ != 0) delete smearFunc_;
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

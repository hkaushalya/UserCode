#ifndef PHOTONJETS_HH
#define PHOTONJETS_HH
//
// Original Author:  samantha hewamanage
// $Id$
//
//


// system include files
#include <memory>
#include <vector>
#include <algorithm>
#include <iostream>

//ROOT
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "assert.h"
#include "TVector3.h"
#include <iomanip>
#include "TPad.h"
#include "TBenchmark.h"
#include "TChain.h"
#include "TFile.h"
#include "TROOT.h"

using namespace std;


//
// class declaration
//

class Factorization {
	public:
		Factorization();
		~Factorization();
		int Run(TChain* chain, int iRunEvents = 0);
		void Init(TChain* outTree);
		int iProgressBy;
		void SetHistFileName(const std::string name) { sHistFileName = name;}
		void BookHistograms();

		static const UInt_t Njet_max = 100;
		Int_t run, event, ls; bool isdata;
		Int_t nVtxPUcut_, vtxSize;
		Int_t nJets;
		double evtWeight_TR;
		double mht_TR, ht_TR, met_TR, mt_TR, meff;
		double mhtphi_TR, metphi_TR;
		double metSgnf_TR, metSgnfProb_TR;
		double jet1pt_TR, jet1eta_TR, jet1phi_TR;
		double jet2pt_TR, jet2eta_TR, jet2phi_TR;
		double jet3pt_TR, jet3eta_TR, jet3phi_TR;
		UInt_t jet_num;
		Float_t alljets_E[Njet_max], alljets_Px[Njet_max], alljets_Py[Njet_max], alljets_Pz[Njet_max];
		Float_t alljets_Pt[Njet_max], alljets_Eta[Njet_max];
		float storedWeight; //for flat MC sample weights
		float PUWeight;
		int mcFlag;   //0=DATA, 1= MC
		int npv;




	private:
		void endJob() ;

		// ----------member data ---------------------------
		//jet collections

		struct RunLumiEvt_t
		{
			unsigned run;
			unsigned lumi;
			unsigned evt;
		};


		int inMinVtx;
		unsigned int uProcessed; // number of processed events
		unsigned int uPassed; //number of events passed the filter
		RunLumiEvt_t kRunLumiEvent;
		int iVerbose; // control print levels
		double dMinHT;
		double dMinMHT;

		struct EventHist_t {
			TH1F* mht;
			TH1F* ht;
			TH1F* njet;
			TH1F* meff;
		};

		struct JetHist_t{
			TH1F* pt;
			TH1F* eta;
			TH1F* phi;
			TH1F* delphi;  //delphi(jet,MHT)
		};

		struct Hist_t {  //hists for each inclusive jet category
			EventHist_t evt;
			JetHist_t jet[5];
			TH1F* delphiMin_jetmet[3];
			TProfile* delphiMin_jetmetVsMHT[3];
			TH1F* delphiMin_jetjet[3];
			TProfile* delphiMin_jetjetVsMHT[3];
		};

		TH1F* hDelPhiMin_mht[8];
		TH1F* hPass[9];
		TH1F* hFail[20];
		TH2F* dphiVsMHT_2D;

		EventHist_t evtHist;
		JetHist_t pf30_jet1Hist, pf30_jet2Hist, pf30_jet3Hist; //for upto 5 jets may be
		JetHist_t pf50eta25_jet1Hist, pf50eta25_jet2Hist, pf50eta25_jet3Hist; //for upto 5 jets may be

		Hist_t Hist; //all hists for various mht ranges
		TH1F* MHT_by_phislice[6];
		void DoDelMinStudy();
		void PrintHeader();
		TLorentzVector vMetVec;
		unsigned uFailMinHTCut, uFailMinPFMHTCut;

		std::vector<float> vMCNvtxDist, vDATANvtxDist;
		double sumLumiWeights; //sum of lumi weights for cross checks
		double Weight; //total weight for an event
		int doLumiWeighing, doEventWeighing;
		TH1F* hSignalRegion[3];

		TFile* histFile;
		std::string sHistFileName;

};

#endif

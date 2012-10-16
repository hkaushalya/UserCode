#ifndef FactorizationBySmearing_H
#define FactorizationBySmearing_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "NtupleSelector.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "SmearFunction.h"
#include <sstream>
#include "IOColors.hh"
#include "TDirectory.h"

class FactorizationBySmearing : public NtupleSelector{

	public:
		FactorizationBySmearing(const TString &inputFileList="foo.txt", const char *outFileName="histo.root"); 
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

		TFile *oFile;
		TH1F *h_Total_HT, *h_Total_MHT, *h_Total_MEff, *h_Total_MHTSig;
		TH1F *h_gen_Njet50eta2p5, *h_gen_Njet30eta5p0, *h_gen_Jet1_pt, *h_gen_Jet1_eta;
		TH1F *h_gen_smeared_Njet50eta2p5, *h_gen_smeared_Njet30eta5p0, *h_gen_smeared_Jet1_pt, *h_gen_smeared_Jet1_eta;
			
		struct JetHist_t {
			TH1F *h_Jet_pt;
			TH1F *h_Jet_eta;
			TH1F *h_Jet_phi;
			TH1F *h_Jet_dphi;
		};
		struct CommonHist_t {
			TH1F *h_Njet50eta2p5;
			TH1F *h_Njet30eta5p0;
			TH1F *h_Mht;      //this Mht/Ht plots are for cut confirmations. Total HT/MHT plots are separately made
 			TH1F *h_Ht;
			TH1F *h_DphiMin;
			TH2F *h_DphiMinVsMht;
			TH1F* pass[6];
			TH1F* passFineBin[6];
			TH1F* fail[6];
			TH1F* failFineBin[6];
			TH1F* signal;
			TH1F* signalFineBin;
			TH1F* sidebandSyst[2];
			TH1F* sidebandSystFineBin[2];
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
		bool debug;
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
				vector<JetHist_t> hist, const int jetcoll, const TLorentzVector& mhtvec, const double& wgt); 
		float DelPhiMin(const vector<TLorentzVector>& jets, 
					const TLorentzVector& mhtVec);
		bool PassRA2dphiCut(const vector<TLorentzVector>& jets, const TLorentzVector& mht);
};
#endif

#ifdef FactorizationBySmearing_cxx

void FactorizationBySmearing::BookHistogram(const char *outFileName) {

	oFile = new TFile(outFileName, "recreate");
	TDirectory *dir = oFile->mkdir("Hist");
	if (debug) cout << __FUNCTION__ << " dir = " << dir->GetName() << endl;

	h_Total_HT   = new TH1F("h_Total_HT",   "h_Total_HT",   100, 0.0, 5000.0);  
	h_Total_MHT  = new TH1F("h_Total_MHT",  "h_Total_MHT",  100, 0.0, 2000.0);
	h_Total_MEff = new TH1F("h_Total_MEff", "h_Total_MEff", 100, 0.0, 5000.0);
	h_Total_MHTSig= new TH1F("h_Total_MHTSig","h_Total_MHTSig",100, 0.0, 100.0);
	h_Total_HT->Sumw2();
	h_Total_MHT->Sumw2();
	h_Total_MEff->Sumw2();
	h_Total_MHTSig->Sumw2();
	//h_Njet50eta2p5->Sumw2();
	//h_Njet30eta5p0->Sumw2();

	//h_Njet50eta2p5   = new TH1F("h_Njet50eta2p5",   "Njets pt>50, |eta|<2.5",   15, 0.0, 15);  
	//h_Njet30eta5p0   = new TH1F("h_Njet30eta5p0",   "Njets pt>30, |eta|<5.0",   15, 0.0, 15);  


	for (unsigned htbin=0; htbin < HtBins_.size() -1; ++htbin)
	{
		vector<Hist_t> h;
		for (unsigned mhtbin=0; mhtbin < MhtBins_.size() -1; ++mhtbin)
		{
			const float htmin = HtBins_.at(htbin);
			const float htmax = HtBins_.at(htbin+1);
			const float mhtmin = MhtBins_.at(mhtbin);
			const float mhtmax = MhtBins_.at(mhtbin+1);
			stringstream folder;
			folder << "HT" << htmin << "to" << htmax
						<<  "MHT" << mhtmin << "to" << mhtmax;
			//oFile->ls();
			//cout << __LINE__ << ":" << folder.str() << endl; 
			dir->cd();
			//gDirectory->pwd();
			TDirectory *subdir = dir->mkdir(folder.str().c_str());
			//cout << "subdir = " << subdir->GetName() << endl;
			//subdir->ls();
			subdir->cd();
			//gDirectory->pwd();

			Hist_t hist;
			GetHist(subdir, hist, htmin, htmax, mhtmin, mhtmax);
			h.push_back(hist);
		}
		Hist.push_back(h);
	}

}


FactorizationBySmearing::FactorizationBySmearing(const TString &inputFileList, const char *outFileName) {

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
	HtBins_.push_back(900);
	HtBins_.push_back(1300);
	HtBins_.push_back(7000);

	MhtBins_.push_back(0);
//	MhtBins_.push_back(200);
//	MhtBins_.push_back(350);
//	MhtBins_.push_back(500);
	MhtBins_.push_back(7000);
	

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



void FactorizationBySmearing::GetHist(TDirectory *dir, Hist_t& hist, 
					const float htmin, const float htmax,
					const float mhtmin, const float mhtmax
					)
{
	stringstream htmhtrange;
	htmhtrange << 	htmin << "<HT<" << htmax << " GeV, " << mhtmin << "<MHT<" << mhtmax << " GeV";
	
	vector<string> jetcoll;
	jetcoll.push_back("reco");
	jetcoll.push_back("gen");
	jetcoll.push_back("smeared");

	dir->cd();
	for (unsigned i = 0;i < jetcoll.size(); ++i)
	{
		stringstream njet50eta2p5title, njet30eta5p0title, httitle, mhttitle;
		njet50eta2p5title << htmhtrange.str().c_str() << ";Njets [ET>50 GeV,| #eta |<2.5];Events;";
		njet30eta5p0title << htmhtrange.str().c_str() << ";Njets [ET>30 GeV,| #eta |<5.0];Events;";
		httitle << htmhtrange.str().c_str() << ";HT [3 Jets, PT>50 GeV, | #eta |<2.5];Events;";
		mhttitle << htmhtrange.str().c_str() << ";MHT [PT>30 GeV, | #eta |<5.0];Events;";

		stringstream njet50eta2p5name, njet30eta5p0name, htname, mhtname;
		njet50eta2p5name << jetcoll.at(i) << "_njet50eta2p5";
		njet30eta5p0name << jetcoll.at(i) << "_njet30eta5p0";
		htname << jetcoll.at(i) << "_ht";
		mhtname << jetcoll.at(i) << "_mht";

		stringstream dphiminname, dphiminvsmhtname;
		stringstream dphimintitle, dphiminvsmhttitle;
		dphimintitle << htmhtrange.str().c_str() << ";#Delta #Phi_{min};Events;";
		dphiminname << jetcoll.at(i) << "_dphimin";
		dphiminvsmhttitle << htmhtrange.str().c_str() << ";MHT;#Delta #Phi_{min};";
		dphiminvsmhtname << jetcoll.at(i) << "_dphiminVsMht";
		

		if (i==0)
		{
			hist.hv_RecoEvt.h_Njet50eta2p5 = new TH1F(njet50eta2p5name.str().c_str(), njet50eta2p5title.str().c_str(), 20, 0, 20);
			hist.hv_RecoEvt.h_Njet30eta5p0 = new TH1F(njet30eta5p0name.str().c_str(), njet30eta5p0title.str().c_str(), 20, 0, 20);
			hist.hv_RecoEvt.h_Mht = new TH1F(mhtname.str().c_str(), mhttitle.str().c_str(), 500, 0, 1000);
			hist.hv_RecoEvt.h_Ht = new TH1F(htname.str().c_str(), httitle.str().c_str(), 300, 0, 3000);
			hist.hv_RecoEvt.h_Mht->Sumw2();
			hist.hv_RecoEvt.h_Ht->Sumw2();
			hist.hv_RecoEvt.h_DphiMin = new TH1F(dphiminname.str().c_str(), dphimintitle.str().c_str(), 125, 0, 2.5);
			hist.hv_RecoEvt.h_DphiMinVsMht = new TH2F(dphiminvsmhtname.str().c_str(), dphiminvsmhttitle.str().c_str(), 1400, 0, 350,250,0,2.5);
			hist.hv_RecoEvt.h_DphiMin->Sumw2();
			hist.hv_RecoEvt.h_DphiMinVsMht->Sumw2();
		} else if (i==1)
		{
			hist.hv_GenEvt.h_Njet50eta2p5 = new TH1F(njet50eta2p5name.str().c_str(), njet50eta2p5title.str().c_str(), 20, 0, 20);
			hist.hv_GenEvt.h_Njet30eta5p0 = new TH1F(njet30eta5p0name.str().c_str(), njet30eta5p0title.str().c_str(), 20, 0, 20);
			hist.hv_GenEvt.h_Mht = new TH1F(mhtname.str().c_str(), mhttitle.str().c_str(), 500, 0, 1000);
			hist.hv_GenEvt.h_Ht = new TH1F(htname.str().c_str(), httitle.str().c_str(), 300, 0, 3000);
			hist.hv_GenEvt.h_Mht->Sumw2();
			hist.hv_GenEvt.h_Ht->Sumw2();
			hist.hv_GenEvt.h_DphiMin = new TH1F(dphiminname.str().c_str(), dphimintitle.str().c_str(), 120, 0, 2.5);
			hist.hv_GenEvt.h_DphiMinVsMht = new TH2F(dphiminvsmhtname.str().c_str(), dphiminvsmhttitle.str().c_str(), 1400, 0, 350,250,0,2.5);
			hist.hv_GenEvt.h_DphiMin->Sumw2();
			hist.hv_GenEvt.h_DphiMinVsMht->Sumw2();
		} else if (i==2)
		{
			hist.hv_SmearedEvt.h_Njet50eta2p5 = new TH1F(njet50eta2p5name.str().c_str(), njet50eta2p5title.str().c_str(), 20, 0, 20);
			hist.hv_SmearedEvt.h_Njet30eta5p0 = new TH1F(njet30eta5p0name.str().c_str(), njet30eta5p0title.str().c_str(), 20, 0, 20);
			hist.hv_SmearedEvt.h_Mht = new TH1F(mhtname.str().c_str(), mhttitle.str().c_str(), 500, 0, 1000);
			hist.hv_SmearedEvt.h_Ht = new TH1F(htname.str().c_str(), httitle.str().c_str(), 300, 0, 3000);
			hist.hv_SmearedEvt.h_Mht->Sumw2();
			hist.hv_SmearedEvt.h_Ht->Sumw2();
			hist.hv_SmearedEvt.h_DphiMin = new TH1F(dphiminname.str().c_str(), dphimintitle.str().c_str(), 125, 0, 2.5);
			hist.hv_SmearedEvt.h_DphiMinVsMht = new TH2F(dphiminvsmhtname.str().c_str(), dphiminvsmhttitle.str().c_str(), 1400, 0, 350,250,0,2.5);
			hist.hv_SmearedEvt.h_DphiMin->Sumw2();
			hist.hv_SmearedEvt.h_DphiMinVsMht->Sumw2();
		}
	}

	GetJetHist(hist.hv_RecoJets,"reco", htmhtrange.str());
	GetJetHist(hist.hv_GenJets,"gen", htmhtrange.str());
	GetJetHist(hist.hv_SmearedJets,"smeared", htmhtrange.str());

	//for ratio plot
	const double evt_mht_max = 1500, evt_mht_bins = 750;
	const double evt_ht_max = 4000, evt_ht_bins = 800;
	const float npassFailHistBins = 14;
	const float passFailHistBins[] = {50,60,70,80,90,100,110,120,140,160,200,350,500,800,1000};
	//pass variations
	hist.hv_SmearedEvt.pass[0] = new TH1F("pass0","PASS from #Delta#Phi_{min} cut >0.15", npassFailHistBins, passFailHistBins); 
	hist.hv_SmearedEvt.pass[1] = new TH1F("pass1","PASS from #Delta#Phi_{min} cut >0.20", npassFailHistBins, passFailHistBins); 
	hist.hv_SmearedEvt.pass[2] = new TH1F("pass2","PASS from #Delta#Phi_{min} cut >0.25", npassFailHistBins, passFailHistBins); 
	hist.hv_SmearedEvt.pass[3] = new TH1F("pass3","PASS from #Delta#Phi_{min} cut >0.30", npassFailHistBins, passFailHistBins); 
	hist.hv_SmearedEvt.pass[4] = new TH1F("pass4","PASS from #Delta#Phi_{min} cut >0.35", npassFailHistBins, passFailHistBins); 
	hist.hv_SmearedEvt.pass[5] = new TH1F("pass5","PASS from #Delta#Phi_{min} cut >0.40", npassFailHistBins, passFailHistBins); 

	hist.hv_SmearedEvt.fail[0] = new TH1F("fail0","FAIL from #Delta#Phi_{min} cut <0.15", npassFailHistBins, passFailHistBins);
	hist.hv_SmearedEvt.fail[1] = new TH1F("fail1","FAIL from #Delta#Phi_{min} cut <0.20", npassFailHistBins, passFailHistBins);
	hist.hv_SmearedEvt.fail[2] = new TH1F("fail2","FAIL from #Delta#Phi_{min} cut <0.25", npassFailHistBins, passFailHistBins);
	hist.hv_SmearedEvt.fail[3] = new TH1F("fail3","FAIL from #Delta#Phi_{min} cut <0.30", npassFailHistBins, passFailHistBins);
	hist.hv_SmearedEvt.fail[4] = new TH1F("fail4","FAIL from #Delta#Phi_{min} cut <0.35", npassFailHistBins, passFailHistBins);
	hist.hv_SmearedEvt.fail[5] = new TH1F("fail5","FAIL from #Delta#Phi_{min} cut <0.40", npassFailHistBins, passFailHistBins);

	hist.hv_SmearedEvt.passFineBin[0] = new TH1F("passFineBin0","PASS from #Delta#Phi_{min} cut >0.15", evt_mht_bins, 0, evt_mht_max); 
	hist.hv_SmearedEvt.passFineBin[1] = new TH1F("passFineBin1","PASS from #Delta#Phi_{min} cut >0.20", evt_mht_bins, 0, evt_mht_max); 
	hist.hv_SmearedEvt.passFineBin[2] = new TH1F("passFineBin2","PASS from #Delta#Phi_{min} cut >0.25", evt_mht_bins, 0, evt_mht_max); 
	hist.hv_SmearedEvt.passFineBin[3] = new TH1F("passFineBin3","PASS from #Delta#Phi_{min} cut >0.30", evt_mht_bins, 0, evt_mht_max); 
	hist.hv_SmearedEvt.passFineBin[4] = new TH1F("passFineBin4","PASS from #Delta#Phi_{min} cut >0.35", evt_mht_bins, 0, evt_mht_max); 
	hist.hv_SmearedEvt.passFineBin[5] = new TH1F("passFineBin5","PASS from #Delta#Phi_{min} cut >0.40", evt_mht_bins, 0, evt_mht_max); 

	hist.hv_SmearedEvt.failFineBin[0] = new TH1F("failFineBin0","FAIL from #Delta#Phi_{min} cut <0.15", evt_mht_bins, 0, evt_mht_max);
	hist.hv_SmearedEvt.failFineBin[1] = new TH1F("failFineBin1","FAIL from #Delta#Phi_{min} cut <0.20", evt_mht_bins, 0, evt_mht_max);
	hist.hv_SmearedEvt.failFineBin[2] = new TH1F("failFineBin2","FAIL from #Delta#Phi_{min} cut <0.25", evt_mht_bins, 0, evt_mht_max);
	hist.hv_SmearedEvt.failFineBin[3] = new TH1F("failFineBin3","FAIL from #Delta#Phi_{min} cut <0.30", evt_mht_bins, 0, evt_mht_max);
	hist.hv_SmearedEvt.failFineBin[4] = new TH1F("failFineBin4","FAIL from #Delta#Phi_{min} cut <0.35", evt_mht_bins, 0, evt_mht_max);
	hist.hv_SmearedEvt.failFineBin[5] = new TH1F("failFineBin5","FAIL from #Delta#Phi_{min} cut <0.40", evt_mht_bins, 0, evt_mht_max);

	for (int i=0; i<6; ++i) 
	{ 
		hist.hv_SmearedEvt.pass[i]->Sumw2(); hist.hv_SmearedEvt.fail[i]->Sumw2(); 
		hist.hv_SmearedEvt.passFineBin[i]->Sumw2(); hist.hv_SmearedEvt.failFineBin[i]->Sumw2(); 
	}

	hist.hv_SmearedEvt.sidebandSyst[0] = new TH1F("sidebandSyst1"," Sideband Syst1: FAIL from dphi (j1 & j2 , mht) >0.5 and dphi (j3, mht)>0.1 ", npassFailHistBins, passFailHistBins);
	hist.hv_SmearedEvt.sidebandSyst[1] = new TH1F("sidebandSyst2"," Sideband Syst2: FAIL from dphi (j1 & j2 , mht) >0.5 and dphi (j3, mht)>0.1 ", npassFailHistBins, passFailHistBins);
	hist.hv_SmearedEvt.sidebandSystFineBin[0] = new TH1F("sidebandSyst1_fineBin","Sideband Syst1 (FineBinned) : FAIL from dphi (j1 & j2 , mht) >0.5 and dphi (j3, mht)>0.1 ", 1500, 0, 1500);
	hist.hv_SmearedEvt.sidebandSystFineBin[1] = new TH1F("sidebandSyst2_fineBin","Sideband Syst2 (FineBinned) : FAIL from dphi (j1 & j2 , mht) >0.4 and dphi (j3, mht)>0.2 ", 1500, 0, 1500);

	hist.hv_SmearedEvt.signal = new TH1F("signal" ,"Signal Region", npassFailHistBins, passFailHistBins);
	hist.hv_SmearedEvt.signalFineBin = new TH1F("signalFineBin" ,"Signal Region", 1500, 0, 1500);

	hist.hv_SmearedEvt.sidebandSyst[0]->Sumw2();
	hist.hv_SmearedEvt.sidebandSyst[1]->Sumw2();
	hist.hv_SmearedEvt.sidebandSystFineBin[0]->Sumw2();
	hist.hv_SmearedEvt.sidebandSystFineBin[1]->Sumw2();

	hist.hv_SmearedEvt.signal->Sumw2();
	hist.hv_SmearedEvt.signalFineBin->Sumw2();
}

void FactorizationBySmearing::GetJetHist(vector<JetHist_t>& Hist, const string jetcoll,
					const string htmhtrange)
{
	//save up to 10 leading jets
	const int njets = 10;
	for (int i=0; i<njets; ++i)
	{
		stringstream ptname, etaname, phiname, dphiname; 
		stringstream pttitle, etatitle, phititle, dphititle; 
		const int index = i + 1;
		ptname << jetcoll << "jet" << index << "_pt";
		etaname << jetcoll << "jet" << index << "_eta";
		phiname << jetcoll << "jet" << index << "_phi";
		dphiname << jetcoll << "jet" << index << "_dphi";
		pttitle << htmhtrange << ";Jet-" << index << " P_{T} [GeV];Events;";
		etatitle << htmhtrange << ";Jet-" << index << " #eta;Events;";
		pttitle << htmhtrange << ";Jet-" << index << " #phi ;Events;";
		pttitle << htmhtrange << ";#delta #phi [Jet-" << index << ", MHT];Events;";

		JetHist_t hist;
		hist.h_Jet_pt   = new TH1F(ptname.str().c_str(), pttitle.str().c_str(),   125, 0.0, 2500.0);  
		hist.h_Jet_eta  = new TH1F(etaname.str().c_str(), etatitle.str().c_str(),  120, -6.0, 6.0);  
		hist.h_Jet_phi  = new TH1F(phiname.str().c_str(), phititle.str().c_str(),  70, -3.5, 3.5);  
		hist.h_Jet_dphi  = new TH1F(dphiname.str().c_str(), dphititle.str().c_str(),  70, 0, 3.5);  
		hist.h_Jet_pt->Sumw2();
		hist.h_Jet_eta->Sumw2();
		hist.h_Jet_phi->Sumw2();
		hist.h_Jet_dphi->Sumw2();
		Hist.push_back(hist);
	}

}



#endif

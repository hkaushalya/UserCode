// -*- C++ -*-
//
// Package:    ZJets
// Class:      ZJets
// 
/**\class ZJets ZJets.cc UserCode/Factorization/src/ZJets.cc

 Description: To study the difference seen in dPhiMin between 
              data and MC. 

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  samantha hewamanage
//         Created:  Tue Jul 26 22:15:44 CDT 2011
// $Id$
//
//


// system include files
#include <memory>
#include <vector>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/print.h"
#include "FWCore/Utilities/interface/Verbosity.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/View.h"

//jet collections
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
//ROOT
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "assert.h"
#include "TVector3.h"
#include <iomanip>
#include "TPad.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

using namespace std;

//
// class declaration
//

class ZJets : public edm::EDFilter {
   public:
      explicit ZJets(const edm::ParameterSet&);
      ~ZJets();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

		struct RunLumiEvt_t
		{
			unsigned run;
			unsigned lumi;
			unsigned evt;
		};

		struct EventHist_t {
			TH1F* nvtx;
			TH1F* mht;
			TH1F* pfmht;
			TH1F* ht;	//calculated by hand
			TH1F* pfht; //ht in PAT
			TH1F* njet30;
			TH1F* njet50;
			TH1F* meff;
			TH1F* dphiMin;
			TH1F* delphiMin_jetmht[10]; // using 3 jets, 4 jets etc
			TH1F* prescaleWeights;
			TH1F* eventWeights;
			//TH2F* trigVsWeight;
		};

		struct JetHist_t{
			TH1F* pt;
			TH1F* eta;
			TH1F* phi;
			TH1F* delphi;  //delphi(jet,MHT)
		};

		struct Hist_t {  //hists for each inclusive jet category
			EventHist_t evt;
			JetHist_t jet[10];
			TH1F* pass[6];
			TH1F* passFineBin[6];
			TH1F* fail[6];
			TH1F* failFineBin[6];
			TH1F* signal;
			TH1F* signalFineBin;
			TH1F* sidebandSyst[2];
			TH1F* sidebandSystFineBin[2];
		};

		struct CutsPassedHists_t
		{
			TH1F *processed;
			TH1F *njet; 
			TH1F *ht;
			TH1F *mht;
		};

		struct matchingHists_t {
			TH1F* etRatio;
			TH1F* delR;
			TH1F* delEta;
			TH1F* etDiff;
		};


   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
		bool storeHLT(const edm::Event& e, const edm::EventSetup& iSetup);
		void PrintHeader();
		void BookHistograms();
		void BookCommonHistograms(TFileDirectory& dir, const float htMin, 
								const float htMax, Hist_t& hist);

		float DelPhiMin(const std::vector<pat::Jet>& jets, const TLorentzVector& mhtVec);
		void FillHistograms(
				edm::Handle<reco::VertexCollection> vertexHandle,
				const std::vector<pat::Jet>&  cleanJets,
				const unsigned& njet50eta25,
				const unsigned& njet30eta50,
				const TLorentzVector mhtVec,
				const float& ht,
				const float dPhiMin, 
				const bool bRA2DphiCut, 
				const bool bSidebandSyst1, 
				const bool bSidebandSyst2,
				const float prescaleWeight,
				const float Weight,
				Hist_t& hist
				);

		bool IsMuon(const pat::Jet& jet);
		void GetCleanJets(std::vector<pat::Jet>& vJets, 
				const std::vector<pat::Muon>&  vMuons, 
				edm::Handle<std::vector<pat::Jet> > jetHandle);
		void GetNjetHt(unsigned& njet50eta25, unsigned& njet30eta50,
				float& ht, const std::vector<pat::Jet>& vJets);
		void GetMht(TLorentzVector& mhtVec, const std::vector<pat::Muon>& vMuon);

		// ----------member data ---------------------------
		//jet collections
		edm::InputTag mhtInputTag_, htInputTag_;
		edm::Handle<edm::View<reco::MET> > mhtHandle;
		edm::Handle<double> htHandle;
		edm::InputTag prescaleWeightInputTag;

		int inMinVtx;
		unsigned int uProcessed; // number of processed events
		unsigned int uPassed; //number of events passed the filter
		RunLumiEvt_t kRunLumiEvent;
		int iVerbose; // control print levels
		double dMinHT;
		double dMinMHT, dMaxMHT;


		CutsPassedHists_t cutsHist;
		matchingHists_t matchedHists;

		EventHist_t evtHist;
		JetHist_t pf30_jet1Hist, pf30_jet2Hist, pf30_jet3Hist; //for upto 5 jets may be
		JetHist_t pf50eta25_jet1Hist, pf50eta25_jet2Hist, pf50eta25_jet3Hist; //for upto 5 jets may be

		edm::InputTag patJetsPFPt50Eta25InputTag_, patJetsPFPt30Eta50InputTag_, muonInputTag_;
		edm::Handle<std::vector<pat::Jet> > pfpt50eta25JetHandle, pfpt30eta50JetHandle;
		edm::Handle<std::vector<pat::Muon> >muonHandle;
		std::vector<pat::Muon> vGoodMuonList;
		std::vector<pat::Jet> vCleanJets;

		TH1F* MHT_by_phislice[6];
		TLorentzVector vMetVec;
		edm::Handle<std::vector<reco::PFMET> >pfMetHandle;
		unsigned uFailMinHTCut, uFailMinPFMHTCut;

		edm::LumiReWeighting LumiWeights_;
		std::vector<float> vMCNvtxDist, vDATANvtxDist;
		double sumLumiWeights; //sum of lumi weights for cross checks
		double Weight; //total weight for an event
		int doLumiWeighing, doEventWeighing;
		edm::Handle<double> prescaleWeightHandle;
		bool usePrescaleWeight;
		std::vector<std::string> triggerPathsToStore_;  // Vector to store list of HLT paths to store results of in ntuple
		edm::InputTag hlTriggerResults_;    // Input tag for TriggerResults
    	std::vector<double> htBins_;
		std::vector<Hist_t> vHist;
		std::vector<float> vDphiVariations;

};


static bool sort_using_less_than(double u, double v)
{
	   return u < v;
}



//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ZJets::ZJets(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
	triggerPathsToStore_ = iConfig.getParameter<std::vector<std::string> >("TriggerPathsToStore");
	hlTriggerResults_    = iConfig.getParameter<edm::InputTag>("HltTriggerResults");	
	patJetsPFPt50Eta25InputTag_ = iConfig.getParameter<edm::InputTag>("patJetsPFPt50Eta25InputTag");
	patJetsPFPt30Eta50InputTag_ = iConfig.getParameter<edm::InputTag>("patJetsPFPt30Eta50InputTag");
	mhtInputTag_    = iConfig.getParameter<edm::InputTag>("mhtInputTag");
	htInputTag_     = iConfig.getParameter<edm::InputTag>("htInputTag");
	//outputFile      = iConfig.getUntrackedParameter<double>("outputFile","Default.root");
	dMinHT          = iConfig.getUntrackedParameter<double>("dMinHT",0.0);
	dMinMHT         = iConfig.getUntrackedParameter<double>("dMinMHT",0.0);
	dMaxMHT         = iConfig.getUntrackedParameter<double>("dMaxMHT",99999.0);
	iVerbose        = iConfig.getUntrackedParameter<int>("verbose",0);
	doLumiWeighing  = iConfig.getUntrackedParameter<int>("ApplyLumiWeighing",0);
	doEventWeighing = iConfig.getUntrackedParameter<int>("ApplyEventWeighing",0);
	prescaleWeightInputTag = iConfig.getParameter<edm::InputTag>("prescaleWeight");
	usePrescaleWeight = iConfig.getUntrackedParameter<int>("usePrescaleWeight",0);
	uProcessed       = 0;
	uPassed          = 0;
	uFailMinHTCut    = 0;
	uFailMinPFMHTCut = 0;
	sumLumiWeights   = 0;
	Weight           = 1;
	htBins_ = iConfig.getParameter<std::vector<double > >("htBins");
	muonInputTag_    = iConfig.getParameter<edm::InputTag>("muonSrc");

	//difference variation of the dPhiMin selections.
	//make sure the book the correct number of histograms 
	//when these are changed!!!
	vDphiVariations.push_back(0.15);
	vDphiVariations.push_back(0.20);
	vDphiVariations.push_back(0.25);
	vDphiVariations.push_back(0.30);
	vDphiVariations.push_back(0.35);
	vDphiVariations.push_back(0.40);

	BookHistograms();

	//prescale weight with respective triggers
	/*const int nBins = (int) triggerPathsToStore_.size();
	hist_trigVsWeight = fs->make<TH2F> ("trigVsprescaleWeights","Trigger Vs Prescale Weights",nBins,0,nBins,2000,0,2000);
	int i = 1;
	for (std::vector<std::string>::const_iterator trigNameTempl = triggerPathsToStore_.begin();
			trigNameTempl != triggerPathsToStore_.end(); ++trigNameTempl) 
	{
		hist_trigVsWeight->GetXaxis()->SetBinLabel(i,(*trigNameTempl).c_str());
		++i;
	}*/

}

void ZJets::BookHistograms()
{
	edm::Service<TFileService> fs;
	assert (htBins_.size()>1 && "ZJets:: htBins_ size must be >1!");
	
	matchedHists.etRatio = fs->make<TH1F> ("matching_etRatio", "Matching Muon-jet;Pt Ratio (jet/muon);Events;", 100, 0, 10);
	matchedHists.delR    = fs->make<TH1F> ("matching_delR"   , "Matching Muon-jet;delta R ;Events;", 50, 0, 5);
	matchedHists.delEta  = fs->make<TH1F> ("matching_delEta" , "Matching Muon-jet;delta Eta ;Events;", 50, 0, 5);
	matchedHists.etDiff  = fs->make<TH1F> ("matching_etDiff" , "Matching Muon-jet;Et difference/Muon Et ;Events;", 50, -1, 1);

	for (unsigned i =0; i < htBins_.size()-1; ++i)
	{
		stringstream dirName;
		dirName << "HT" << htBins_.at(i) << "to" << htBins_.at(i+1);
		TFileDirectory subDir = fs->mkdir(dirName.str());
		Hist_t hist;
		vHist.push_back(hist);
		BookCommonHistograms(subDir, htBins_.at(i), htBins_.at(i+1), vHist.at(i));
	}

}
void ZJets::BookCommonHistograms(TFileDirectory& dir, const float htMin, 
								const float htMax, Hist_t& hist)
{

	//ht/mht/jet et,eta,phi,dphi-jet-mht, jet mass

	const double evt_mht_max = 1500, evt_mht_bins = 750;
	const double evt_ht_max = 4000, evt_ht_bins = 800;
	hist.evt.nvtx   = dir.make<TH1F> ("nvtx"  ,"RA2 Good Vertices; N Vtx;Events;", 40, 0, 40);
	hist.evt.pfmht  = dir.make<TH1F> ("pfmht","MHT from PFmetHandle;MHT [GeV];Events;", evt_mht_bins, 0, evt_mht_max);
	hist.evt.mht    = dir.make<TH1F> ("mht"  ,"Calculated MHT;MHT [GeV];Events;", evt_mht_bins, 0, evt_mht_max);
	hist.evt.ht     = dir.make<TH1F> ("ht"   ,"HT from pfHT Handle ;pfHT [GeV];Events;", evt_ht_bins, 0, evt_ht_max);
	hist.evt.pfht   = dir.make<TH1F> ("pfht" ,"HT from Jets ET>50 GeV && |#Eta|<2.4;HT [GeV];Events;", evt_ht_bins, 0, evt_ht_max);
	hist.evt.njet30 = dir.make<TH1F> ("njet30" ,"Njets (Et>30 GeV && | #Eta |<5.0;NJETS;Events;", 20, 0, 20);
	hist.evt.njet50 = dir.make<TH1F> ("njet50" ,"Njets (Et>50 GeV && | #Eta |<2.4;NJETS;Events;", 20, 0, 20);
	hist.evt.meff   = dir.make<TH1F> ("meff" ,"RA2:;MEff;Events;", 50, 0, 5000);
	hist.evt.prescaleWeights = dir.make<TH1F> ("prescaleWeights","Prescale Weights",2000,0,2000);
	hist.evt.eventWeights    = dir.make<TH1F> ("totalEventWeights","Total Event Weights",2000,0,2000);
	hist.evt.dphiMin = dir.make<TH1F> ("dphiMin"  ,"RA2: #Delta#Phi_{min}" , 150, 0, 3);

	hist.evt.pfmht->Sumw2();
	hist.evt.mht->Sumw2();
	hist.evt.ht->Sumw2();
	hist.evt.pfht->Sumw2();
	hist.evt.njet30->Sumw2();
	hist.evt.njet50->Sumw2();
	hist.evt.meff->Sumw2();
	hist.evt.prescaleWeights->Sumw2();
	hist.evt.eventWeights->Sumw2();

	const double pt_bins = 200, pt_max = 2000;
	for (int i =0; i < 10; ++i)
	{
		const int jetnum = i+1;
		stringstream pt_name, pt_title, eta_name, eta_title, phi_name, phi_title, dphi_name, dphi_title;
		pt_name << "jet"<< jetnum << "_pt";
		eta_name << "jet"<< jetnum << "_eta";
		phi_name << "jet"<< jetnum << "_phi";
		dphi_name << "jet"<< jetnum << "_dphi";
		pt_title << "Jet-"<< jetnum <<" pt;pt [GeV];Events";
		eta_title << "Jet-"<< jetnum <<" eta;eta;Events;";
		phi_title << "Jet-"<< jetnum <<" phi;phi;Events;";
		dphi_title << "Jet-"<< jetnum <<" #Delta#Phi; #Delta#Phi (jet, MHT) ;Events;";

		hist.jet[i].pt  = dir.make<TH1F> (pt_name.str().c_str(),pt_title.str().c_str(), pt_bins, 0, pt_max);
		hist.jet[i].eta = dir.make<TH1F> (eta_name.str().c_str(),eta_title.str().c_str(), 100, -5, 5);
		hist.jet[i].phi = dir.make<TH1F> (phi_name.str().c_str(),phi_title.str().c_str(), 160, -8, 8);
		hist.jet[i].delphi = dir.make<TH1F> (dphi_name.str().c_str(),dphi_title.str().c_str(), 40, 0, 4);

		hist.jet[i].pt->Sumw2();
		hist.jet[i].eta->Sumw2();
		hist.jet[i].phi->Sumw2();
		hist.jet[i].delphi->Sumw2();
	}
	
	//for ratio plot
	const float npassFailHistBins = 14;
	const float passFailHistBins[] = {50,60,70,80,90,100,110,120,140,160,200,350,500,800,1000};
	//pass variations
	hist.pass[0] = dir.make<TH1F> ("pass0","PASS from #Delta#Phi_{min} cut >0.15", npassFailHistBins, passFailHistBins); 
	hist.pass[1] = dir.make<TH1F> ("pass1","PASS from #Delta#Phi_{min} cut >0.20", npassFailHistBins, passFailHistBins); 
	hist.pass[2] = dir.make<TH1F> ("pass2","PASS from #Delta#Phi_{min} cut >0.25", npassFailHistBins, passFailHistBins); 
	hist.pass[3] = dir.make<TH1F> ("pass3","PASS from #Delta#Phi_{min} cut >0.30", npassFailHistBins, passFailHistBins); 
	hist.pass[4] = dir.make<TH1F> ("pass4","PASS from #Delta#Phi_{min} cut >0.35", npassFailHistBins, passFailHistBins); 
	hist.pass[5] = dir.make<TH1F> ("pass5","PASS from #Delta#Phi_{min} cut >0.40", npassFailHistBins, passFailHistBins); 

	hist.fail[0] = dir.make<TH1F> ("fail0","FAIL from #Delta#Phi_{min} cut <0.15", npassFailHistBins, passFailHistBins);
	hist.fail[1] = dir.make<TH1F> ("fail1","FAIL from #Delta#Phi_{min} cut <0.20", npassFailHistBins, passFailHistBins);
	hist.fail[2] = dir.make<TH1F> ("fail2","FAIL from #Delta#Phi_{min} cut <0.25", npassFailHistBins, passFailHistBins);
	hist.fail[3] = dir.make<TH1F> ("fail3","FAIL from #Delta#Phi_{min} cut <0.30", npassFailHistBins, passFailHistBins);
	hist.fail[4] = dir.make<TH1F> ("fail4","FAIL from #Delta#Phi_{min} cut <0.35", npassFailHistBins, passFailHistBins);
	hist.fail[5] = dir.make<TH1F> ("fail5","FAIL from #Delta#Phi_{min} cut <0.40", npassFailHistBins, passFailHistBins);

	hist.passFineBin[0] = dir.make<TH1F> ("passFineBin0","PASS from #Delta#Phi_{min} cut >0.15", evt_mht_bins, 0, evt_mht_max); 
	hist.passFineBin[1] = dir.make<TH1F> ("passFineBin1","PASS from #Delta#Phi_{min} cut >0.20", evt_mht_bins, 0, evt_mht_max); 
	hist.passFineBin[2] = dir.make<TH1F> ("passFineBin2","PASS from #Delta#Phi_{min} cut >0.25", evt_mht_bins, 0, evt_mht_max); 
	hist.passFineBin[3] = dir.make<TH1F> ("passFineBin3","PASS from #Delta#Phi_{min} cut >0.30", evt_mht_bins, 0, evt_mht_max); 
	hist.passFineBin[4] = dir.make<TH1F> ("passFineBin4","PASS from #Delta#Phi_{min} cut >0.35", evt_mht_bins, 0, evt_mht_max); 
	hist.passFineBin[5] = dir.make<TH1F> ("passFineBin5","PASS from #Delta#Phi_{min} cut >0.40", evt_mht_bins, 0, evt_mht_max); 

	hist.failFineBin[0] = dir.make<TH1F> ("failFineBin0","FAIL from #Delta#Phi_{min} cut <0.15", evt_mht_bins, 0, evt_mht_max);
	hist.failFineBin[1] = dir.make<TH1F> ("failFineBin1","FAIL from #Delta#Phi_{min} cut <0.20", evt_mht_bins, 0, evt_mht_max);
	hist.failFineBin[2] = dir.make<TH1F> ("failFineBin2","FAIL from #Delta#Phi_{min} cut <0.25", evt_mht_bins, 0, evt_mht_max);
	hist.failFineBin[3] = dir.make<TH1F> ("failFineBin3","FAIL from #Delta#Phi_{min} cut <0.30", evt_mht_bins, 0, evt_mht_max);
	hist.failFineBin[4] = dir.make<TH1F> ("failFineBin4","FAIL from #Delta#Phi_{min} cut <0.35", evt_mht_bins, 0, evt_mht_max);
	hist.failFineBin[5] = dir.make<TH1F> ("failFineBin5","FAIL from #Delta#Phi_{min} cut <0.40", evt_mht_bins, 0, evt_mht_max);

	for (int i=0; i<6; ++i) 
	{ 
		hist.pass[i]->Sumw2(); hist.fail[i]->Sumw2(); 
		hist.passFineBin[i]->Sumw2(); hist.failFineBin[i]->Sumw2(); 
	}

	hist.sidebandSyst[0] = dir.make<TH1F> ("sidebandSyst1"," Sideband Syst1: FAIL from dphi (j1 & j2 , mht) >0.5 and dphi (j3, mht)>0.1 ", npassFailHistBins, passFailHistBins);
	hist.sidebandSyst[1] = dir.make<TH1F> ("sidebandSyst2"," Sideband Syst2: FAIL from dphi (j1 & j2 , mht) >0.5 and dphi (j3, mht)>0.1 ", npassFailHistBins, passFailHistBins);
	hist.sidebandSystFineBin[0] = dir.make<TH1F> ("sidebandSyst1_fineBin","Sideband Syst1 (FineBinned) : FAIL from dphi (j1 & j2 , mht) >0.5 and dphi (j3, mht)>0.1 ", 1500, 0, 1500);
	hist.sidebandSystFineBin[1] = dir.make<TH1F> ("sidebandSyst2_fineBin","Sideband Syst2 (FineBinned) : FAIL from dphi (j1 & j2 , mht) >0.4 and dphi (j3, mht)>0.2 ", 1500, 0, 1500);

	hist.signal = dir.make<TH1F> ("signal" ,"Signal Region", npassFailHistBins, passFailHistBins);
	hist.signalFineBin = dir.make<TH1F> ("signalFineBin" ,"Signal Region", 1500, 0, 1500);

	hist.sidebandSyst[0]->Sumw2();
	hist.sidebandSyst[1]->Sumw2();
	hist.sidebandSystFineBin[0]->Sumw2();
	hist.sidebandSystFineBin[1]->Sumw2();

	hist.signal->Sumw2();
	hist.signalFineBin->Sumw2();
}


ZJets::~ZJets()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
ZJets::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
	++uProcessed;
	//cutsHist.processed->Fill(1);
	kRunLumiEvent.run  = iEvent.id().run();
	kRunLumiEvent.lumi = iEvent.id().luminosityBlock(); 
	kRunLumiEvent.evt  = iEvent.id().event();

	if (iVerbose) std::cout << __LINE__<< ":: Processing event: "
			<< kRunLumiEvent.run << ":" << kRunLumiEvent.lumi 
			<< ":" << kRunLumiEvent.evt << std::endl;

	//vertex cut is applied in the RA2Cleanning process.
	Handle<reco::VertexCollection> vertexHandle;
	iEvent.getByLabel("goodVerticesRA2", vertexHandle);
   if (! vertexHandle.isValid())
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":vertexHandle handle not found!" << std::endl;
		assert(false);
	}
	

	iEvent.getByLabel(muonInputTag_, muonHandle);
   if (! muonHandle.isValid())
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":muonHandle handle not found!" << std::endl;
		assert(false);
	}

	//Z->mumu selection for DATA.
	//require exactly 2 muons within kinematic regions
	//in anycae, extra muons should be vetoed as done in Znunu MC
	//then reconstruct Z and apply tight mass window cut
	//to reduce bakgrounds.
	
	//read in the objects using View
	//edm::Handle<edm::View<pat::Electron> > electrons;
	//iEvent.getByLabel(electronSrc_, electrons);
	// // check which ones to keep
	//  84   std::auto_ptr<std::vector<pat::Electron> > prod(new std::vector<pat::Electron>()); 
	//   85   if (vertices->size() > 0) {
	//    86     for (edm::View<pat::Electron>::const_iterator e = electrons->begin(); e != electrons->end(); ++e) {
	//     87       // acceptance cuts
	//      88       if (e->pt() < minElePt_) continue;
	//
	//prod->push_back(pat::Electron(*e));
	//// determine result before losing ownership of the pointer
	//140   bool result = (doEleVeto_ ? (prod->size() == 0) : true);
	//:w
	//


	//apply lepton veto
	
	//look for any 2 muons
	if ( muonHandle->size() < 2) return 0;
	
	vGoodMuonList.clear();
	for (std::vector<pat::Muon>::const_iterator m = muonHandle->begin();
			m != muonHandle->end(); ++m)
	{
		if (m->pt()<20 || fabs(m->eta())>5.0) continue;
//		cout << "muon = " << m->pt() << endl;
		vGoodMuonList.push_back(*m);	
	}

	if (vGoodMuonList.size()<2) return 0;

//	cout << "muons / goodmuons = " << muonHandle->size() << " / " << vGoodMuonList.size() << endl;

	/**********************************************************
	 * detemined different weights and generate total weight
	 **********************************************************/
	Weight = 1;

	/* PU S3 reweighting for QCD MC sample
	*/
	double lumiWeight = 1;
	if ( doLumiWeighing && ! iEvent.isRealData() )
	{
		lumiWeight = LumiWeights_.weight( iEvent );
		//std::cout << "lum wgt = " << lumiWeight << std::endl;
		sumLumiWeights += lumiWeight;
		Weight *= lumiWeight;
	}

	//event weights for flat QCD samples
	double storedWeight = 1;
	if ( doEventWeighing  && ! iEvent.isRealData() )
	{
		edm::Handle<GenEventInfoProduct> genEvtInfoHandle;
		iEvent.getByLabel("generator", genEvtInfoHandle);
		storedWeight = genEvtInfoHandle->weight();
		//std::cout << "storedWeight = " << storedWeight << std::endl;
		Weight *= storedWeight;
	}

	/* 
	 * prescaled trigger weights
	 */
	double prescaleWeight = 1;
	if ( iEvent.isRealData() && usePrescaleWeight)
	{
		iEvent.getByLabel(prescaleWeightInputTag, prescaleWeightHandle);
		if ( ! prescaleWeightHandle.isValid()) 
		{ 
			std::cout << "prescaleWeightHandle found!" << std::endl;	
			assert (false);
		}
		prescaleWeight = (*prescaleWeightHandle);
		Weight *= prescaleWeight;
		//std::cout << "Weight = " << Weight << std::endl;
	}

	iEvent.getByLabel(patJetsPFPt30Eta50InputTag_, pfpt30eta50JetHandle);
	if (! pfpt30eta50JetHandle.isValid()) 
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":pfpt30eta50JetHandle handle not found!" << std::endl;
		assert(false);
	}

	GetCleanJets(vCleanJets, vGoodMuonList, pfpt30eta50JetHandle);

	unsigned njet50eta25 = 0, njet30eta50 = 0;
	float ht = 0;
	TLorentzVector mhtVec(0,0,0,0);

	GetMht(mhtVec, vGoodMuonList);
	GetNjetHt(njet50eta25, njet30eta50, ht, vCleanJets);

	//cout << ">>>> njet50/njet30/mht/ht = " << njet50eta25 << "/" << njet30eta50 << "/" << mhtVec.Pt() << "/" << ht << endl;
	cout << ">>>> njet50/ht/mht = " << njet50eta25 << " / " << ht << " / " << mhtVec.Pt() << endl;
	//APPLY RA2 cuts
	if (njet50eta25 < 3) 
	{ 
		cout << "fail njet " << endl; 
		return 0; 
	}
	if (ht < dMinHT) 
	{ 
		cout << "fail ht " << endl; 
		++uFailMinHTCut; return 0; 
	}
	if (mhtVec.Pt() < dMinMHT && mhtVec.Pt() > dMaxMHT) 
	{ 
		cout << "fail mht " << endl; 
		++uFailMinPFMHTCut; 
		return 0; 
	}

	cout << ">>>>>> pass all cuts! " << endl;
	const float dPhiMin = DelPhiMin(vCleanJets,	mhtVec);

	//Check if RA2 dphi cuts are satisfied
	bool bPassRA2dphiCut = true;
	bPassRA2dphiCut = bPassRA2dphiCut && (std::abs(reco::deltaPhi(vCleanJets.at(0).phi(), mhtVec.Phi())) > 0.5);
	bPassRA2dphiCut = bPassRA2dphiCut && (std::abs(reco::deltaPhi(vCleanJets.at(1).phi(), mhtVec.Phi())) > 0.5);
	bPassRA2dphiCut = bPassRA2dphiCut && (std::abs(reco::deltaPhi(vCleanJets.at(2).phi(), mhtVec.Phi())) > 0.3);


	/**************************************************/
	//                 for systematics
	/**************************************************/
	bool bSidebandSyst1 = true;
	bSidebandSyst1 = bSidebandSyst1 && (std::abs(reco::deltaPhi(vCleanJets.at(0).phi(), mhtVec.Phi())) > 0.5);
	bSidebandSyst1 = bSidebandSyst1 && (std::abs(reco::deltaPhi(vCleanJets.at(1).phi(), mhtVec.Phi())) > 0.5);
	bSidebandSyst1 = bSidebandSyst1 && (std::abs(reco::deltaPhi(vCleanJets.at(2).phi(), mhtVec.Phi())) > 0.1);

	bool bSidebandSyst2 = true;
	bSidebandSyst2 = bSidebandSyst2 && (std::abs(reco::deltaPhi(vCleanJets.at(0).phi(), mhtVec.Phi())) > 0.4);
	bSidebandSyst2 = bSidebandSyst2 && (std::abs(reco::deltaPhi(vCleanJets.at(1).phi(), mhtVec.Phi())) > 0.4);
	bSidebandSyst2 = bSidebandSyst2 && (std::abs(reco::deltaPhi(vCleanJets.at(2).phi(), mhtVec.Phi())) > 0.2);


	//loop over all ht bins
	for (unsigned i =0; i < htBins_.size()-1; ++i)
	{
		//stringstream dirName;
		//dirName << "HT" << htBins_.at(i) << "to" << htBins_.at(i+1);
	
		if (ht >= htBins_.at(i) && ht < htBins_.at(i+1)) 
		{
			FillHistograms(
				vertexHandle,
				vCleanJets,
				njet50eta25,
				njet30eta50,
				mhtVec,
				ht, 
				dPhiMin, 
				bPassRA2dphiCut, 
				bSidebandSyst1, 
				bSidebandSyst2,
				prescaleWeight,
				Weight,
				vHist.at(i)
				);
		}
	}

	++uPassed;
   return true;
}

// ------------ method called once each job just before starting event loop  ------------
void 
ZJets::beginJob()
{
	//2011 Pileup Scenarios https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupInformation#2011_Pileup_Scenarios
	//"Flat to 10 plus tail" scenario shown above, the relative normalization of each bin is: 
	vMCNvtxDist.push_back(0.069286816); //1
	vMCNvtxDist.push_back(0.069286816); //2
	vMCNvtxDist.push_back(0.069286816); //3
	vMCNvtxDist.push_back(0.069286816); //4
	vMCNvtxDist.push_back(0.069286816); //5
	vMCNvtxDist.push_back(0.069286816); //6
	vMCNvtxDist.push_back(0.069286816); //7
	vMCNvtxDist.push_back(0.069286816); //8
	vMCNvtxDist.push_back(0.069286816); //9
	vMCNvtxDist.push_back(0.069286816); //10
	vMCNvtxDist.push_back(0.069286816); //11
	vMCNvtxDist.push_back(0.06518604 ); //12
	vMCNvtxDist.push_back(0.053861878); //13
	vMCNvtxDist.push_back(0.040782032); //14
	vMCNvtxDist.push_back(0.030135062); //15
	vMCNvtxDist.push_back(0.019550796); //16
	vMCNvtxDist.push_back(0.012264707); //17
	vMCNvtxDist.push_back(0.007449117); //18
	vMCNvtxDist.push_back(0.004502075); //19
	vMCNvtxDist.push_back(0.002194605); //20
	vMCNvtxDist.push_back(0.001166276); //21
	vMCNvtxDist.push_back(0.000476543); //22
	vMCNvtxDist.push_back(0.000188109); //23
	vMCNvtxDist.push_back(7.52436E-05); //24
	vMCNvtxDist.push_back(1.25406E-05); //25

	//DATA Nvtx distribution from Pileup_2011_EPS_8_jul.root
	vDATANvtxDist.push_back(1.45417e+07); //1
	vDATANvtxDist.push_back(3.47743e+07); //2
	vDATANvtxDist.push_back(7.89247e+07); //3
	vDATANvtxDist.push_back(1.26467e+08); //4
	vDATANvtxDist.push_back(1.59329e+08); //5
	vDATANvtxDist.push_back(1.67603e+08); //6
	vDATANvtxDist.push_back(1.52684e+08); //7
	vDATANvtxDist.push_back(1.23794e+08); //8
	vDATANvtxDist.push_back(9.09462e+07); //9
	vDATANvtxDist.push_back(6.13973e+07); //10
	vDATANvtxDist.push_back(3.8505e+07); //11
	vDATANvtxDist.push_back(2.2628e+07); //12
	vDATANvtxDist.push_back(1.25503e+07); //13
	vDATANvtxDist.push_back(6.61051e+06); //14
	vDATANvtxDist.push_back(3.32403e+06); //15
	vDATANvtxDist.push_back(1.60286e+06); //16
	vDATANvtxDist.push_back(743920); //17
	vDATANvtxDist.push_back(333477); //18
	vDATANvtxDist.push_back(144861); //19
	vDATANvtxDist.push_back(61112.7); //20
	vDATANvtxDist.push_back(25110.2); //21
	vDATANvtxDist.push_back(10065.1); //22
	vDATANvtxDist.push_back(3943.98); //23
	vDATANvtxDist.push_back(1513.54); //24
	vDATANvtxDist.push_back(896.161); //25

	LumiWeights_ = edm::LumiReWeighting(vMCNvtxDist, vDATANvtxDist);


}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZJets::endJob() {

	std::cout << __FILE__ << "::" << __FUNCTION__  << "\n" << std::endl;
	std::cout << "[ZJT:00] Job settings " << std::endl;
	std::cout << "[ZJT:01] Events Processed --- = " << uProcessed << std::endl;
	std::cout << "[ZJT:02] Events Passed ------ = " << uPassed << std::endl;
	std::cout << "[ZJT:03] Lumi Weighing? ----- = " << doLumiWeighing << std::endl;
	std::cout << "[ZJT:04] Event Weighing? ---- = " << doEventWeighing << std::endl;
	std::cout << "[ZJT:05] Prescale Weighing? - = " << usePrescaleWeight << std::endl;
	std::cout << "[ZJT:31] Minimum HT --------- = " << dMinHT << std::endl;
	std::cout << "[ZJT:32] MHT min/max -------- = " << dMinMHT << " - " << dMaxMHT << std::endl;
	std::cout << "[ZJT:40] PASS Summary --------- " << std::endl;
	std::cout << "[ZJT:43] Pass ht ------------ = " << (uProcessed - uFailMinHTCut) << std::endl;
	std::cout << "[ZJT:44] Pass mht ----------- = " << (uProcessed - uFailMinHTCut - uFailMinPFMHTCut)  << std::endl;
	std::cout << "[ZJT:50] LumiWeights Avg ---- = ";
	if ( uPassed>0) std::cout << sumLumiWeights/(double)uPassed << std::endl;
	else std::cout << "0" << std::endl;

	for (std::vector<std::string>::const_iterator trigNameTempl = triggerPathsToStore_.begin();
			trigNameTempl != triggerPathsToStore_.end(); ++trigNameTempl) 
	{
		std::cout << "[" << *trigNameTempl << "]";
	}
	std::cout << std::endl;

}

// ------------ method called when starting to processes a run  ------------
bool 
ZJets::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
ZJets::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
ZJets::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
ZJets::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZJets::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

float ZJets::DelPhiMin(const std::vector<pat::Jet>& cleanJets, 
				const TLorentzVector& mhtVec
				)
{
	PrintHeader();
	std::vector<float> vDelPhi_jetmht;
	const float mht = mhtVec.Pt();

	for (unsigned i = 0 ; i < cleanJets.size() ; ++i)
	{
		if (i>2) break; //use only three leading jets
		//std::cout << "i = " <<  i << std::endl;
		const float delphi_jetmht = fabs(TVector2::Phi_mpi_pi(cleanJets.at(i).phi() - mhtVec.Phi()));
		vDelPhi_jetmht.push_back(delphi_jetmht);
	}

	std::sort(vDelPhi_jetmht.begin(), vDelPhi_jetmht.end(), sort_using_less_than);	

	assert (vDelPhi_jetmht.size() == 3 && "ERROR: More than 3 dPhiMin calculations found!");
	std::cout << "vDelPhi_jetmht size = " << vDelPhi_jetmht.size() << std::endl;
	
	const float dPhiMin = vDelPhi_jetmht.at(0);
	return dPhiMin;
}

void ZJets::FillHistograms(
				edm::Handle<reco::VertexCollection> vertexHandle,
				const std::vector<pat::Jet>&  cleanJets,
				const unsigned& njet50eta25,
				const unsigned& njet30eta50,
				const TLorentzVector mhtVec,
				const float& ht,
				const float dPhiMin, 
				const bool bRA2DphiCut, 
				const bool bSidebandSyst1, 
				const bool bSidebandSyst2,
				const float prescaleWeight,
				const float Weight,
				Hist_t& hist
				)
{

	hist.evt.prescaleWeights->Fill(prescaleWeight);
	hist.evt.eventWeights->Fill(Weight);

	hist.evt.nvtx->Fill(vertexHandle->size(), Weight);

	const float mht = mhtVec.Pt();

	hist.evt.njet50->Fill(njet50eta25, Weight);
	hist.evt.njet30->Fill(njet30eta50, Weight);

	hist.evt.pfmht->Fill(mht, Weight);
	hist.evt.pfht->Fill(ht, Weight);
	hist.evt.meff->Fill(mht+ht, Weight);
	hist.evt.mht->Fill(mht, Weight);
	hist.evt.ht->Fill(ht, Weight);
	
	const int maxJets =  (njet30eta50>10) ? 10 : njet30eta50;
	for (int i=0; i < maxJets; ++i)
	{
		const float pt  = cleanJets.at(i).pt();
		const float phi = cleanJets.at(i).phi();
		const float eta = cleanJets.at(i).eta();
		const float dphi= fabs(TVector2::Phi_mpi_pi(phi - mhtVec.Phi()));

		hist.jet[i].pt->Fill(pt, Weight);
		hist.jet[i].phi->Fill(phi, Weight);
		hist.jet[i].eta->Fill(eta, Weight);
		hist.jet[i].delphi->Fill(dphi, Weight);
	}

	hist.evt.dphiMin->Fill(dPhiMin, Weight);
	if (bRA2DphiCut) 
	{
		hist.signal->Fill(mht,Weight);
		hist.signalFineBin->Fill(mht,Weight);
	}

	for (unsigned i=0; i < vDphiVariations.size(); ++i)
	{
		if (dPhiMin > vDphiVariations.at(i)) 
		{
			hist.pass[i]->Fill(mht, Weight);
			hist.passFineBin[i]->Fill(mht, Weight);
		} else 
		{
			hist.fail[i]->Fill(mht, Weight);
			hist.failFineBin[i]->Fill(mht, Weight);
		}
	}


	if (! bSidebandSyst1) //we want the denominator (or failed events)
	{
		hist.sidebandSyst[0]->Fill(mht, Weight);
		hist.sidebandSystFineBin[0]->Fill(mht, Weight);
	}
	if (! bSidebandSyst2)
	{
		hist.sidebandSyst[1]->Fill(mht, Weight);
		hist.sidebandSystFineBin[1]->Fill(mht, Weight);
	}
}


void ZJets::PrintHeader()
{
	std::cout  << " >>>>>>>>>>>>>>>  Run:Lumi:Event: "
			<< kRunLumiEvent.run << ":" << kRunLumiEvent.lumi 
			<< ":" << kRunLumiEvent.evt 
			<< "<<<<<<<<<<<<<<" << std::endl;
}

bool ZJets::storeHLT(const edm::Event& e, const edm::EventSetup& iSetup){
	//////////////////////////////////////////////////////////////////////
	////  Trigger Section: Analyzing HLT Trigger Results (TriggerResults) // 
	////////////////////////////////////////////////////////////////////////
	using namespace edm;

	bool accept(kFALSE);      
	int ntrigs(0);
	// get hold of TriggerResults
	Handle<TriggerResults> TrgResultsHandle;
	e.getByLabel(hlTriggerResults_, TrgResultsHandle);

	if (TrgResultsHandle.isValid()) {
		const TriggerResults *hltResults = TrgResultsHandle.product();
		const TriggerNames TrgNames = e.triggerNames(*hltResults);
		ntrigs=TrgNames.size();
		//std::cout << "%HLTInfo --  Number of HLT Triggers: " << ntrigs << std::endl;

		for (int i=0;i < ntrigs; ++i)
		{
			const std::string trigName = TrgNames.triggerName(i);
			std::cout << "Triger [" << i << "] " << trigName << std::endl;
		}

		for ( std::vector<std::string>::const_iterator trigNameTempl = triggerPathsToStore_.begin();  trigNameTempl != triggerPathsToStore_.end(); ++trigNameTempl) {
			int index = TrgNames.triggerIndex(*trigNameTempl);
			std::cout << "my trig: " << *trigNameTempl << " & index = " << index << std::endl;
			if (index !=  ntrigs )
			{
				accept = TrgResultsHandle->accept(index);
				if (accept) 
				{
					std::cout << "my trig: " << *trigNameTempl << " & accept = " << accept << std::endl;
					//return accept;
				}
			}
		}


		// add prescale information
/*		for ( std::vector<std::string>::const_iterator trigNameTempl = triggerPrescaleToStore_.begin();  trigNameTempl != triggerPrescaleToStore_.end(); ++trigNameTempl) {
			int index = TrgNames.triggerIndex(*trigNameTempl);
			int hlt_prescale = 0;
			if (index !=  ntrigs )
			{
				hlt_prescale = hltConfig_.prescaleValue(e, iSetup, *trigNameTempl);
				std::cout << __LINE__ << "hlt_prescale = " << hlt_prescale << std::endl;
			}
		}
*/		

	} else { std::cout << "%HLTInfo -- No Trigger Result" << std::endl;}

	return accept;
}

bool ZJets::IsMuon(const pat::Jet& jet)
{
	const TLorentzVector jVec(jet.px(), jet.py(), jet.pz(), jet.energy());
	int nMatch = 0;
	for (vector<pat::Muon>::const_iterator m = vGoodMuonList.begin();
			m != vGoodMuonList.end(); ++m)
	{
		TLorentzVector mVec(m->px(), m->py(),m->pz(), m->energy());
		if (mVec.DeltaR(jVec) >0.6) continue;
		else 
		{
			//cout << "Match found for jet with muon!" << endl;
			++nMatch;
			const float jetPt = jVec.Pt();
			const float muPt = mVec.Pt();
			matchedHists.etRatio->Fill(jetPt/muPt);
			matchedHists.etDiff->Fill((jetPt-muPt)/muPt);
			matchedHists.delR->Fill(mVec.DeltaR(jVec));
			matchedHists.delEta->Fill(fabs(jVec.Eta()-mVec.Eta()));
		}
	}

	if (nMatch>1) 
	{
		cout << __FILE__ << ":" << __FUNCTION__ << ":" 
			<< "Warning !! more than one jet matched to a muon!";
		PrintHeader();
	} //else cout << "No Match found for jet with muon!" << endl;


	return (bool) nMatch;
}

void ZJets::GetCleanJets(std::vector<pat::Jet>& vJets, 
				const std::vector<pat::Muon>&  vMuons, 
				edm::Handle<std::vector<pat::Jet> > jetHandle)
{
	vJets.clear();

	for (unsigned i=0; i<jetHandle->size(); ++i)
	{
		if (IsMuon((*jetHandle)[i])) continue;
		vJets.push_back((*jetHandle)[i]);
	}

}

void ZJets::GetNjetHt(unsigned& njet50eta25, unsigned& njet30eta50, 
				float& ht, const std::vector<pat::Jet>& vJets)
{

	for (unsigned i = 0; i < vJets.size(); ++i)
	{
		const TLorentzVector iJetVec(
				vJets.at(i).px(), vJets.at(i).py(),
				vJets.at(i).pz(), vJets.at(i).energy()
				);

		if (vJets.at(i).pt()>30.0 && fabs(vJets.at(i).eta() < 5.0) )
		{
			++njet30eta50;
		}
		if (vJets.at(i).pt()>50.0 && fabs(vJets.at(i).eta() < 2.5) )
		{
			++njet50eta25;
			ht += vJets.at(i).pt();
		}
	}
}
void ZJets::GetMht(TLorentzVector& mhtVec, const std::vector<pat::Muon>& vMuon)
{
	for (unsigned i = 0; i < vMuon.size(); ++i)
	{
		const TLorentzVector iMuonVec(
				vMuon.at(i).px(), vMuon.at(i).py(),
				vMuon.at(i).pz(), vMuon.at(i).energy()
				);

			mhtVec += iMuonVec;
	}
}



//define this as a plug-in
DEFINE_FWK_MODULE(ZJets);

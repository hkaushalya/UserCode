// -*- C++ -*-
//
// Package:    RA2bQCDvetoAna
// Class:      RA2bQCDvetoAna
// 
/**\class RA2bQCDvetoAna RA2bQCDvetoAna.cc UserCode/RA2bQCDvetoAna/src/RA2bQCDvetoAna.cc

 Description: Studies the DelPhi(jet,MET) for QCD veto.

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  samantha hewamanage
//         Created:  Tue Jul 26 22:15:44 CDT 2011
// $Id: RA2bQCDvetoAna.cc,v 1.2 2011/10/11 15:19:36 samantha Exp $
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

//jet collections
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "DataFormats/Math/interface/deltaPhi.h"

//ROOT
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "assert.h"
#include "TVector3.h"
#include <iomanip>
#include "TPad.h"

//
// class declaration
//

class RA2bQCDvetoAna : public edm::EDFilter {
   public:
      explicit RA2bQCDvetoAna(const edm::ParameterSet&);
      ~RA2bQCDvetoAna();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------
		//jet collections
		edm::InputTag caloJetInputTag_;
		edm::Handle<reco::CaloJetCollection> caloJetcollection;
		edm::InputTag pfJetInputTag_;
		edm::Handle<reco::PFJetCollection> pfJetcollection;
		edm::InputTag mhtInputTag_, htInputTag_;
		edm::Handle<edm::View<reco::MET> > mhtHandle;
		edm::Handle<double> htHandle;

		struct RunLumiEvt_t
		{
			unsigned run;
			unsigned lumi;
			unsigned evt;
		};

		
		int inMinVtx;
		double dMinMet;
		unsigned int uProcessed; // number of processed events
		unsigned int uPassed; //number of events passed the filter
		int inMinNdofVtx; //minimum number of dof for the vtx
		RunLumiEvt_t kRunLumiEvent;
		int iVerbose; // control print levels
		double dMaxPrimVtxZ;
		double dMaxPrimVtxRho;
		double dMinHT;
		double dMinMHT;

		struct EventHist_t {
			TH1F* nvtx;
			TH1F* vtxz;
			TH1F* mht;
			TH1F* ht;
			TH1F* njet;
			TH1F* metphi;
			TH1F* met;
			TH1F* meff;
		};

		struct JetHist_t{
			TH1F* pt;
			TH1F* eta;
			TH1F* phi;
			TH1F* delphi;  //delphi(jet,MHT)
			//TH1F* delT;
			TH1F* delTDevidedByJetPt; //delT/Jet Pt
			TProfile* delTvsJetPt;
			TProfile* delTDevidedByJetPtvsJetPt; //delT/Jet Pt
		};

		struct Hist_t {  //hists for each inclusive jet category
			EventHist_t evt;
			JetHist_t jet[5];
			TH1F* delphiMin_jetmet[3];
			TProfile* delphiMin_jetmetVsMHT[3];
			TH1F* delphiMin_jetjet[3];
			TProfile* delphiMin_jetjetVsMHT[3];
		};

		TH1F* hDelPhiMin[5];
		TH1F* hDelPhiMinNorm[5];
		TH1F* hDelPhiMinNorm_mht[8];
		TProfile* hDelPhiMinVsMET;
	//	TH2F* hDelPhiMinNormVshDelPhiMin;
		TProfile* hDelPhiMinNormVsMET;
		TH1F* hPass;
		TH1F* hFail;
		TH1F* hPass_Norm[10];
		TH1F* hFail_Norm[10];
		TH1F* hPass_Norm_mht[10];
		TH1F* hFail_Norm_mht[10];

		//TH1F* hist_delphiMin_jetmet[5];
		//TH2F* hist_delphiMin_jetmetVsMHT[5]; //0=all jets, 1=only 3 jets, 2=upto 4 jets 
		//TH1F* hist_delphiMin_jetmetVsHT[5];

		EventHist_t evtHist;
		JetHist_t pf30_jet1Hist, pf30_jet2Hist, pf30_jet3Hist; //for upto 5 jets may be
		JetHist_t pf50eta25_jet1Hist, pf50eta25_jet2Hist, pf50eta25_jet3Hist; //for upto 5 jets may be

		edm::InputTag patJetsPFPt30InputTag_;
		edm::InputTag patJetsPFPt50Eta25InputTag_;
		edm::Handle<std::vector<pat::Jet> > jetHandle;
		edm::Handle<std::vector<pat::Jet> > pfpt50eta25JetHandle;
		edm::Handle<std::vector<pat::Jet> > pfpt30JetHandle;
		edm::Handle<std::vector<pat::Jet> > alljetHandle;

		Hist_t Hist[8]; //all hists for various mht ranges
		TH1F* MHT_by_phislice[6];
		void BookHistograms(edm::Service<TFileService>& fs, Hist_t& hist
						, const float mHt_min, const float mHt_max);
		void FillHistograms(Hist_t& hist
					, edm::Handle<edm::View<reco::MET> > mhtHandle 
				   , edm::Handle<std::vector<pat::Jet> > jetHandle);
		unsigned minDelPhi(std::vector<float>& vDelPhi_jetjet
				, std::vector<float>& vDelPhi_jetmet
				, const unsigned njets
				, edm::Handle<std::vector<pat::Jet> > jetHandle
				, edm::Handle<edm::View<reco::MET> > mhtHandle); 
			//void DoDelMinStudy(edm::Handle<edm::View<reco::MET> > mhtHandle 
			//	, edm::Handle<std::vector<pat::Jet> > pt30jetHandle);
			void DoDelMinStudy(edm::Handle<std::vector<reco::PFMET> >pfMetHandle
				, edm::Handle<std::vector<pat::Jet> > pt30jetHandle
				, const std::vector<unsigned> vLead3JetIndices
				, edm::Handle<edm::View<reco::MET> > mhtHandle
				, edm::Event& iEvent);
		void PrintHeader();
		TLorentzVector vMetVec;
		edm::Handle<std::vector<reco::PFMET> >pfMetHandle;
		unsigned uFailNvtxCut, uFailNjet50Eta24Cut, uFailMinHTCut, uFailMinPFMetCut;
		unsigned uFailMinPFMHTCut;

		edm::LumiReWeighting LumiWeights_;
		std::vector<float> vMCNvtxDist, vDATANvtxDist;
		double sumLumiWeights; //sum of lumi weights for cross checks
		float lumiWeight;

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
RA2bQCDvetoAna::RA2bQCDvetoAna(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
	//triggerPathsToStore_ = iConfig.getParameter<std::vector<std::string> >("TriggerPathsToStore");
	//hlTriggerResults_    = iConfig.getParameter<edm::InputTag>("HltTriggerResults");	
	//req_trigger = iConfig.getUntrackedParameter<bool>("req_trigger",true);
	patJetsPFPt30InputTag_ = iConfig.getParameter<edm::InputTag>("patJetsPFPt30InputTag");
	patJetsPFPt50Eta25InputTag_ = iConfig.getParameter<edm::InputTag>("patJetsPFPt50Eta25InputTag");
	mhtInputTag_ = iConfig.getParameter<edm::InputTag>("mhtInputTag");
	htInputTag_ = iConfig.getParameter<edm::InputTag>("htInputTag");
	inMinVtx = iConfig.getUntrackedParameter<int>("nMinVtx",1);
	dMinMet = iConfig.getUntrackedParameter<double>("minMet",0.0);
	dMinHT = iConfig.getUntrackedParameter<double>("dMinHT",350.0);
	dMinMHT = iConfig.getUntrackedParameter<double>("dMinMHT",0.0);
	inMinNdofVtx = iConfig.getUntrackedParameter<int>("ndofVtx",0);
	dMaxPrimVtxZ = iConfig.getUntrackedParameter<double>("maxDelzVtx",24.0);
	dMaxPrimVtxRho = iConfig.getUntrackedParameter<double>("maxDelRho",2.0);
	//caloJetInputTag_ = iConfig.getParameter<edm::InputTag>("caloJetInputTag_");
	//pfJetInputTag_ = iConfig.getParameter<edm::InputTag>("pfJetInputTag_");
	iVerbose = iConfig.getUntrackedParameter<int>("verbose",0);
	uProcessed = 0;
	uPassed = 0;
	uFailNvtxCut = 0;
	uFailNjet50Eta24Cut = 0;
	uFailMinHTCut = 0;
	uFailMinPFMetCut = 0;
	uFailMinPFMHTCut = 0;
	lumiWeight = 1; //event-by-event lumi weight
	sumLumiWeights = 0;


	//generate hists
	edm::Service<TFileService> fs;
	/*BookHistograms(fs, Hist[0], 60, 80); 
	BookHistograms(fs, Hist[1], 80,100); 
	BookHistograms(fs, Hist[2],100,120); 
	BookHistograms(fs, Hist[3],120,140); 
	BookHistograms(fs, Hist[4],140,170); 
	BookHistograms(fs, Hist[5],170,200); 
	BookHistograms(fs, Hist[6],200,7000); 
*/
	const float met_min = 0, met_max=500, met_bins=100;
	MHT_by_phislice[0] = fs->make<TH1F> ("mht_phislice_lt0.1" ,"MHT (#Delta#Phi_{min}<0.1);MHT [GeV];Events;", met_bins, 0, met_max);
	MHT_by_phislice[1] = fs->make<TH1F> ("mht_phislice_lt0.2" ,"MHT (0.1<#Delta#Phi_{min}<0.2);MHT [GeV];Events;", met_bins, 0, met_max);
	MHT_by_phislice[2] = fs->make<TH1F> ("mht_phislice_lt0.3" ,"MHT (0.2<#Delta#Phi_{min}<0.3);MHT [GeV];Events;", met_bins, 0, met_max);
	MHT_by_phislice[3] = fs->make<TH1F> ("mht_phislice_lt0.5" ,"MHT (0.3<#Delta#Phi_{min}<0.5);MHT [GeV];Events;", met_bins, 0, met_max);
	MHT_by_phislice[4] = fs->make<TH1F> ("mht_phislice_lt0.8" ,"MHT (0.5<#Delta#Phi_{min}<0.8);MHT [GeV];Events;", met_bins, 0, met_max);
	MHT_by_phislice[5] = fs->make<TH1F> ("mht_phislice_gt0.8" ,"MHT (#Delta#Phi_{min}>0.8);MHT [GeV];Events;", met_bins, 0, met_max);

	//these are general event hist to check the PAT tuples cuts
	const double evt_met_max = 800, evt_met_bins = 400;
	const double evt_ht_max = 4000, evt_ht_bins = 80;
	/*
	 * All plots with ra2b_* prefix should remain as it is. They are all made with MET
	 * and using std. RA2b selection.
	 */
	evtHist.nvtx = fs->make<TH1F> ("ra2b_nvtx" ,"RA2b: Nvtx;Nvtx;Events;", 15, 0, 15);
	evtHist.vtxz = fs->make<TH1F> ("ra2b_vtxz","RA2b:;Primary Vertex z [cm];Arbitrary;",30,0,30);
	evtHist.met = fs->make<TH1F> ("ra2b_met" ,"RA2b: (MET from PFmetHandle);MET [GeV];Events;", evt_met_bins, 0, evt_met_max);
	evtHist.mht = fs->make<TH1F> ("ra2b_mht" ,"RA2b: (MHT from PFmetHandle);MHT [GeV];Events;", evt_met_bins, 0, evt_met_max);
	evtHist.ht = fs->make<TH1F> ("ra2b_ht" ,"RA2b: HT from Jets ET>50 && |#Eta|<2.4;HT [GeV];Events;", evt_ht_bins, 0, evt_ht_max);
	evtHist.njet = fs->make<TH1F> ("ra2b_njet_et50eta24" ,"RA2b: Njets (Et>50 && |#Eta|<2.4;NJETS;Events;", 10, 0, 10);
	evtHist.metphi = fs->make<TH1F> ("metphi" ,"PAT-tuple quantities for debugging;met phi;Events;", 160, -8, 8);
	evtHist.meff = fs->make<TH1F> ("ra2b_meteff" ,"RA2b:;MEff;Events;", 50, 0, 5000);


	const double pt_bins = 200, pt_max = 2000;
	pf30_jet1Hist.pt = fs->make<TH1F> ("ra2b_pf30_jet1_pt" ,"RA2b: PF30-Jet1 pt;pt [GeV];Events;", pt_bins, 0, pt_max);
	pf30_jet1Hist.eta = fs->make<TH1F> ("ra2b_pf30_jet1_eta" ,"RA2b: PF30-Jet1 eta;eta;Events;", 100, -5, 5);
	pf30_jet1Hist.phi = fs->make<TH1F> ("ra2b_pf30_jet1_phi" ,"RA2b: PF30-Jet1 phi;phi;Events;", 160, -8, 8);
//	pf30_jet1Hist.delphi = fs->make<TH1F> ("pf30_jet1_delphi" ,"RA2b: PF30-Jet1: delphi;delphi;Events;", 160, -8, 8);

	pf30_jet2Hist.pt = fs->make<TH1F> ("ra2b_pf30_jet2_pt" ,"RA2b: PF30-Jet2 pt;pt [GeV];Events;", pt_bins, 0, pt_max);
	pf30_jet2Hist.eta = fs->make<TH1F> ("ra2b_pf30_jet2_eta" ,"RA2b: PF30-Jet2 eta;eta;Events;", 100, -5, 5);
	pf30_jet2Hist.phi = fs->make<TH1F> ("ra2b_pf30_jet2_phi" ,"RA2b: PF30-Jet2 phi;phi;Events;", 160, -8, 8);
//	pf30_jet2Hist.delphi = fs->make<TH1F> ("ra2b_pf30_jet2_delphi" ,"RA2b: PF30-Jet2: delphi;delphi;Events;", 160, -8, 8);

	pf30_jet3Hist.pt = fs->make<TH1F> ("ra2b_pf30_jet3_pt" ,"RA2b: PF30-Jet3 pt;pt [GeV];Events;", pt_bins, 0, pt_max);
	pf30_jet3Hist.eta = fs->make<TH1F> ("ra2b_pf30_jet3_eta" ,"RA2b: PF30-Jet3 eta;eta;Events;", 100, -5, 5);
	pf30_jet3Hist.phi = fs->make<TH1F> ("ra2b_pf30_jet3_phi" ,"RA2b: PF30-Jet3 phi;phi;Events;", 160, -8, 8);
//	pf30_jet3Hist.delphi = fs->make<TH1F> ("ra2b_pf30_jet3_delphi" ,"RA2b: PF30-Jet3: delphi;delphi;Events;", 160, -8, 8);


	//pf30_jet1Hist.delT = fs->make<TH1F> ("ra2b_pf30_jet1_delT" ,"RA2b: PF30-Jet1: #Delta T; #Delta T;Events;", 200, 0, 1000); 
	pf30_jet1Hist.delTDevidedByJetPt = fs->make<TH1F> ("ra2b_pf30_jet1_delTJetPt" ,"RA2b: PF30-Jet1; #Delta T/P_{T}^{Jet1}; Events;", 200, 0, 200); 
	//pf30_jet2Hist.delT = fs->make<TH1F> ("ra2b_pf30_jet2_delT" ,"RA2b: PF30-Jet2: #Delta T; #Delta T;Events;", 200, 0, 1000); 
	pf30_jet2Hist.delTDevidedByJetPt = fs->make<TH1F> ("ra2b_pf30_jet2_delTJetPt" ,"RA2b: PF30-Jet2; #Delta T/P_{T}^{Jet2}; Events;", 200, 0, 200); 
	//pf30_jet3Hist.delT = fs->make<TH1F> ("ra2b_pf30_jet3_delT" ,"RA2b: PF30-Jet3: #Delta T; #Delta T;Events;", 200, 0, 1000); 
	pf30_jet3Hist.delTDevidedByJetPt = fs->make<TH1F> ("ra2b_pf30_jet3_delTJetPt" ,"RA2b: PF30-Jet3; #Delta T/P_{T}^{Jet3}; Events;", 200, 0, 200); 

	pf30_jet1Hist.delTvsJetPt = fs->make<TProfile> ("ra2b_pf30_jet1_delTvsJetPt" ,"RA2b: PF30-Jet1: #Delta T vs Jet1 Pt; Jet1 Pt;#Delta T;", 50, 0, 2000); 
	pf30_jet1Hist.delTDevidedByJetPtvsJetPt = fs->make<TProfile> ("ra2b_pf30_jet1_delTDevidedByJetPtvsJetPt" ,"RA2b: PF30-Jet1; Jet1 Pt;#Delta T/P_{T}^{Jet1};", 50, 0, 2000); 
	pf30_jet2Hist.delTvsJetPt = fs->make<TProfile> ("ra2b_pf30_jet2_delTvsJetPt" ,"RA2b: PF30-Jet1: #Delta T vs Jet2 Pt; Jet2 Pt;#Delta T;", 50, 0, 2000); 
	pf30_jet2Hist.delTDevidedByJetPtvsJetPt = fs->make<TProfile> ("ra2b_pf30_jet2_delTDevidedByJetPtvsJetPt" ,"RA2b: PF30-Jet2; Jet2 Pt;#Delta T/P_{T}^{Jet2};", 50, 0, 2000); 
	pf30_jet3Hist.delTvsJetPt = fs->make<TProfile> ("ra2b_pf30_jet3_delTvsJetPt" ,"RA2b: PF30-Jet1: #Delta T vs Jet3 Pt; Jet3 Pt;#Delta T;", 50, 0, 2000); 
	pf30_jet3Hist.delTDevidedByJetPtvsJetPt = fs->make<TProfile> ("ra2b_pf30_jet3_delTDevidedByJetPtvsJetPt" ,"RA2b: PF30-Jet3; Jet3 Pt;#Delta T/P_{T}^{Jet3};", 50, 0, 2000); 

	hDelPhiMin[0] = fs->make<TH1F> ("ra2b_delPhiMin","RA2b: #Delta#Phi_{min} distribution", 40, 0, 4);
	hDelPhiMin[1] = fs->make<TH1F> ("ra2b_delPhiMin_MET_0to50","RA2b: #Delta#Phi_{min} distribution for 0<#slash{E}_{T}<50 GeV", 40, 0, 4);
	hDelPhiMin[2] = fs->make<TH1F> ("ra2b_delPhiMin_MET_50to100","RA2b: #Delta#Phi_{min} distribution for 50<#slash{E}_{T}<100 GeV", 40, 0, 4);
	hDelPhiMin[3] = fs->make<TH1F> ("ra2b_delPhiMin_MET_100to150","RA2b: #Delta#Phi_{min} distribution for 100<#slash{E}_{T}<150 GeV", 40, 0, 4);
	hDelPhiMin[4] = fs->make<TH1F> ("ra2b_delPhiMin_MET_150up","RA2b: #Delta#Phi_{min} distribution for #slash{E}_{T}>150 GeV", 40, 0, 4);
	hDelPhiMinNorm[0] = fs->make<TH1F> ("ra2b_delPhiMinNorm","RA2b: #Delta#Phi_{min}^{norm} distribution", 100, 0, 20);
	hDelPhiMinNorm[1] = fs->make<TH1F> ("ra2b_delPhiMinNorm_MET_0to50","RA2b: #Delta#Phi_{min}^{norm} distribution for 0<#slash{E}_{T}<50 GeV", 100, 0, 20);
	hDelPhiMinNorm[2] = fs->make<TH1F> ("ra2b_delPhiMinNorm_MET_50to100","RA2b: #Delta#Phi_{min}^{norm} distribution for 50<#slash{E}_{T}<100 GeV", 100, 0, 20);
	hDelPhiMinNorm[3] = fs->make<TH1F> ("ra2b_delPhiMinNorm_MET_100to150","RA2b: #Delta#Phi_{min}^{norm} distribution for 100<#slash{E}_{T}<150 GeV", 100, 0, 20);
	hDelPhiMinNorm[4] = fs->make<TH1F> ("ra2b_delPhiMinNorm_MET_150up","RA2b: #Delta#Phi_{min}^{norm} distribution for #slash{E}_{T}>150 GeV", 100, 0, 20);
	//hDelPhiMinNormVshDelPhiMin = fs->make<TH2F> ("delPhiMinNormVsdelPhiMin","delPhiMinNormalized VS delPhiMin", 500,0,5, 2000, 0, 20);
	hDelPhiMinVsMET = fs->make<TProfile> ("ra2b_delPhiMinVsMET","RA2b: #Delta#Phi_{min} vs #slash{E}_{T}", 30,0,600, 0, 4);
	hDelPhiMinNormVsMET = fs->make<TProfile> ("ra2b_delPhiMinNormVsMET","RA2b: #Delta#Phi_{min}^{norm} vs #slash{E}_{T}",30,0,600, 0, 20);

	//const float npassFailHistBins = 16;
	//const float passFailHistBins[] = {0,20,40,60,80,100,120,140,160,180,200,250,300,350,400,500,600};
	const float npassFailHistBins = 15;
	const float passFailHistBins[] = {0,25,50,75,100,125,150,175,200,225,250,275,300,350,400,500};

	hPass = fs->make<TH1F> ("ra2b_Pass","RA2b: PASS from #Delta#Phi_{min} cut", npassFailHistBins, passFailHistBins);
	hFail = fs->make<TH1F> ("ra2b_Fail","RA2b: FAIL from #Delta#Phi_{min} cut", npassFailHistBins, passFailHistBins);
	hPass->Sumw2();
	hFail->Sumw2();

	for (int i=0; i< 10; ++i)
	{
		std::stringstream name1, title1, name2, title2;
		name1 << "ra2b_Pass_Norm_" << i+1;
		title1 << "RA2b: PASS from #Delta#Phi_{min}^{norm} > "<< i+1 << ";#slash{E}_{T};Events;";
		name2 << "ra2b_Fail_Norm_" << i+1;
		title2 << "RA2b: FAIL from #Delta#Phi_{min}^{norm} < "<< i+1 << ";#slash{E}_{T};Events;";

		hPass_Norm[i] = fs->make<TH1F> (name1.str().c_str(), title1.str().c_str(), npassFailHistBins, passFailHistBins);
		hPass_Norm[i]->Sumw2();
		hFail_Norm[i] = fs->make<TH1F> (name2.str().c_str(), title2.str().c_str(), npassFailHistBins, passFailHistBins);
		hFail_Norm[i]->Sumw2();

		//RA2 plots
		std::stringstream name11, title11, name22, title22;
		name11 << "ra2_Pass_Norm_" << i+1;
		title11 << "RA2: PASS from #Delta#Phi_{min}^{norm} > "<< i+1 << ";#slash{H}_{T};Events;";
		name22 << "ra2_Fail_Norm_" << i+1;
		title22 << "RA2: FAIL from #Delta#Phi_{min}^{norm} < "<< i+1 << ";#slash{H}_{T};Events;";

		hPass_Norm_mht[i] = fs->make<TH1F> (name11.str().c_str(), title11.str().c_str(), npassFailHistBins, passFailHistBins);
		hPass_Norm_mht[i]->Sumw2();
		hFail_Norm_mht[i] = fs->make<TH1F> (name22.str().c_str(), title22.str().c_str(), npassFailHistBins, passFailHistBins);
		hFail_Norm_mht[i]->Sumw2();

	}


	hDelPhiMinNorm_mht[0] = fs->make<TH1F> ("ra2_delPhiMinNorm","RA2: #Delta#Phi_{min}^{norm} distribution", 100, 0, 20);
	hDelPhiMinNorm_mht[1] = fs->make<TH1F> ("ra2_delPhiMinNorm_MHT_0to50","RA2: #Delta#Phi_{min}^{norm} distribution for 0<#slash{H}_{T}<50 GeV", 100, 0, 20);
	hDelPhiMinNorm_mht[2] = fs->make<TH1F> ("ra2_delPhiMinNorm_MHT_50to100","RA2: #Delta#Phi_{min}^{norm} distribution for 50<#slash{H}_{T}<100 GeV", 100, 0, 20);
	hDelPhiMinNorm_mht[3] = fs->make<TH1F> ("ra2_delPhiMinNorm_MHT_100to200","RA2: #Delta#Phi_{min}^{norm} distribution for 100<#slash{H}_{T}<200 GeV", 100, 0, 20);
	hDelPhiMinNorm_mht[4] = fs->make<TH1F> ("ra2_delPhiMinNorm_MHT_200to300","RA2: #Delta#Phi_{min}^{norm} distribution for 200<#slash{H}_{T}<300 GeV", 100, 0, 20);
	hDelPhiMinNorm_mht[5] = fs->make<TH1F> ("ra2_delPhiMinNorm_MHT_300to400","RA2: #Delta#Phi_{min}^{norm} distribution for 300<#slash{H}_{T}<400 GeV", 100, 0, 20);
	hDelPhiMinNorm_mht[6] = fs->make<TH1F> ("ra2_delPhiMinNorm_MHT_400to500","RA2: #Delta#Phi_{min}^{norm} distribution for 400<#slash{H}_{T}<500 GeV", 100, 0, 20);
	hDelPhiMinNorm_mht[7] = fs->make<TH1F> ("ra2_delPhiMinNorm_MHT_500up","RA2: #Delta#Phi_{min}^{norm} distribution for #slash{H}_{T}>500 GeV", 100, 0, 20);

	//RA2 plots
	pf50eta25_jet1Hist.pt = fs->make<TH1F> ("pat_pf50eta25_jet1_pt" ,"PAT-tuple quantities for debugging: pf50eta25-Jet1 pt;pt [GeV];Events;", pt_bins, 0, pt_max);
	pf50eta25_jet1Hist.eta = fs->make<TH1F> ("pat_pf50eta25_jet1_eta" ,"PAT-tuple quantities for debugging: pf50eta25-Jet1 eta;eta;Events;", 100, -5, 5);
	pf50eta25_jet1Hist.phi = fs->make<TH1F> ("pat_pf50eta25_jet1_phi" ,"PAT-tuple quantities for debugging: pf50eta25-Jet1 phi;phi;Events;", 160, -8, 8);
	pf50eta25_jet1Hist.delphi = fs->make<TH1F> ("pat_pf50eta25_jet1_delphi" ,"PAT-tuple quantities for debugging: pf50eta25-Jet1: delphi;delphi;Events;", 160, -8, 8);

	pf50eta25_jet2Hist.pt = fs->make<TH1F> ("pat_pf50eta25_jet2_pt" ,"PAT-tuple quantities for debugging: pf50eta25-Jet2 pt;pt [GeV];Events;", pt_bins, 0, pt_max);
	pf50eta25_jet2Hist.eta = fs->make<TH1F> ("pat_pf50eta25_jet2_eta" ,"PAT-tuple quantities for debugging: pf50eta25-Jet2 eta;eta;Events;", 100, -5, 5);
	pf50eta25_jet2Hist.phi = fs->make<TH1F> ("pat_pf50eta25_jet2_phi" ,"PAT-tuple quantities for debugging: pf50eta25-Jet2 phi;phi;Events;", 160, -8, 8);
	pf50eta25_jet2Hist.delphi = fs->make<TH1F> ("pat_pf50eta25_jet2_delphi" ,"PAT-tuple quantities for debugging: pf50eta25-Jet2: delphi;delphi;Events;", 160, -8, 8);

	pf50eta25_jet3Hist.pt = fs->make<TH1F> ("pat_pf50eta25_jet3_pt" ,"PAT-tuple quantities for debugging: pf50eta25-Jet3 pt;pt [GeV];Events;", pt_bins, 0, pt_max);
	pf50eta25_jet3Hist.eta = fs->make<TH1F> ("pat_pf50eta25_jet3_eta" ,"PAT-tuple quantities for debugging: pf50eta25-Jet3 eta;eta;Events;", 100, -5, 5);
	pf50eta25_jet3Hist.phi = fs->make<TH1F> ("pat_pf50eta25_jet3_phi" ,"PAT-tuple quantities for debugging: pf50eta25-Jet3 phi;phi;Events;", 160, -8, 8);
	pf50eta25_jet3Hist.delphi = fs->make<TH1F> ("pat_pf50eta25_jet3_delphi" ,"PAT-tuple quantities for debugging: pf50eta25-Jet3: delphi;delphi;Events;", 160, -8, 8);


//	hist_delphiMin_jetmet[0] = fs->make<TH1F> ("delphiMin_jetmet" ,"delphi Min.;delphi min;Events;", 160, -8, 8);
//	hist_delphiMin_jetmetVsMHT[0] = fs->make<TH2F> ("delphiMin_jetmetVsMHT" ,"delphi min Vs MHT;MHT;delphi min;", 1000, 0, 1000,400,0,4);



}


RA2bQCDvetoAna::~RA2bQCDvetoAna()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
RA2bQCDvetoAna::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
	++uProcessed;
	bool bSaveEvent = false;
	kRunLumiEvent.run  = iEvent.id().run();
	kRunLumiEvent.lumi = iEvent.id().luminosityBlock(); 
	kRunLumiEvent.evt  = iEvent.id().event();

//    std::cout << __LINE__<< ":: Processing event: "
//		 	<<  kRunLumiEvent.run << ":" << kRunLumiEvent.lumi 
//			<< ":" << kRunLumiEvent.evt << std::endl;
//			return 0;

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   //Handle<ExampleData> pIn;
   //iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   //ESHandle<SetupData> pSetup;
   //iSetup.get<SetupRecord>().get(pSetup);
#endif
	
    if (iVerbose) std::cout << __LINE__<< ":: Processing event: "
		 	<<  kRunLumiEvent.run << ":" << kRunLumiEvent.lumi 
			<< ":" << kRunLumiEvent.evt << std::endl;

	/* PU S3 reweighting for QCD MC sample
	 *
	 */
	lumiWeight = LumiWeights_.weight( iEvent );
	//std::cout << "lum wgt = " << lumiWeight << std::endl;
	sumLumiWeights += lumiWeight;

	/* NO need. I am pplying STD RA2 Selection
	 *
	 */
	Handle<reco::VertexCollection> vertexHandle;
	iEvent.getByLabel("goodVerticesRA2", vertexHandle);
	double dNvtx = 0;
   if (vertexHandle.isValid())
	{
     reco::VertexCollection vertexCollection = *(vertexHandle.product());
     if (iVerbose) std::cout << __LINE__<< "::" << __FUNCTION__ << "::VtxSize = "<<  vertexCollection.size() << std::endl;
     reco::VertexCollection::const_iterator v = vertexCollection.begin();
	  if (v->isValid())
	  {
		  for (; v != vertexCollection.end(); ++v)
		  {
			  if (v->isFake()) continue;
			  double dPrimVtx_z = v->z();
			  if (fabs(dPrimVtx_z) > dMaxPrimVtxZ) continue;  // do not use the event. 
			  if (v->ndof() < inMinNdofVtx) continue; 
			  //if (fabs(v->rho()) < dMaxPrimVtxRho) continue; 
			  if (iVerbose) std::cout << "Valid Vertex Found " << std::endl;
			  ++dNvtx;
			  evtHist.vtxz->Fill(dPrimVtx_z);
		  }
	  }
	}
	

	evtHist.nvtx->Fill(dNvtx);
	
	/*  Get all the handles needed and check their
	 *  validity.
	 */
	iEvent.getByLabel(htInputTag_, htHandle);
	if (! htHandle.isValid()) 
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":htHandle handle not found!" << std::endl;
		assert(false);
	} else
	{
		if (iVerbose) std::cout << __FUNCTION__ << ":" << __LINE__ << ":" << htInputTag_<< " found!" << std::endl;  
	}

	iEvent.getByLabel("pfMet",pfMetHandle);
	if (! pfMetHandle.isValid())
	{
		std::cout << "Valid PFMET Collection Not Found!" << std::endl;
		assert(false);
	}

	iEvent.getByLabel(patJetsPFPt50Eta25InputTag_, pfpt50eta25JetHandle);
	if (! pfpt50eta25JetHandle.isValid()) 
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":pfpt50eta25JetHandle handle not found!" << std::endl;
		assert(false);
	} else
	{
		if (iVerbose) std::cout << __FUNCTION__ << ":" << __LINE__ << ":" << patJetsPFPt50Eta25InputTag_ << " found!" << std::endl;  
	}
	iEvent.getByLabel(patJetsPFPt30InputTag_, pfpt30JetHandle);
	if (! pfpt30JetHandle.isValid()) 
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":pfpt30JetHandle handle not found!" << std::endl;
		assert(false);
	}

	iEvent.getByLabel(mhtInputTag_, mhtHandle);
	if (! mhtHandle.isValid()) 
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":MHT handle not found!" << std::endl;
		assert(false);
	} else
	{
		if (iVerbose) std::cout << __FUNCTION__ << ":" << __LINE__ << ":" 
			<< mhtInputTag_ << " found!" << std::endl;  
	}

	/*replicate RA2b HT which uses only jet with et>50and eta<2.4
	 */
	//std::cout << "====================================" << std::endl;
	unsigned njet50eta24 = 0;
	TLorentzVector vSumPt30JetVec(0,0,0,0), vSumJetEt50Eta24(0,0,0,0);
	double dHt_et50eta24 = 0;
	std::vector<unsigned> vLead3JetIndices;
	TLorentzVector vJet1Vec(0,0,0,0), vJet2Vec(0,0,0,0), vJet3Vec(0,0,0,0);
	for (unsigned i = 0 ; i < pfpt30JetHandle->size() ; ++i)
	{
		const TLorentzVector iJetVec((*pfpt30JetHandle)[i].px(),
										(*pfpt30JetHandle)[i].py(),
										(*pfpt30JetHandle)[i].pz(),
										(*pfpt30JetHandle)[i].energy());
		//std::cout << "--- i, et =" << i << ", " << iJetVec.Pt() << std::endl;
		vSumPt30JetVec += iJetVec;
		if (iJetVec.Pt()<50. || fabs(iJetVec.Eta())>2.4) continue;
		//loose jet id stuff
		if ((*pfpt30JetHandle)[i].neutralHadronEnergyFraction() >=0.99) continue;
		if ((*pfpt30JetHandle)[i].neutralEmEnergyFraction() >=0.99) continue;
		if (((*pfpt30JetHandle)[i].getPFConstituents()).size() <=1) continue;
		if ((*pfpt30JetHandle)[i].chargedHadronEnergyFraction() <=0) continue;
		if ((*pfpt30JetHandle)[i].chargedMultiplicity() <=0) continue;
		if ((*pfpt30JetHandle)[i].chargedEmEnergyFraction() >=0.99) continue;
		
		if (vLead3JetIndices.size()<3) vLead3JetIndices.push_back(i);
		vSumJetEt50Eta24 += iJetVec;
		dHt_et50eta24 += iJetVec.Pt();
		++njet50eta24;
	}


	//APPLY RA2b cuts
	if (dNvtx<1) { ++uFailNvtxCut; return 0;}
	if (njet50eta24<3) { ++uFailNjet50Eta24Cut; return 0;}
	if (dHt_et50eta24<dMinHT) { ++uFailMinHTCut; return 0; }
	if ((*pfMetHandle)[0].pt() < dMinMet) { ++uFailMinPFMetCut; return 0;}

	//NOTE: htHandle gives the same value as the sumPt of all jets with Et>30 in pfpt30JetHandle
	//std::cout << "jet ht/ht handle = " << vSumJetEt50Eta24.Pt() << "/ " << (*htHandle)  
	//				<< " ( " << dHt_et50eta24 << std::endl;


	/*
	 * Fill hists after RA2b base selection
	 */
	if ( (*mhtHandle)[0].pt() < dMinMHT ) {++uFailMinPFMHTCut; return 0; }


	evtHist.nvtx->Fill(dNvtx);
	evtHist.njet->Fill(njet50eta24);
	evtHist.metphi->Fill((*pfMetHandle)[0].phi());
	evtHist.met->Fill( (*pfMetHandle)[0].pt());
	evtHist.mht->Fill( (*mhtHandle)[0].pt());
	evtHist.ht->Fill(dHt_et50eta24);
	evtHist.meff->Fill( (*mhtHandle)[0].pt() + (*htHandle) );

	pf30_jet1Hist.pt->Fill((*pfpt30JetHandle)[vLead3JetIndices.at(0)].pt());
	pf30_jet2Hist.pt->Fill((*pfpt30JetHandle)[vLead3JetIndices.at(1)].pt());
	pf30_jet3Hist.pt->Fill((*pfpt30JetHandle)[vLead3JetIndices.at(2)].pt());
	pf30_jet1Hist.eta->Fill((*pfpt30JetHandle)[vLead3JetIndices.at(0)].eta());
	pf30_jet2Hist.eta->Fill((*pfpt30JetHandle)[vLead3JetIndices.at(1)].eta());
	pf30_jet3Hist.eta->Fill((*pfpt30JetHandle)[vLead3JetIndices.at(2)].eta());

	DoDelMinStudy(pfMetHandle, pfpt30JetHandle, vLead3JetIndices, 
					mhtHandle, iEvent);


	//	RA2 Stuff from here

	//met calcualted from all jets with et>30

//	const float mht = (*mhtHandle)[0].pt();


/*
	if (mht>60 && mht<=80) FillHistograms(Hist[0], mhtHandle, jetHandle); 
	else if (mht>80 && mht<=100) FillHistograms(Hist[1], mhtHandle, jetHandle); 
	else if (mht>100&& mht<=120) FillHistograms(Hist[2], mhtHandle, jetHandle); 
	else if (mht>120&& mht<=140) FillHistograms(Hist[3], mhtHandle, jetHandle); 
	else if (mht>140&& mht<=170) FillHistograms(Hist[4], mhtHandle, jetHandle); 
	else if (mht>170&& mht<=200) FillHistograms(Hist[5], mhtHandle, jetHandle); 
	else if (mht>200) FillHistograms(Hist[6], mhtHandle, jetHandle); 
*/	


/*	unsigned numpt30Jets = pfpt30JetHandle->size();
	//std::cout << "pfpt30JetHandle numjets = " << numpt30Jets << std::endl;
	TLorentzVector oJetsVec(0,0,0,0);
	for (unsigned i = 0 ; i < numpt30Jets ; ++i)
	{
		//if ((*numpt50eta25Jets)[i].pt() > minRa2JetPt && (*numpt50eta25Jets)[i].eta()<2.5) ++Njet50;  
		//if ((*numpt50eta25Jets)[i].pt() > minJetPt) ++Njet30;  
			const TLorentzVector iJetVec((*pfpt30JetHandle)[i].px(),(*pfpt30JetHandle)[i].py(),(*pfpt30JetHandle)[i].pz(),(*pfpt30JetHandle)[i].energy());
			oJetsVec += iJetVec;
	}


	if (iVerbose)
	{
		std::cout << "mht/my mht = " << (*mhtHandle)[0].pt() << "/" << oJetsVec.Pt() << std::endl;
		std::cout << "ht/my ht = " << (*htHandle) << "/" << myHt << std::endl;
	}

*/

	++uPassed;
   return true;
}

// ------------ method called once each job just before starting event loop  ------------
void 
RA2bQCDvetoAna::beginJob()
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
RA2bQCDvetoAna::endJob() {

	std::cout << __FILE__ << "::" << __FUNCTION__  << "\n" << std::endl;
	std::cout << "[ATM:00] Job settings " << std::endl;
	std::cout << "[ATM:01] Events Processed --- = " << uProcessed << std::endl;
	std::cout << "[ATM:02] Events Passed ------ = " << uPassed << std::endl;
	std::cout << "[ATM:10] Primary Vtx Min ---- = " << inMinVtx << std::endl; 
	std::cout << "[ATM:11] Primary Vtx Min ndof = " << inMinNdofVtx << std::endl; 
/*	std::cout << "[ATM:20] Triggers Used ------ = ";
	if (req_trigger)
	{
		for (std::vector<std::string>::const_iterator trigNameTempl = triggerPathsToStore_.begin();
				trigNameTempl != triggerPathsToStore_.end(); ++trigNameTempl) 
		{
			std::cout << "[" << *trigNameTempl << "]";
		}
		std::cout << std::endl;
	} else {
		std::cout << "NONE" << std::endl; 
	}
	*/
	std::cout << "[ATM:30] Minimum MET -------- = " << dMinMet << std::endl;
	std::cout << "[ATM:31] Minimum HT --------- = " << dMinHT << std::endl;
	std::cout << "[ATM:32] Minimum MHT -------- = " << dMinMHT << std::endl;
	std::cout << "[ATM:40] PASS Summary --------- " << std::endl;
	std::cout << "[ATM:41] Pass nvtx ---------- = " << (uProcessed - uFailNvtxCut) << " (" << uFailNvtxCut << ")" << std::endl;
	std::cout << "[ATM:42] Pass njet ---------- = " << (uProcessed - uFailNvtxCut - uFailNjet50Eta24Cut) 
																	<< " (" << uFailNjet50Eta24Cut << ")" << std::endl;
	std::cout << "[ATM:43] Pass ht ------------ = " << (uProcessed - uFailNvtxCut - uFailNjet50Eta24Cut - uFailMinHTCut) 
																	<< " (" << uFailMinHTCut << ")" << std::endl;
	std::cout << "[ATM:44] Pass met ----------- = " << (uProcessed - uFailNvtxCut - uFailNjet50Eta24Cut - uFailMinHTCut - uFailMinPFMetCut) 
																	<< " (" << uFailMinPFMetCut << ")" << std::endl;
	std::cout << "[ATM:45] Pass mht ----------- = " << (uProcessed - uFailNvtxCut - uFailNjet50Eta24Cut - uFailMinHTCut - uFailMinPFMetCut - uFailMinPFMHTCut) 
																	<< " (" << uFailMinPFMHTCut << ")" << std::endl;
	std::cout << "[ATM:50] LumiWeights Avg ---- = " << sumLumiWeights/(double)uPassed << std::endl;

}

// ------------ method called when starting to processes a run  ------------
bool 
RA2bQCDvetoAna::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
RA2bQCDvetoAna::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
RA2bQCDvetoAna::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
RA2bQCDvetoAna::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
RA2bQCDvetoAna::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
void
RA2bQCDvetoAna::BookHistograms(edm::Service<TFileService>& fs, Hist_t& hist , const float mHt_min, const float mHt_max)
{
	std::stringstream mht_range;
	mht_range<< mHt_min << "_MET_" << mHt_max; 

	std::stringstream title_nvtx, title_vtxz, title_mht, title_ht, title_njet, title_metphi;
	std::stringstream name_nvtx, name_vtxz, name_mht, name_ht, name_njet, name_metphi;
	title_nvtx << mht_range.str() << ";Nvtx;Events;";
	title_vtxz << mht_range.str() << ";Primary Vertex z [cm];Arbitrary;";
	title_mht << mht_range.str() << ";MHT [GeV];Events;";
	title_ht << mht_range.str() << ";H_{T} [GeV];Events;";
	title_njet << mht_range.str() << ";N-JET(E_{T}>15) [GeV];Events;";
	title_metphi << mht_range.str() << ";#phi(#slash{E}_{T});Events;";

	name_nvtx << mht_range.str() << "_nvtx";
	name_vtxz << mht_range.str() << "_vtxz";
	name_mht << mht_range.str() << "_mht";
	name_ht << mht_range.str() << "_ht";
	name_njet << mht_range.str() << "_njet";
	name_metphi << mht_range.str() << "_metphi";

	const double met_max = 4000, met_bins = 400;
	hist.evt.nvtx = fs->make<TH1F> (name_nvtx.str().c_str() ,title_nvtx.str().c_str(), 15, 0, 15);
	hist.evt.vtxz = fs->make<TH1F> (name_vtxz.str().c_str(),title_vtxz.str().c_str(),30,0,30);
	hist.evt.mht = fs->make<TH1F> (name_mht.str().c_str() ,title_mht.str().c_str(), met_bins, 0, met_max);
	hist.evt.ht = fs->make<TH1F> (name_ht.str().c_str() ,title_ht.str().c_str(), met_bins, 0, met_max);
	hist.evt.njet = fs->make<TH1F> (name_njet.str().c_str() ,title_njet.str().c_str(), 10, 0, 10);
	hist.evt.metphi = fs->make<TH1F> (name_metphi.str().c_str() ,title_metphi.str().c_str(), 160, -8, 8);

	const float pt_bins = 200, pt_max = 2000;
	const float eta_bins = 100, eta_min = -5, eta_max = 5;
	const float phi_bins = 160, phi_min = -8, phi_max = 8;
	
	for (int i=0; i <5; ++i)
	{
		std::stringstream name_pt, name_eta, name_phi, name_delphi, title_pt, title_eta, title_phi, title_delphi;
		std::stringstream jetname;
		jetname << "jet" << (++i) << "_";
		name_pt << mht_range.str() << "_"<< jetname.str() << "pt";
		name_eta << mht_range.str() << "_" << jetname.str() << "eta";
		name_phi << mht_range.str() << "_" << jetname.str() << "phi";
		name_delphi << mht_range.str() << "_" << jetname.str() << "delphi";
		title_pt << mht_range.str() << ": " << jetname.str() << ": " << ";pt [GeV];Events;";
		title_eta << mht_range.str() << ";#eta;Events;";
		title_phi << mht_range.str() << ";#Phi;Events;";
		title_delphi << mht_range.str() << ";#Delta#Phi(jet,#slash{E}_{T});Events;";

		hist.jet[i].pt = fs->make<TH1F> (name_pt.str().c_str(),title_pt.str().c_str(), pt_bins, 0, pt_max);
		hist.jet[i].eta = fs->make<TH1F> (name_eta.str().c_str(),title_eta.str().c_str(), eta_bins, eta_min, eta_max);
		hist.jet[i].phi = fs->make<TH1F> (name_phi.str().c_str(),title_phi.str().c_str(), phi_bins, phi_min, phi_max);
		hist.jet[i].delphi = fs->make<TH1F> (name_delphi.str().c_str(),title_delphi.str().c_str(), phi_bins, phi_min, phi_max);
	}

	std::stringstream dphimin_jetmet_j123_name, dphimin_jetmet_j123prof_name;
	std::stringstream dphimin_jetmet_j123_title, dphimin_jetmet_j123prof_title;
	dphimin_jetmet_j123_name << mht_range.str() << "_delphiMin_jetmet_j123";
	dphimin_jetmet_j123prof_name << mht_range.str() << "_delphiMin_jetmetVsMHT_j123";
	dphimin_jetmet_j123_title << mht_range.str() << ";#Delta #phi_{min}(j_{1,2,3},#slash{E}_{T});Events;";
	dphimin_jetmet_j123prof_title << mht_range.str() << ";MHT;#Delta #phi_{min}(j_{1,2,3},#slash{E}_{T});";
	std::stringstream dphimin_jetmet_j1234_name, dphimin_jetmet_j1234prof_name;
	std::stringstream dphimin_jetmet_j1234_title, dphimin_jetmet_j1234prof_title;
	dphimin_jetmet_j1234_name << mht_range.str() << "_delphiMin_jetmet_j1234";
	dphimin_jetmet_j1234prof_name << mht_range.str() << "_delphiMin_jetmetVsMHT_j1234";
	dphimin_jetmet_j1234_title << mht_range.str() << ";#Delta #phi_{min}(j_{1,2,3,4},#slash{E}_{T});Events;";
	dphimin_jetmet_j1234prof_title << mht_range.str() << ";MHT;#Delta #phi_{min}(j_{1,2,3,4},#slash{E}_{T});";
	std::stringstream dphimin_jetmet_j12345_name, dphimin_jetmet_j12345prof_name;
	std::stringstream dphimin_jetmet_j12345_title, dphimin_jetmet_j12345prof_title;
	dphimin_jetmet_j12345_name << mht_range.str() << "_delphiMin_jetmet_j12345";
	dphimin_jetmet_j12345prof_name << mht_range.str() << "_delphiMin_jetmetVsMHT_j12345";
	dphimin_jetmet_j12345_title << mht_range.str() << ";#Delta #phi_{min}(j_{1,2,3,4,5},#slash{E}_{T});Events;";
	dphimin_jetmet_j12345prof_title << mht_range.str() << ";MHT;#Delta #phi_{min}(j_{1,2,3,4,5},#slash{E}_{T});";

	hist.delphiMin_jetmet[0] = fs->make<TH1F> (dphimin_jetmet_j123_name.str().c_str(), dphimin_jetmet_j123_title.str().c_str(), 400, 0, 4);
	hist.delphiMin_jetmetVsMHT[0] = fs->make<TProfile> (dphimin_jetmet_j123prof_name.str().c_str(), dphimin_jetmet_j123prof_title.str().c_str(), 160, 0, 800,0,4);
	hist.delphiMin_jetmet[1] = fs->make<TH1F> (dphimin_jetmet_j1234_name.str().c_str(), dphimin_jetmet_j1234_title.str().c_str(), 400, 0, 4);
	hist.delphiMin_jetmetVsMHT[1] = fs->make<TProfile> (dphimin_jetmet_j1234prof_name.str().c_str(), dphimin_jetmet_j1234prof_title.str().c_str(), 160, 0, 800,0,4);
	hist.delphiMin_jetmet[2] = fs->make<TH1F> (dphimin_jetmet_j12345_name.str().c_str(), dphimin_jetmet_j12345_title.str().c_str(), 400, 0, 4);
	hist.delphiMin_jetmetVsMHT[2] = fs->make<TProfile> (dphimin_jetmet_j12345prof_name.str().c_str(), dphimin_jetmet_j12345prof_title.str().c_str(), 160, 0, 800,0,4);

	std::stringstream dphimin_jetjet_j123_name, dphimin_jetjet_j123prof_name;
	std::stringstream dphimin_jetjet_j123_title, dphimin_jetjet_j123prof_title;
	dphimin_jetjet_j123_name << mht_range.str() << "_delphiMin_jetjet_j123";
	dphimin_jetjet_j123prof_name << mht_range.str() << "_delphiMin_jetjetVsMHT_j123";
	dphimin_jetjet_j123_title << mht_range.str() << ";#Delta #phi_{min}(j_{1,2,3},#slash{E}_{T});Events;";
	dphimin_jetjet_j123prof_title << mht_range.str() << ";MHT;#Delta #phi_{min}(j_{1,2,3},#slash{E}_{T});";
	std::stringstream dphimin_jetjet_j1234_name, dphimin_jetjet_j1234prof_name;
	std::stringstream dphimin_jetjet_j1234_title, dphimin_jetjet_j1234prof_title;
	dphimin_jetjet_j1234_name << mht_range.str() << "_delphiMin_jetjet_j1234";
	dphimin_jetjet_j1234prof_name << mht_range.str() << "_delphiMin_jetjetVsMHT_j1234";
	dphimin_jetjet_j1234_title << mht_range.str() << ";#Delta #phi_{min}(j_{1,2,3,4},#slash{E}_{T});Events;";
	dphimin_jetjet_j1234prof_title << mht_range.str() << ";MHT;#Delta #phi_{min}(j_{1,2,3,4},#slash{E}_{T});";
	std::stringstream dphimin_jetjet_j12345_name, dphimin_jetjet_j12345prof_name;
	std::stringstream dphimin_jetjet_j12345_title, dphimin_jetjet_j12345prof_title;
	dphimin_jetjet_j12345_name << mht_range.str() << "_delphiMin_jetjet_j12345";
	dphimin_jetjet_j12345prof_name << mht_range.str() << "_delphiMin_jetjetVsMHT_j12345";
	dphimin_jetjet_j12345_title << mht_range.str() << ";#Delta #phi_{min}(j_{1,2,3,4,5},#slash{E}_{T});Events;";
	dphimin_jetjet_j12345prof_title << mht_range.str() << ";MHT;#Delta #phi_{min}(j_{1,2,3,4,5},#slash{E}_{T});";

	hist.delphiMin_jetjet[0] = fs->make<TH1F> (dphimin_jetjet_j123_name.str().c_str(), dphimin_jetjet_j123_title.str().c_str(), 400, 0, 4);
	hist.delphiMin_jetjetVsMHT[0] = fs->make<TProfile> (dphimin_jetjet_j123prof_name.str().c_str(), dphimin_jetjet_j123prof_title.str().c_str(), 160, 0, 800,0,4);
	hist.delphiMin_jetjet[1] = fs->make<TH1F> (dphimin_jetjet_j1234_name.str().c_str(), dphimin_jetjet_j1234_title.str().c_str(), 400, 0, 4);
	hist.delphiMin_jetjetVsMHT[1] = fs->make<TProfile> (dphimin_jetjet_j1234prof_name.str().c_str(), dphimin_jetjet_j1234prof_title.str().c_str(), 160, 0, 800,0,4);
	hist.delphiMin_jetjet[2] = fs->make<TH1F> (dphimin_jetjet_j12345_name.str().c_str(), dphimin_jetjet_j12345_title.str().c_str(), 400, 0, 4);
	hist.delphiMin_jetjetVsMHT[2] = fs->make<TProfile> (dphimin_jetjet_j12345prof_name.str().c_str(), dphimin_jetjet_j12345prof_title.str().c_str(), 160, 0, 800,0,4);

	const float phislices[]={0.1,0.3,0.5,0.7,0.8};
	for (int i=0; i <5; ++i)
	{
		std::stringstream mht_by_phislice_name, mht_by_phislice_title;
		mht_by_phislice_name <<mht_range.str() << "_MHT_phiSlice";
		mht_by_phislice_title << mht_range.str() << ": ";
		if (i==0)
		{ 
			mht_by_phislice_name << "_lt" << phislices[i];
			mht_by_phislice_title << "|#Delta#Phi|<" << phislices[i];
		} else if (i>0 && i<4) 
		{
			mht_by_phislice_name << "_" << phislices[i] << "to"<<phislices[i+1];
			mht_by_phislice_title << "_" << phislices[i] << "<|#Delta#Phi|<" << phislices[i+1];
		} else if (i==4) 
		{
			mht_by_phislice_name << "_gt" << phislices[i];
			mht_by_phislice_title << "|#Delta#Phi|>" << phislices[i];
		}
		//hist.MHT_by_phislice[i] = fs->make<TH1F> (mht_by_phislice_name.str().c_str() ,mht_by_phislice_title.str().c_str(), met_bins, 0, met_max);
	}

}

void
RA2bQCDvetoAna::FillHistograms(Hist_t& hist
				, edm::Handle<edm::View<reco::MET> > mhtHandle 
				, edm::Handle<std::vector<pat::Jet> > jetHandle)
{
	const float minRa2JetPt= 50.;
	const float minJetPt= 30.;
	int Njet50 = 0;
	int Njet30 = 0;
	unsigned numJets = jetHandle->size();
	std::cout << "numjets = " << numJets << std::endl;
	for (unsigned i = 0 ; i < numJets ; ++i)
	{
		if ((*jetHandle)[i].pt() > minRa2JetPt && (*jetHandle)[i].eta()<2.5) ++Njet50;  
		if ((*jetHandle)[i].pt() > minJetPt) ++Njet30;  
	}
	if (Njet50 < 3) return;

	hist.evt.mht->Fill((*mhtHandle)[0].pt());

	evtHist.njet->Fill(Njet30);
	std::vector<float> vDelPhi_jetjet_3jet, vDelPhi_jetjet_4jet, vDelPhi_jetjet_5jet;
	std::vector<float> vDelPhi_jetmet_3jet, vDelPhi_jetmet_4jet, vDelPhi_jetmet_5jet;
	std::cout << __LINE__ << "3 jet min delphi" << std::endl;
	minDelPhi(vDelPhi_jetjet_3jet, vDelPhi_jetmet_3jet, 3, jetHandle, mhtHandle);

	hist.delphiMin_jetmet[0]->Fill(vDelPhi_jetmet_3jet.at(0));
	hist.delphiMin_jetjet[0]->Fill(vDelPhi_jetjet_3jet.at(0));
	hist.delphiMin_jetmetVsMHT[0]->Fill((*mhtHandle)[0].pt(), vDelPhi_jetmet_3jet.at(0));
	hist.delphiMin_jetjetVsMHT[0]->Fill((*mhtHandle)[0].pt(), vDelPhi_jetjet_3jet.at(0));
	if (Njet30>=4)
	{
		std::cout << __LINE__ << "4 jet min delphi" << std::endl;
		minDelPhi(vDelPhi_jetjet_4jet, vDelPhi_jetmet_4jet, 4, jetHandle, mhtHandle);
		hist.delphiMin_jetmet[1]->Fill(vDelPhi_jetmet_4jet.at(0));
		hist.delphiMin_jetjet[1]->Fill(vDelPhi_jetjet_4jet.at(0));
		hist.delphiMin_jetmetVsMHT[1]->Fill((*mhtHandle)[0].pt(), vDelPhi_jetmet_4jet.at(0));
		hist.delphiMin_jetjetVsMHT[1]->Fill((*mhtHandle)[0].pt(), vDelPhi_jetjet_4jet.at(0));
	}
	if (Njet30>=5) 
	{
		std::cout << __LINE__ << "5 jet min delphi" << std::endl;
		minDelPhi(vDelPhi_jetjet_5jet, vDelPhi_jetmet_5jet, 5, jetHandle, mhtHandle);
		hist.delphiMin_jetmet[2]->Fill(vDelPhi_jetmet_5jet.at(0));
		hist.delphiMin_jetjet[2]->Fill(vDelPhi_jetjet_5jet.at(0));
		hist.delphiMin_jetmetVsMHT[2]->Fill((*mhtHandle)[0].pt(), vDelPhi_jetmet_5jet.at(0));
		hist.delphiMin_jetjetVsMHT[2]->Fill((*mhtHandle)[0].pt(), vDelPhi_jetjet_5jet.at(0));
	}

	if (vDelPhi_jetmet_3jet.at(0)<0.1) MHT_by_phislice[0]->Fill((*mhtHandle)[0].pt());
	else if (vDelPhi_jetmet_3jet.at(0)>0.1 && vDelPhi_jetmet_3jet.at(0)<0.2) MHT_by_phislice[0]->Fill((*mhtHandle)[0].pt());
	else if (vDelPhi_jetmet_3jet.at(0)>0.2 && vDelPhi_jetmet_3jet.at(0)<0.3) MHT_by_phislice[0]->Fill((*mhtHandle)[0].pt());
	else if (vDelPhi_jetmet_3jet.at(0)>0.3 && vDelPhi_jetmet_3jet.at(0)<0.4) MHT_by_phislice[0]->Fill((*mhtHandle)[0].pt());
	else if (vDelPhi_jetmet_3jet.at(0)>0.5 && vDelPhi_jetmet_3jet.at(0)<0.8) MHT_by_phislice[0]->Fill((*mhtHandle)[0].pt());
	else if (vDelPhi_jetmet_3jet.at(0)>0.8) MHT_by_phislice[0]->Fill((*mhtHandle)[0].pt());

}

unsigned RA2bQCDvetoAna::minDelPhi(std::vector<float>& vDelPhi_jetjet
				, std::vector<float>& vDelPhi_jetmet
				, const unsigned njets
				, edm::Handle<std::vector<pat::Jet> > jetHandle
				, edm::Handle<edm::View<reco::MET> > mhtHandle 
				)
{
/********
 * finds the jet closes to MET using speficied number of jets.
 * for e.g. if njets ==3, this will find the jet closest to MET
 * using only the 3leading jets
 * ******/

	unsigned numJets = jetHandle->size();
	if (njets>0 && njets<numJets) numJets = njets;
	unsigned jetClose2Met = -1;
	unsigned minDelPhi = 99.99;
	for (unsigned i = 0 ; i < numJets ; ++i)
	{
		const float delphi_jetmet = fabs(TVector2::Phi_mpi_pi((*jetHandle)[i].phi() - (*mhtHandle)[0].phi()));
		vDelPhi_jetmet.push_back(delphi_jetmet);
		if (delphi_jetmet<minDelPhi)
		{
			minDelPhi = delphi_jetmet;
			jetClose2Met = i;
		}

		TLorentzVector oJetsVec(0,0,0,0);
		for (unsigned j = 0 ; j < numJets ; ++j)
		{
			if ( i == j ) continue;
			//std::cout << " other jets [" << j << "][" << (*jetHandle)[j].pt() << "][" << (*jetHandle)[j].phi() << std::endl;  
			const TLorentzVector jJetVec((*jetHandle)[j].px(),(*jetHandle)[j].py(),(*jetHandle)[j].pz(),(*jetHandle)[j].energy());
			oJetsVec += jJetVec;
//			if (j==njets) break; //uses only njets for the caculation
		}
		const float delphi = fabs(TVector2::Phi_mpi_pi((*jetHandle)[i].phi() - oJetsVec.Phi()));
		//std::cout << "delphi min for ijet = " << i << "[" << delphi << "]" << std::endl;
		vDelPhi_jetjet.push_back(delphi);
//		if (i==njets) break;  //uses only njets for the caculation
	}

	std::sort(vDelPhi_jetjet.begin(), vDelPhi_jetjet.end(), sort_using_less_than);	
	std::sort(vDelPhi_jetmet.begin(), vDelPhi_jetmet.end(), sort_using_less_than);	
	
	if (iVerbose)
	{
		std::cout << "delphi ordered (jetjet)= ";
		for (std::vector<float>::const_iterator it = vDelPhi_jetjet.begin(); it != vDelPhi_jetjet.end(); ++it)
		{
			std::cout << "  " << (*it);
		}
		std::cout << std::endl;
		std::cout << "delphi ordered (jetmet)= ";
		for (std::vector<float>::const_iterator it = vDelPhi_jetmet.begin(); it != vDelPhi_jetmet.end(); ++it)
		{
			std::cout << "  " << (*it);
		}
		std::cout << std::endl;
	}
	std::cout << __FUNCTION__ << ": jet close to met = " << jetClose2Met << std::endl;
	return jetClose2Met;

}


void RA2bQCDvetoAna::DoDelMinStudy(edm::Handle<std::vector<reco::PFMET> >pfMetHandle
				, edm::Handle<std::vector<pat::Jet> > pt30jetHandle
				, const std::vector<unsigned> vLead3JetIndices
				, edm::Handle<edm::View<reco::MET> > mhtHandle
				, edm::Event& iEvent
				)
{

	//PrintHeader();
	std::vector<float> vDelPhi_jetmet, vDelPhiNorm_jetmet;
	std::vector<float> vDelPhi_jetmht, vDelPhiNorm_jetmht;

	const float met = (*pfMetHandle)[0].pt();

	//const float mht = (*mhtHandle)[0].pt();
	iEvent.getByLabel("patJetsAK5PF", alljetHandle);
	TLorentzVector vSumAllJets(0,0,0,0);
	for (unsigned i = 0 ; i < alljetHandle->size() ; ++i)
	{
		const TLorentzVector iJetVec((*alljetHandle)[i].px(),
										(*alljetHandle)[i].py(),
										(*alljetHandle)[i].pz(),
										(*alljetHandle)[i].energy());
		//std::cout << "i, pt , eta = " << i << " / " <<  iJetVec.Pt() << " / " << iJetVec.Eta() << std::endl;
		if (iJetVec.Pt()<10) continue;
		vSumAllJets += iJetVec;
	}

	const float mht = vSumAllJets.Pt();

	unsigned njet = 0; //consider cross product only from lead 3 jets
	TLorentzVector vSumJetPt30Vec(0,0,0,0);

	for (unsigned i = 0 ; i < vLead3JetIndices.size() ; ++i)
	{
		const TLorentzVector iJetVec((*pt30jetHandle)[vLead3JetIndices.at(i)].px(),
												(*pt30jetHandle)[vLead3JetIndices.at(i)].py(),
												(*pt30jetHandle)[vLead3JetIndices.at(i)].pz(),
												(*pt30jetHandle)[vLead3JetIndices.at(i)].energy());
		vSumJetPt30Vec += iJetVec;
		const float delphi_jetmet = fabs(TVector2::Phi_mpi_pi((*pt30jetHandle)[i].phi() - (*pfMetHandle)[0].phi()));
		const float delphi_jetmht = fabs(TVector2::Phi_mpi_pi((*pt30jetHandle)[i].phi() - (*mhtHandle)[0].phi()));
		//const float delphi_jetmht_opt = std::abs(reco::deltaPhi( (*pt30jetHandle)[i].phi() , (*mhtHandle)[0].phi() ));
		//std::cout << "my jet-mht dphi" << delphi_jetmht <<  "(" << delphi_jetmht_opt << std::endl;

		vDelPhi_jetmet.push_back(delphi_jetmet);
		vDelPhi_jetmht.push_back(delphi_jetmht);

		/*
		 * loop over rest of the jets to calcualte deltaT for i th jet
		 */
		float sumCP2 = 0;  //sqaure of sum of cross products
		for (unsigned j = 0 ; j < pt30jetHandle->size() ; ++j)
		{
			//Skip cross product with itself
			if ( vLead3JetIndices.at(i) == j ) continue;
			//loose jet id stuff
			if ((*pt30jetHandle)[j].neutralHadronEnergyFraction() >=0.99) continue;
			if ((*pt30jetHandle)[j].neutralEmEnergyFraction() >=0.99) continue;
			if (((*pt30jetHandle)[j].getPFConstituents()).size() <=1) continue;
			if ((*pt30jetHandle)[j].chargedHadronEnergyFraction() <=0) continue;
			if ((*pt30jetHandle)[j].chargedMultiplicity() <=0) continue;
			if ((*pt30jetHandle)[j].chargedEmEnergyFraction() >=0.99) continue;

			const TLorentzVector jJetVec((*pt30jetHandle)[j].px(),(*pt30jetHandle)[j].py(),
										(*pt30jetHandle)[j].pz(),(*pt30jetHandle)[j].energy());
			sumCP2 += pow(iJetVec.Px() * jJetVec.Py() - iJetVec.Py() * jJetVec.Px() ,2);
		} //END of pt30jet loop

		const float delT_i = (0.1*sqrt(sumCP2))/iJetVec.Pt();
		const float delPhi_i = delphi_jetmet/atan2(delT_i,met);
		const float delPhi_i_mht = delphi_jetmht/atan2(delT_i,mht);
		vDelPhiNorm_jetmet.push_back(delPhi_i);
		vDelPhiNorm_jetmht.push_back(delPhi_i_mht);

		if (i==0) //first jet
		{
			//pf30_jet1Hist.delT->Fill(sqrt(sumCP2));
			pf30_jet1Hist.delTDevidedByJetPt->Fill(delT_i);
			pf30_jet1Hist.delTvsJetPt->Fill(iJetVec.Pt(),sqrt(sumCP2));
			pf30_jet1Hist.delTDevidedByJetPtvsJetPt->Fill(iJetVec.Pt(), delT_i);
		} else if (i==1) //2nd jet
		{
			//pf30_jet2Hist.delT->Fill(sqrt(sumCP2));
			pf30_jet2Hist.delTDevidedByJetPt->Fill(delT_i);
			pf30_jet2Hist.delTvsJetPt->Fill(iJetVec.Pt(),sqrt(sumCP2));
			pf30_jet2Hist.delTDevidedByJetPtvsJetPt->Fill(iJetVec.Pt(), delT_i);
		} else if (i==2) //3rd jet
		{
			//pf30_jet3Hist.delT->Fill(sqrt(sumCP2));
			pf30_jet3Hist.delTDevidedByJetPt->Fill(delT_i);
			pf30_jet3Hist.delTvsJetPt->Fill(iJetVec.Pt(),sqrt(sumCP2));
			pf30_jet3Hist.delTDevidedByJetPtvsJetPt->Fill(iJetVec.Pt(), delT_i);
		}
		
		//std::cout << "jet/dphi/dphi_norm = " << i << " | " << delphi_jetmet
		//			<< " | " << delPhi_i << std::endl;

	} //END of 3 lead jet loop

	std::sort(vDelPhi_jetmet.begin(), vDelPhi_jetmet.end(), sort_using_less_than);	
	std::sort(vDelPhiNorm_jetmet.begin(), vDelPhiNorm_jetmet.end(), sort_using_less_than);	
	std::sort(vDelPhiNorm_jetmht.begin(), vDelPhiNorm_jetmht.end(), sort_using_less_than);	
	
	hDelPhiMin[0]->Fill(vDelPhi_jetmet.at(0));
	hDelPhiMinNorm[0]->Fill(vDelPhiNorm_jetmet.at(0));
	hDelPhiMinNorm_mht[0]->Fill(vDelPhiNorm_jetmht.at(0));

	//Fill dPhi and dPhiN in slices of MET
	if ( met < 50 )
	{
		hDelPhiMin[1]->Fill(vDelPhi_jetmet.at(0));
		hDelPhiMinNorm[1]->Fill(vDelPhiNorm_jetmet.at(0));

	} else if ( met >= 50 && met < 100 )
	{
		hDelPhiMin[2]->Fill(vDelPhi_jetmet.at(0));
		hDelPhiMinNorm[2]->Fill(vDelPhiNorm_jetmet.at(0));
	} else if ( met >= 100 && met < 150 )
	{
		hDelPhiMin[3]->Fill(vDelPhi_jetmet.at(0));
		hDelPhiMinNorm[3]->Fill(vDelPhiNorm_jetmet.at(0));
	} else if ( met >= 150)
	{
		hDelPhiMin[4]->Fill(vDelPhi_jetmet.at(0));
		hDelPhiMinNorm[4]->Fill(vDelPhiNorm_jetmet.at(0));
	}


	//Fill dPhiN in slices of MHT
	
	if ( mht<50 ) hDelPhiMinNorm_mht[1]->Fill(vDelPhiNorm_jetmht.at(0));
	else if ( mht >= 50 && mht< 100 ) hDelPhiMinNorm_mht[2]->Fill(vDelPhiNorm_jetmht.at(0));
	else if ( mht >= 100 && mht< 200 ) hDelPhiMinNorm_mht[3]->Fill(vDelPhiNorm_jetmht.at(0));
	else if ( mht >= 200 && mht< 300 ) hDelPhiMinNorm_mht[4]->Fill(vDelPhiNorm_jetmht.at(0));
	else if ( mht >= 300 && mht< 400 ) hDelPhiMinNorm_mht[5]->Fill(vDelPhiNorm_jetmht.at(0));
	else if ( mht >= 400 && mht< 500 ) hDelPhiMinNorm_mht[6]->Fill(vDelPhiNorm_jetmht.at(0));
	else if ( mht >= 500 ) hDelPhiMinNorm_mht[7]->Fill(vDelPhiNorm_jetmht.at(0));


	//hDelPhiMinNormVshDelPhiMin->Fill(vDelPhi_jetmet.at(0),vDelPhiNorm_jetmet.at(0));
	hDelPhiMinVsMET->Fill(met, vDelPhi_jetmet.at(0));
	hDelPhiMinNormVsMET->Fill(met, vDelPhiNorm_jetmet.at(0));

	//make PASS/FAIL plots with RA2b cuts
	
	for (int cut=1; cut <= 10; ++cut)
	{
		if (vDelPhiNorm_jetmet.at(0)> cut) hPass_Norm[cut-1]->Fill(met);
		else 	hFail_Norm[cut-1]->Fill(met);


		if (vDelPhiNorm_jetmht.at(0)> cut) hPass_Norm_mht[cut-1]->Fill(mht);
		else 	hFail_Norm_mht[cut-1]->Fill(mht);
	}

	if (vDelPhi_jetmet.at(0)>0.3) 
	{
		hPass->Fill(met);
	} else 
	{
		hFail->Fill(met);
	}

/*	std::cout << "delphi ordered (jetmet)= ";
	for (std::vector<float>::const_iterator it = vDelPhi_jetmet.begin(); it != vDelPhi_jetmet.end(); ++it)
	{
		std::cout << "  " << (*it);
	}
	std::cout << std::endl;
	std::cout << "delphiNorm ordered (jetmet)= ";
	for (std::vector<float>::const_iterator it = vDelPhiNorm_jetmet.begin(); it != vDelPhiNorm_jetmet.end(); ++it)
	{
		std::cout << "  " << (*it);
	}
	std::cout << std::endl;
*/
}

void RA2bQCDvetoAna::PrintHeader()
{
	std::cout  << "======= Run/Event: " << kRunLumiEvent.run
			<< ":" << kRunLumiEvent.evt << std::endl;
}
//define this as a plug-in
DEFINE_FWK_MODULE(RA2bQCDvetoAna);

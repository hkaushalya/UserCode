// -*- C++ -*-
//
// Package:    RA2QCDvetoAna
// Class:      RA2QCDvetoAna
// 
/**\class RA2QCDvetoAna RA2QCDvetoAna.cc UserCode/RA2QCDvetoAna/src/RA2QCDvetoAna.cc

 Description: Studies the DelPhi(jet,MET) for QCD veto.

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  samantha hewamanage
//         Created:  Tue Jul 26 22:15:44 CDT 2011
// $Id: RA2QCDvetoAna.cc,v 1.5 2012/02/20 20:51:23 samantha Exp $
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

class RA2QCDvetoAna : public edm::EDFilter {
   public:
      explicit RA2QCDvetoAna(const edm::ParameterSet&);
      ~RA2QCDvetoAna();

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
		edm::InputTag mhtInputTag_, htInputTag_, metInputTag_;
		edm::Handle<edm::View<reco::MET> > mhtHandle;
		edm::Handle<double> htHandle;
		edm::Handle<std::vector<reco::PFMET> >pfMetHandle;

		struct RunLumiEvt_t
		{
			unsigned run;
			unsigned lumi;
			unsigned evt;
		};

		
		int inMinVtx;
		unsigned int iProcessed; // number of processed events
		unsigned int iPassed; //number of events passed the filter
		int inMinNdofVtx; //minimum number of dof for the vtx
		RunLumiEvt_t kRunLumiEvent;
		int iVerbose; // control print levels
		double dMaxPrimVtxZ;
		double dMaxPrimVtxRho;
		double dMinHT, dMinMHt;
		double dMinJetEt4MHt, dMaxJetEta4MHt;

		struct EventHist_t {
			TH1F* nvtx;
			TH1F* vtxz;
			TH1F* mht;
			TH1F* ht;
			TH1F* njet;
			TH1F* njet_et30eta5;
			TH1F* njet_et50eta25;
			TH1F* metphi;
			TH1F* meff;
		};

		struct JetHist_t{
			TH1F* pt;
			TH1F* eta;
			TH1F* phi;
			TH1F* delphi;  //delphi(jet,MHT)
		};

		TH1F* hDelPhiMin_mht;
		TH1F* hDelPhiMinNorm_mht;
		TH2F* hDelPhiMinVsMHT;
		TH2F* hDelPhiMinNormVshDelPhiMin_mht;
		TH2F* hDelPhiMinNormVsMHT;
		TH1F* hPassFail_mht;
		TH1F* hFail_mht;
		TH1F* hPassFail_Norm_mht;
		TH1F* hPass_Norm_mht;
		TH1F* hFail_Norm_mht;

		TH2F* hDelPhiMinNormVshDelPhiMin_met;
		TH2F* hDelPhiMinVsMET;
		TH2F* hDelPhiMinNormVsMET;
		TH1F* hPassFail_met;
		TH1F* hFail_met;
		TH1F* hPassFail_Norm_met;
		TH1F* hPass_Norm_met;
		TH1F* hFail_Norm_met;


		TH1F* hDelPhiMin_bymhtSlice[5];
		TH1F* hDelPhiMinNorm_bymhtSlice[5];
		TH1F* hDelPhiMin_bymetSlice[5];
		TH1F* hDelPhiMinNorm_bymetSlice[5];


		TH1F* myMht;
		TH1F* myMet;
		TH1F* metsig_delphi_pass;
		TH1F* metsig_delphi_fail;
		TH1F* metsig_delphinorm_pass;
		TH1F* metsig_delphinorm_fail;


		EventHist_t evtHist, evtHist_b4;
		JetHist_t pf30_jet1Hist, pf30_jet2Hist, pf30_jet3Hist; //for upto 5 jets may be

		JetHist_t pf50eta25_jet1Hist, pf50eta25_jet2Hist, pf50eta25_jet3Hist; //for upto 5 jets may be

		edm::InputTag patJetsPFPt30InputTag_;
		edm::InputTag patJetsPFPt50Eta25InputTag_;
		edm::Handle<std::vector<pat::Jet> > pfpt30JetHandle;
		edm::Handle<std::vector<pat::Jet> > pfpt50Eta25JetHandle;

		void DoDelMinStudy(edm::Event& iEvent, edm::Handle<std::vector<pat::Jet> > pt30jetHandle);
		bool hasMinNjetHtMhtFromMyJets(edm::Handle<std::vector<pat::Jet> > jetHandle);
		bool hasValidVertex(edm::Handle<reco::VertexCollection> vertexHandle);
		void PrintHeader();
		void FillJetHists(JetHist_t& hist, const TLorentzVector& jetVec, const float met_phi);
		TLorentzVector vMHtVec;
		double dHt, dMEff;
		unsigned int fail_njet, fail_ht, fail_mht, fail_nvtx;

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
RA2QCDvetoAna::RA2QCDvetoAna(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
	patJetsPFPt30InputTag_ = iConfig.getParameter<edm::InputTag>("patJetsPFPt30InputTag");
	patJetsPFPt50Eta25InputTag_ = iConfig.getParameter<edm::InputTag>("patJetsPFPt50Eta25InputTag");
	mhtInputTag_ = iConfig.getParameter<edm::InputTag>("mhtInputTag");
	htInputTag_ = iConfig.getParameter<edm::InputTag>("htInputTag");
	inMinVtx = iConfig.getUntrackedParameter<int>("nMinVtx",1);
	dMinJetEt4MHt = iConfig.getUntrackedParameter<double>("minJetEt4MHt",50.0);
	dMaxJetEta4MHt = iConfig.getUntrackedParameter<double>("maxJetEta4MHt",5.0);
	dMinHT = iConfig.getUntrackedParameter<double>("dMinHT",350.0);
	dMinMHt = iConfig.getUntrackedParameter<double>("dMinMHt",0.0);
	inMinNdofVtx = iConfig.getUntrackedParameter<int>("ndofVtx",0);
	dMaxPrimVtxZ = iConfig.getUntrackedParameter<double>("maxDelzVtx",24.0);
	dMaxPrimVtxRho = iConfig.getUntrackedParameter<double>("maxDelRho",15.0);
	iVerbose = iConfig.getUntrackedParameter<int>("verbose",0);
	iProcessed = 0;
	iPassed = 0;
	vMHtVec.SetPxPyPzE(0,0,0,0);
	dHt = 0;
	dMEff = 0;
	fail_njet = 0;
	fail_ht = 0;
	fail_mht = 0;
	fail_nvtx = 0;

	//generate hists
	edm::Service<TFileService> fs;
	const double evt_met_max = 4000, evt_met_bins = 400;
	const double evt_ht_max = 6000, evt_ht_bins = 120;

	//hists before any cuts
	evtHist_b4.nvtx = fs->make<TH1F> ("pat_nvtx_b4" ,"PAT-tuple quantities for debugging: before any cut;Nvtx;Events;", 15, 0, 15);
	evtHist_b4.vtxz = fs->make<TH1F> ("pat_vtxz_b4","PAT-tuple quantities for debugging: before any cut;Primary Vertex z [cm];Arbitrary;",30,0,30);
	evtHist_b4.mht = fs->make<TH1F> ("pat_mht_b4" ,"PAT-tuple quantities for debugging: before any cut;MHT [GeV];Events;", evt_met_bins, 0, evt_met_max);
	evtHist_b4.ht = fs->make<TH1F> ("pat_ht_b4" ,"PAT-tuple quantities for debugging: before any cut;HT [GeV];Events;", evt_ht_bins, 0, evt_ht_max);
	evtHist_b4.njet = fs->make<TH1F> ("pat_njet_b4" ,"PAT-tuple quantities for debugging: before any cut;NJET [GeV];Events;", 10, 0, 10);
	evtHist_b4.njet_et50eta25 = fs->make<TH1F> ("pat_njet_et50eta25_b4" ,"PAT-tuple quantities for debugging: before any cut;NJET (Et50Eta25);Events;", 10, 0, 10);
	evtHist_b4.metphi = fs->make<TH1F> ("pat_metphi_b4" ,"PAT-tuple quantities for debugging: before any cut;met phi;Events;", 160, -8, 8);
	evtHist_b4.meff = fs->make<TH1F> ("pat_meff_b4" ,"PAT-tuple quantities for debugging: before any cut;MEff;Events;", evt_ht_bins, 0, evt_ht_max);

	//these are general event hist to check the PAT tuples cuts
	evtHist.nvtx = fs->make<TH1F> ("pat_nvtx" ,"PAT-tuple quantities for debugging;Nvtx;Events;", 15, 0, 15);
	evtHist.vtxz = fs->make<TH1F> ("pat_vtxz","PAT-tuple quantities for debugging;Primary Vertex z [cm];Arbitrary;",30,0,30);
	evtHist.mht = fs->make<TH1F> ("pat_mht" ,"PAT-tuple quantities for debugging;MHT [GeV];Events;", evt_met_bins, 0, evt_met_max);
	evtHist.ht = fs->make<TH1F> ("pat_ht" ,"PAT-tuple quantities for debugging;HT [GeV];Events;", evt_ht_bins, 0, evt_ht_max);
	evtHist.njet = fs->make<TH1F> ("pat_njet" ,"PAT-tuple quantities for debugging;NJET [GeV];Events;", 10, 0, 10);
	evtHist.njet_et50eta25 = fs->make<TH1F> ("pat_njet_et50eta25" ,"PAT-tuple quantities for debugging: before any cut;NJET (Et50Eta25);Events;", 10, 0, 10);
	evtHist.metphi = fs->make<TH1F> ("pat_metphi" ,"PAT-tuple quantities for debugging;met phi;Events;", 160, -8, 8);
	evtHist.meff = fs->make<TH1F> ("pat_meff" ,"PAT-tuple quantities for debugging;MEff;Events;", evt_ht_bins, 0, evt_ht_max);




	const double pt_bins = 800, pt_max = 4000;
	pf30_jet1Hist.pt = fs->make<TH1F> ("pat_pf30_jet1_pt" ,"PAT-tuple quantities for debugging: PF30-Jet1 pt;pt [GeV];Events;", pt_bins, 0, pt_max);
	pf30_jet1Hist.eta = fs->make<TH1F> ("pat_pf30_jet1_eta" ,"PAT-tuple quantities for debugging: PF30-Jet1 eta;eta;Events;", 100, -5, 5);
	pf30_jet1Hist.phi = fs->make<TH1F> ("pat_pf30_jet1_phi" ,"PAT-tuple quantities for debugging: PF30-Jet1 phi;phi;Events;", 160, -8, 8);
	pf30_jet1Hist.delphi = fs->make<TH1F> ("pat_pf30_jet1_delphi" ,"PAT-tuple quantities for debugging: PF30-Jet1: delphi;delphi;Events;", 160, -8, 8);

	pf30_jet2Hist.pt = fs->make<TH1F> ("pat_pf30_jet2_pt" ,"PAT-tuple quantities for debugging: PF30-Jet2 pt;pt [GeV];Events;", pt_bins, 0, pt_max);
	pf30_jet2Hist.eta = fs->make<TH1F> ("pat_pf30_jet2_eta" ,"PAT-tuple quantities for debugging: PF30-Jet2 eta;eta;Events;", 100, -5, 5);
	pf30_jet2Hist.phi = fs->make<TH1F> ("pat_pf30_jet2_phi" ,"PAT-tuple quantities for debugging: PF30-Jet2 phi;phi;Events;", 160, -8, 8);
	pf30_jet2Hist.delphi = fs->make<TH1F> ("pat_pf30_jet2_delphi" ,"PAT-tuple quantities for debugging: PF30-Jet2: delphi;delphi;Events;", 160, -8, 8);

	pf30_jet3Hist.pt = fs->make<TH1F> ("pat_pf30_jet3_pt" ,"PAT-tuple quantities for debugging: PF30-Jet3 pt;pt [GeV];Events;", pt_bins, 0, pt_max);
	pf30_jet3Hist.eta = fs->make<TH1F> ("pat_pf30_jet3_eta" ,"PAT-tuple quantities for debugging: PF30-Jet3 eta;eta;Events;", 100, -5, 5);
	pf30_jet3Hist.phi = fs->make<TH1F> ("pat_pf30_jet3_phi" ,"PAT-tuple quantities for debugging: PF30-Jet3 phi;phi;Events;", 160, -8, 8);
	pf30_jet3Hist.delphi = fs->make<TH1F> ("pat_pf30_jet3_delphi" ,"PAT-tuple quantities for debugging: PF30-Jet3: delphi;delphi;Events;", 160, -8, 8);


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


	const float npassFailHistBins = 16;
	const float passFailHistBins[] = {0,20,40,60,80,100,120,140,160,180,200,250,300,350,400,500,600};

//	hist_delphiMin_jetmet[0] = fs->make<TH1F> ("delphiMin_jetmet" ,"delphi Min.;delphi min;Events;", 160, -8, 8);
//	hist_delphiMin_jetmetVsMHT[0] = fs->make<TH2F> ("delphiMin_jetmetVsMHT" ,"delphi min Vs MHT;MHT;delphi min;", 1000, 0, 1000,400,0,4);

	hDelPhiMin_mht = fs->make<TH1F> ("delPhiMin_mht","delPhiMin", 400, 0, 4);
	hDelPhiMinNorm_mht = fs->make<TH1F> ("delPhiMinNorm_mht","delPhiMinNorm", 2000, 0, 20);
	hDelPhiMinNormVshDelPhiMin_mht = fs->make<TH2F> ("delPhiMinNormVsdelPhiMin_mht","delPhiMinNormalized VS delPhiMin", 500,0,5, 2000, 0, 20);
	hDelPhiMinVsMHT = fs->make<TH2F> ("delPhiMinVsMHT","delPhiMin VS MHT", 150,0,1500, 400, 0, 10);
	hDelPhiMinNormVsMHT = fs->make<TH2F> ("delPhiMinNormVsMHT","delPhiMinNorm VS MET",150,0,1500, 400, 0, 20);
	//hPassFail_mht = fs->make<TH1F> ("hPassFail_mht","PASS/FAIL",30, 0, 600);
	hPassFail_mht = fs->make<TH1F> ("hPassFail_mht","PASS/FAIL", npassFailHistBins, passFailHistBins);
	hPassFail_mht->Sumw2();
	hFail_mht = fs->make<TH1F> ("hFail_mht","FAIL",npassFailHistBins, passFailHistBins);
	hFail_mht->Sumw2();
	hPass_Norm_mht = fs->make<TH1F> ("hPass_Norm_mht","Normalized: PASS",npassFailHistBins, passFailHistBins);
	hPass_Norm_mht->Sumw2();
	hPassFail_Norm_mht = fs->make<TH1F> ("hPassFail_Norm_mht","Normalized: PASS/FAIL;MHT;Events;",npassFailHistBins, passFailHistBins);
	hPassFail_Norm_mht->Sumw2();
	hFail_Norm_mht = fs->make<TH1F> ("hFail_Norm_mht","Normalized: FAIL;MHT;Events;",npassFailHistBins, passFailHistBins);
	hFail_Norm_mht->Sumw2();



	hDelPhiMinNormVshDelPhiMin_met = fs->make<TH2F> ("delPhiMinNormVsdelPhiMin_met","delPhiMinNormalized VS delPhiMin", 500,0,5, 2000, 0, 20);
	hDelPhiMinVsMET = fs->make<TH2F> ("delPhiMinVsMET","delPhiMin VS MET", 150,0,1500, 400, 0, 10);
	hDelPhiMinNormVsMET = fs->make<TH2F> ("delPhiMinNormVsMET","delPhiMinNorm VS MET",150,0,1500, 400, 0, 20);
	hPassFail_met = fs->make<TH1F> ("hPassFail_met","PASS/FAIL", npassFailHistBins, passFailHistBins);
	hPassFail_met->Sumw2();
	hFail_met = fs->make<TH1F> ("hFail_met","FAIL", npassFailHistBins, passFailHistBins);
	hFail_met->Sumw2();
	hPass_Norm_met = fs->make<TH1F> ("hPass_Norm_met","Normalized: PASS", npassFailHistBins, passFailHistBins);
	hPass_Norm_met->Sumw2();
	hPassFail_Norm_met = fs->make<TH1F> ("hPassFail_Norm_met","Normalized: PASS/FAIL;MET;Events;", npassFailHistBins, passFailHistBins);
	hPassFail_Norm_met->Sumw2();
	hFail_Norm_met = fs->make<TH1F> ("hFail_Norm_met","Normalized: FAIL;MET;Events;", npassFailHistBins, passFailHistBins);
	hFail_Norm_met->Sumw2();



	hDelPhiMin_bymhtSlice[0] = fs->make<TH1F> ("delPhiMin_bymhtSlice0","delPhiMin (50<MHT<100)", 400, 0, 4);
	hDelPhiMin_bymhtSlice[1] = fs->make<TH1F> ("delPhiMin_bymhtSlice1","delPhiMin (100<MHT<150)", 400, 0, 4);
	hDelPhiMin_bymhtSlice[2] = fs->make<TH1F> ("delPhiMin_bymhtSlice2","delPhiMin (150<MHT<200)", 400, 0, 4);
	hDelPhiMin_bymhtSlice[3] = fs->make<TH1F> ("delPhiMin_bymhtSlice3","delPhiMin (200<MHT<250)", 400, 0, 4);
	hDelPhiMin_bymhtSlice[4] = fs->make<TH1F> ("delPhiMin_bymhtSlice4","delPhiMin (MHT>250)", 400, 0, 4);
	hDelPhiMinNorm_bymhtSlice[0] = fs->make<TH1F> ("delPhiMinNorm_bymhtSlice0","delPhiMinNorm (50<MHT<100)", 2000, 0, 20);
	hDelPhiMinNorm_bymhtSlice[1] = fs->make<TH1F> ("delPhiMinNorm_bymhtSlice1","delPhiMinNorm (100<MHT<150)", 2000, 0, 20);
	hDelPhiMinNorm_bymhtSlice[2] = fs->make<TH1F> ("delPhiMinNorm_bymhtSlice2","delPhiMinNorm (150<MHT<200)", 2000, 0, 20);
	hDelPhiMinNorm_bymhtSlice[3] = fs->make<TH1F> ("delPhiMinNorm_bymhtSlice3","delPhiMinNorm (200<MHT<250)", 2000, 0, 20);
	hDelPhiMinNorm_bymhtSlice[4] = fs->make<TH1F> ("delPhiMinNorm_bymhtSlice4","delPhiMinNorm (MHT>250)", 2000, 0, 20);

	hDelPhiMin_bymetSlice[0] = fs->make<TH1F> ("delPhiMin_bymetSlice0","delPhiMin (50<MET<100)", 400, 0, 4);
	hDelPhiMin_bymetSlice[1] = fs->make<TH1F> ("delPhiMin_bymetSlice1","delPhiMin (100<MET<150)", 400, 0, 4);
	hDelPhiMin_bymetSlice[2] = fs->make<TH1F> ("delPhiMin_bymetSlice2","delPhiMin (150<MET<200)", 400, 0, 4);
	hDelPhiMin_bymetSlice[3] = fs->make<TH1F> ("delPhiMin_bymetSlice3","delPhiMin (200<MET<250)", 400, 0, 4);
	hDelPhiMin_bymetSlice[4] = fs->make<TH1F> ("delPhiMin_bymetSlice4","delPhiMin (MET>250)", 400, 0, 4);
	hDelPhiMinNorm_bymetSlice[0] = fs->make<TH1F> ("delPhiMinNorm_bymetSlice0","delPhiMinNorm (50<MET<100)", 2000, 0, 20);
	hDelPhiMinNorm_bymetSlice[1] = fs->make<TH1F> ("delPhiMinNorm_bymetSlice1","delPhiMinNorm (100<MET<150)", 2000, 0, 20);
	hDelPhiMinNorm_bymetSlice[2] = fs->make<TH1F> ("delPhiMinNorm_bymetSlice2","delPhiMinNorm (150<MET<200)", 2000, 0, 20);
	hDelPhiMinNorm_bymetSlice[3] = fs->make<TH1F> ("delPhiMinNorm_bymetSlice3","delPhiMinNorm (200<MET<250)", 2000, 0, 20);
	hDelPhiMinNorm_bymetSlice[4] = fs->make<TH1F> ("delPhiMinNorm_bymetSlice4","delPhiMinNorm (MET>250)", 2000, 0, 20);


	myMht = fs->make<TH1F> ("my_mht" ,"My MHT;MHT [GeV];Events;", evt_met_bins, 0, evt_met_max);
	myMet = fs->make<TH1F> ("my_met" ,"My MET;MET [GeV];Events;", evt_met_bins, 0, evt_met_max);

	metsig_delphi_pass = fs->make<TH1F> ("metsig_delphi_pass" ,"METsig of events pass from #Delta#Phi_{min};METsig;Events;", 100, 0, 10);
	metsig_delphi_fail = fs->make<TH1F> ("metsig_delphi_fail" ,"METsig of events fail from #Delta#Phi_{min};METsig;Events;", 100, 0, 10);
	metsig_delphinorm_pass = fs->make<TH1F> ("metsig_delphinorm_pass" ,"METsig of events pass from Normalized #Delta#Phi_{min};METsig;Events;", 100, 0, 10);
	metsig_delphinorm_fail = fs->make<TH1F> ("metsig_delphinorm_fail" ,"METsig of events fail from Normalized #Delta#Phi_{min};METsig;Events;", 100, 0, 10);
}


RA2QCDvetoAna::~RA2QCDvetoAna()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
RA2QCDvetoAna::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
	++iProcessed;
	kRunLumiEvent.run  = iEvent.id().run();
	kRunLumiEvent.lumi = iEvent.id().luminosityBlock(); 
	kRunLumiEvent.evt  = iEvent.id().event();

	
	if (iVerbose) PrintHeader();

	Handle<reco::VertexCollection> vertexHandle;
	iEvent.getByLabel("goodVerticesRA2", vertexHandle);
   if (! vertexHandle.isValid())
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":vertexHandle handle not found!" << std::endl;
		assert(false);
	}

	iEvent.getByLabel(htInputTag_, htHandle);
	if (! htHandle.isValid()) 
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":htHandle handle not found!" << std::endl;
		assert(false);
	}
	
 	// ==========================================================
	// MHT Information 
	iEvent.getByLabel(mhtInputTag_, mhtHandle);
	if (! mhtHandle.isValid()) 
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":MHT handle not found!" << std::endl;
		assert(false);
	}


	iEvent.getByLabel(patJetsPFPt30InputTag_, pfpt30JetHandle);
	if (! pfpt30JetHandle.isValid()) 
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":pfpt30JetHandle handle not found!" << std::endl;
		assert(false);
	}

	iEvent.getByLabel(patJetsPFPt50Eta25InputTag_, pfpt50Eta25JetHandle);
	if (! pfpt50Eta25JetHandle.isValid()) 
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":pfpt50Eta25JetHandle handle not found!" << std::endl;
		assert(false);
	}
	
	/*************
	 *before any cuts fill the hists to check
	 * PAT contents
	 *************/
	evtHist_b4.metphi->Fill((*mhtHandle)[0].phi());
	evtHist_b4.mht->Fill( (*mhtHandle)[0].pt());
	evtHist_b4.ht->Fill( (*htHandle) );
	evtHist_b4.njet_et50eta25->Fill(pfpt50Eta25JetHandle->size());
	evtHist_b4.nvtx->Fill(vertexHandle->size());
	evtHist_b4.meff->Fill((*mhtHandle)[0].pt() + (*htHandle));

	float mht_phi = (*mhtHandle)[0].phi();	
	//int nvtx = vertexHandle->size();
	//int njets_et50eta25 = pfpt50Eta25JetHandle->size();
	for (unsigned i = 0 ; i < pfpt50Eta25JetHandle->size() ; ++i)
	{
		const TLorentzVector iJetVec((*pfpt50Eta25JetHandle)[i].px(),(*pfpt50Eta25JetHandle)[i].py(),(*pfpt50Eta25JetHandle)[i].pz(),(*pfpt50Eta25JetHandle)[i].energy());
		if (i==0) FillJetHists(pf50eta25_jet1Hist, iJetVec, mht_phi);
		else if (i==1) FillJetHists(pf50eta25_jet2Hist, iJetVec, mht_phi);
		else if (i==2) FillJetHists(pf50eta25_jet3Hist, iJetVec, mht_phi);
	}


	//basic event selection
	if (! hasValidVertex(vertexHandle)) { ++fail_nvtx; return 0; }
	if (! hasMinNjetHtMhtFromMyJets(pfpt30JetHandle))
	{
		if (iVerbose) std::cout << __FUNCTION__ << ": Fail NjetMht!" << std::endl;
		return 0;
	}

	evtHist.metphi->Fill((*mhtHandle)[0].phi());
	evtHist.mht->Fill( (*mhtHandle)[0].pt());
	evtHist.ht->Fill( (*htHandle) );
	evtHist.meff->Fill((*mhtHandle)[0].pt() + (*htHandle));

	DoDelMinStudy(iEvent, pfpt30JetHandle);

	++iPassed;
   return true;
}

// ------------ method called once each job just before starting event loop  ------------
void 
RA2QCDvetoAna::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
RA2QCDvetoAna::endJob() {
	hPassFail_mht->Divide(hFail_mht);
	hPassFail_Norm_mht->Divide(hFail_Norm_mht);
	hPassFail_met->Divide(hFail_met);
	hPassFail_Norm_met->Divide(hFail_Norm_met);

	std::cout << __FILE__ << "::" << __FUNCTION__  << "\n" << std::endl;
	std::cout << "[ATM:00] Job settings " << std::endl;
	std::cout << "[ATM:01] Events Processed --- = " << iProcessed << std::endl;
	std::cout << "[ATM:02] Events Passed ------ = " << iPassed << std::endl;
	std::cout << "[ATM:10] Primary Vtx Min ---- = " << inMinVtx << std::endl; 
	std::cout << "[ATM:11] Primary Vtx Min ndof = " << inMinNdofVtx << std::endl; 
	std::cout << "[ATM:12] Primary Vtx Max z -- = " << dMaxPrimVtxZ << std::endl; 
	std::cout << "[ATM:13] Primary Vtx Max rho  = " << dMaxPrimVtxRho << std::endl; 

	std::cout << "[ATM:14] MHT: Min Jet Et ---- = " << dMinJetEt4MHt << std::endl; 
	std::cout << "[ATM:15] MHT: Max Jet Eta --- = " << dMaxJetEta4MHt << std::endl; 
	std::cout << "[ATM:15] Primary Vtx Max z -- = " << dMaxPrimVtxZ << std::endl; 
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
	std::cout << "[ATM:31] Minimum HT --------- = " << dMinHT << std::endl;
	std::cout << "[ATM:32] Minimum MHT -------- = " << dMinMHt << std::endl;
	std::cout << "[ATM:33] Filter Summary------------- " << std::endl;
	std::cout << "[ATM:34] Pass nvtx            = " << iProcessed - fail_nvtx << " (" << fail_nvtx << ")" << std::endl;
	std::cout << "[ATM:35] Pass njet            = " << iProcessed - fail_nvtx - fail_njet 
									<< " (" << fail_njet << ")" << std::endl;
	std::cout << "[ATM:36] Pass mht             = " << iProcessed - fail_nvtx - fail_njet - fail_mht
									<< " (" << fail_mht << ")" << std::endl;
	std::cout << "[ATM:37] Pass ht             = " << iProcessed - fail_nvtx - fail_njet - fail_mht - fail_ht
									<< " (" << fail_ht << ")" << std::endl;
}

// ------------ method called when starting to processes a run  ------------
bool 
RA2QCDvetoAna::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
RA2QCDvetoAna::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
RA2QCDvetoAna::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
RA2QCDvetoAna::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
RA2QCDvetoAna::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void RA2QCDvetoAna::DoDelMinStudy(edm::Event& iEvent,
		edm::Handle<std::vector<pat::Jet> > jetHandle)
{
	//PrintHeader();
	if (iVerbose) std::cout <<  "In "<< __FUNCTION__ << std::endl;

	iEvent.getByLabel("pfMet",pfMetHandle);
	if (! pfMetHandle.isValid())
	{
		std::cout << "Valid PFMET Collection Not Found!" << std::endl;
		assert(false);
	}

	std::vector<float> vDelPhi_jetmht, vDelPhiNorm_jetmht;
	std::vector<float> vDelPhi_jetmet, vDelPhiNorm_jetmet;

	const float mht = vMHtVec.Pt();
	const float met = (*pfMetHandle)[0].pt();
	const float	metsig = (*pfMetHandle)[0].mEtSig();      

	unsigned njet = 0; //consider cross product only from lead 3 jets

	for (unsigned i = 0 ; i < jetHandle->size() ; ++i)
	{
		const TLorentzVector iJetVec((*jetHandle)[i].px(),(*jetHandle)[i].py(),(*jetHandle)[i].pz(),(*jetHandle)[i].energy());
		if (iJetVec.Pt()<50. || fabs(iJetVec.Eta())>2.5) continue;
		if (iJetVec.Pt()<dMinJetEt4MHt || fabs(iJetVec.Eta())>dMaxJetEta4MHt) continue;
		if (iVerbose) std::cout << __LINE__ << ": pass jet et/eta cuts [" << i << std::endl; 
		//loose jet id stuff
/*		if ((*jetHandle)[i].neutralHadronEnergyFraction() >=0.99) continue;
		if ((*jetHandle)[i].neutralEmEnergyFraction() >=0.99) continue;
		if (((*jetHandle)[i].getPFConstituents()).size() <=1) continue;
		if ((*jetHandle)[i].chargedHadronEnergyFraction() <=0) continue;
		if ((*jetHandle)[i].chargedMultiplicity() <=0) continue;
		if ((*jetHandle)[i].chargedEmEnergyFraction() >=0.99) continue;
*/

		++njet;
		if (njet <=3)
		{
			const float delphi_jetmet = fabs(TVector2::Phi_mpi_pi((*jetHandle)[i].phi() - (*pfMetHandle)[0].phi()));
			const float delphi_jetmht = fabs(TVector2::Phi_mpi_pi((*jetHandle)[i].phi() - vMHtVec.Phi()));
			vDelPhi_jetmht.push_back(delphi_jetmht);
			vDelPhi_jetmet.push_back(delphi_jetmet);
			//loop over rest of the jets to calcualte deltaT for i th jet
			float sumCP2 = 0;  //sqaure of sum of cross products
			for (unsigned j = 0 ; j < jetHandle->size() ; ++j)
			{
				if ( i == j ) continue;
				//if (iJetVec.Pt()<30. || fabs(iJetVec.Eta())>2.5) continue;
				if (iJetVec.Pt()<30) continue;
				const TLorentzVector jJetVec((*jetHandle)[j].px(),(*jetHandle)[j].py(),(*jetHandle)[j].pz(),(*jetHandle)[j].energy());
				sumCP2 += pow(iJetVec.Px() * jJetVec.Py() + iJetVec.Py() * jJetVec.Px() ,2);
			}
			const float delT_i = (0.1*sqrt(sumCP2))/iJetVec.Pt();
			const float delPhi_mht_i = delphi_jetmht/atan2(delT_i,mht);
			const float delPhi_met_i = delphi_jetmet/atan2(delT_i,met);
			vDelPhiNorm_jetmht.push_back(delPhi_mht_i);
			vDelPhiNorm_jetmet.push_back(delPhi_met_i);
			if (iVerbose) std::cout << "jet/dphi/dphi_norm = " << i << " | " << delphi_jetmht
						<< " | " << delPhi_mht_i << std::endl;
		}

	}

	std::sort(vDelPhi_jetmht.begin(), vDelPhi_jetmht.end(), sort_using_less_than);	
	std::sort(vDelPhiNorm_jetmht.begin(), vDelPhiNorm_jetmht.end(), sort_using_less_than);	
	std::sort(vDelPhi_jetmet.begin(), vDelPhi_jetmet.end(), sort_using_less_than);	
	std::sort(vDelPhiNorm_jetmet.begin(), vDelPhiNorm_jetmet.end(), sort_using_less_than);	
	
	if (vDelPhi_jetmht.size())
	{
		if (iVerbose) std::cout << __FUNCTION__ << ":" << __LINE__ 
				<< ": Filling Histograms." << std::endl;

		myMht->Fill(mht);
		myMet->Fill(met);

		if (mht>50 && mht<100)
		{
			hDelPhiMin_bymhtSlice[0]->Fill(vDelPhi_jetmht.at(0));
			hDelPhiMinNorm_bymhtSlice[0]->Fill(vDelPhiNorm_jetmht.at(0));
		} else if (mht>100 && mht < 150)
		{
			hDelPhiMin_bymhtSlice[1]->Fill(vDelPhi_jetmht.at(0));
			hDelPhiMinNorm_bymhtSlice[1]->Fill(vDelPhiNorm_jetmht.at(0));
		} else if (mht>150 && mht < 200)
		{
			hDelPhiMin_bymhtSlice[2]->Fill(vDelPhi_jetmht.at(0));
			hDelPhiMinNorm_bymhtSlice[2]->Fill(vDelPhiNorm_jetmht.at(0));
		} else if (mht>200 && mht < 250)
		{
			hDelPhiMin_bymhtSlice[3]->Fill(vDelPhi_jetmht.at(0));
			hDelPhiMinNorm_bymhtSlice[3]->Fill(vDelPhiNorm_jetmht.at(0));
		} else if (mht>250)
		{
			hDelPhiMin_bymhtSlice[4]->Fill(vDelPhi_jetmht.at(0));
			hDelPhiMinNorm_bymhtSlice[4]->Fill(vDelPhiNorm_jetmht.at(0));
		}

		if (met>50 && met<100)
		{
			hDelPhiMin_bymetSlice[0]->Fill(vDelPhi_jetmet.at(0));
			hDelPhiMinNorm_bymetSlice[0]->Fill(vDelPhiNorm_jetmet.at(0));
		} else if (met>100 && met < 150)
		{
			hDelPhiMin_bymetSlice[1]->Fill(vDelPhi_jetmet.at(0));
			hDelPhiMinNorm_bymetSlice[1]->Fill(vDelPhiNorm_jetmet.at(0));
		} else if (met>150 && met < 200)
		{
			hDelPhiMin_bymetSlice[2]->Fill(vDelPhi_jetmet.at(0));
			hDelPhiMinNorm_bymetSlice[2]->Fill(vDelPhiNorm_jetmet.at(0));
		} else if (met>200 && met < 250)
		{
			hDelPhiMin_bymetSlice[3]->Fill(vDelPhi_jetmet.at(0));
			hDelPhiMinNorm_bymetSlice[3]->Fill(vDelPhiNorm_jetmet.at(0));
		} else if (met>250)
		{
			hDelPhiMin_bymetSlice[4]->Fill(vDelPhi_jetmet.at(0));
			hDelPhiMinNorm_bymetSlice[4]->Fill(vDelPhiNorm_jetmet.at(0));
		}


		hDelPhiMin_mht->Fill(vDelPhi_jetmht.at(0));
		hDelPhiMinNorm_mht->Fill(vDelPhiNorm_jetmht.at(0));
		hDelPhiMinNormVshDelPhiMin_mht->Fill(vDelPhi_jetmht.at(0),vDelPhiNorm_jetmht.at(0));
		hDelPhiMinVsMHT->Fill(mht, vDelPhi_jetmht.at(0));
		hDelPhiMinNormVsMHT->Fill(mht, vDelPhiNorm_jetmht.at(0));


		hDelPhiMinNormVshDelPhiMin_met->Fill(vDelPhi_jetmet.at(0),vDelPhiNorm_jetmet.at(0));
		hDelPhiMinVsMET->Fill(mht, vDelPhi_jetmet.at(0));
		hDelPhiMinNormVsMET->Fill(met, vDelPhiNorm_jetmet.at(0));

		if (vDelPhiNorm_jetmht.at(0)>2) 
		{
			hPassFail_Norm_mht->Fill(mht);
			hPass_Norm_mht->Fill(mht);
		} else 
		{
			hFail_Norm_mht->Fill(mht);
		}

		if (vDelPhi_jetmht.at(0)>0.3) 
		{
			hPassFail_mht->Fill(mht);
		} else 
		{
			hFail_mht->Fill(mht);
		}


		if (vDelPhiNorm_jetmet.at(0)>3) 
		{
			hPassFail_Norm_met->Fill(met);
			hPass_Norm_met->Fill(met);
			metsig_delphinorm_pass->Fill(metsig);
		} else 
		{
			hFail_Norm_met->Fill(mht);
			metsig_delphinorm_fail->Fill(metsig);
		}

		if (vDelPhi_jetmet.at(0)>0.3) 
		{
			hPassFail_met->Fill(mht);
			metsig_delphi_pass->Fill(metsig);
		} else 
		{
			hFail_met->Fill(mht);
			metsig_delphi_fail->Fill(metsig);
		}

	}
/*	std::cout << "delphi ordered (jetmet)= ";
	for (std::vector<float>::const_iterator it = vDelPhi_jetmht.begin(); it != vDelPhi_jetmht.end(); ++it)
	{
		std::cout << "  " << (*it);
	}
	std::cout << std::endl;
	std::cout << "delphiNorm ordered (jetmet)= ";
	for (std::vector<float>::const_iterator it = vDelPhiNorm_jetmht.begin(); it != vDelPhiNorm_jetmht.end(); ++it)
	{
		std::cout << "  " << (*it);
	}
	std::cout << std::endl;
*/
}

void RA2QCDvetoAna::PrintHeader()
{
	std::cout  << "======= Run/Event: " << kRunLumiEvent.run
			<< ":" << kRunLumiEvent.evt << std::endl;
}

bool RA2QCDvetoAna::hasValidVertex(edm::Handle<reco::VertexCollection> vertexHandle)
{
	double dNvtx = 0;
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
			const math::XYZPoint vPoint = v->position();
			if (vPoint.Rho() >= dMaxPrimVtxRho) continue;
			if (iVerbose) std::cout << "Valid Vertex Found " << std::endl;
			++dNvtx;
			if (dNvtx == 1) evtHist.vtxz->Fill(dPrimVtx_z);
		}
	}

	if (dNvtx >= inMinVtx)
	{
		evtHist.nvtx->Fill(dNvtx);
		return 1;
	} else return 0;
}

bool RA2QCDvetoAna::hasMinNjetHtMhtFromMyJets(edm::Handle<std::vector<pat::Jet> > jetHandle)
{
	unsigned njets = 0;
	//unsigned nextrajets=0;
	TLorentzVector vSumJetVecMHt(0,0,0,0);
	for (unsigned i = 0 ; i < jetHandle->size() ; ++i)
	{
		const TLorentzVector iJetVec((*jetHandle)[i].px(),(*jetHandle)[i].py(),(*jetHandle)[i].pz(),(*jetHandle)[i].energy());
		if (iJetVec.Pt()>dMinJetEt4MHt && fabs(iJetVec.Eta())<dMaxJetEta4MHt)
		//if (iJetVec.Pt()>50 && fabs(iJetVec.Eta())<2.5)
		{
			dHt += iJetVec.Pt();
			++njets;
		}
		if (iJetVec.Pt()>30 && fabs(iJetVec.Eta())<5)
		{	
			vSumJetVecMHt += iJetVec;
		}
		//loose jet id stuff
		/*if ((*jetHandle)[i].neutralHadronEnergyFraction() >=0.99) continue;
		if ((*jetHandle)[i].neutralEmEnergyFraction() >=0.99) continue;
		if (((*jetHandle)[i].getPFConstituents()).size() <=1) continue;
		if ((*jetHandle)[i].chargedHadronEnergyFraction() <=0) continue;
		if ((*jetHandle)[i].chargedMultiplicity() <=0) continue;
		if ((*jetHandle)[i].chargedEmEnergyFraction() >=0.99) continue;
		*/
	}

	if (njets<3) { ++fail_njet; return 0; }
	if (vSumJetVecMHt.Pt() < dMinMHt) { ++fail_mht; return 0; }
	if (dHt < dMinHT) { ++fail_ht; return 0; }
	
	dMEff = dHt + vSumJetVecMHt.Pt();
	vMHtVec = vSumJetVecMHt;
	return 1;
}


void RA2QCDvetoAna::FillJetHists(JetHist_t& hist, const TLorentzVector& jetVec, const float met_phi)
{
	hist.pt->Fill(jetVec.Pt());
	hist.eta->Fill(jetVec.Eta());
	hist.phi->Fill(jetVec.Phi());
	hist.delphi->Fill(fabs(TVector2::Phi_mpi_pi(jetVec.Phi() - met_phi)));
}

//define this as a plug-in
DEFINE_FWK_MODULE(RA2QCDvetoAna);

// -*- C++ -*-
//
// Package:    AnomTrkMETFilter
// Class:      AnomTrkMETFilter
// 
/**\class AnomTrkMETFilter AnomTrkMETFilter.cc AnomTrkFilter/AnomTrkMETFilter/src/AnomTrkMETFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  samantha hewamanage
//         Created:  Tue Jul  5 17:28:00 CDT 2011
// $Id: AnomTrkMETFilter.cc,v 1.4 2011/07/19 16:22:26 samantha Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
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

//ROOT
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "assert.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "TVector3.h"
#include <iomanip>
#include "TPad.h"
//
// class declaration
//

class AnomTrkMETFilter : public edm::EDFilter {
   public:
      explicit AnomTrkMETFilter(const edm::ParameterSet&);
      ~AnomTrkMETFilter();

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
		struct Hist_t{
			TProfile *hAtanVsMet;
			TH1F *hPfMet;
			TH1F *hPfMetSig;
			TH1F *hPfMet_highPureTrks_atanGt7;
			TH1F *hPfMetSig_highPureTrks_atanGt7;
			TH1F *hPfMet_highPureTrks_atanLt7;
			TH1F *hPfMetSig_highPureTrks_atanLt7;
		};

		struct TrackHist_t{
			TH1F *hQuality;
			TH1F *hRatioIt4toIt01;
			TH1F *hRatioIt5toIt01;
			TH1F *hRatioPixlesstoIt01;
			TH1F *hRatioIt4toIt01pure;
			TH1F *hRatioIt5toIt01pure;
			TH1F *hRatioPixlesstoIt01pure;
		};
		TH1F* hPrimVtxz; //z cdt of the primary vertex
		TH1F* hNVtx; //primary vertices
		struct FinalHist_t {
			TH1F* hPFMet;
			TH1F* hTCMet;
			TH1F* hCaloMet;
			TH1F* hCaloPFMet;
			TH1F* hTCPFMet;
			TH1F* hCaloTCMet;
			TH1F* hCalo2PFMet;
			TH1F* hTC2PFMet;
			TH1F* hCalo2TCMet;
		};
		struct RunLumiEvt_t
		{
			unsigned run;
			unsigned lumi;
			unsigned evt;
		};
		struct TrkRatio_t
		{
			double run;
			double event;
			double ratio;
			double ratio2;
			double ratio3;
		};
		struct EvtInfo_t
		{
			double run;
			double lumi;
			double event;
			double pfmet;
			double tcmet;
			double calomet;
		};
		struct CutHists_t {
			TH1F* nTrksAssoWithVtx2ntrks;
			TH1F* ntrks_nopixhits2ntrks;
			TH1F* nTrksNotAssoWithVtx2ntrks;
			TH1F* AfterNtrkCut_ratio; 
			TH1F* AfterNtrkCut_ratio2; 
			TH1F* AfterNtrkCut_ratio3;
			TH1F* AfterNtrkRatioCut_ratio2; 
			TH1F* AfterNtrkRatioCut_ratio3; 
			TH1F* AfterNtrkRatioRatio2Cut_ratio3; 
		};
		struct TrkInfo_t
		{
			TH1F *ntrks;
			TH1F *algo;
			TH1F *validPixHits;
			TH1F *purity;
			TH1F *pt;
			TH1F *wgt;
			TH1F *ptErr;
			TH1F* dxyErr;
			TH1F* d0Err;
			TH1F* exptInnerHits;
			TH1F* exptOuterHits;
			TProfile *ptVsptErr;
			TH1F* validHitFraction;
		};
		struct PFHists_t {
			TH1F* chargeHadFraction;
		//	TH1F* ;
		};
		void AnalyseTracks(const edm::Handle<reco::TrackCollection> tracks, 
							const reco::VertexCollection::const_iterator vtx,
							std::vector<unsigned>& vMatchedTracks);
		void FillUnmatchTrackHists(edm::Handle<reco::TrackCollection> tracks,
					const std::vector<unsigned> vMatchedTracks);
		bool storeHLT(const edm::Event& e, const edm::EventSetup& iSetup);
		void PFlowAnalyse(edm::Event&, const edm::EventSetup&);
		void FillTrackInfoHists(TrkInfo_t trkHist, reco::TrackCollection::const_iterator it);
		edm::InputTag m_trackCollection;
		edm::InputTag m_vtxCollection;
		//const double m_vtxzcut;
		std::vector<std::string> triggerPathsToStore_;  // Vector to store list of HLT paths to store results of in ntuple
		edm::InputTag hlTriggerResults_;    // Input tag for TriggerResults
		bool req_trigger;
		int inMinVtx;
		double dminMet;
		unsigned int iProcessed; // number of processed events
		unsigned int iPassed; //number of events passed the filter
		int inMinNdofVtx; //minimum number of dof for the vtx
		double dMaxPrimVtxZ; //maximum seperation between a track and the primary vertex.

		FinalHist_t hFinal_all, hFinal_my, hFinal_Andrea;

		TrackHist_t hTrkIter0, hTrkIter1, hTrkIter2, hTrkIter3, hTrkIter4, hTrkIter5;
		Hist_t histsIter4, histsIter5, histsPixless;

		std::vector<RunLumiEvt_t> it4Rej, it5Rej, pxlessRej;
		std::vector<RunLumiEvt_t>::iterator rejIt;

		//jet collections
		edm::InputTag caloJetInputTag_;
		edm::Handle<reco::CaloJetCollection> caloJetcollection;

		edm::InputTag pfJetInputTag_;
		edm::Handle<reco::PFJetCollection> pfJetcollection;
		
		std::vector<TrkRatio_t> vTrkRatio;
		std::vector<std::pair<edm::RunNumber_t, edm::EventNumber_t> > vBadEvents;
		bool AnomEvent(const RunLumiEvt_t runlumevt);
		bool processBadOnly;

		std::vector<EvtInfo_t> vBadEvents_fromMyCuts, vBadEvents_fromAndreasCuts;
		int iVerbose; // control print levels

		CutHists_t hCutHist;
		bool print_hists;
		std::string GetEPSname(const std::string name) 
		{
			std::stringstream s;
			s << name << ".png";
			return s.str();
		}
		void PrintHist(TH1 *hist);

		TH1F* trkPurityFraction;
		TH1F* nTrksVtxWgtlessThan1; //0.1
		TH1F* nTrksVtxWgtlessThan2; //0.2
		TH1F* nTrksVtxWgtlessThan3; //0.3
		TH1F* nTrksVtxWgtlessThan5; //0.5
		TH1F* nTrksVtxWgtmoreThan5; //0.5

		TH1F* nTrksVtxWgtlessThan1_pureTrkFraction;
		TH1F* nTrksVtxWgtlessThan2_pureTrkFraction;
		TH1F* nTrksVtxWgtlessThan3_pureTrkFraction;
		TH1F* nTrksVtxWgtlessThan5_pureTrkFraction;
		TH1F* nTrksVtxWgtmoreThan5_pureTrkFraction;

		TH1F* nTrksVtxWgtlessThan1_algo;
		TH1F* nTrksVtxWgtlessThan2_algo;
		TH1F* nTrksVtxWgtlessThan3_algo;
		TH1F* nTrksVtxWgtlessThan5_algo;
		TH1F* nTrksVtxWgtmoreThan5_algo;

		TH1F* nTrksVtxWgtlessThan1_validPixHits;
		TH1F* nTrksVtxWgtlessThan2_validPixHits;
		TH1F* nTrksVtxWgtlessThan3_validPixHits;
		TH1F* nTrksVtxWgtlessThan5_validPixHits;
		TH1F* nTrksVtxWgtmoreThan5_validPixHits;

		TH1F* nTrksVtxWgtlessThan1_pt;
		TH1F* nTrksVtxWgtlessThan2_pt;
		TH1F* nTrksVtxWgtlessThan3_pt;
		TH1F* nTrksVtxWgtlessThan5_pt;
		TH1F* nTrksVtxWgtmoreThan5_pt;

		TrkInfo_t hTrksAssoWithVtx, hTrksNotAssoWithVtx;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
AnomTrkMETFilter::AnomTrkMETFilter(const edm::ParameterSet& iConfig):
	m_trackCollection(iConfig.getParameter<edm::InputTag>("trackCollection"))
{
   //now do what ever initialization is needed
	triggerPathsToStore_ = iConfig.getParameter<std::vector<std::string> >("TriggerPathsToStore");
	hlTriggerResults_    = iConfig.getParameter<edm::InputTag>("HltTriggerResults");	
	req_trigger = iConfig.getUntrackedParameter<bool>("req_trigger",true);
	inMinVtx = iConfig.getUntrackedParameter<int>("nMinVtx",1);
	dminMet = iConfig.getUntrackedParameter<double>("minMet",0.0);
	inMinNdofVtx = iConfig.getUntrackedParameter<int>("ndofVtx",0);
	dMaxPrimVtxZ = iConfig.getUntrackedParameter<double>("maxDelzTrkVtx",200.0);
	caloJetInputTag_ = iConfig.getParameter<edm::InputTag>("caloJetInputTag_");
	pfJetInputTag_ = iConfig.getParameter<edm::InputTag>("pfJetInputTag_");
	iVerbose = iConfig.getUntrackedParameter<int>("verbose",0);
	print_hists = iConfig.getUntrackedParameter<bool>("printHists", false);
	processBadOnly = iConfig.getUntrackedParameter<bool>("processBadOnly", false);
	iProcessed = 0;
	iPassed = 0;

	//generate hists
	edm::Service<TFileService> fs;

	const double met_max = 4000, met_bins = 400;

	histsIter4.hAtanVsMet = fs->make<TProfile> ("Atan2VsMet_iter4" ,"atan2(iter4/iter0+1) Vs PF#slash{E}_{T}", met_bins,0,met_max);
	histsIter4.hPfMet  = fs->make<TH1F> ("PfMet_iter4" ," Raw #slash{E}_{T} (PFMETCollection::pfMet)", met_bins,0,met_max);
	histsIter4.hPfMetSig  = fs->make<TH1F> ("PfMetSig_iter4" ," #slash{E}_{T}-Sig (PFMETCollection::metSig)", 200,0,20);

	//high purity tracks
	histsIter4.hPfMet_highPureTrks_atanGt7  = fs->make<TH1F> ("PfMet_highPureTrks_atanGt7_iter4" ,"High Purity Tracks atan2(iter4/iter0+1)>0.7:  PF#slash{E}_{T} (PFMETCollection::pfMet)", met_bins,0,met_max);
	histsIter4.hPfMetSig_highPureTrks_atanGt7  = fs->make<TH1F> ("PfMetSig_highPureTrks_atanGt7_iter4" ,"High Purity Tracks atan2(iter4/iter0+1)>0.7: PF#slash{E}_{T}-Sig (PFMETCollection::metSig)", 200,0,20);
	histsIter4.hPfMet_highPureTrks_atanLt7  = fs->make<TH1F> ("PfMet_highPureTrks_atanLt7_iter4" ,"High Purity Tracks atan2(iter4/iter0+1)<0.7:  PF#slash{E}_{T} (PFMETCollection::pfMet)", met_bins,0,met_max);
	histsIter4.hPfMetSig_highPureTrks_atanLt7  = fs->make<TH1F> ("PfMetSig_highPureTrks_atanLt7_iter4" ,"High Purity Tracks atan2(iter4/iter0+1)<0.7: PF#slash{E}_{T}-Sig (PFMETCollection::metSig)", 200,0,20);

	histsIter5.hAtanVsMet = fs->make<TProfile> ("Atan2VsMet_iter5" ,"atan2(iter5/iter0+1) Vs PF#slash{E}_{T}", 200,0,1000);

	//high purity tracks
	histsIter5.hPfMet_highPureTrks_atanGt7  = fs->make<TH1F> ("PfMet_highPureTrks_atanGt7_iter5" ,"High Purity Tracks atan>(iter5/iter0+1)0.7:  PF#slash{E}_{T} (PFMETCollection::pfMet)", met_bins,0,met_max);
	histsIter5.hPfMetSig_highPureTrks_atanGt7  = fs->make<TH1F> ("PfMetSig_highPureTrks_atanGt7_iter5" ,"High Purity Tracks atan2(iter5/iter0+1)>0.7: PF#slash{E}_{T}-Sig (PFMETCollection::metSig)", 200,0,20);

	histsIter5.hPfMet_highPureTrks_atanLt7  = fs->make<TH1F> ("PfMet_highPureTrks_atanLt7_iter5" ,"High Purity Tracks atan<(iter5/iter0+1)0.7:  PF#slash{E}_{T} (PFMETCollection::pfMet)", met_bins,0,met_max);
	histsIter5.hPfMetSig_highPureTrks_atanLt7  = fs->make<TH1F> ("PfMetSig_highPureTrks_atanLt7_iter5" ,"High Purity Tracks atan2(iter5/iter0+1)<0.7: PF#slash{E}_{T}-Sig (PFMETCollection::metSig)", 200,0,20);


	//pixelless stuff from pure tracks
	histsPixless.hAtanVsMet = fs->make<TProfile> ("Atan2VsMet_Pixless" ,"atan2(Pixelless/iter0+1) Vs PF#slash{E}_{T}", met_bins,0,met_max);

	histsPixless.hPfMet_highPureTrks_atanGt7  = fs->make<TH1F> ("PfMet_highPureTrks_atanGt7_Pixless" ,"High Purity Tracks atan>(Pixless/iter0+1)0.7:  PF#slash{E}_{T} (PFMETCollection::pfMet)", met_bins,0,met_max);
	histsPixless.hPfMetSig_highPureTrks_atanGt7  = fs->make<TH1F> ("PfMetSig_highPureTrks_atanGt7_Pixless" ,"High Purity Tracks atan2(Pixless/iter0+1)>0.7: PF#slash{E}_{T}-Sig (PFMETCollection::metSig)", 200,0,20);

	histsPixless.hPfMet_highPureTrks_atanLt7  = fs->make<TH1F> ("PfMet_highPureTrks_atanLt7_Pixless" ,"High Purity Tracks atan2(Pixless/iter0+1)<0.7:  PF#slash{E}_{T} (PFMETCollection::pfMet)", met_bins,0,met_max);
	histsPixless.hPfMetSig_highPureTrks_atanLt7  = fs->make<TH1F> ("PfMetSig_highPureTrks_atanLt7_Pixless" ,"High Purity Tracks atan2(Pixless/iter0+1)<0.7: PF#slash{E}_{T}-Sig (PFMETCollection::metSig)", 200,0,20);

	hTrkIter0.hQuality = fs->make<TH1F> ("TrkQuality","Trk Quality Flag",10,0,10);
	hTrkIter0.hRatioIt4toIt01 = fs->make<TH1F> ("hRatioIt4toIt01","General Tracks: atan2(iter4/Iter0+1);atan2(ratio);Events;",200,0,2);
	hTrkIter0.hRatioIt5toIt01 = fs->make<TH1F> ("hRatioIt5toIt01","General Tracks: atan2(iter5/Iter0+1);atan2(ratio);Events;",200,0,2);
	hTrkIter0.hRatioPixlesstoIt01 = fs->make<TH1F> ("hRatioPixlesstoIt01","General Tracks: atan2(Pixless/Iter0+1);atan2(ratio);Events;",200,0,2);
	hTrkIter0.hRatioIt4toIt01pure = fs->make<TH1F> ("hRatioIt4toIt01pure","High Purity Tracks: atan2(iter4/Iter0+1);atan2(ratio);Events;",200,0,2);
	hTrkIter0.hRatioIt5toIt01pure = fs->make<TH1F> ("hRatioIt5toIt01pure","High Purity Tracks: atan2(iter5/Iter0+1);atan2(ratio);Events;",200,0,2);
	hTrkIter0.hRatioPixlesstoIt01pure = fs->make<TH1F> ("hRatioPixlesstoIt01pure","High Purity Tracks: atan2(Pixless/Iter0+1);atan2(ratio);Events;",200,0,2);

	hPrimVtxz = fs->make<TH1F> ("PrimVtxZ","Primary Vertex z-position;z [cm];Events;",300,0,300);
	hNVtx = fs->make<TH1F> ("nVtx","Number of Primary Vertices;N Vertices;Events;",30,0,30);

	const double ratio_bins = 300, ratio_min = -2., ratio_max = 1.;
	const double ratio2_bins = 1000, ratio2_min = -10., ratio2_max = 10.;
	hFinal_all.hPFMet = fs->make<TH1F> ("all_PFMet" ,"All Events: PFMET;PFMET [GeV];Events;", met_bins, 0,met_max);
	hFinal_all.hTCMet = fs->make<TH1F> ("all_TCMet" ,"All Events: TCMET;TCMET [GeV];Events;", met_bins, 0,met_max);
	hFinal_all.hCaloMet = fs->make<TH1F> ("all_CaloMet" ,"All Events: CALOMET;CALOMET [GeV];Events;", met_bins, 0,met_max);
	hFinal_all.hCaloPFMet = fs->make<TH1F> ("all_CaloPFMet" ,"All Events: #frac{CALOMET - PFMET}{CALOMET+PFMET};Ratio;Events;", ratio_bins, ratio_min, ratio_max);
	hFinal_all.hTCPFMet = fs->make<TH1F> ("all_TCPFMet" ,"All Events: #frac{TCMET - PFMET}{TCMET+PFMET};Ratio;Events;", ratio_bins, ratio_min, ratio_max);
	hFinal_all.hCaloTCMet = fs->make<TH1F> ("all_CaloTCMet" ,"All Events: #frac{CALOMET - TCMET}{CALOMET+TCMET};Ratio;Events;", ratio_bins, ratio_min, ratio_max);
	hFinal_all.hCalo2PFMet = fs->make<TH1F> ("all_Calo2PFMet" ,"All Events: #frac{CALOMET - PFMET}{PFMET};Ratio;Events;", ratio2_bins, ratio2_min, ratio2_max);
	hFinal_all.hTC2PFMet = fs->make<TH1F> ("all_TC2PFMet" ,"All Events: #frac{TCMET - PFMET}{PFMET};Ratio;Events;", ratio2_bins, ratio2_min, ratio2_max);
	hFinal_all.hCalo2TCMet = fs->make<TH1F> ("all_Calo2TCMet" ,"All Events: #frac{CALOMET - TCMET}{TCMET};Ratio;Events;", ratio2_bins, ratio2_min, ratio2_max);

	hFinal_my.hPFMet = fs->make<TH1F> ("my_PFMet" ,"My Events: PFMET;PFMET [GeV];Events;", met_bins, 0,met_max);
	hFinal_my.hTCMet = fs->make<TH1F> ("my_TCMet" ,"My Events: TCMET;TCMET [GeV];Events;", met_bins, 0,met_max);
	hFinal_my.hCaloMet = fs->make<TH1F> ("my_CaloMet" ,"My Events: CALOMET;CALOMET [GeV];Events;", met_bins, 0,met_max);
	hFinal_my.hCaloPFMet = fs->make<TH1F> ("my_CaloPFMet" ,"My Events: #frac{CALOMET - PFMET}{CALOMET+PFMET};Ratio;Events;", ratio_bins, ratio_min, ratio_max);
	hFinal_my.hTCPFMet = fs->make<TH1F> ("my_TCPFMet" ,"My Events: #frac{TCMET - PFMET}{TCMET+PFMET};Ratio;Events;", ratio_bins, ratio_min, ratio_max);
	hFinal_my.hCaloTCMet = fs->make<TH1F> ("my_CaloTCMet" ,"My Events: #frac{CALOMET - TCMET}{CALOMET+TCMET};Ratio;Events;", ratio_bins, ratio_min, ratio_max);
	hFinal_my.hCalo2PFMet = fs->make<TH1F> ("my_Calo2PFMet" ,"My Events: #frac{CALOMET - PFMET}{PFMET};Ratio;Events;", ratio2_bins, ratio2_min, ratio2_max);
	hFinal_my.hTC2PFMet = fs->make<TH1F> ("my_TC2PFMet" ,"My Events: #frac{TCMET - PFMET}{PFMET};Ratio;Events;", ratio2_bins, ratio2_min, ratio2_max);
	hFinal_my.hCalo2TCMet = fs->make<TH1F> ("my_Calo2TCMet" ,"My Events: #frac{CALOMET - TCMET}{TCMET};Ratio;Events;", ratio2_bins, ratio2_min, ratio2_max);

	hFinal_Andrea.hPFMet = fs->make<TH1F> ("andrea_PFMet" ,"Andrea Events: PFMET;PFMET [GeV];Events;", met_bins, 0,met_max);
	hFinal_Andrea.hTCMet = fs->make<TH1F> ("andrea_TCMet" ,"Andrea Events: TCMET;TCMET [GeV];Events;", met_bins, 0,met_max);
	hFinal_Andrea.hCaloMet = fs->make<TH1F> ("andrea_CaloMet" ,"Andrea Events: CALOMET;CALOMET [GeV];Events;", met_bins, 0,met_max);
	hFinal_Andrea.hCaloPFMet = fs->make<TH1F> ("andrea_CaloPFMet" ,"Andrea Events: #frac{CALOMET - PFMET}{CALOMET+PFMET};Ratio;Events;", ratio_bins, ratio_min, ratio_max);
	hFinal_Andrea.hTCPFMet = fs->make<TH1F> ("andrea_TCPFMet" ,"Andrea Events: #frac{TCMET - PFMET}{TCMET+PFMET};Ratio;Events;", ratio_bins, ratio_min, ratio_max);
	hFinal_Andrea.hCaloTCMet = fs->make<TH1F> ("andrea_CaloTCMet" ,"Andrea Events: #frac{CALOMET - TCMET}{CALOMET+TCMET};Ratio;Events;", ratio_bins, ratio_min, ratio_max);
	hFinal_Andrea.hCalo2PFMet = fs->make<TH1F> ("andrea_Calo2PFMet" ,"Andrea Events: #frac{CALOMET - PFMET}{PFMET};Ratio;Events;", ratio2_bins, ratio2_min, ratio2_max);
	hFinal_Andrea.hTC2PFMet = fs->make<TH1F> ("andrea_TC2PFMet" ,"Andrea Events: #frac{TCMET - PFMET}{PFMET};Ratio;Events;", ratio2_bins, ratio2_min, ratio2_max);
	hFinal_Andrea.hCalo2TCMet = fs->make<TH1F> ("andrea_Calo2TCMet" ,"Andrea Events: #frac{CALOMET - TCMET}{TCMET};Ratio;Events;", ratio2_bins, ratio2_min, ratio2_max);


	const float ratio3_bins = 100, ratio3_min = 0, ratio3_max = 1;
	hCutHist.nTrksAssoWithVtx2ntrks = fs->make<TH1F> ("nTrksAssoWithVtx2ntrks" ,"nTrksAssoWithVtx2ntrks;Ratio;Events;", ratio3_bins, ratio3_min, ratio3_max);
	hCutHist.ntrks_nopixhits2ntrks = fs->make<TH1F> ("ntrks_nopixhits2ntrks" ,"ntrks_nopixhits2ntrks;Ratio;Events;", ratio3_bins, ratio3_min, ratio3_max);
	hCutHist.nTrksNotAssoWithVtx2ntrks = fs->make<TH1F> ("nTrksNotAssoWithVtx2ntrks" ,"nTrksNotAssoWithVtx2ntrks;Ratio;Events;", ratio3_bins, ratio3_min, ratio3_max);
	hCutHist.AfterNtrkCut_ratio = fs->make<TH1F> ("AfterNtrkCut_ratio" ,"AfterNtrkCut_ratio;Ratio;Events;", ratio3_bins, ratio3_min, ratio3_max);
	hCutHist.AfterNtrkCut_ratio2 = fs->make<TH1F> ("AfterNtrkCut_ratio2" ,"AfterNtrkCut_ratio2;Ratio;Events;", ratio3_bins, ratio3_min, ratio3_max);
	hCutHist.AfterNtrkCut_ratio3 = fs->make<TH1F> ("AfterNtrkCut_ratio3" ,"AfterNtrkCut_ratio3;Ratio;Events;", ratio3_bins, ratio3_min, ratio3_max);
	hCutHist.AfterNtrkRatioCut_ratio2 = fs->make<TH1F> ("AfterNtrkRatioCut_ratio2" ,"AfterNtrkRatioCut_ratio2;Ratio;Events;", ratio3_bins, ratio3_min, ratio3_max);
	hCutHist.AfterNtrkRatioCut_ratio3 = fs->make<TH1F> ("AfterNtrkRatioCut_ratio3" ,"AfterNtrkRatioCut_ratio3;Ratio;Events;", ratio3_bins, ratio3_min, ratio3_max);
	hCutHist.AfterNtrkRatioRatio2Cut_ratio3 = fs->make<TH1F> ("AfterNtrkRatioRatio2Cut_ratio3" ,"AfterNtrkRatioRatio2Cut_ratio3;Ratio;Events;", ratio3_bins, ratio3_min, ratio3_max);


	//track quality check
	hTrksAssoWithVtx.ntrks = fs->make<TH1F> ("nTrksAssoWithVtx_ntrks" ,"nTrksAssoWithVtx;N tracks;Events;", 700,0 ,700);
	hTrksAssoWithVtx.algo = fs->make<TH1F> ("nTrksAssoWithVtx_algo" ,"nTrksAssoWithVtx algo;Algo;;", 11, -.05, 10.5);
	hTrksAssoWithVtx.validPixHits = fs->make<TH1F> ("nTrksAssoWithVtx_pixHits" ,"nTrksAssoWithVtx validPixHist;PixHist;;", 20, 0, 20);
	hTrksAssoWithVtx.purity = fs->make<TH1F> ("nTrksAssoWithVtx_purity" ,"nTrksAssoWithVtx purity;Algo;Events;", 11, -.05, 10.5);
	hTrksAssoWithVtx.pt = fs->make<TH1F> ("nTrksAssoWithVtx_pt" ,"nTrksAssoWithVtx pt;pt;;", 100, 0, 500);
	hTrksAssoWithVtx.wgt = fs->make<TH1F> ("nTrksAssoWithVtx_wgt" ,"nTrksAssoWithVtx weight;Weight;;", 100, 0, 1);
	hTrksAssoWithVtx.ptErr = fs->make<TH1F> ("nTrksAssoWithVtx_ptErr" ,"nTrksAssoWithVtx ptErr;#Delta pt/pt;;", 400, 0, 20);
	hTrksAssoWithVtx.dxyErr = fs->make<TH1F> ("nTrksAssoWithVtx_dxyErr" ,"nTrksAssoWithVtx dxyErr;#Delta dxy;;", 300, -30, 30);
	hTrksAssoWithVtx.d0Err = fs->make<TH1F> ("nTrksAssoWithVtx_d0Err" ,"nTrksAssoWithVtx d0Err;#Delta d0;;", 300, -30, 30);
	hTrksAssoWithVtx.exptInnerHits = fs->make<TH1F> ("nTrksAssoWithVtx_exptInnerHits" ,"nTrksAssoWithVtx exptInnerHits;#Delta exptInnerHits;;", 20, 0, 20);
	hTrksAssoWithVtx.exptOuterHits = fs->make<TH1F> ("nTrksAssoWithVtx_exptOuterHits" ,"nTrksAssoWithVtx exptOuterHits;#Delta exptOuterHits;;", 20, 0, 20);
	hTrksAssoWithVtx.ptVsptErr = fs->make<TProfile> ("nTrksAssoWithVtx_ptVsptErr" ,"nTrksAssoWithVtx ptVsptErr;pt;#Delta pt/pt;", 100, 0, 500);
	hTrksAssoWithVtx.validHitFraction = fs->make<TH1F> ("nTrksAssoWithVtx_validHitFraction" ,"nTrksAssoWithVtx validHitFractio;Valid Hit Fraction;;", 100, 0, 1);

	hTrksNotAssoWithVtx.ntrks = fs->make<TH1F> ("nTrksNotAssoWithVtx_ntrks" ,"nTrksNotAssoWithVtx;N tracks;Events;", 700,0 ,700);
	hTrksNotAssoWithVtx.algo = fs->make<TH1F> ("nTrksNotAssoWithVtx_algo" ,"nTrksNotAssoWithVtx algo;Algo;;", 11, -.05, 10.5);
	hTrksNotAssoWithVtx.validPixHits = fs->make<TH1F> ("nTrksNotAssoWithVtx_pixHits" ,"nTrksNotAssoWithVtx validPixHist;PixHist;;", 20, 0, 20);
	hTrksNotAssoWithVtx.purity = fs->make<TH1F> ("nTrksNotAssoWithVtx_purity" ,"nTrksNotAssoWithVtx purity;Algo;Events;", 11, -.05, 10.5);
	hTrksNotAssoWithVtx.pt = fs->make<TH1F> ("nTrksNotAssoWithVtx_pt" ,"nTrksNotAssoWithVtx pt;pt;;", 100, 0, 500);
	hTrksNotAssoWithVtx.ptErr = fs->make<TH1F> ("nTrksNotAssoWithVtx_ptErr" ,"nTrksNotAssoWithVtx ptErr;#Delta pt/pt;;", 400, 0, 20);
	hTrksNotAssoWithVtx.dxyErr = fs->make<TH1F> ("nTrksNotAssoWithVtx_dxyErr" ,"nTrksNotAssoWithVtx dxyErr;#Delta dxy;;", 300, -30, 30);
	hTrksNotAssoWithVtx.d0Err = fs->make<TH1F> ("nTrksNotAssoWithVtx_d0Err" ,"nTrksNotAssoWithVtx d0Err;#Delta d0;;", 300, -30, 30);
	hTrksNotAssoWithVtx.exptInnerHits = fs->make<TH1F> ("nTrksNotAssoWithVtx_exptInnerHits" ,"nTrksNotAssoWithVtx exptInnerHits;#Delta exptInnerHits;;", 20, 0, 20);
	hTrksNotAssoWithVtx.exptOuterHits = fs->make<TH1F> ("nTrksNotAssoWithVtx_exptOuterHits" ,"nTrksNotAssoWithVtx exptOuterHits;#Delta exptOuterHits;;", 20, 0, 20);
	hTrksNotAssoWithVtx.ptVsptErr = fs->make<TProfile> ("nTrksNotAssoWithVtx_ptVsptErr" ,"nTrksNotAssoWithVtx ptVsptErr;pt;#Delta pt/pt;", 100, 0, 500);
	hTrksNotAssoWithVtx.validHitFraction = fs->make<TH1F> ("nTrksNotAssoWithVtx_validHitFraction" ,"nTrksNotAssoWithVtx validHitFractio;Valid Hit Fraction;;", 100, 0, 1);


	trkPurityFraction  = fs->make<TH1F> ("trkPurityFraction" ,"trkPurityFraction (high purity tracks/all tracks);ratio;;", 100, 0, 1);
	nTrksVtxWgtlessThan1  = fs->make<TH1F> ("nTrksWgt_lt1" ,"Tracks Associated with a vertex: NTrks with Vtx wgt <0.1 ;ratio;;", 500, 0, 500);
	nTrksVtxWgtlessThan2  = fs->make<TH1F> ("nTrksWgt_lt2" ,"Tracks Associated with a vertex: NTrks with Vtx wgt <0.2 ;ratio;;", 500, 0, 500);
	nTrksVtxWgtlessThan3  = fs->make<TH1F> ("nTrksWgt_lt3" ,"Tracks Associated with a vertex: NTrks with Vtx wgt <0.3 ;ratio;;", 500, 0, 500);
	nTrksVtxWgtlessThan5  = fs->make<TH1F> ("nTrksWgt_lt5" ,"Tracks Associated with a vertex: NTrks with Vtx wgt <0.5 ;ratio;;", 500, 0, 500);
	nTrksVtxWgtmoreThan5  = fs->make<TH1F> ("nTrksWgt_gt5" ,"Tracks Associated with a vertex: NTrks with Vtx wgt >=0.5 ;ratio;;", 500, 0, 500);

	nTrksVtxWgtlessThan1_pureTrkFraction  = fs->make<TH1F> ("nTrksWgt_lt1_pureTrkFrac" ,"Tracks Associated with a vertex: NTrks with Vtx wgt <0.1: High purity track fraction;ratio;;", 100, 0, 1);
	nTrksVtxWgtlessThan2_pureTrkFraction  = fs->make<TH1F> ("nTrksWgt_lt2_pureTrkFrac" ,"Tracks Associated with a vertex: NTrks with Vtx wgt <0.2: High purity track fraction;ratio;;", 100, 0, 1);
	nTrksVtxWgtlessThan3_pureTrkFraction  = fs->make<TH1F> ("nTrksWgt_lt3_pureTrkFrac" ,"Tracks Associated with a vertex: NTrks with Vtx wgt <0.3: High purity track fraction;ratio;;", 100, 0, 1);
	nTrksVtxWgtlessThan5_pureTrkFraction  = fs->make<TH1F> ("nTrksWgt_lt5_pureTrkFrac" ,"Tracks Associated with a vertex: NTrks with Vtx wgt <0.5: High purity track fraction;ratio;;", 100, 0, 1);
	nTrksVtxWgtmoreThan5_pureTrkFraction  = fs->make<TH1F> ("nTrksWgt_gt5_pureTrkFrac" ,"Tracks Associated with a vertex: NTrks with Vtx wgt>=0.5: High purity track fraction;ratio;;", 100, 0, 1);

	nTrksVtxWgtlessThan1_algo  = fs->make<TH1F> ("nTrksWgt_lt1_algo" ,"Tracks Associated with a vertex: NTrks with Vtx wgt <0.1: algo;ratio;;", 10, 0, 10);
	nTrksVtxWgtlessThan2_algo  = fs->make<TH1F> ("nTrksWgt_lt2_algo" ,"Tracks Associated with a vertex: NTrks with Vtx wgt <0.2: algo;ratio;;", 10, 0, 10);
	nTrksVtxWgtlessThan3_algo  = fs->make<TH1F> ("nTrksWgt_lt3_algo" ,"Tracks Associated with a vertex: NTrks with Vtx wgt <0.3: algo;ratio;;", 10, 0, 10);
	nTrksVtxWgtlessThan5_algo  = fs->make<TH1F> ("nTrksWgt_lt5_algo" ,"Tracks Associated with a vertex: NTrks with Vtx wgt <0.5: algo;ratio;;", 10, 0, 10);
	nTrksVtxWgtmoreThan5_algo  = fs->make<TH1F> ("nTrksWgt_gt5_algo" ,"Tracks Associated with a vertex: NTrks with Vtx wgt>=0.5: algo;ratio;;", 10, 0, 10);

	nTrksVtxWgtlessThan1_validPixHits  = fs->make<TH1F> ("nTrksWgt_lt1_validPixHits" ,"Tracks Associated with a vertex: NTrks with Vtx wgt <0.1: validPixHits;ratio;;", 20, 0, 20);
	nTrksVtxWgtlessThan2_validPixHits  = fs->make<TH1F> ("nTrksWgt_lt2_validPixHits" ,"Tracks Associated with a vertex: NTrks with Vtx wgt <0.2: validPixHits;ratio;;", 20, 0, 20);
	nTrksVtxWgtlessThan3_validPixHits  = fs->make<TH1F> ("nTrksWgt_lt3_validPixHits" ,"Tracks Associated with a vertex: NTrks with Vtx wgt <0.3: validPixHits;ratio;;", 20, 0, 20);
	nTrksVtxWgtlessThan5_validPixHits  = fs->make<TH1F> ("nTrksWgt_lt5_validPixHits" ,"Tracks Associated with a vertex: NTrks with Vtx wgt <0.5: validPixHits;ratio;;", 20, 0, 20);
	nTrksVtxWgtmoreThan5_validPixHits  = fs->make<TH1F> ("nTrksWgt_gt5_validPixHits" ,"Tracks Associated with a vertex: NTrks with Vtx wgt>=0.5: validPixHits;ratio;;", 20, 0, 20);

	nTrksVtxWgtlessThan1_pt  = fs->make<TH1F> ("Tracks Associated with a vertex: nTrksWgt_lt1_pt" ,"NTrks with Vtx wgt <0.1: pt;ratio;;", 500, 0, 500);
	nTrksVtxWgtlessThan2_pt  = fs->make<TH1F> ("Tracks Associated with a vertex: nTrksWgt_lt2_pt" ,"NTrks with Vtx wgt <0.2: pt;ratio;;", 500, 0, 500);
	nTrksVtxWgtlessThan3_pt  = fs->make<TH1F> ("Tracks Associated with a vertex: nTrksWgt_lt3_pt" ,"NTrks with Vtx wgt <0.3: pt;ratio;;", 500, 0, 500);
	nTrksVtxWgtlessThan5_pt  = fs->make<TH1F> ("Tracks Associated with a vertex: nTrksWgt_lt5_pt" ,"NTrks with Vtx wgt <0.5: pt;ratio;;", 500, 0, 500);
	nTrksVtxWgtmoreThan5_pt  = fs->make<TH1F> ("Tracks Associated with a vertex: nTrksWgt_gt5_pt" ,"NTrks with Vtx wgt>=0.5: pt;ratio;;", 500, 0, 500);
}


AnomTrkMETFilter::~AnomTrkMETFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
AnomTrkMETFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
	++iProcessed;
	const double fMinTrkPt = 0.0;
	bool bSaveEvent = false;

	RunNumber_t kRun   = iEvent.id().run();
	EventNumber_t kEvent = iEvent.id().event();
	LuminosityBlockNumber_t kLumi  = iEvent.id().luminosityBlock(); 
	RunLumiEvt_t kRunLumiEvent;
	kRunLumiEvent.run = kRun;
	kRunLumiEvent.lumi = kLumi;
	kRunLumiEvent.evt = kEvent;
	if (processBadOnly)
	{
		std::cout << "Processing only bad events." << std::endl;
		if (! AnomEvent(kRunLumiEvent)) return 0;
	}


	if (iVerbose) std::cout << "====== processing event " << kRun << ", " << kEvent << std::endl;

	//========== Trigger Selections
	if (req_trigger)
	{
		if (! storeHLT(iEvent, iSetup)) return 0;
	}

	//************** PF Collections analyser
	//PFlowAnalyse(iEvent, iSetup);

 	// ==========================================================
	// MET Information 
	Handle<reco::PFMETCollection> pfMet;
	iEvent.getByLabel("pfMet",pfMet);
	double dPFMet = 0.0, dPFMetSig = 0.;
	if (pfMet.isValid())
	{
		for (reco::PFMETCollection::const_iterator it = pfMet->begin(); it != pfMet->end(); ++it)
		{
			dPFMet = it->pt();
			dPFMetSig = it->mEtSig();
			if (dPFMet<dminMet) return 0;  //minimum met requirement
		}
	} else
	{
		std::cout << "Valid PFMETCollection Not Found!" << std::endl;
		assert(false);
	}
	
	Handle<reco::CaloMETCollection> caloMet;
	iEvent.getByLabel("met",caloMet);
	double dCaloMet = 0.0, dCaloMetSig = 0.;
	if (caloMet.isValid())
	{
		for (reco::CaloMETCollection::const_iterator it = caloMet->begin(); it != caloMet->end(); ++it)
		{
			dCaloMet = it->pt();
			dCaloMetSig = it->mEtSig();
			if (dCaloMet<dminMet) return 0;  //minimum met requirement
		}
	} else
	{
		std::cout << "Valid CaloMETCollection Not Found!" << std::endl;
		assert(false);
	}

	Handle<reco::METCollection> tcMet;
	iEvent.getByLabel("tcMet",tcMet);
	double dTcMet = 0.0, dTcMetSig = 0.;
	if (tcMet.isValid())
	{
		for (reco::METCollection::const_iterator it = tcMet->begin(); it != tcMet->end(); ++it)
		{
			dTcMet = it->pt();
			dTcMetSig = it->mEtSig();
			if (dTcMet<dminMet) return 0;  //minimum met requirement
		}
	} else
	{
		std::cout << "Valid TCMETCollection Not Found!" << std::endl;
		assert(false);
	}

	//NEED TO DO THIS AFTER REMOVING EVENTS WITH MUONS AS CALO MET DOES NOT ACCOUNT FOR MUONS?
	const double absdiff_CaloPFMet = (dCaloMet - dPFMet)/(dCaloMet + dPFMet);
	const double absdiff_TCPFMet = (dTcMet - dPFMet)/(dTcMet + dTcMet);
	const double absdiff_CaloTCMet = (dCaloMet - dTcMet)/(dCaloMet + dTcMet);
	const double rat_CaloMet2PFMet = (dCaloMet-dPFMet)/dPFMet;
	const double rat_TCMet2PFMet = (dTcMet-dPFMet)/dPFMet;
	const double rat_CaloMet2TCMet = (dCaloMet - dTcMet)/dTcMet;

	EvtInfo_t eventInfo;
	eventInfo.run = kRun;
	eventInfo.event = kEvent;
	eventInfo.lumi = kLumi;
	eventInfo.pfmet = dPFMet;
	eventInfo.tcmet = dTcMet;
	eventInfo.calomet = dCaloMet;

	/**********************************************************/
	/* DO MY SELECTION METHOD *********************************/
	/**********************************************************/
	
	Handle<reco::TrackCollection> tracks;
	iEvent.getByLabel("generalTracks", tracks);
	double ntrks = 0;
	double ntrks_above9 = 0;  //number of tracks with pt<0.9 GeV to remove event with bunch of very low pt tracks
	double ntrks_nopixhits = 0;
	int jj = 0 ;
	for (reco::TrackCollection::const_iterator it = tracks->begin(); it != tracks->end(); ++it)
	{
			const double trkpt = it->pt();
			if (trkpt< fMinTrkPt) continue;
			jj++;
			//std::cout << "trk [" << jj << "]"<< std::setw(10) << it->pt() << "/"<< std::setw(10) << it->phi() << "/"<< std::setw(10)<< it->eta() << std::endl;
			ntrks++;
			if (trkpt>0.9) ++ntrks_above9;
			//if (it->hitPattern().numberOfValidPixelHit() ==0) ++ntrks_nopixhits; 
			if (it->hitPattern().numberOfValidPixelHits() <3) ++ntrks_nopixhits; 
	}

 // ==========================================================
 //Vertex information

	Handle<reco::VertexCollection> vertexHandle;
  	iEvent.getByLabel("offlinePrimaryVertices", vertexHandle);
	double dNvtx = 0;
	double dPrimVtx_z = 0;
	double nTrksAssoWithVtx = 0;
	std::vector<unsigned> vMatchedTracks;

   if (vertexHandle.isValid())
	{
     reco::VertexCollection vertexCollection = *(vertexHandle.product());
     if (iVerbose) std::cout << __LINE__<< "::" << __FUNCTION__ << "::VtxSize = "<<  vertexCollection.size() << std::endl;
     reco::VertexCollection::const_iterator v = vertexCollection.begin();
	  if (v->isValid())
	  {
		  int i=0;
		  for (; v != vertexCollection.end(); ++v)
		  {
			  i++;
			  if (v->isFake()) return 0;
			  dPrimVtx_z = v->z();
			  if (fabs(dPrimVtx_z) > dMaxPrimVtxZ) return 0;  // do not use the event. 
			  if (v->ndof() < inMinNdofVtx) return 0; 
			  
			  ++dNvtx;

			  hPrimVtxz->Fill(dPrimVtx_z);
				//std::cout << "vertex["<< i << "] ndof/chi2/z = " <<v->ndof() << "/ " << v->chi2() << v->z() 
			  	//		<< "  ][ntrks associated = " << v->tracksSize() << std::endl; 
				AnalyseTracks(tracks, v, vMatchedTracks);
				//loop over tracks associated with the vertex
				reco::Vertex::trackRef_iterator trackIter = v->tracks_begin();
				double nTrk = 0;
				//hTrksPerVtx->Fill(v->tracksSize());

				if (iVerbose) std::cout << "Vtx[" << i << "] nTrks[" << v->tracksSize() << "," << v->nTracks(-.01) <<"]<0.1, 0.3, 0.5,>0.5["<< v->nTracks(0.1) << ", "<<
								v->nTracks(0.3) << ", "<< v->nTracks(0.5) << std::endl;
				double pureTrkFrac[] = {0,0,0,0,0};
				for (;trackIter != v->tracks_end(); trackIter++)
				{
					const double trkpt = (*trackIter)->pt();
					if (trkpt<fMinTrkPt) continue;
					nTrk++;
					const double trk_z = (*trackIter)->innerPosition().Z();
					//hTrkVtxSeparation->Fill(abs(dPrimVtx_z - trk_z));
					if ( ( ( (float) (v->tracksSize() - v->nTracks(0.1)))/(float)v->tracksSize() ) > 0.1)
					{
						nTrksVtxWgtlessThan1->Fill(v->nTracks(0.1));
						nTrksVtxWgtlessThan1_algo->Fill((*trackIter)->algo());
						nTrksVtxWgtlessThan1_validPixHits->Fill((*trackIter)->hitPattern().numberOfValidPixelHits());
						nTrksVtxWgtlessThan1_pt->Fill((*trackIter)->pt());
						if ((*trackIter)->quality(reco::TrackBase::highPurity)) ++pureTrkFrac[0];
					}
					if ( ( ( (float)(v->tracksSize() - v->nTracks(0.2)))/v->tracksSize() ) > 0.2)
					{
						nTrksVtxWgtlessThan2->Fill(v->nTracks(0.2));
						nTrksVtxWgtlessThan2_algo->Fill((*trackIter)->algo());
						nTrksVtxWgtlessThan2_validPixHits->Fill((*trackIter)->hitPattern().numberOfValidPixelHits());
						nTrksVtxWgtlessThan2_pt->Fill((*trackIter)->pt());
						if ((*trackIter)->quality(reco::TrackBase::highPurity)) ++pureTrkFrac[1];
					}
					if ( ( ( (float)(v->tracksSize() - v->nTracks(0.3)))/v->tracksSize() ) > 0.3)
					{
						nTrksVtxWgtlessThan3->Fill(v->nTracks(0.3));
						nTrksVtxWgtlessThan3_algo->Fill((*trackIter)->algo());
						nTrksVtxWgtlessThan3_validPixHits->Fill((*trackIter)->hitPattern().numberOfValidPixelHits());
						nTrksVtxWgtlessThan3_pt->Fill((*trackIter)->pt());
						if ((*trackIter)->quality(reco::TrackBase::highPurity)) ++pureTrkFrac[2];
					}
					if ( ( ( (float)(v->tracksSize() - v->nTracks(0.5)))/v->tracksSize() ) > 0.5)
					{
						nTrksVtxWgtlessThan5->Fill(v->nTracks(0.5));
						nTrksVtxWgtlessThan5_algo->Fill((*trackIter)->algo());
						nTrksVtxWgtlessThan5_validPixHits->Fill((*trackIter)->hitPattern().numberOfValidPixelHits());
						nTrksVtxWgtlessThan5_pt->Fill((*trackIter)->pt());
						if ((*trackIter)->quality(reco::TrackBase::highPurity)) ++pureTrkFrac[3];
					} 
					if ( ( ( (float)(v->tracksSize() - v->nTracks(0.5)))/v->tracksSize() ) < 0.5)
					{
						nTrksVtxWgtmoreThan5->Fill(v->nTracks(0.5));
						nTrksVtxWgtmoreThan5_algo->Fill((*trackIter)->algo());
						nTrksVtxWgtmoreThan5_validPixHits->Fill((*trackIter)->hitPattern().numberOfValidPixelHits());
						nTrksVtxWgtmoreThan5_pt->Fill((*trackIter)->pt());
						if ((*trackIter)->quality(reco::TrackBase::highPurity)) ++pureTrkFrac[4];
					}

				}
				nTrksAssoWithVtx += nTrk;

				if ( ( (float)(v->tracksSize() - v->nTracks(0.1))/v->tracksSize() ) > 0.1)
				{
					nTrksVtxWgtlessThan1_pureTrkFraction->Fill(pureTrkFrac[0]/v->tracksSize());
				}
				if ( ( (float)(v->tracksSize() - v->nTracks(0.2))/v->tracksSize() ) > 0.2)
				{
					nTrksVtxWgtlessThan2_pureTrkFraction->Fill(pureTrkFrac[1]/v->tracksSize());
				}
				if ( ( (float)(v->tracksSize() - v->nTracks(0.3))/v->tracksSize() ) > 0.3)
				{
					nTrksVtxWgtlessThan3_pureTrkFraction->Fill(pureTrkFrac[2]/v->tracksSize());
				}
				if ( ( (float)(v->tracksSize() - v->nTracks(0.5))/v->tracksSize() ) > 0.5)
				{
					nTrksVtxWgtlessThan5_pureTrkFraction->Fill(pureTrkFrac[3]/v->tracksSize());
				}
				if ( ( (float)(v->tracksSize() - v->nTracks(0.5))/v->tracksSize() ) < 0.5)
				{
					nTrksVtxWgtmoreThan5_pureTrkFraction->Fill(pureTrkFrac[4]/v->tracksSize());
				}

				//if (ntrks>0) hTrksPerVtxNtrksRatio->Fill(nTrk/ntrks);

//				std::cout << "for vtx[" << i << "] trks for bins [0-5,5-10,15-20,20,25] =" << 
//							trkPer5GevBin[0] << ", " << trkPer5GevBin[1] << ", " << trkPer5GevBin[2]
//							<< ", " << trkPer5GevBin[3] << ", " << trkPer5GevBin[4] << ", " << trkPer5GevBin[5] << std::endl;
//				std::cout << "for vtx[" << i << "] ratio to ntrks: "<< nTrk << "/" << ntrks <<" = " << nTrk/ntrks << std::endl;

			  //break;
		  }

		 hNVtx->Fill(dNvtx);
		 //std::cout << __LINE__ << "::ntrks = " << vMatchedTracks.size() << std::endl;
		 hTrksAssoWithVtx.ntrks->Fill(vMatchedTracks.size());
	  }

		//fill all events hist after the general good vertex selection
	  hFinal_all.hPFMet->Fill(dPFMet);
	  hFinal_all.hTCMet->Fill(dTcMet);
	  hFinal_all.hCaloMet->Fill(dCaloMet);
	  hFinal_all.hCaloPFMet->Fill(absdiff_CaloPFMet);
	  hFinal_all.hTCPFMet->Fill(absdiff_TCPFMet);
	  hFinal_all.hCaloTCMet->Fill(absdiff_CaloTCMet);
	  hFinal_all.hCalo2PFMet->Fill(rat_CaloMet2PFMet);
	  hFinal_all.hTC2PFMet->Fill(rat_TCMet2PFMet);
	  hFinal_all.hCalo2TCMet->Fill(rat_CaloMet2TCMet);


	  //std::cout << " {{ All trks assoc with vtx to all trks ratio =  nTrksAssoWithVtx / ntrks =" << nTrksAssoWithVtx / ntrks << "}}" << std::endl; 
	  TrkRatio_t tr;
	  tr.run = kRun;
	  tr.event = kEvent;
	  tr.ratio  = ( (ntrks > 0) ? nTrksAssoWithVtx / ntrks : 0);
	  tr.ratio2 = ( (ntrks > 0) ? ntrks_nopixhits / ntrks : 0);
	  tr.ratio3 = ( (ntrks - nTrksAssoWithVtx >0) ? ntrks_nopixhits/(ntrks - nTrksAssoWithVtx) : 0);
	  //vTrkRatio.push_back(tr);
	  //if  (ntrks>10 && tr.ratio<0.50 && tr.ratio2>0.4 && tr.ratio3>0.3) bSaveEvent = true;
	 // if  (ntrks_above9>3 && tr.ratio<0.30 && tr.ratio2>0.3 && tr.ratio3>0.3) bSaveEvent = true;

		hCutHist.nTrksAssoWithVtx2ntrks->Fill(tr.ratio);
		hCutHist.ntrks_nopixhits2ntrks->Fill(tr.ratio2);
		hCutHist.nTrksNotAssoWithVtx2ntrks->Fill(tr.ratio3);

		if (ntrks>10)
		{
			hCutHist.AfterNtrkCut_ratio->Fill(tr.ratio);
			hCutHist.AfterNtrkCut_ratio2->Fill(tr.ratio2);
			hCutHist.AfterNtrkCut_ratio3->Fill(tr.ratio3);
			if (tr.ratio<0.55)
			{
				hCutHist.AfterNtrkRatioCut_ratio2->Fill(tr.ratio2);
				hCutHist.AfterNtrkRatioCut_ratio3->Fill(tr.ratio3);

				if (tr.ratio2>0.29)
				{
					hCutHist.AfterNtrkRatioRatio2Cut_ratio3->Fill(tr.ratio3);
				}
			}
		}
		//Tai's event with low threshold
		//19822 ====== FOUND EVENT event 165415:103:91591777 , /Jet/Run2011A-PromptReco-v4/RECO, /store/data/Run2011A/Jet/RECO/PromptReco-v4/000/165/415/A44BE14A-BD85-E011-964F-001D09F24024.root
		//19823 ntrks/ratio/ratio2/ratio3 = 79/0.109649/0.298246/0.334975
		//163270:409:252238770 (/METBTag/Run2011A-PromptReco-v2/RECO, /store/data/Run2011A/METBTag/RECO/PromptReco-v2/000/163/270/52B6A9A1-A66E-E011-A45D-003048F024DE.root
		//167151:20:22174344 (found in "/MET/Run2011A-PromptReco-v4/RECO" /store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/151/94223088-409B-E011-8FB5-003048F11C58.root
		//
//	 	std::cout << "ntrks/ratio/ratio2/ratio3 = "<< ntrks_above9 << "/" <<
//					tr.ratio << "/"<< tr.ratio2 <<"/"<< tr.ratio3 << std::endl;
	  if  (ntrks>10 &&  tr.ratio<0.55 && tr.ratio2>0.35 && tr.ratio3>0.4) 
	  //if  (ntrks>10 &&  tr.ratio<0.55 && tr.ratio2>0.29 && tr.ratio3>0.3) 
	  {
		  if (dPFMet>dminMet || dTcMet>dminMet) bSaveEvent = true;
		  vBadEvents_fromMyCuts.push_back(eventInfo);
		  hFinal_my.hPFMet->Fill(dPFMet);
		  hFinal_my.hTCMet->Fill(dTcMet);
		  hFinal_my.hCaloMet->Fill(dCaloMet);
		  hFinal_my.hCaloPFMet->Fill(absdiff_CaloPFMet);
		  hFinal_my.hTCPFMet->Fill(absdiff_TCPFMet);
		  hFinal_my.hCaloTCMet->Fill(absdiff_CaloTCMet);
		  hFinal_my.hCalo2PFMet->Fill(rat_CaloMet2PFMet);
		  hFinal_my.hTC2PFMet->Fill(rat_TCMet2PFMet);
		  hFinal_my.hCalo2TCMet->Fill(rat_CaloMet2TCMet);
	  }
	} else 
	{
		std::cout << "AnomTrkAna: Could not find vertex collection" << std::endl;
		assert(false);	
   }
   
	//Fill info on tracks that has no association with vertices.
	FillUnmatchTrackHists(tracks, vMatchedTracks);
	//std::cout << __LINE__ << "::ntrks = " << vMatchedTracks.size() << std::endl;
	hTrksNotAssoWithVtx.ntrks->Fill(tracks->size()-vMatchedTracks.size());


	/**********************************************************/
	/* DO ANDREAS METHOD **************************************/
	/**********************************************************/

	if (tracks.isValid())
	{
		//histsIter4.hNtrks->Fill(tracks->size());
		double iNiter01trks = 0, iNiter4trks = 0, iNiter5trks = 0, iNpxlLess = 0;
		double iNiter01trkspure = 0, iNiter4trkspure = 0, iNiter5trkspure = 0, iNpxlLesspure = 0;
				
		for (reco::TrackCollection::const_iterator it = tracks->begin(); it != tracks->end(); ++it)
		{
			const double trk_z = it->innerPosition().Z();
			const double trkpt = it->pt();
			//	++i;
			if (it->quality(reco::TrackBase::highPurity))
			{
				hTrkIter0.hQuality->Fill(reco::TrackBase::highPurity);
				if (it->algo() == 4 || it->algo() == 5) ++iNiter01trkspure;
				if (it->algo() == 8) ++iNiter4trkspure;
				if (it->algo() == 9) ++iNiter5trkspure;
				if (it->hitPattern().numberOfValidPixelHits() == 0) ++iNpxlLesspure;
			} else 
			{
				if (it->algo() == 4 || it->algo() == 5) ++iNiter01trks;
				if (it->algo() == 8) ++iNiter4trks;
				if (it->algo() == 9) ++iNiter5trks;
				if (it->hitPattern().numberOfValidPixelHits() == 0) ++iNpxlLess;
			}
		}

		//atan2 is defined except for x==y==0 
		const double dRatIter4 = atan2(iNiter4trks, iNiter01trks);
		const double dRatIter5 = atan2(iNiter5trks, iNiter01trks);
		const double dRatPixless = atan2(iNpxlLess, iNiter01trks);
		const double dRatIter4pure = atan2(iNiter4trkspure, iNiter01trkspure);
		const double dRatIter5pure = atan2(iNiter5trkspure, iNiter01trkspure);
		const double dRatPixlesspure = atan2(iNpxlLesspure, iNiter01trkspure);

		//std::cout << "Ratio = Nit4/Nit0+1 = " << dRatIter4 << std::endl; 
		//std::cout << "Ratio = Nit5/Nit0+1 = " << dRatIter5 << std::endl; 
		hTrkIter0.hRatioIt4toIt01->Fill(dRatIter4);
		hTrkIter0.hRatioIt5toIt01->Fill(dRatIter5);
		hTrkIter0.hRatioPixlesstoIt01->Fill(dRatPixless);

		hTrkIter0.hRatioIt4toIt01pure->Fill(dRatIter4pure);
		hTrkIter0.hRatioIt5toIt01pure->Fill(dRatIter5pure);
		hTrkIter0.hRatioPixlesstoIt01pure->Fill(dRatPixlesspure);

		//make MET plots for different atan regions
		
		histsIter4.hPfMet->Fill(dPFMet);
		histsIter4.hPfMetSig->Fill(dPFMetSig);
		histsIter4.hAtanVsMet->Fill(dPFMet, dRatIter4);
		//high purity tracks: atan(Iter4/iter0+1)>0.7
		if (dRatIter4pure > 0.7)
		{
			histsIter4.hPfMet_highPureTrks_atanGt7->Fill(dPFMet);
			//std::cout << "dRatIter4pure == " << dRatIter4pure << "(>0.7) event:: " << kRun << ", " << kEvent << std::endl;
			it4Rej.push_back(kRunLumiEvent);
		  if (dPFMet>dminMet || dTcMet>dminMet) bSaveEvent = true;
		} else 
		{
			//histsIter4.hNvtx_highPureTrks_atanLt7->Fill(dNvtx);
			histsIter4.hPfMet_highPureTrks_atanLt7->Fill(dPFMet);
			histsIter4.hPfMetSig_highPureTrks_atanLt7->Fill(dPFMetSig);
		}


		histsIter5.hAtanVsMet->Fill(dPFMet, dRatIter5);
		//high purity tracks: atan(Iter5/iter0+1)>0.7
		if (dRatIter5pure > 0.7)
		{
			histsIter5.hPfMet_highPureTrks_atanGt7->Fill(dPFMet);
			//std::cout << "dRatIter5pure == " << dRatIter5pure << "(>0.7) event:: " << kRun << ", " << kEvent << std::endl;
			it5Rej.push_back(kRunLumiEvent);
		  if (dPFMet>dminMet || dTcMet>dminMet) bSaveEvent = true;
		} else 
		{
			//histsIter5.hNvtx_highPureTrks_atanLt7->Fill(dNvtx);
			histsIter5.hPfMet_highPureTrks_atanLt7->Fill(dPFMet);
			histsIter5.hPfMetSig_highPureTrks_atanLt7->Fill(dPFMetSig);
		}
	
		//pixless stuff
		histsPixless.hAtanVsMet->Fill(dPFMet, dRatPixlesspure);
		//high purity tracks: atan(PIXLESS/iter0+1)>0.7
		if (dRatPixlesspure > 0.7)
		{
			//std::cout << "dRatPixless == " << dRatPixlesspure << "(>0.7) event:: " << kRun << ", " << kEvent << std::endl;
			histsPixless.hPfMet_highPureTrks_atanGt7->Fill(dPFMet);
			pxlessRej.push_back(kRunLumiEvent);
		  if (dPFMet>dminMet || dTcMet>dminMet) bSaveEvent = true;
		} else 
		{
			//histsPixless.hNvtx_highPureTrks_atanLt7->Fill(dNvtx);
			histsPixless.hPfMet_highPureTrks_atanLt7->Fill(dPFMet);
			histsPixless.hPfMetSig_highPureTrks_atanLt7->Fill(dPFMetSig);
		}

		assert (histsIter4.hPfMet->GetEntries() == 
				(histsIter4.hPfMet_highPureTrks_atanLt7->GetEntries()
				 + histsIter4.hPfMet_highPureTrks_atanGt7->GetEntries())
					&& "PFMET entries do not match PFMET_HP_atanLT7!");
		
		//FILL HISTS with Andreas bad events
		if (dRatIter4pure > 0.7 || dRatIter5pure > 0.7 || dRatPixlesspure > 0.7)
		{
		  	vBadEvents_fromAndreasCuts.push_back(eventInfo);
			hFinal_Andrea.hPFMet->Fill(dPFMet);
			hFinal_Andrea.hTCMet->Fill(dTcMet);
			hFinal_Andrea.hCaloMet->Fill(dCaloMet);
			hFinal_Andrea.hCaloPFMet->Fill(absdiff_CaloPFMet);
			hFinal_Andrea.hTCPFMet->Fill(absdiff_TCPFMet);
			hFinal_Andrea.hCaloTCMet->Fill(absdiff_CaloTCMet);
			hFinal_Andrea.hCalo2PFMet->Fill(rat_CaloMet2PFMet);
			hFinal_Andrea.hTC2PFMet->Fill(rat_TCMet2PFMet);
			hFinal_Andrea.hCalo2TCMet->Fill(rat_CaloMet2TCMet);
		}


	
	} else
	{
		std::cout << "Valid Track Collection Not Found!" << std::endl;
		assert(false);
	}
	
	//std::cout << "std::cout>> number of tracks = " << tracks->size() << std::endl;

	//bSaveEvent = true;
	if (bSaveEvent) 
	{
		++iPassed;
		std::cout << "$$$$$$$$$$$$$$$$$$$$$ Saved event: " << kRun << ", " << kEvent << std::endl;
	}
	return bSaveEvent;
}
bool AnomTrkMETFilter::storeHLT(const edm::Event& e, const edm::EventSetup& iSetup){
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

/*		for( int itrig=0; itrig < ntrigs; itrig++){
			bool accept = TrgResultsHandle->accept(itrig);
			std::string trigName=TrgNames.triggerName(itrig);      
			if (accept){
				//std::cout << "%HLTInfo --  HLTTrigger(" << itrig << "): " << trigName << " = " << accept << std::endl;
			}
		}
*/
		for ( std::vector<std::string>::const_iterator trigNameTempl = triggerPathsToStore_.begin();  trigNameTempl != triggerPathsToStore_.end(); ++trigNameTempl) {
			int index = TrgNames.triggerIndex(*trigNameTempl);
			if (index !=  ntrigs )
				accept = TrgResultsHandle->accept(index);      
				if (accept) 
				{
					//std::cout << "my trig: " << *trigNameTempl << " & accept = " << accept << std::endl;
					return accept;
				}
		}


		// add prescale information
		/*for ( vector<string>::const_iterator trigNameTempl = triggerPrescaleToStore_.begin();  trigNameTempl != triggerPrescaleToStore_.end(); ++trigNameTempl) {
			int index = TrgNames.triggerIndex(*trigNameTempl);
			int hlt_prescale = 0;
			if (index !=  ntrigs )
				hlt_prescale = hltConfig_.prescaleValue(e, iSetup, *trigNameTempl);

		}
		*/

	} else { std::cout << "%HLTInfo -- No Trigger Result" << std::endl;}

	return accept;
}


// ------------ method called once each job just before starting event loop  ------------
void 
AnomTrkMETFilter::beginJob()
{
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent1(167913,151287912);  //xxxxxx 
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent2(167913,275441701);  //
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent3(167913,286153055);
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent4(167898,966619531);    //xxxx     x
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent5(167898,1278287198);   //xxxx      x
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent6(167898,123612500);
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent7(167898,1905907041);
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent8(167830,907176040); 	  //xxxxxx		 x
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent9(167830,1222131814);			 
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent10(167830,296322641);	  //xxxxxxx		 x
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent11(167830,406445608);
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent12(167830,72090649);    //xxxxxx 	 x
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent13(167830,235429927);	//xxxxxx		 x
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent14(167830,944492575);    //xxxxxx     x

	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent15(167830,1249641785);
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent16(167830,262946685);
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent17(167830,241148247);
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent18(167830,529803472);

	//tai's events
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent19(167151,22174344);
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent20(165415,91591777);
	//tai's event from METB sample
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent21(163270,252238770);

	//3 of my tagged events from metpd with pfmet>1000
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent22(165548,727893578);
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent23(166049,451805356);
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent24(167281,357116007);


	//from jetpd
	/*std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent(167754,38808855);
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent(167754,59289073);
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent(167754,80598975);
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent(167754,88336663);
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent(167740,162687353);
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent(167740,165964901);
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent(167740,84270631);
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent(167740,84724785);
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent(167740,96611397);
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent(167676,113022480);
	
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent(167676,172195782);
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent(167676,177483455);
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent(167676,180490166);
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent(167676,184654252);
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent(167676,193516492);
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent(167676,195189761);
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent(167676,207422354);
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent(,);
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent(,);
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent(,);
	*/
	//
	vBadEvents.push_back(kRunEvent1);
	vBadEvents.push_back(kRunEvent2);
	vBadEvents.push_back(kRunEvent3);
	vBadEvents.push_back(kRunEvent4);
	vBadEvents.push_back(kRunEvent5);
	vBadEvents.push_back(kRunEvent6);
	vBadEvents.push_back(kRunEvent7);
	vBadEvents.push_back(kRunEvent8);
	vBadEvents.push_back(kRunEvent9);
	vBadEvents.push_back(kRunEvent10);
	vBadEvents.push_back(kRunEvent11);
	vBadEvents.push_back(kRunEvent12);
	vBadEvents.push_back(kRunEvent13);
	vBadEvents.push_back(kRunEvent14);
	vBadEvents.push_back(kRunEvent15);
	vBadEvents.push_back(kRunEvent16);
	vBadEvents.push_back(kRunEvent17);
	vBadEvents.push_back(kRunEvent18);

	vBadEvents.push_back(kRunEvent19);
	vBadEvents.push_back(kRunEvent20);
	vBadEvents.push_back(kRunEvent21);

	vBadEvents.push_back(kRunEvent22);
	vBadEvents.push_back(kRunEvent23);
	vBadEvents.push_back(kRunEvent24);
}

bool AnomTrkMETFilter::AnomEvent(const RunLumiEvt_t runevt)
{
//	std::cout << "****** bad evt check! " << vBadEvents.size() << std::endl;
	for (std::vector<std::pair<edm::RunNumber_t, edm::EventNumber_t> >::const_iterator it = vBadEvents.begin(); it != vBadEvents.end(); ++it)
	{
//			std::cout << "****** bad evt = " << runevt.first << " [" << runevt.second << "]" << it->first << ", " << it->second << std::endl;
			if (it->first == runevt.run && it->second == runevt.evt)
			{
				std::cout << "########################### MATCH FOUND: " 
					<< runevt.run << ":" << runevt.lumi << ":" 
					<< runevt.evt << std::endl;
				return true;
			}
	}
	return false;
}


// ------------ method called once each job just after ending the event loop  ------------
void 
AnomTrkMETFilter::endJob() {


	std::cout << __FILE__ << "::" << __FUNCTION__  << "\n" << std::endl;
	std::cout << "[ATM:00] Job settings " << std::endl;
	std::cout << "[ATM:01] Events Processed --- = " << iProcessed << std::endl;
	std::cout << "[ATM:02] Events Passed ------ = " << iPassed << std::endl;
	std::cout << "[ATM:10] Primary Vtx Min ---- = " << inMinVtx << std::endl; 
	std::cout << "[ATM:11] Primary Vtx Min ndof = " << inMinNdofVtx << std::endl; 
	std::cout << "[ATM:12] Max Vtx-Trk Seper. - = " << dMaxPrimVtxZ << std::endl; 
	std::cout << "[ATM:20] Triggers Used ------ = ";
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
	std::cout << "[ATM:30] Minimum MET -------- = " << dminMet << std::endl; 
	std::cout << "[ATM:31] Process only bad? -- = " << processBadOnly << std::endl; 


	/*	std::cout << histsIter4.hPfMet->GetName() << ": Entries = " 
		<< histsIter4.hPfMet->GetEntries() << std::endl;
		std::cout << histsIter4.hPfMet_highPureTrks_atanLt7->GetName() << ": Entries = " 
		<< histsIter4.hPfMet_highPureTrks_atanLt7->GetEntries() << std::endl;
		std::cout << ">>>>>>>>>>> Events with atan>0.7:: iter4 <<<<<<<<<<<<<" << std::endl;
		for (rejIt = it4Rej.begin(); rejIt != it4Rej.end(); rejIt++)
		{
		std::cout << rejIt->first << ", " << rejIt->second << std::endl;
		}
		std::cout << ">>>>>>>>>>> Events with atan>0.7:: iter5 <<<<<<<<<<<<<" << std::endl;
		for (rejIt = it5Rej.begin(); rejIt != it5Rej.end(); rejIt++)
		{
		std::cout << rejIt->first << ", " << rejIt->second << std::endl;
		}
		std::cout << ">>>>>>>>>>> Events with atan>0.7:: pixless <<<<<<<<<<<<<" << std::endl;
		for (rejIt = pxlessRej.begin(); rejIt != pxlessRej.end(); rejIt++)
		{
		std::cout << rejIt->first << ", " << rejIt->second << std::endl;
		}
		*/
	std::cout << " <<<<<<<<<<<<<<<< track ratios "<< std::endl;

	for (std::vector<TrkRatio_t>::const_iterator it = vTrkRatio.begin(); it != vTrkRatio.end(); ++it)
	{
		std::cout << it->run << " , "<< std::setw(15) << std::setprecision(4) << it->event << " = " << it->ratio  << " / " <<std::setw(10) << it->ratio2 << " / " <<std::setw(10) << it->ratio3 << std::endl;

	}

	std::cout << "<<<<<<<<<<<<<< MY SELECTION OF BAD EVENTS " << std::endl;
	std::cout << std::setw(10) << "PassBoth?" << std::setw(10) << "Run"<< std::setw(6) << "Lumi" << std::setw(15) << "Event#" 
		<< std::setw(10) << "PFMET" << std::setw(10) << "TCMET"
		<< std::setw(10) << "CALOMET" << std::setw(12) << "CALO-PF/PF" << std::endl;

	for (std::vector<EvtInfo_t>::const_iterator it = vBadEvents_fromMyCuts.begin(); it != vBadEvents_fromMyCuts.end(); ++it)
	{
		int found = 0;
		for (std::vector<EvtInfo_t>::const_iterator it2 = vBadEvents_fromAndreasCuts.begin(); it2 != vBadEvents_fromAndreasCuts.end(); ++it2)
		{
			if (it->run == it2->run && it->event == it2->event) 
			{
				found = 1;
				break;
			}
		}
		std::cout << std::setw(5) << found << std::setprecision(14) << std::setw(15) << it->run
			<< std::setw(10) << it->lumi << std::setw(15) << it->event 
			<< std::setprecision(5) << std::setw(10) << it->pfmet << std::setw(10) << it->tcmet
			<< std::setw(10) << it->calomet << std::setw(12) << (it->calomet -it->pfmet)/it->pfmet << std::endl;
	}

	std::cout << "<<<<<<<<<<<<<< ANDREAS SELECTION OF BAD EVENTS " << std::endl;
	for (std::vector<EvtInfo_t>::const_iterator it2 = vBadEvents_fromAndreasCuts.begin(); it2 != vBadEvents_fromAndreasCuts.end(); ++it2)
	{
		int found = 0;
		for (std::vector<EvtInfo_t>::const_iterator it = vBadEvents_fromMyCuts.begin(); it != vBadEvents_fromMyCuts.end(); ++it)
		{
			if (it->run == it2->run && it->event == it2->event) 
			{
				found = 1;
				break;
			}
		}
		std::cout << std::setw(5) << found << std::setprecision(14) << std::setw(15) << it2->run 
			<< std::setw(10) << it2->lumi << std::setw(15) << it2->event 
			<< std::setprecision(5) << std::setw(10) << it2->pfmet << std::setw(10) << it2->tcmet
			<< std::setw(10) << it2->calomet << std::setw(12) << (it2->calomet -it2->pfmet)/it2->pfmet << std::endl;
	}

	if (print_hists)
	{
		PrintHist(hCutHist.nTrksAssoWithVtx2ntrks);
		PrintHist(hCutHist.ntrks_nopixhits2ntrks);
		PrintHist(hCutHist.nTrksNotAssoWithVtx2ntrks);
		PrintHist(hCutHist.AfterNtrkCut_ratio);
		PrintHist(hCutHist.AfterNtrkCut_ratio2);
		PrintHist(hCutHist.AfterNtrkCut_ratio3);
		PrintHist(hCutHist.AfterNtrkRatioCut_ratio2);
		PrintHist(hCutHist.AfterNtrkRatioCut_ratio3);
		PrintHist(hCutHist.AfterNtrkRatioRatio2Cut_ratio3);
	}
}

void AnomTrkMETFilter::PrintHist(TH1 *hist)
{
	hist->Draw();
	gPad->Print(GetEPSname(hist->GetName()).c_str());
}
void AnomTrkMETFilter::AnalyseTracks(const edm::Handle<reco::TrackCollection> tracks, 
				const reco::VertexCollection::const_iterator vtx,
				std::vector<unsigned>& vMatchedTracks)

{
/* Dump addition info about the traks associated with a vertex.
 *
 *
 */
	int j = 0;
	int matched = 0;
	if (iVerbose) std::cout << std::setw(10) << "Match?" << std::setw(12) << "highPure?" << std::setw(5) << "algo" << std::setw(12) << "pt" <<  std::setw(15) << "numberOfValidPixelHits" <<  std::endl;   
	for (reco::TrackCollection::const_iterator it = tracks->begin(); it != tracks->end(); ++it)
	{
		++j;
		int i = 0;
		bool matchFound = false;
		for (reco::Vertex::trackRef_iterator vtxTrkIter = vtx->tracks_begin();vtxTrkIter != vtx->tracks_end(); vtxTrkIter++)
		{
			++i;

			if ( (*vtxTrkIter)->pt() == it->pt()
					&& (*vtxTrkIter)->eta() == it->eta()
				   && (*vtxTrkIter)->phi() == it->phi() ) 
			{
				matchFound = true;
				++matched;
				vMatchedTracks.push_back(j);
				//std::cout << "Matching [vtx trk=trklbk][dz] =" << i << " = " << j << "][" << fabs(vtx->z() - it->innerPosition().Z())  << std::endl;
				//hTrkVtxSeparationGt->Fill(fabs(vtx->z() - it->innerPosition().Z()));
				if (iVerbose) std::cout << std::setw(10) << "YES ["<<j << "," << i << "]" 
						<< std::setw(12) << it->quality(reco::TrackBase::highPurity) 
						<< std::setw(5) << it->algo() << std::setw(12) << it->pt() 
						<<  std::setw(15) << it->hitPattern().numberOfValidPixelHits() <<  std::endl;   

				hTrksAssoWithVtx.wgt->Fill(vtx->trackWeight((*vtxTrkIter)));
				FillTrackInfoHists(hTrksAssoWithVtx,it);
			}

		}
		if (! matchFound)
		{
			//std::cout << " NOT Matching trk : innerPosition= " << it->innerPosition().Z()  << std::endl;
			//std::cout << std::setw(12) << "highPure?" << std::setw(5) << "algo" << std::setw(12) << "pt" <<  std::setw(15) << "numberOfValidPixelHits" <<  std::endl;   
			//std::cout << std::setw(12) << it->quality(reco::TrackBase::highPurity) << std::setw(5) << it->algo() << std::setw(12) << it->pt() <<  std::setw(15) << it->hitPattern().numberOfValidPixelHits() <<  std::endl;   
			if (iVerbose) std::cout << std::setw(10) << "NO ["<< j << "," << i<< "]" 
				<< std::setw(12) << it->quality(reco::TrackBase::highPurity) 
				<< std::setw(5) << it->algo() << std::setw(12) << it->pt() 
				<<  std::setw(15) << it->hitPattern().numberOfValidPixelHits() <<  std::endl;   
		}
	}
	if (iVerbose) 
	{
		std::cout << "Found " << matched << " matches from " << j << " tracks." << std::endl;
		std::cout << " >>>>>>>>>>> Total matched tracks = " << vMatchedTracks.size() << std::endl;
	}
	

}

void AnomTrkMETFilter::FillUnmatchTrackHists(edm::Handle<reco::TrackCollection> tracks,
					const std::vector<unsigned> vMatchedTracks)
{

	unsigned k =0;
   double nPureTrks = 0;
	for (reco::TrackCollection::const_iterator it = tracks->begin(); it != tracks->end(); ++it)
	{
		++k;
		if (it->quality(reco::TrackBase::highPurity)) ++nPureTrks;
		bool match = false;
		for (std::vector<unsigned>::const_iterator it2 = vMatchedTracks.begin(); it2 != vMatchedTracks.end(); ++it2)
		{
			if (k == (*it2))
			{
				match = true;
				break;
			}
		}
		if (! match)
		{
			FillTrackInfoHists(hTrksNotAssoWithVtx,it);
		}
	}

	trkPurityFraction->Fill(nPureTrks/tracks->size());

}

// ------------ method called when starting to processes a run  ------------
bool 
AnomTrkMETFilter::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
AnomTrkMETFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
AnomTrkMETFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
AnomTrkMETFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
AnomTrkMETFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


void AnomTrkMETFilter::PFlowAnalyse(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	edm::Handle<reco::PFJetCollection> pfJetcollection;
	iEvent.getByLabel(pfJetInputTag_, pfJetcollection);

	//std::cout << __FILE__ << ":" << __FUNCTION__ << std::endl;
	if (! pfJetcollection.isValid()) 
	{
		std::cout << __FILE__ << ":" << __FUNCTION__ << ": Valid PFJetCollection not found!" << std::endl;
		assert(false);
	}  
	for (reco::PFJetCollection::const_iterator it = (*pfJetcollection).begin(); it != (*pfJetcollection).end(); ++it)
	{
		//if (iVerbose)
			//std::cout << __LINE__ << ":" << it->chargedMultiplicity() << std::endl;

	}


//	Handle<reco::PFCandidateCollection> pfCandidates;
//	iEvent.getByLabel("particleFlow",pfCandidates);
/*

 typedef math::XYZTLorentzVector LorentzVector;
 typedef math::XYZPoint Point;
 double sum_et = 0.0;
 double sum_ex = 0.0;
 double sum_ey = 0.0;
 double sum_ez = 0.0;
 double NeutralEMEt = 0.0;
 double NeutralHadEt = 0.0;
 double ChargedEMEt = 0.0;
 double ChargedHadEt = 0.0;
 double MuonEt = 0.0;
 double type6Et = 0.0;
 double type7Et = 0.0;
 for (unsigned int pfc=0;pfc<pfCandidates.size();++pfc) {
    double phi   = pfCandidates[pfc].phi();
    double theta = pfCandidates[pfc].theta();
    double e     = pfCandidates[pfc].energy();
    double et    = e*sin(theta);
    sum_ez += e*cos(theta);
    sum_et += et;
    sum_ex += et*cos(phi);
    sum_ey += et*sin(phi);

    // compute met specific data:
    if (pfCandidates[pfc].particleId() == 1) ChargedHadEt += et;
    if (pfCandidates[pfc].particleId() == 2) ChargedEMEt += et;
    if (pfCandidates[pfc].particleId() == 3) MuonEt += et;
    if (pfCandidates[pfc].particleId() == 4) NeutralEMEt += et;
    if (pfCandidates[pfc].particleId() == 5) NeutralHadEt += et;
    if (pfCandidates[pfc].particleId() == 6) type6Et += et;
    if (pfCandidates[pfc].particleId() == 7) type7Et += et;

  }

  const double Et_total=NeutralEMEt+NeutralHadEt+ChargedEMEt+ChargedHadEt+MuonEt+type6Et+type7Et;

  double met = sqrt( sum_ex*sum_ex + sum_ey*sum_ey );
  const LorentzVector p4( -sum_ex, -sum_ey, 0.0, met);
  const Point vtx(0.0,0.0,0.0);
 
  SpecificPFMETData specific;
  specific.NeutralEMFraction = NeutralEMEt/Et_total;
  specific.NeutralHadFraction = NeutralHadEt/Et_total;
  specific.ChargedEMFraction = ChargedEMEt/Et_total;
  specific.ChargedHadFraction = ChargedHadEt/Et_total;
  specific.MuonFraction = MuonEt/Et_total;
  specific.Type6Fraction = type6Et/Et_total;
  specific.Type7Fraction = type7Et/Et_total;

  reco::PFMET specificPFMET( specific, sum_et, p4, vtx );
  return specificPFMET;
*/
}

void AnomTrkMETFilter::FillTrackInfoHists(TrkInfo_t trkHist, reco::TrackCollection::const_iterator it)
{
	trkHist.algo->Fill(it->algo());
	trkHist.validPixHits->Fill(it->hitPattern().numberOfValidPixelHits());
	trkHist.purity->Fill(it->quality(reco::TrackBase::highPurity));
	trkHist.pt->Fill(it->pt());
	trkHist.ptErr->Fill(it->ptError()/it->pt());
	trkHist.dxyErr->Fill(it->dxyError()/it->dxy());
	trkHist.d0Err->Fill(it->d0Error()/it->d0());
	trkHist.exptInnerHits->Fill(it->trackerExpectedHitsInner().numberOfLostTrackerHits());
	trkHist.exptOuterHits->Fill(it->trackerExpectedHitsOuter().numberOfLostTrackerHits());
	trkHist.ptVsptErr->Fill(it->pt(),it->ptError()/it->pt());
	trkHist.validHitFraction->Fill(it->validFraction());
}

//define this as a plug-in
DEFINE_FWK_MODULE(AnomTrkMETFilter);

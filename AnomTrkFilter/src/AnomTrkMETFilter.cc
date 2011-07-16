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
// $Id: AnomTrkMETFilter.cc,v 1.2 2011/07/12 03:49:00 samantha Exp $
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
		void AnalyseTracks(const edm::Handle<reco::TrackCollection> tracks, const reco::VertexCollection::const_iterator vtx);
		bool storeHLT(const edm::Event& e, const edm::EventSetup& iSetup);
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

		//	 std::vector<edm::RunNumber_t> startrun;
		//	 std::vector<edm::EventNumber_t> endrun;

		struct Hist_t{
			TProfile *hAtanVsMet;
			//TH1F *hNtrks;
			//TH1F *hTrksPerVtx;
			//TH1F *hNvtx;
			//TH1F *hNpets;
			TH1F *hPfMet;
			TH1F *hPfMetSig;
			//High Purity Tracks atan>0.7
			//TH1F *hNtrks_highPureTrks_atanGt7;
			//TH1F *hTrksPerVtx_highPureTrks_atanGt7;
			//TH1F *hNvtx_highPureTrks_atanGt7;
			TH1F *hPfMet_highPureTrks_atanGt7;
			TH1F *hPfMetSig_highPureTrks_atanGt7;
			//High Purity Tracks atan<0.7
			//TH1F *hNtrks_highPureTrks_atanLt7;
			//TH1F *hTrksPerVtx_highPureTrks_atanLt7;
			//TH1F *hNvtx_highPureTrks_atanLt7;
			TH1F *hPfMet_highPureTrks_atanLt7;
			TH1F *hPfMetSig_highPureTrks_atanLt7;
		};

		struct TrackHist_t{
			TH1F *hQuality;
			TH1F *hNPixelHits;
			TH1F *hNStripHits;
			TH1F *hTrkDelZfromPrimVtx;
			TH1F *hRatioIt4toIt01;
			TH1F *hRatioIt5toIt01;
			TH1F *hRatioPixlesstoIt01;
			TH1F *hRatioIt4toIt01pure;
			TH1F *hRatioIt5toIt01pure;
			TH1F *hRatioPixlesstoIt01pure;
		};
		TH1F* hPrimVtxz; //z cdt of the primary vertex
		//TH1F* hTrksPerVtx;
		//TH2F* hTrksPerVtxGevBins;
		//TH1F* hTrkVtxSeparation;
		//TH1F* hTrksPerVtxNtrksRatio;
		//TH1F* hTrkVtxSeparationGt;
		//TH1F* hTrkVtxSeparationLt;
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
		FinalHist_t hFinal_all, hFinal_my, hFinal_Andrea;

		TrackHist_t hTrkIter0, hTrkIter1, hTrkIter2, hTrkIter3, hTrkIter4, hTrkIter5;
		Hist_t histsIter4, histsIter5, histsPixless;

		struct RunLumEvt_t
		{
			unsigned run;
			unsigned lumi;
			unsigned evt;
		};

		std::vector<RunLumEvt_t> it4Rej, it5Rej, pxlessRej;
		std::vector<RunLumEvt_t>::iterator rejIt;

		//jet collections
		edm::InputTag caloJetInputTag_;
		edm::Handle<reco::CaloJetCollection> caloJetcollection;

		edm::InputTag pfJetInputTag_;
		edm::Handle<reco::PFJetCollection> pfJetcollection;
		
		struct TrkRatio_t
		{
			double run;
			double event;
			double ratio;
			double ratio2;
			double ratio3;
		};
		std::vector<TrkRatio_t> vTrkRatio;
		std::vector<std::pair<edm::RunNumber_t, edm::EventNumber_t> > vBadEvents;
		bool AnomEvent(const std::pair<edm::RunNumber_t, edm::EventNumber_t> runevt);
		struct EvtInfo_t
		{
			double run;
			double lumi;
			double event;
			double pfmet;
			double tcmet;
			double calomet;
		};

		std::vector<EvtInfo_t> vBadEvents_fromMyCuts, vBadEvents_fromAndreasCuts;
		int iVerbose; // control print levels
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
		CutHists_t hCutHist;
		bool print_hists;
		std::string GetEPSname(const std::string name) 
		{
			std::stringstream s;
			s << name << ".eps";
			return s.str();
		}
		void PrintHist(TH1 *hist);
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
	print_hists = iConfig.getUntrackedParameter<bool>("print_hists", false);
	iProcessed = 0;
	iPassed = 0;

	//generate hists
	edm::Service<TFileService> fs;

	const double met_max = 4000, met_bins = 400;

	histsIter4.hAtanVsMet = fs->make<TProfile> ("Atan2VsMet_iter4" ,"atan2(iter4/iter0+1) Vs PF#slash{E}_{T}", met_bins,0,met_max);
//	histsIter4.hNvtx  = fs->make<TH1F> ("NVtx_iter4" ," N Vertices", 50,0,50);
	histsIter4.hPfMet  = fs->make<TH1F> ("PfMet_iter4" ," Raw #slash{E}_{T} (PFMETCollection::pfMet)", met_bins,0,met_max);
	histsIter4.hPfMetSig  = fs->make<TH1F> ("PfMetSig_iter4" ," #slash{E}_{T}-Sig (PFMETCollection::metSig)", 200,0,20);
//	histsIter4.hTrksPerVtx = fs->make<TH1F> ("NtrksPerVtx_iter4" ," N Tracks per Vertex", 250,0,500);

	//high purity tracks
//	histsIter4.hNtrks_highPureTrks_atanGt7 = fs->make<TH1F> ("Ntrks_highPureTrks_atanGt7_iter4" ,"High Purity Tracks atan2(iter4/iter0+1)>0.7: N Tracks", 250,0,500);
//	histsIter4.hNvtx_highPureTrks_atanGt7  = fs->make<TH1F> ("NVtx_highPureTrks_atanGt7_iter4" ,"High Purity Tracks atan2(iter4/iter0+1)>0.7:  N Vertices", 50,0,50);
	histsIter4.hPfMet_highPureTrks_atanGt7  = fs->make<TH1F> ("PfMet_highPureTrks_atanGt7_iter4" ,"High Purity Tracks atan2(iter4/iter0+1)>0.7:  PF#slash{E}_{T} (PFMETCollection::pfMet)", met_bins,0,met_max);
	histsIter4.hPfMetSig_highPureTrks_atanGt7  = fs->make<TH1F> ("PfMetSig_highPureTrks_atanGt7_iter4" ,"High Purity Tracks atan2(iter4/iter0+1)>0.7: PF#slash{E}_{T}-Sig (PFMETCollection::metSig)", 200,0,20);
//	histsIter4.hTrksPerVtx_highPureTrks_atanGt7 = fs->make<TH1F> ("NtrksPerVtx_highPureTrks_atanGt7_iter4" ,"High Purity Tracks atan2(iter4/iter0+1)>0.7:  N Tracks per Vertex", 250,0,500);

//	histsIter4.hNtrks_highPureTrks_atanLt7 = fs->make<TH1F> ("Ntrks_highPureTrks_atanLt7_iter4" ,"High Purity Tracks atan2(iter4/iter0+1)<0.7: N Tracks", 250,0,500);
//	histsIter4.hNvtx_highPureTrks_atanLt7  = fs->make<TH1F> ("NVtx_highPureTrks_atanLt7_iter4" ,"High Purity Tracks atan2(iter4/iter0+1)<0.7:  N Vertices", 50,0,50);
	histsIter4.hPfMet_highPureTrks_atanLt7  = fs->make<TH1F> ("PfMet_highPureTrks_atanLt7_iter4" ,"High Purity Tracks atan2(iter4/iter0+1)<0.7:  PF#slash{E}_{T} (PFMETCollection::pfMet)", met_bins,0,met_max);
	histsIter4.hPfMetSig_highPureTrks_atanLt7  = fs->make<TH1F> ("PfMetSig_highPureTrks_atanLt7_iter4" ,"High Purity Tracks atan2(iter4/iter0+1)<0.7: PF#slash{E}_{T}-Sig (PFMETCollection::metSig)", 200,0,20);
//	histsIter4.hTrksPerVtx_highPureTrks_atanLt7 = fs->make<TH1F> ("NtrksPerVtx_highPureTrks_atanLt7_iter4" ,"High Purity Tracks atan2(iter4/iter0+1)<0.7:  N Tracks per Vertex", 250,0,500);


	histsIter5.hAtanVsMet = fs->make<TProfile> ("Atan2VsMet_iter5" ,"atan2(iter5/iter0+1) Vs PF#slash{E}_{T}", 200,0,1000);

	//high purity tracks
//	histsIter5.hNtrks_highPureTrks_atanGt7 = fs->make<TH1F> ("Ntrks_highPureTrks_atanGt7_iter5" ,"High Purity Tracks atan>(iter5/iter0+1)0.7: N Tracks", 250,0,500);
//	histsIter5.hNvtx_highPureTrks_atanGt7  = fs->make<TH1F> ("NVtx_highPureTrks_atanGt7_iter5" ,"High Purity Tracks atan>(iter5/iter0+1)0.7:  N Vertices", 50,0,50);
	histsIter5.hPfMet_highPureTrks_atanGt7  = fs->make<TH1F> ("PfMet_highPureTrks_atanGt7_iter5" ,"High Purity Tracks atan>(iter5/iter0+1)0.7:  PF#slash{E}_{T} (PFMETCollection::pfMet)", met_bins,0,met_max);
	histsIter5.hPfMetSig_highPureTrks_atanGt7  = fs->make<TH1F> ("PfMetSig_highPureTrks_atanGt7_iter5" ,"High Purity Tracks atan2(iter5/iter0+1)>0.7: PF#slash{E}_{T}-Sig (PFMETCollection::metSig)", 200,0,20);
//	histsIter5.hTrksPerVtx_highPureTrks_atanGt7 = fs->make<TH1F> ("NtrksPerVtx_highPureTrks_atanGt7_iter5" ,"High Purity Tracks atan>(iter5/iter0+1)0.7:  N Tracks per Vertex", 250,0,500);

//	histsIter5.hNtrks_highPureTrks_atanLt7 = fs->make<TH1F> ("Ntrks_highPureTrks_atanLt7_iter5" ,"High Purity Tracks atan<(iter5/iter0+1)0.7: N Tracks", 250,0,500);
//	histsIter5.hNvtx_highPureTrks_atanLt7  = fs->make<TH1F> ("NVtx_highPureTrks_atanLt7_iter5" ,"High Purity Tracks atan<(iter5/iter0+1)0.7:  N Vertices", 50,0,50);
	histsIter5.hPfMet_highPureTrks_atanLt7  = fs->make<TH1F> ("PfMet_highPureTrks_atanLt7_iter5" ,"High Purity Tracks atan<(iter5/iter0+1)0.7:  PF#slash{E}_{T} (PFMETCollection::pfMet)", met_bins,0,met_max);
	histsIter5.hPfMetSig_highPureTrks_atanLt7  = fs->make<TH1F> ("PfMetSig_highPureTrks_atanLt7_iter5" ,"High Purity Tracks atan2(iter5/iter0+1)<0.7: PF#slash{E}_{T}-Sig (PFMETCollection::metSig)", 200,0,20);
//	histsIter5.hTrksPerVtx_highPureTrks_atanLt7 = fs->make<TH1F> ("NtrksPerVtx_highPureTrks_atanLt7_iter5" ,"High Purity Tracks atan<(iter5/iter0+1)0.7:  N Tracks per Vertex", 250,0,500);


	//pixelless stuff from pure tracks
	histsPixless.hAtanVsMet = fs->make<TProfile> ("Atan2VsMet_Pixless" ,"atan2(Pixelless/iter0+1) Vs PF#slash{E}_{T}", met_bins,0,met_max);

//	histsPixless.hNtrks_highPureTrks_atanGt7 = fs->make<TH1F> ("Ntrks_highPureTrks_atanGt7_Pixless" ,"High Purity Tracks atan>(Pixless/iter0+1)0.7: N Tracks", 250,0,500);
//	histsPixless.hNvtx_highPureTrks_atanGt7  = fs->make<TH1F> ("NVtx_highPureTrks_atanGt7_Pixless" ,"High Purity Tracks atan>(Pixless/iter0+1)0.7:  N Vertices", 50,0,50);
	histsPixless.hPfMet_highPureTrks_atanGt7  = fs->make<TH1F> ("PfMet_highPureTrks_atanGt7_Pixless" ,"High Purity Tracks atan>(Pixless/iter0+1)0.7:  PF#slash{E}_{T} (PFMETCollection::pfMet)", met_bins,0,met_max);
	histsPixless.hPfMetSig_highPureTrks_atanGt7  = fs->make<TH1F> ("PfMetSig_highPureTrks_atanGt7_Pixless" ,"High Purity Tracks atan2(Pixless/iter0+1)>0.7: PF#slash{E}_{T}-Sig (PFMETCollection::metSig)", 200,0,20);
//	histsPixless.hTrksPerVtx_highPureTrks_atanGt7 = fs->make<TH1F> ("NtrksPerVtx_highPureTrks_atanGt7_Pixless" ,"High Purity Tracks atan>(Pixless/iter0+1)0.7:  N Tracks per Vertex", 250,0,500);

//	histsPixless.hNtrks_highPureTrks_atanLt7 = fs->make<TH1F> ("Ntrks_highPureTrks_atanLt7_Pixless" ,"High Purity Tracks atan2(Pixless/iter0+1)<0.7: N Tracks", 250,0,500);
//	histsPixless.hNvtx_highPureTrks_atanLt7  = fs->make<TH1F> ("NVtx_highPureTrks_atanLt7_Pixless" ,"High Purity Tracks atan2(Pixless/iter0+1)<0.7:  N Vertices", 50,0,50);
	histsPixless.hPfMet_highPureTrks_atanLt7  = fs->make<TH1F> ("PfMet_highPureTrks_atanLt7_Pixless" ,"High Purity Tracks atan2(Pixless/iter0+1)<0.7:  PF#slash{E}_{T} (PFMETCollection::pfMet)", met_bins,0,met_max);
	histsPixless.hPfMetSig_highPureTrks_atanLt7  = fs->make<TH1F> ("PfMetSig_highPureTrks_atanLt7_Pixless" ,"High Purity Tracks atan2(Pixless/iter0+1)<0.7: PF#slash{E}_{T}-Sig (PFMETCollection::metSig)", 200,0,20);

	hTrkIter0.hQuality = fs->make<TH1F> ("TrkQuality","Trk Quality Flag",10,0,10);
	//hTrkIter0.hNPixelHits = fs->make<TH1F> ("NPixelHits","Trk N Pixel Hits",100,0,100);
	//hTrkIter0.hNStripHits = fs->make<TH1F> ("NStripHits","Trk N Strip Hits",100,0,100);
	//hTrkIter0.hTrkDelZfromPrimVtx = fs->make<TH1F> ("TrkDelZfromPrimVtx","#Delta z = abs(z^{trk} - z^{prim. vtx})",200,0,200);
	hTrkIter0.hRatioIt4toIt01 = fs->make<TH1F> ("hRatioIt4toIt01","General Tracks: atan2(iter4/Iter0+1);atan2(ratio);Events;",200,0,2);
	hTrkIter0.hRatioIt5toIt01 = fs->make<TH1F> ("hRatioIt5toIt01","General Tracks: atan2(iter5/Iter0+1);atan2(ratio);Events;",200,0,2);
	hTrkIter0.hRatioPixlesstoIt01 = fs->make<TH1F> ("hRatioPixlesstoIt01","General Tracks: atan2(Pixless/Iter0+1);atan2(ratio);Events;",200,0,2);
	hTrkIter0.hRatioIt4toIt01pure = fs->make<TH1F> ("hRatioIt4toIt01pure","High Purity Tracks: atan2(iter4/Iter0+1);atan2(ratio);Events;",200,0,2);
	hTrkIter0.hRatioIt5toIt01pure = fs->make<TH1F> ("hRatioIt5toIt01pure","High Purity Tracks: atan2(iter5/Iter0+1);atan2(ratio);Events;",200,0,2);
	hTrkIter0.hRatioPixlesstoIt01pure = fs->make<TH1F> ("hRatioPixlesstoIt01pure","High Purity Tracks: atan2(Pixless/Iter0+1);atan2(ratio);Events;",200,0,2);

	hPrimVtxz = fs->make<TH1F> ("PrimVtxZ","Primvary Vertex z-position;z [cm];Events;",300,0,300);
	//hTrksPerVtx = fs->make<TH1F> ("trksPerVtx","Tracks per vtx;Tracks per vtx;Events;",300,0,300);
	//hTrksPerVtxGevBins = fs->make<TH2F> ("trksPerVtxGeVBins","Tracks per vtx in GeVBin;Tracks per vtx in 5GeV bins;Events;",20,0,100,50,0,50);
	//hTrkVtxSeparation = fs->make<TH1F> ("trkVtxSeparation","Track-vtx z separation;#Delta z;Events",100,0,100);
	//hTrksPerVtxNtrksRatio	 = fs->make<TH1F> ("trksPetrVtxNtrksRation","Traks per vertex/ntracks;ratio;Events",100,0,1);
	//hTrkVtxSeparationGt = fs->make<TH1F> ("trkVtxSeparationGt","Track-vtx z separation for trk ratio > ;#Delta z;Events",100,0,100);
	//hTrkVtxSeparationLt = fs->make<TH1F> ("trkVtxSeparationLt","Track-vtx z separation for trk ratio < ;#Delta z;Events",100,0,100);

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
/*#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
*/
	++iProcessed;
	const double fMinTrkPt = 0.9;
	bool bSaveEvent = false;

	RunNumber_t kRun   = iEvent.id().run();
	EventNumber_t kEvent = iEvent.id().event();
	LuminosityBlockNumber_t kLumi  = iEvent.id().luminosityBlock(); 
	RunLumEvt_t kRunLumiEvent;
	kRunLumiEvent.run = kRun;
	kRunLumiEvent.lumi = kLumi;
	kRunLumiEvent.evt = kEvent;

	//if (! AnomEvent(kRunLumiEvent)) return 0;
	//std::cout << "====== processing event " << kRun << ", " << kEvent << std::endl;
	//if (kEvent != 22174344 || kEvent != 91591777) return 0;
	//if (kEvent != 91591777) return 0;
	//std::cout << "====== FOUND EVENT event " << kRun << ", " << kEvent << std::endl;
	if (iVerbose) std::cout << "====== processing event " << kRun << ", " << kEvent << std::endl;

	//========== Trigger Selections
	if (req_trigger)
	{
		if (! storeHLT(iEvent, iSetup)) return 0;
	}

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
			//if (trkpt< fMinTrkPt) continue;
			jj++;
			//std::cout << "trk [" << jj << "]"<< std::setw(10) << it->pt() << "/"<< std::setw(10) << it->phi() << "/"<< std::setw(10)<< it->eta() << std::endl;
			ntrks++;
			if (trkpt>0.9) ++ntrks_above9;
			if (it->hitPattern().numberOfValidPixelHits() ==0) ++ntrks_nopixhits; 
	}

 // ==========================================================
 //Vertex information

	Handle<reco::VertexCollection> vertexHandle;
  	iEvent.getByLabel("offlinePrimaryVertices", vertexHandle);
	double dNvtx = 0;
	double dPrimVtx_z = 0;
	double nTrksAssoWithVtx = 0;
   if (vertexHandle.isValid())
	{
     reco::VertexCollection vertexCollection = *(vertexHandle.product());
     //std::cout << "VtxSize = "<<  vertexCollection.size() << std::endl;
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
				//AnalyseTracks(tracks, v);
				//loop over tracks associated with the vertex
				reco::Vertex::trackRef_iterator trackIter = v->tracks_begin();
				double nTrk = 0;
				//hTrksPerVtx->Fill(v->tracksSize());

				for (;trackIter != v->tracks_end(); trackIter++)
				{
					const double trkpt = (*trackIter)->pt();
					if (trkpt<fMinTrkPt) continue;
					nTrk++;
					const double trk_z = (*trackIter)->innerPosition().Z();
					//hTrkVtxSeparation->Fill(abs(dPrimVtx_z - trk_z));

				}
				nTrksAssoWithVtx += nTrk;
				//if (ntrks>0) hTrksPerVtxNtrksRatio->Fill(nTrk/ntrks);

//				std::cout << "for vtx[" << i << "] trks for bins [0-5,5-10,15-20,20,25] =" << 
//							trkPer5GevBin[0] << ", " << trkPer5GevBin[1] << ", " << trkPer5GevBin[2]
//							<< ", " << trkPer5GevBin[3] << ", " << trkPer5GevBin[4] << ", " << trkPer5GevBin[5] << std::endl;
//				std::cout << "for vtx[" << i << "] ratio to ntrks: "<< nTrk << "/" << ntrks <<" = " << nTrk/ntrks << std::endl;

			  //break;
		  }
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
		//19822 ====== FOUND EVENT event 165415, 91591777
		//19823 ntrks/ratio/ratio2/ratio3 = 79/0.109649/0.298246/0.334975
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

}

bool AnomTrkMETFilter::AnomEvent(const std::pair<edm::RunNumber_t, edm::EventNumber_t> runevt)
{
//	std::cout << "****** bad evt check! " << vBadEvents.size() << std::endl;
	for (std::vector<std::pair<edm::RunNumber_t, edm::EventNumber_t> >::const_iterator it = vBadEvents.begin(); it != vBadEvents.end(); ++it)
	{
//			std::cout << "****** bad evt = " << runevt.first << " [" << runevt.second << "]" << it->first << ", " << it->second << std::endl;
			if (it->first == runevt.first && it->second == runevt.second)
			{
				std::cout << "########################### MATCH FOUND " << std::endl;
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
	std::cout << "[ATM:01] Events Passed ------ = " << iPassed << std::endl;
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
void AnomTrkMETFilter::AnalyseTracks(const edm::Handle<reco::TrackCollection> tracks, const reco::VertexCollection::const_iterator vtx)
{
	int j = 0;
	int matched = 0;
	for (reco::TrackCollection::const_iterator it = tracks->begin(); it != tracks->end(); ++it)
	{
		++j;
		int i = 0;
		bool matchFound = false;
		for (reco::Vertex::trackRef_iterator vtxTrkIter = vtx->tracks_begin();vtxTrkIter != vtx->tracks_end(); vtxTrkIter++)
		{
			++i;

			if ((*vtxTrkIter)->pt() == it->pt()
					&& (*vtxTrkIter)->eta() == it->eta()
				 && (*vtxTrkIter)->phi() == it->phi()) 
			{
				matchFound = true;
				++matched;
				//std::cout << "Matching [vtx trk=trklbk][dz] =" << i << " = " << j << "][" << fabs(vtx->z() - it->innerPosition().Z())  << std::endl;
				//hTrkVtxSeparationGt->Fill(fabs(vtx->z() - it->innerPosition().Z()));
				//std::cout << std::setw(12) << "highPure?" << std::setw(5) << "algo" << std::setw(12) << "pt" <<  std::setw(15) << "numberOfValidPixelHits" <<  std::endl;   
				//std::cout << std::setw(12) << it->quality(reco::TrackBase::highPurity) << std::setw(5) << it->algo() << std::setw(12) << it->pt() <<  std::setw(15) << it->hitPattern().numberOfValidPixelHits() <<  std::endl;   
			}

		}
		if (! matchFound)
		{
			//std::cout << " NOT Matching trk : innerPosition= " << it->innerPosition().Z()  << std::endl;
			//std::cout << std::setw(12) << "highPure?" << std::setw(5) << "algo" << std::setw(12) << "pt" <<  std::setw(15) << "numberOfValidPixelHits" <<  std::endl;   
			//std::cout << std::setw(12) << it->quality(reco::TrackBase::highPurity) << std::setw(5) << it->algo() << std::setw(12) << it->pt() <<  std::setw(15) << it->hitPattern().numberOfValidPixelHits() <<  std::endl;   
		}
	}
	std::cout << "Found " << matched << " matches from " << j << " tracks." << std::endl;

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
//define this as a plug-in
DEFINE_FWK_MODULE(AnomTrkMETFilter);

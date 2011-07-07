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
// $Id$
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
#include "TH1F.h"
#include "assert.h"
#include "TProfile.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
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

		//	 std::vector<edm::RunNumber_t> startrun;
		//	 std::vector<edm::EventNumber_t> endrun;

		struct Hist_t{
			TProfile *hAtanVsMet;
			//TH1F *hNtrks;
			TH1F *hTrksPerVtx;
			TH1F *hNvtx;
			//TH1F *hNJets;
			TH1F *hRawHt;
			TH1F *hRawMet;
			TH1F *hPfHt;
			TH1F *hPfMet;
			TH1F *hPfMetSig;
			//General Tracks atan>0.7
			TH1F *hNtrks_genTrks_atanGt7;
			TH1F *hTrksPerVtx_genTrks_atanGt7;
			TH1F *hNvtx_genTrks_atanGt7;
			TH1F *hRawHt_genTrks_atanGt7;
			TH1F *hRawMet_genTrks_atanGt7;
			TH1F *hPfHt_genTrks_atanGt7;
			TH1F *hPfMet_genTrks_atanGt7;
			TH1F *hPfMetSig_genTrks_atanGt7;
			//General Tracks atan<0.7
			TH1F *hNtrks_genTrks_atanLt7;
			TH1F *hTrksPerVtx_genTrks_atanLt7;
			TH1F *hNvtx_genTrks_atanLt7;
			TH1F *hRawHt_genTrks_atanLt7;
			TH1F *hRawMet_genTrks_atanLt7;
			TH1F *hPfHt_genTrks_atanLt7;
			TH1F *hPfMet_genTrks_atanLt7;
			TH1F *hPfMetSig_genTrks_atanLt7;
			//High Purity Tracks atan>0.7
			TH1F *hNtrks_highPureTrks_atanGt7;
			TH1F *hTrksPerVtx_highPureTrks_atanGt7;
			TH1F *hNvtx_highPureTrks_atanGt7;
			TH1F *hRawHt_highPureTrks_atanGt7;
			TH1F *hRawMet_highPureTrks_atanGt7;
			TH1F *hPfHt_highPureTrks_atanGt7;
			TH1F *hPfMet_highPureTrks_atanGt7;
			TH1F *hPfMetSig_highPureTrks_atanGt7;
			//High Purity Tracks atan<0.7
			TH1F *hNtrks_highPureTrks_atanLt7;
			TH1F *hTrksPerVtx_highPureTrks_atanLt7;
			TH1F *hNvtx_highPureTrks_atanLt7;
			TH1F *hRawHt_highPureTrks_atanLt7;
			TH1F *hRawMet_highPureTrks_atanLt7;
			TH1F *hPfHt_highPureTrks_atanLt7;
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

		TrackHist_t hTrkIter0, hTrkIter1, hTrkIter2, hTrkIter3, hTrkIter4, hTrkIter5;
		Hist_t histsIter4, histsIter5, histsPixless;

		std::set<std::pair<edm::EventNumber_t, edm::RunNumber_t> > it4Rej, it5Rej, pxlessRej;
		std::set<std::pair<edm::EventNumber_t, edm::RunNumber_t> >::iterator rejIt;

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
	iProcessed = 0;
	iPassed = 0;

	//generate hists
	edm::Service<TFileService> fs;

	const double met_max = 4000, met_bins = 400;
	const double ht_max = 4000, ht_bins = 400;

	histsIter4.hAtanVsMet = fs->make<TProfile> ("Atan2VsMet_iter4" ,"atan2(iter4/iter0+1) Vs PF#slash{E}_{T}", met_bins,0,met_max);
	//histsIter4.hNtrks = fs->make<TH1F> ("Ntrks_iter4" ," N Tracks", 250,0,500);
	histsIter4.hNvtx  = fs->make<TH1F> ("NVtx_iter4" ," N Vertices", 50,0,50);
//	histsIter4.hNJets  = fs->make<TH1F> ("NJets_iter4" ," NJets", 50,0,50);
	histsIter4.hRawMet  = fs->make<TH1F> ("CaloRawMet_iter4" ," Raw #slash{E}_{T} (CaloMETCollection::met)", met_bins,0,met_max);
	histsIter4.hRawHt  = fs->make<TH1F> ("CaloRawHt_iter4" ," H_{T} (CaloMETCollection::H_{T})", ht_bins,0,ht_max);
	histsIter4.hPfMet  = fs->make<TH1F> ("PfMet_iter4" ," Raw #slash{E}_{T} (PFMETCollection::pfMet)", met_bins,0,met_max);
	histsIter4.hPfHt  = fs->make<TH1F> ("PfHt_iter4" ," H_{T} (PFMETCollection::H_{T})", ht_bins,0,ht_max);
	histsIter4.hPfMetSig  = fs->make<TH1F> ("PfMetSig_iter4" ," #slash{E}_{T}-Sig (PFMETCollection::metSig)", 200,0,20);
	histsIter4.hTrksPerVtx = fs->make<TH1F> ("NtrksPerVtx_iter4" ," N Tracks per Vertex", 250,0,500);


	//general tracks
	histsIter4.hNtrks_genTrks_atanGt7 = fs->make<TH1F> ("Ntrks_genTrks_atanGt7_iter4" ,"General Tracks atan2(iter4/iter0+1)>0.7: N Tracks", 250,0,500);
	histsIter4.hNvtx_genTrks_atanGt7  = fs->make<TH1F> ("NVtx_genTrks_atanGt7_iter4" ,"General Tracks atan2(iter4/iter0+1)>0.7:  N Vertices", 50,0,50);
	histsIter4.hRawMet_genTrks_atanGt7  = fs->make<TH1F> ("CaloRawMet_genTrks_atanGt7_iter4" ,"General Tracks atan2(iter4/iter0+1)>0.7:  Raw #slash{E}_{T} (CaloMETCollection::met)", met_bins,0,met_max);
	histsIter4.hRawHt_genTrks_atanGt7  = fs->make<TH1F> ("CaloRawHt_genTrks_atanGt7_iter4" ,"General Tracks atan2(iter4/iter0+1)>0.7:  H_{T} (CaloMETCollection::H_{T})", ht_bins,0,ht_max);
	histsIter4.hPfMet_genTrks_atanGt7  = fs->make<TH1F> ("PfMet_genTrks_atanGt7_iter4" ,"General Tracks atan2(iter4/iter0+1)>0.7:  PF#slash{E}_{T} (PFMETCollection::pfMet)", met_bins,0,met_max);
	histsIter4.hPfMetSig_genTrks_atanGt7  = fs->make<TH1F> ("PfMetSig_genTrks_atanGt7_iter4" ,"General Tracks atan2(iter4/iter0+1)<0.7: PF#slash{E}_{T}-Sig (PFMETCollection::metSig)", 200,0,20);
	histsIter4.hPfHt_genTrks_atanGt7  = fs->make<TH1F> ("PfHt_genTrks_atanGt7_iter4" ,"General Tracks atan2(iter4/iter0+1)>0.7:  PFH_{T} (PFMETCollection::H_{T})", ht_bins,0,ht_max);
	histsIter4.hTrksPerVtx_genTrks_atanGt7 = fs->make<TH1F> ("NtrksPerVtx_genTrks_atanGt7_iter4" ,"General Tracks atan2(iter4/iter0+1)>0.7:  N Tracks per Vertex", 250,0,500);

	histsIter4.hNtrks_genTrks_atanLt7 = fs->make<TH1F> ("Ntrks_genTrks_atanLt7_iter4" ,"General Tracks atan2(iter4/iter0+1)<0.7: N Tracks", 250,0,500);
	histsIter4.hNvtx_genTrks_atanLt7  = fs->make<TH1F> ("NVtx_genTrks_atanLt7_iter4" ,"General Tracks atan2(iter4/iter0+1)<0.7:  N Vertices", 50,0,50);
	histsIter4.hRawMet_genTrks_atanLt7  = fs->make<TH1F> ("CaloRawMet_genTrks_atanLt7_iter4" ,"General Tracks atan2(iter4/iter0+1)<0.7:  Raw #slash{E}_{T} (CaloMETCollection::met)", met_bins,0,met_max);
	histsIter4.hRawHt_genTrks_atanLt7  = fs->make<TH1F> ("CaloRawHt_genTrks_atanLt7_iter4" ,"General Tracks atan2(iter4/iter0+1)<0.7:  H_{T} (CaloMETCollection::H_{T})", ht_bins,0,ht_max);
	histsIter4.hPfMet_genTrks_atanLt7  = fs->make<TH1F> ("PfMet_genTrks_atanLt7_iter4" ,"General Tracks atan2(iter4/iter0+1)<0.7:  PF#slash{E}_{T} (PFMETCollection::pfMet)", met_bins,0,met_max);
	histsIter4.hPfMetSig_genTrks_atanLt7  = fs->make<TH1F> ("PfMetSig_genTrks_atanLt7_iter4" ,"General Tracks atan2(iter4/iter0+1)<0.7: PF#slash{E}_{T}-Sig (PFMETCollection::metSig)", 200,0,20);
	histsIter4.hPfHt_genTrks_atanLt7  = fs->make<TH1F> ("PfHt_genTrks_atanLt7_iter4" ,"General Tracks atan2(iter4/iter0+1)<0.7:  PFH_{T} (PFMETCollection::H_{T})", ht_bins,0,ht_max);
	histsIter4.hTrksPerVtx_genTrks_atanLt7 = fs->make<TH1F> ("NtrksPerVtx_genTrks_atanLt7_iter4" ,"General Tracks atan2(iter4/iter0+1)<0.7:  N Tracks per Vertex", 250,0,500);


	//high purity tracks
	histsIter4.hNtrks_highPureTrks_atanGt7 = fs->make<TH1F> ("Ntrks_highPureTrks_atanGt7_iter4" ,"High Purity Tracks atan2(iter4/iter0+1)>0.7: N Tracks", 250,0,500);
	histsIter4.hNvtx_highPureTrks_atanGt7  = fs->make<TH1F> ("NVtx_highPureTrks_atanGt7_iter4" ,"High Purity Tracks atan2(iter4/iter0+1)>0.7:  N Vertices", 50,0,50);
	histsIter4.hRawMet_highPureTrks_atanGt7  = fs->make<TH1F> ("CaloRawMet_highPureTrks_atanGt7_iter4" ,"High Purity Tracks atan2(iter4/iter0+1)>0.7:  Raw #slash{E}_{T} (CaloMETCollection::met)", met_bins,0,met_max);
	histsIter4.hRawHt_highPureTrks_atanGt7  = fs->make<TH1F> ("CaloRawHt_highPureTrks_atanGt7_iter4" ,"High Purity Tracks atan2(iter4/iter0+1)>0.7:  H_{T} (CaloMETCollection::H_{T})", ht_bins,0,ht_max);
	histsIter4.hPfMet_highPureTrks_atanGt7  = fs->make<TH1F> ("PfMet_highPureTrks_atanGt7_iter4" ,"High Purity Tracks atan2(iter4/iter0+1)>0.7:  PF#slash{E}_{T} (PFMETCollection::pfMet)", met_bins,0,met_max);
	histsIter4.hPfMetSig_highPureTrks_atanGt7  = fs->make<TH1F> ("PfMetSig_highPureTrks_atanGt7_iter4" ,"High Purity Tracks atan2(iter4/iter0+1)>0.7: PF#slash{E}_{T}-Sig (PFMETCollection::metSig)", 200,0,20);
	histsIter4.hPfHt_highPureTrks_atanGt7  = fs->make<TH1F> ("PfHt_highPureTrks_atanGt7_iter4" ,"High Purity Tracks atan2(iter4/iter0+1)>0.7:  PFH_{T} (PFMETCollection::H_{T})", ht_bins,0,ht_max);
	histsIter4.hTrksPerVtx_highPureTrks_atanGt7 = fs->make<TH1F> ("NtrksPerVtx_highPureTrks_atanGt7_iter4" ,"High Purity Tracks atan2(iter4/iter0+1)>0.7:  N Tracks per Vertex", 250,0,500);

	histsIter4.hNtrks_highPureTrks_atanLt7 = fs->make<TH1F> ("Ntrks_highPureTrks_atanLt7_iter4" ,"High Purity Tracks atan2(iter4/iter0+1)<0.7: N Tracks", 250,0,500);
	histsIter4.hNvtx_highPureTrks_atanLt7  = fs->make<TH1F> ("NVtx_highPureTrks_atanLt7_iter4" ,"High Purity Tracks atan2(iter4/iter0+1)<0.7:  N Vertices", 50,0,50);
	histsIter4.hRawMet_highPureTrks_atanLt7  = fs->make<TH1F> ("CaloRawMet_highPureTrks_atanLt7_iter4" ,"High Purity Tracks atan2(iter4/iter0+1)<0.7:  Raw #slash{E}_{T} (CaloMETCollection::met)", met_bins,0,met_max);
	histsIter4.hRawHt_highPureTrks_atanLt7  = fs->make<TH1F> ("CaloRawHt_highPureTrks_atanLt7_iter4" ,"High Purity Tracks atan2(iter4/iter0+1)<0.7:  H_{T} (CaloMETCollection::H_{T})", ht_bins,0,ht_max);
	histsIter4.hPfMet_highPureTrks_atanLt7  = fs->make<TH1F> ("PfMet_highPureTrks_atanLt7_iter4" ,"High Purity Tracks atan2(iter4/iter0+1)<0.7:  PF#slash{E}_{T} (PFMETCollection::pfMet)", met_bins,0,met_max);
	histsIter4.hPfMetSig_highPureTrks_atanLt7  = fs->make<TH1F> ("PfMetSig_highPureTrks_atanLt7_iter4" ,"High Purity Tracks atan2(iter4/iter0+1)<0.7: PF#slash{E}_{T}-Sig (PFMETCollection::metSig)", 200,0,20);
	histsIter4.hPfHt_highPureTrks_atanLt7  = fs->make<TH1F> ("PfHt_highPureTrks_atanLt7_iter4" ,"High Purity Tracks atan2(iter4/iter0+1)<0.7:  PFH_{T} (PFMETCollection::H_{T})", ht_bins,0,ht_max);
	histsIter4.hTrksPerVtx_highPureTrks_atanLt7 = fs->make<TH1F> ("NtrksPerVtx_highPureTrks_atanLt7_iter4" ,"High Purity Tracks atan2(iter4/iter0+1)<0.7:  N Tracks per Vertex", 250,0,500);


	histsIter5.hAtanVsMet = fs->make<TProfile> ("Atan2VsMet_iter5" ,"atan2(iter5/iter0+1) Vs PF#slash{E}_{T}", 200,0,1000);
	//general tracks
	histsIter5.hNtrks_genTrks_atanGt7 = fs->make<TH1F> ("Ntrks_genTrks_atanGt7_iter5" ,"General Tracks atan2(iter5/iter0+1)>0.7: N Tracks", 250,0,500);
	histsIter5.hNvtx_genTrks_atanGt7  = fs->make<TH1F> ("NVtx_genTrks_atanGt7_iter5" ,"General Tracks atan2(iter5/iter0+1)>0.7:  N Vertices", 50,0,50);
	histsIter5.hRawMet_genTrks_atanGt7  = fs->make<TH1F> ("CaloRawMet_genTrks_atanGt7_iter5" ,"General Tracks atan2(iter5/iter0+1)>0.7:  Raw #slash{E}_{T} (CaloMETCollection::met)", met_bins,0,met_max);
	histsIter5.hRawHt_genTrks_atanGt7  = fs->make<TH1F> ("CaloRawHt_genTrks_atanGt7_iter5" ,"General Tracks atan2(iter5/iter0+1)>0.7:  H_{T} (CaloMETCollection::H_{T})", ht_bins,0,ht_max);
	histsIter5.hPfMet_genTrks_atanGt7  = fs->make<TH1F> ("PfMet_genTrks_atanGt7_iter5" ,"General Tracks atan2(iter5/iter0+1)>0.7:  PF#slash{E}_{T} (PFMETCollection::pfMet)", met_bins,0,met_max);
	histsIter5.hPfMetSig_genTrks_atanGt7  = fs->make<TH1F> ("PfMetSig_genTrks_atanGt7_iter5" ,"General Tracks atan2(iter5/iter0+1)>0.7: PF#slash{E}_{T}-Sig (PFMETCollection::metSig)", 200,0,20);
	histsIter5.hPfHt_genTrks_atanGt7  = fs->make<TH1F> ("PfHt_genTrks_atanGt7_iter5" ,"General Tracks atan2(iter5/iter0+1)>0.7:  PFH_{T} (PFMETCollection::H_{T})", ht_bins,0,ht_max);
	histsIter5.hTrksPerVtx_genTrks_atanGt7 = fs->make<TH1F> ("NtrksPerVtx_genTrks_atanGt7_iter5" ,"General Tracks atan2(iter5/iter0+1)>0.7:  N Tracks per Vertex", 250,0,500);

	histsIter5.hNtrks_genTrks_atanLt7 = fs->make<TH1F> ("Ntrks_genTrks_atanLt7_iter5" ,"General Tracks atan2(iter5/iter0+1)<0.7: N Tracks", 250,0,500);
	histsIter5.hNvtx_genTrks_atanLt7  = fs->make<TH1F> ("NVtx_genTrks_atanLt7_iter5" ,"General Tracks atan2(iter5/iter0+1)<0.7:  N Vertices", 50,0,50);
	histsIter5.hRawMet_genTrks_atanLt7  = fs->make<TH1F> ("CaloRawMet_genTrks_atanLt7_iter5" ,"General Tracks atan2(iter5/iter0+1)<0.7:  Raw #slash{E}_{T} (CaloMETCollection::met)", met_bins,0,met_max);
	histsIter5.hRawHt_genTrks_atanLt7  = fs->make<TH1F> ("CaloRawHt_genTrks_atanLt7_iter5" ,"General Tracks atan2(iter5/iter0+1)<0.7:  H_{T} (CaloMETCollection::H_{T})", ht_bins,0,ht_max);
	histsIter5.hPfMet_genTrks_atanLt7  = fs->make<TH1F> ("PfMet_genTrks_atanLt7_iter5" ,"General Tracks atan2(iter5/iter0+1)<0.7:  PF#slash{E}_{T} (PFMETCollection::pfMet)", met_bins,0,met_max);
	histsIter5.hPfMetSig_genTrks_atanLt7  = fs->make<TH1F> ("PfMetSig_genTrks_atanLt7_iter5" ,"General Tracks atan2(iter5/iter0+1)<0.7: PF#slash{E}_{T}-Sig (PFMETCollection::metSig)", 200,0,20);
	histsIter5.hPfHt_genTrks_atanLt7  = fs->make<TH1F> ("PfHt_genTrks_atanLt7_iter5" ,"General Tracks atan2(iter5/iter0+1)<0.7:  PFH_{T} (PFMETCollection::H_{T})", ht_bins,0,ht_max);
	histsIter5.hTrksPerVtx_genTrks_atanLt7 = fs->make<TH1F> ("NtrksPerVtx_genTrks_atanLt7_iter5" ,"General Tracks atan2(iter5/iter0+1)<0.7:  N Tracks per Vertex", 250,0,500);


	//high purity tracks
	histsIter5.hNtrks_highPureTrks_atanGt7 = fs->make<TH1F> ("Ntrks_highPureTrks_atanGt7_iter5" ,"High Purity Tracks atan>(iter5/iter0+1)0.7: N Tracks", 250,0,500);
	histsIter5.hNvtx_highPureTrks_atanGt7  = fs->make<TH1F> ("NVtx_highPureTrks_atanGt7_iter5" ,"High Purity Tracks atan>(iter5/iter0+1)0.7:  N Vertices", 50,0,50);
	histsIter5.hRawMet_highPureTrks_atanGt7  = fs->make<TH1F> ("CaloRawMet_highPureTrks_atanGt7_iter5" ,"High Purity Tracks atan>(iter5/iter0+1)0.7:  Raw #slash{E}_{T} (CaloMETCollection::met)", met_bins,0,met_max);
	histsIter5.hRawHt_highPureTrks_atanGt7  = fs->make<TH1F> ("CaloRawHt_highPureTrks_atanGt7_iter5" ,"High Purity Tracks atan>(iter5/iter0+1)0.7:  H_{T} (CaloMETCollection::H_{T})", ht_bins,0,ht_max);
	histsIter5.hPfMet_highPureTrks_atanGt7  = fs->make<TH1F> ("PfMet_highPureTrks_atanGt7_iter5" ,"High Purity Tracks atan>(iter5/iter0+1)0.7:  PF#slash{E}_{T} (PFMETCollection::pfMet)", met_bins,0,met_max);
	histsIter5.hPfMetSig_highPureTrks_atanGt7  = fs->make<TH1F> ("PfMetSig_highPureTrks_atanGt7_iter5" ,"High Purity Tracks atan2(iter5/iter0+1)>0.7: PF#slash{E}_{T}-Sig (PFMETCollection::metSig)", 200,0,20);
	histsIter5.hPfHt_highPureTrks_atanGt7  = fs->make<TH1F> ("PfHt_highPureTrks_atanGt7_iter5" ,"High Purity Tracks atan>(iter5/iter0+1)0.7:  PFH_{T} (PFMETCollection::H_{T})", ht_bins,0,ht_max);
	histsIter5.hTrksPerVtx_highPureTrks_atanGt7 = fs->make<TH1F> ("NtrksPerVtx_highPureTrks_atanGt7_iter5" ,"High Purity Tracks atan>(iter5/iter0+1)0.7:  N Tracks per Vertex", 250,0,500);

	histsIter5.hNtrks_highPureTrks_atanLt7 = fs->make<TH1F> ("Ntrks_highPureTrks_atanLt7_iter5" ,"High Purity Tracks atan<(iter5/iter0+1)0.7: N Tracks", 250,0,500);
	histsIter5.hNvtx_highPureTrks_atanLt7  = fs->make<TH1F> ("NVtx_highPureTrks_atanLt7_iter5" ,"High Purity Tracks atan<(iter5/iter0+1)0.7:  N Vertices", 50,0,50);
	histsIter5.hRawMet_highPureTrks_atanLt7  = fs->make<TH1F> ("CaloRawMet_highPureTrks_atanLt7_iter5" ,"High Purity Tracks atan<(iter5/iter0+1)0.7:  Raw #slash{E}_{T} (CaloMETCollection::met)", met_bins,0,met_max);
	histsIter5.hRawHt_highPureTrks_atanLt7  = fs->make<TH1F> ("CaloRawHt_highPureTrks_atanLt7_iter5" ,"High Purity Tracks atan<(iter5/iter0+1)0.7:  H_{T} (CaloMETCollection::H_{T})", ht_bins,0,ht_max);
	histsIter5.hPfMet_highPureTrks_atanLt7  = fs->make<TH1F> ("PfMet_highPureTrks_atanLt7_iter5" ,"High Purity Tracks atan<(iter5/iter0+1)0.7:  PF#slash{E}_{T} (PFMETCollection::pfMet)", met_bins,0,met_max);
	histsIter5.hPfMetSig_highPureTrks_atanLt7  = fs->make<TH1F> ("PfMetSig_highPureTrks_atanLt7_iter5" ,"High Purity Tracks atan2(iter5/iter0+1)<0.7: PF#slash{E}_{T}-Sig (PFMETCollection::metSig)", 200,0,20);
	histsIter5.hPfHt_highPureTrks_atanLt7  = fs->make<TH1F> ("PfHt_highPureTrks_atanLt7_iter5" ,"High Purity Tracks atan<(iter5/iter0+1)0.7:  PFH_{T} (PFMETCollection::H_{T})", ht_bins,0,ht_max);
	histsIter5.hTrksPerVtx_highPureTrks_atanLt7 = fs->make<TH1F> ("NtrksPerVtx_highPureTrks_atanLt7_iter5" ,"High Purity Tracks atan<(iter5/iter0+1)0.7:  N Tracks per Vertex", 250,0,500);


	//pixelless stuff from pure tracks
	histsPixless.hAtanVsMet = fs->make<TProfile> ("Atan2VsMet_Pixless" ,"atan2(Pixelless/iter0+1) Vs PF#slash{E}_{T}", met_bins,0,met_max);
	histsPixless.hNtrks_highPureTrks_atanLt7 = fs->make<TH1F> ("Ntrks_highPureTrks_atanLt7_Pixless" ,"High Purity Tracks atan2(Pixless/iter0+1)<0.7: N Tracks", 250,0,500);
	histsPixless.hNvtx_highPureTrks_atanLt7  = fs->make<TH1F> ("NVtx_highPureTrks_atanLt7_Pixless" ,"High Purity Tracks atan2(Pixless/iter0+1)<0.7:  N Vertices", 50,0,50);
	histsPixless.hPfMet_highPureTrks_atanLt7  = fs->make<TH1F> ("PfMet_highPureTrks_atanLt7_Pixless" ,"High Purity Tracks atan2(Pixless/iter0+1)<0.7:  PF#slash{E}_{T} (PFMETCollection::pfMet)", met_bins,0,met_max);
	histsPixless.hPfMetSig_highPureTrks_atanLt7  = fs->make<TH1F> ("PfMetSig_highPureTrks_atanLt7_Pixless" ,"High Purity Tracks atan2(Pixless/iter0+1)<0.7: PF#slash{E}_{T}-Sig (PFMETCollection::metSig)", 200,0,20);
	histsPixless.hPfHt_highPureTrks_atanLt7  = fs->make<TH1F> ("PfHt_highPureTrks_atanLt7_Pixless" ,"High Purity Tracks atan2(Pixless/iter0+1)<0.7:  PFH_{T} (PFMETCollection::H_{T})", ht_bins,0,ht_max);






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
	bool bSaveEvent = false;

	RunNumber_t kRun   = iEvent.id().run();
	EventNumber_t kEvent = iEvent.id().event();
	std::pair<edm::RunNumber_t, edm::EventNumber_t> kRunEvent (kRun, kEvent);

	//========== Trigger Selections
	if (req_trigger)
	{
		if (! storeHLT(iEvent, iSetup)) return 0;
	}


 // ==========================================================
 //Vertex information

	Handle<reco::VertexCollection> vertexHandle;
  	iEvent.getByLabel("offlinePrimaryVertices", vertexHandle);
	double dNvtx = 0;
   if (vertexHandle.isValid())
	{
     reco::VertexCollection vertexCollection = *(vertexHandle.product());
     //std::cout << "VtxSize = "<<  vertexCollection.size() << std::endl;
     reco::VertexCollection::const_iterator v = vertexCollection.begin();
	  if (v->isValid())
	  {
		  for (; v != vertexCollection.end(); ++v)
		  {
			  if (v->isFake()) continue;
			  ++dNvtx;
		  }
	  }
	} else 
	{
		std::cout << "AnomTrkAna: Could not find vertex collection" << std::endl;
		assert(false);	
   }
   
	//std::cout << "nvtx " << dNvtx << std::endl;
	if (dNvtx < (double)inMinVtx) return 0;  // require a primary vertex 
   
 	// ==========================================================
	// MET Information 
	Handle<reco::PFMETCollection> pfMet;
	iEvent.getByLabel("pfMet",pfMet);
	double dPFMet = 0.0, dPFHt = 0. , dPFMetSig = 0.;
	if (pfMet.isValid())
	{
		for (reco::PFMETCollection::const_iterator it = pfMet->begin(); it != pfMet->end(); ++it)
		{
			dPFMet = it->pt();
			dPFHt = it->sumEt();
			dPFMetSig = it->mEtSig();
			if (dPFMet<dminMet) return 0;  //minimum met requirement

		}
	} else
	{
		std::cout << "Valid PFMETCollection Not Found!" << std::endl;
		assert(false);
	}


	Handle<reco::TrackCollection> tracks;
	iEvent.getByLabel("generalTracks", tracks);

	if (tracks.isValid())
	{
		//histsIter4.hNtrks->Fill(tracks->size());
		double iNiter01trks = 0, iNiter4trks = 0, iNiter5trks = 0, iNpxlLess = 0;
		double iNiter01trkspure = 0, iNiter4trkspure = 0, iNiter5trkspure = 0, iNpxlLesspure = 0;

		for (reco::TrackCollection::const_iterator it = tracks->begin(); it != tracks->end(); ++it)
		{
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

		//std::cout << "Ratio = Nit4/Nit0+1 = " 
		//	<< dRatIter4
		//	<< std::endl; 
		//std::cout << "Ratio = Nit5/Nit0+1 = " 
		//	<< dRatIter5 
		//	<< std::endl; 
		hTrkIter0.hRatioIt4toIt01->Fill(dRatIter4);
		hTrkIter0.hRatioIt5toIt01->Fill(dRatIter5);
		hTrkIter0.hRatioPixlesstoIt01->Fill(dRatPixless);

		hTrkIter0.hRatioIt4toIt01pure->Fill(dRatIter4pure);
		hTrkIter0.hRatioIt5toIt01pure->Fill(dRatIter5pure);
		hTrkIter0.hRatioPixlesstoIt01pure->Fill(dRatPixlesspure);


		//make MET plots for different atan regions


		histsIter4.hPfMet->Fill(dPFMet);
		histsIter4.hPfHt->Fill(dPFHt);
		histsIter4.hPfMetSig->Fill(dPFMetSig);
		histsIter4.hAtanVsMet->Fill(dPFMet, dRatIter4);
		//high purity tracks: atan(Iter4/iter0+1)>0.7
		if (dRatIter4pure > 0.7)
		{
			histsIter4.hPfMet_highPureTrks_atanGt7->Fill(dPFMet);
			std::cout << "dRatIter4pure == " << dRatIter4pure << "(>0.7) event:: " << kRun << ", " << kEvent << std::endl;
			it4Rej.insert(kRunEvent);
			bSaveEvent = true;
		} else 
		{
			histsIter4.hNvtx_highPureTrks_atanLt7->Fill(dNvtx);
			//histsIter4.hRawHt_highPureTrks_atanLt7->Fill(dRawHt);
			//histsIter4.hRawMet_highPureTrks_atanLt7->Fill(dRawMet);
			histsIter4.hPfHt_highPureTrks_atanLt7->Fill(dPFHt);
			histsIter4.hPfMet_highPureTrks_atanLt7->Fill(dPFMet);
			histsIter4.hPfMetSig_highPureTrks_atanLt7->Fill(dPFMetSig);
		}


		histsIter5.hAtanVsMet->Fill(dPFMet, dRatIter5);
		//high purity tracks: atan(Iter5/iter0+1)>0.7
		if (dRatIter5pure > 0.7)
		{
			histsIter5.hPfMet_highPureTrks_atanGt7->Fill(dPFMet);
			std::cout << "dRatIter5pure == " << dRatIter5pure << "(>0.7) event:: " << kRun << ", " << kEvent << std::endl;
			it5Rej.insert(kRunEvent);
			bSaveEvent = true;
		} else 
		{
			histsIter5.hNvtx_highPureTrks_atanLt7->Fill(dNvtx);
			//histsIter5.hRawHt_highPureTrks_atanLt7->Fill(dRawHt);
			//histsIter5.hRawMet_highPureTrks_atanLt7->Fill(dRawMet);
			histsIter5.hPfHt_highPureTrks_atanLt7->Fill(dPFHt);
			histsIter5.hPfMet_highPureTrks_atanLt7->Fill(dPFMet);
			histsIter5.hPfMetSig_highPureTrks_atanLt7->Fill(dPFMetSig);
		}
	
		//pixless stuff
		histsPixless.hAtanVsMet->Fill(dPFMet, dRatPixlesspure);
		//high purity tracks: atan(PIXLESS/iter0+1)>0.7
		if (dRatPixlesspure > 0.7)
		{
			std::cout << "dRatPixless == " << dRatPixlesspure << "(>0.7) event:: " << kRun << ", " << kEvent << std::endl;
			histsPixless.hPfMet_highPureTrks_atanLt7->Fill(dPFMet);
			pxlessRej.insert(kRunEvent);
			bSaveEvent = true;
		} else 
		{
			histsPixless.hNvtx_highPureTrks_atanLt7->Fill(dNvtx);
			histsPixless.hPfHt_highPureTrks_atanLt7->Fill(dPFHt);
			histsPixless.hPfMet_highPureTrks_atanLt7->Fill(dPFMet);
			histsPixless.hPfMetSig_highPureTrks_atanLt7->Fill(dPFMetSig);
		}

		assert (histsIter4.hPfMet->GetEntries() == 
				(histsIter4.hPfMet_highPureTrks_atanLt7->GetEntries()
				 + histsIter4.hPfMet_highPureTrks_atanGt7->GetEntries())
					&& "PFMET entries do not match PFMET_HP_atanLT7!");
	
	} else
	{
		std::cout << "Valid Track Collection Not Found!" << std::endl;
		assert(false);
	}
	//std::cout << "std::cout>> number of tracks = " << tracks->size() << std::endl;

	if (bSaveEvent) ++iPassed;
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
}

// ------------ method called once each job just after ending the event loop  ------------
void 
AnomTrkMETFilter::endJob() {

	
	std::cout << __FILE__ << "::" << __FUNCTION__  << "\n" << std::endl;
	std::cout << "[ATM:00] Job settings " << std::endl;
	std::cout << "[ATM:01] Events Processed --- = " << iProcessed << std::endl;
	std::cout << "[ATM:01] Events Passed ------ = " << iPassed << std::endl;
	std::cout << "[ATM:10] Primary Vtx Min ---- = " << inMinVtx << std::endl; 
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
	std::cout << histsIter4.hPfMet->GetName() << ": Entries = " 
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

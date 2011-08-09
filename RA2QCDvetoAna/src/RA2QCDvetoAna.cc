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
		double dminMet;
		unsigned int iProcessed; // number of processed events
		unsigned int iPassed; //number of events passed the filter
		int inMinNdofVtx; //minimum number of dof for the vtx
		RunLumiEvt_t kRunLumiEvent;
		int iVerbose; // control print levels
		double dMaxPrimVtxZ;

		struct EventHist_t {
			TH1F* nvtx;
			TH1F* vtxz;
			TH1F* mht;
			TH1F* ht;
			TH1F* njet;
			TH1F* metphi;
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


		//TH1F* hist_delphiMin_jetmet[5];
		//TH2F* hist_delphiMin_jetmetVsMHT[5]; //0=all jets, 1=only 3 jets, 2=upto 4 jets 
		//TH1F* hist_delphiMin_jetmetVsHT[5];

		EventHist_t evtHist;
		JetHist_t jet1Hist, jet2Hist, jet3Hist; //for upto 5 jets may be

		edm::InputTag patJetsPFPt30InputTag_;
		edm::InputTag patJetsPFPt50Eta25InputTag_;
		edm::Handle<std::vector<pat::Jet> > jetHandle;
		edm::Handle<std::vector<pat::Jet> > pfpt50eta25JetHandle;

		Hist_t Hist[8]; //all hists for various mht ranges
		TH1F* MHT_by_phislice[6];
		void BookHistograms(edm::Service<TFileService>& fs, Hist_t& hist
						, const float mHt_min, const float mHt_max);
		void FillHistograms(Hist_t& hist
					, edm::Handle<edm::View<reco::MET> > mhtHandle 
				   , edm::Handle<std::vector<pat::Jet> > jetHandle);
		void minDelPhi(std::vector<float>& vDelPhi_jetjet
				, std::vector<float>& vDelPhi_jetmet
				, const unsigned njets
				, edm::Handle<std::vector<pat::Jet> > jetHandle
				, edm::Handle<edm::View<reco::MET> > mhtHandle); 
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
	//triggerPathsToStore_ = iConfig.getParameter<std::vector<std::string> >("TriggerPathsToStore");
	//hlTriggerResults_    = iConfig.getParameter<edm::InputTag>("HltTriggerResults");	
	//req_trigger = iConfig.getUntrackedParameter<bool>("req_trigger",true);
	patJetsPFPt30InputTag_ = iConfig.getParameter<edm::InputTag>("patJetsPFPt30InputTag");
	patJetsPFPt50Eta25InputTag_ = iConfig.getParameter<edm::InputTag>("patJetsPFPt50Eta25InputTag");
	mhtInputTag_ = iConfig.getParameter<edm::InputTag>("mhtInputTag");
	htInputTag_ = iConfig.getParameter<edm::InputTag>("htInputTag");
	inMinVtx = iConfig.getUntrackedParameter<int>("nMinVtx",1);
	dminMet = iConfig.getUntrackedParameter<double>("minMet",0.0);
	inMinNdofVtx = iConfig.getUntrackedParameter<int>("ndofVtx",0);
	dMaxPrimVtxZ = iConfig.getUntrackedParameter<double>("maxDelzVtx",24.0);
	//caloJetInputTag_ = iConfig.getParameter<edm::InputTag>("caloJetInputTag_");
	//pfJetInputTag_ = iConfig.getParameter<edm::InputTag>("pfJetInputTag_");
	iVerbose = iConfig.getUntrackedParameter<int>("verbose",0);
	iProcessed = 0;
	iPassed = 0;

	//generate hists
	edm::Service<TFileService> fs;
	BookHistograms(fs, Hist[0], 60, 80); 
	BookHistograms(fs, Hist[1], 80,100); 
	BookHistograms(fs, Hist[2],100,120); 
	BookHistograms(fs, Hist[3],120,140); 
	BookHistograms(fs, Hist[4],140,170); 
	BookHistograms(fs, Hist[5],170,200); 
	BookHistograms(fs, Hist[6],200,7000); 

	const float met_min = 0, met_max=500, met_bins=100;
	MHT_by_phislice[0] = fs->make<TH1F> ("mht_phislice_lt0.1" ,"MHT (#Delta#Phi_{min}<0.1);MHT [GeV];Events;", met_bins, 0, met_max);
	MHT_by_phislice[1] = fs->make<TH1F> ("mht_phislice_lt0.2" ,"MHT (0.1<#Delta#Phi_{min}<0.2);MHT [GeV];Events;", met_bins, 0, met_max);
	MHT_by_phislice[2] = fs->make<TH1F> ("mht_phislice_lt0.3" ,"MHT (0.2<#Delta#Phi_{min}<0.3);MHT [GeV];Events;", met_bins, 0, met_max);
	MHT_by_phislice[3] = fs->make<TH1F> ("mht_phislice_lt0.5" ,"MHT (0.3<#Delta#Phi_{min}<0.5);MHT [GeV];Events;", met_bins, 0, met_max);
	MHT_by_phislice[4] = fs->make<TH1F> ("mht_phislice_lt0.8" ,"MHT (0.5<#Delta#Phi_{min}<0.8);MHT [GeV];Events;", met_bins, 0, met_max);
	MHT_by_phislice[5] = fs->make<TH1F> ("mht_phislice_gt0.8" ,"MHT (#Delta#Phi_{min}>0.8);MHT [GeV];Events;", met_bins, 0, met_max);

/*	const double met_max = 4000, met_bins = 400;
	evtHist.nvtx = fs->make<TH1F> ("nvtx" ,"Nvtx;Nvtx;Events;", 15, 0, 15);
	evtHist.vtxz = fs->make<TH1F> ("vtxz","Primary Vertex z;z [cm];Arbitrary;",30,0,30);
	evtHist.mht = fs->make<TH1F> ("mht" ,"MHT;MHT [GeV];Events;", met_bins, 0, met_max);
	evtHist.ht = fs->make<TH1F> ("ht" ,"HT;HT [GeV];Events;", met_bins, 0, met_max);
	evtHist.njet = fs->make<TH1F> ("njet" ,"NJET;NJET [GeV];Events;", 10, 0, 10);
	evtHist.metphi = fs->make<TH1F> ("metphi" ,"met phi;met phi;Events;", 160, -8, 8);

	const double pt_bins = 200, pt_max = 2000;
	jet1Hist.pt = fs->make<TH1F> ("jet1_pt" ,"Jet1 pt;pt [GeV];Events;", pt_bins, 0, pt_max);
	jet1Hist.eta = fs->make<TH1F> ("jet1_eta" ,"Jet1 eta;eta;Events;", 100, -5, 5);
	jet1Hist.phi = fs->make<TH1F> ("jet1_phi" ,"Jet1 phi;phi;Events;", 160, -8, 8);
	jet1Hist.delphi = fs->make<TH1F> ("jet1_delphi" ,"Jet1: delphi;delphi;Events;", 160, -8, 8);

	jet2Hist.pt = fs->make<TH1F> ("jet2_pt" ,"Jet2 pt;pt [GeV];Events;", pt_bins, 0, pt_max);
	jet2Hist.eta = fs->make<TH1F> ("jet2_eta" ,"Jet2 eta;eta;Events;", 100, -5, 5);
	jet2Hist.phi = fs->make<TH1F> ("jet2_phi" ,"Jet2 phi;phi;Events;", 160, -8, 8);
	jet2Hist.delphi = fs->make<TH1F> ("jet2_delphi" ,"Jet2: delphi;delphi;Events;", 160, -8, 8);

	jet3Hist.pt = fs->make<TH1F> ("jet3_pt" ,"Jet3 pt;pt [GeV];Events;", pt_bins, 0, pt_max);
	jet3Hist.eta = fs->make<TH1F> ("jet3_eta" ,"Jet3 eta;eta;Events;", 100, -5, 5);
	jet3Hist.phi = fs->make<TH1F> ("jet3_phi" ,"Jet3 phi;phi;Events;", 160, -8, 8);
	jet3Hist.delphi = fs->make<TH1F> ("jet3_delphi" ,"Jet3: delphi;delphi;Events;", 160, -8, 8);

	hist_delphiMin_jetmet[0] = fs->make<TH1F> ("delphiMin_jetmet" ,"delphi Min.;delphi min;Events;", 160, -8, 8);
	hist_delphiMin_jetmetVsMHT[0] = fs->make<TH2F> ("delphiMin_jetmetVsMHT" ,"delphi min Vs MHT;MHT;delphi min;", 1000, 0, 1000,400,0,4);
*/
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
	bool bSaveEvent = false;
	kRunLumiEvent.run  = iEvent.id().run();
	kRunLumiEvent.lumi = iEvent.id().luminosityBlock(); 
	kRunLumiEvent.evt  = iEvent.id().event();


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
			  ++dNvtx;
			  if (v->isFake()) continue;
			  double dPrimVtx_z = v->z();
			  if (fabs(dPrimVtx_z) > dMaxPrimVtxZ) continue;  // do not use the event. 
			  if (v->ndof() < inMinNdofVtx) continue; 
			  if (iVerbose) std::cout << "Valid Vertex Found " << std::endl;
			  //evtHist.vtxz->Fill(dPrimVtx_z);
		  }
	  }
	}

	//evtHist.nvtx->Fill(dNvtx);
	
 	// ==========================================================
	// MET Information 
	iEvent.getByLabel(mhtInputTag_, mhtHandle);
	if (! mhtHandle.isValid()) 
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":MHT handle not found!" << std::endl;
		assert(false);
	}
	const float mht = (*mhtHandle)[0].pt();

	iEvent.getByLabel(patJetsPFPt30InputTag_, jetHandle);
	if (! jetHandle.isValid()) 
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":JET handle not found!" << std::endl;
		assert(false);
	}

	if (mht>60 && mht<=80) FillHistograms(Hist[0], mhtHandle, jetHandle); 
	else if (mht>80 && mht<=100) FillHistograms(Hist[1], mhtHandle, jetHandle); 
	else if (mht>100&& mht<=120) FillHistograms(Hist[2], mhtHandle, jetHandle); 
	else if (mht>120&& mht<=140) FillHistograms(Hist[3], mhtHandle, jetHandle); 
	else if (mht>140&& mht<=170) FillHistograms(Hist[4], mhtHandle, jetHandle); 
	else if (mht>170&& mht<=200) FillHistograms(Hist[5], mhtHandle, jetHandle); 
	else if (mht>200) FillHistograms(Hist[6], mhtHandle, jetHandle); 
	
	//mhtphi->Fill( (*mht)[0].phi());
//	evtHist.metphi->Fill((*mhtHandle)[0].phi());
//	evtHist.mht->Fill( (*mhtHandle)[0].pt());
//	evtHist.ht->Fill( (*htHandle) );

	iEvent.getByLabel(patJetsPFPt50Eta25InputTag_, pfpt50eta25JetHandle);
	if (! pfpt50eta25JetHandle.isValid()) 
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":pfpt50eta25JetHandle handle not found!" << std::endl;
		assert(false);
	}
	unsigned numpt50eta25Jets = pfpt50eta25JetHandle->size();
	std::cout << "pfpt50eta25JetHandle numjets = " << numpt50eta25Jets << std::endl;
	TLorentzVector oJetsVec(0,0,0,0);
	for (unsigned i = 0 ; i < numpt50eta25Jets ; ++i)
	{
		//if ((*numpt50eta25Jets)[i].pt() > minRa2JetPt && (*numpt50eta25Jets)[i].eta()<2.5) ++Njet50;  
		//if ((*numpt50eta25Jets)[i].pt() > minJetPt) ++Njet30;  
			const TLorentzVector iJetVec((*pfpt50eta25JetHandle)[i].px(),(*pfpt50eta25JetHandle)[i].py(),(*pfpt50eta25JetHandle)[i].pz(),(*pfpt50eta25JetHandle)[i].energy());
			oJetsVec += iJetVec;
	}

	iEvent.getByLabel(htInputTag_, htHandle);

	std::cout << "mht/my mht = " << (*mhtHandle)[0].pt() << "/" << oJetsVec.Pt() << std::endl;

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
	std::cout << __FILE__ << "::" << __FUNCTION__  << "\n" << std::endl;
	std::cout << "[ATM:00] Job settings " << std::endl;
	std::cout << "[ATM:01] Events Processed --- = " << iProcessed << std::endl;
	std::cout << "[ATM:02] Events Passed ------ = " << iPassed << std::endl;
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
	std::cout << "[ATM:30] Minimum MET -------- = " << dminMet << std::endl;

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
void
RA2QCDvetoAna::BookHistograms(edm::Service<TFileService>& fs, Hist_t& hist , const float mHt_min, const float mHt_max)
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
RA2QCDvetoAna::FillHistograms(Hist_t& hist
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

//	evtHist.njet->Fill(numJets);
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

void RA2QCDvetoAna::minDelPhi(std::vector<float>& vDelPhi_jetjet
				, std::vector<float>& vDelPhi_jetmet
				, const unsigned njets
				, edm::Handle<std::vector<pat::Jet> > jetHandle
				, edm::Handle<edm::View<reco::MET> > mhtHandle 
				)
{
	unsigned numJets = jetHandle->size();
	for (unsigned i = 0 ; i < numJets ; ++i)
	{
		const float delphi_jetmet = fabs(TVector2::Phi_mpi_pi((*jetHandle)[i].phi() - (*mhtHandle)[0].phi()));
		vDelPhi_jetmet.push_back(delphi_jetmet);

		TLorentzVector oJetsVec(0,0,0,0);
		for (unsigned j = 0 ; j < numJets ; ++j)
		{
			if ( i == j ) continue;
			//std::cout << " other jets [" << j << "][" << (*jetHandle)[j].pt() << "][" << (*jetHandle)[j].phi() << std::endl;  
			const TLorentzVector jJetVec((*jetHandle)[j].px(),(*jetHandle)[j].py(),(*jetHandle)[j].pz(),(*jetHandle)[j].energy());
			oJetsVec += jJetVec;
			if (j==njets) break; //uses only njets for the caculation
		}
		const float delphi = fabs(TVector2::Phi_mpi_pi((*jetHandle)[i].phi() - oJetsVec.Phi()));
		//std::cout << "delphi min for ijet = " << i << "[" << delphi << "]" << std::endl;
		vDelPhi_jetjet.push_back(delphi);
		if (i==njets) break;  //uses only njets for the caculation
	}

	std::sort(vDelPhi_jetjet.begin(), vDelPhi_jetjet.end(), sort_using_less_than);	
	std::sort(vDelPhi_jetmet.begin(), vDelPhi_jetmet.end(), sort_using_less_than);	
	
	//if (iVerbose)
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

}
//define this as a plug-in
DEFINE_FWK_MODULE(RA2QCDvetoAna);
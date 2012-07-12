// -*- C++ -*-
//
// Package:    MetSigForQcd
// Class:      MetSigForQcd
// 
/**\class MetSigForQcd MetSigForQcd.cc UserCode/MetSigForQcd/src/MetSigForQcd.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  samantha hewamanage
//         Created:  Thu Mar 15 13:34:25 CDT 2012
// $Id: MetSigForQcd.cc,v 1.2 2012/05/02 18:21:10 samantha Exp $
//
//


// system include files
#include <memory>
#include <sstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TLorentzVector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/Math/interface/deltaPhi.h"
//jet collections
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "TH2F.h"
#include "TProfile.h"
//
// class declaration
//

using namespace std;
const static float met_ndof = 2;

class MetSigForQcd : public edm::EDAnalyzer {
   public:
      explicit MetSigForQcd(const edm::ParameterSet&);
      ~MetSigForQcd();

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
			TH1F* ht;
			TH1F* njet30;
			TH1F* njet50;
			TH1F* meff;
			TH1F* metsig;
			TH1F* metprob;
			TH2F* metsigVsSumEt;
			TH2F* metprobVsSumEt;
			TProfile* metsigVsSumEt_prof;
			TProfile* metprobVsSumEt_prof;
			TProfile* metsigVsNjet30_prof;
			TProfile* metprobVsNjet30_prof;
			TProfile* metsigVsNjet50_prof;
			TProfile* metprobVsNjet50_prof;
			TProfile* metsigVsHT;
			TProfile* metprobVsHT;
			TProfile* metsigVsMHT;
			TProfile* metprobVsMHT;
			TH1F* metsigerr;
			TH1F* prescaleWeights;
			TH1F* eventWeights;
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
		};

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
		void PrintHeader();
		void BookHistograms();
		void BookCommonHistograms(TFileDirectory& dir, const float htMin, 
								const float htMax, Hist_t& hist);
		void FillHistograms(
				edm::Handle<reco::PFMETCollection> pfmet,
				edm::Handle<reco::VertexCollection> vertexHandle,
				edm::Handle<std::vector<pat::Jet> > pf30eta50JetHandle,
				edm::Handle<std::vector<pat::Jet> > pf50eta25JetHandle,
				const float ht, const float& mht, const float mhtphi,
				const float prescaleWeight,
				const float Weight,
				Hist_t& hist
				);

      // ----------member data ---------------------------
		//jet collections
		edm::InputTag mhtInputTag_, htInputTag_;
		edm::Handle<edm::View<reco::MET> > mhtHandle;
		edm::Handle<double> htHandle;
		edm::InputTag prescaleWeightInputTag;

		RunLumiEvt_t kRunLumiEvent;
		int iVerbose; // control print levels
		double dMinHT_, dMinMHT_, dMinNjet_, dMaxNjet_, dMinMetSig_;
		edm::InputTag patJetsPFPt50Eta25InputTag_, patJetsPFPt30Eta50InputTag_;
		edm::Handle<std::vector<pat::Jet> > pfpt50eta25JetHandle, pfpt30eta50JetHandle;
		unsigned uFailMinHTCut, uFailMinPFMHTCut, uFailNjetCut, uFailDphiCut, uFailMetSigCut;

		std::vector<float> vMCNvtxDist, vDATANvtxDist;
		double sumLumiWeights; //sum of lumi weights for cross checks
		double Weight; //total weight for an event
		int doLumiWeighing, doEventWeighing;
		edm::Handle<double> prescaleWeightHandle;
		bool usePrescaleWeight;
		bool bApplyNjetCut, bApplyHtCut, bApplyMhtCut, bApplyDphiCut;

		edm::Handle<std::vector<reco::PFMET> >pfMetHandle;
    	std::vector<double> htBins_;
		std::vector<Hist_t> vHist;
		float myHt, myMht;
		TH1F* processed;
		TH1F* passed;
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
MetSigForQcd::MetSigForQcd(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
	patJetsPFPt50Eta25InputTag_ = iConfig.getParameter<edm::InputTag>("patJetsPFPt50Eta25InputTag");
	patJetsPFPt30Eta50InputTag_ = iConfig.getParameter<edm::InputTag>("patJetsPFPt30Eta50InputTag");
	mhtInputTag_    = iConfig.getParameter<edm::InputTag>("mhtInputTag");
	htInputTag_     = iConfig.getParameter<edm::InputTag>("htInputTag");
	dMinNjet_       = iConfig.getUntrackedParameter<double>("dMinNjet", 3.0);
	dMaxNjet_       = iConfig.getUntrackedParameter<double>("dMaxNjet", 3.0);
	dMinHT_         = iConfig.getUntrackedParameter<double>("dMinHT",0.0);
	dMinMHT_        = iConfig.getUntrackedParameter<double>("dMinMHT",0.0);
	iVerbose        = iConfig.getUntrackedParameter<int>("verbose",0);
	doLumiWeighing  = iConfig.getUntrackedParameter<int>("ApplyLumiWeighing",0);
	doEventWeighing = iConfig.getUntrackedParameter<int>("ApplyEventWeighing",0);
	prescaleWeightInputTag = iConfig.getParameter<edm::InputTag>("prescaleWeight");
	usePrescaleWeight = iConfig.getUntrackedParameter<int>("usePrescaleWeight",0);
	sumLumiWeights   = 0;
	Weight           = 1;
	htBins_ = iConfig.getParameter<std::vector<double > >("htBins");
	bApplyNjetCut = iConfig.getUntrackedParameter<int>("applyNjetCut",1);
	bApplyHtCut   = iConfig.getUntrackedParameter<int>("applyHtCut",1);
	bApplyMhtCut  = iConfig.getUntrackedParameter<int>("applyMhtCut",1);
	bApplyDphiCut = iConfig.getUntrackedParameter<int>("applyDphiCut",1);
	dMinMetSig_   = iConfig.getUntrackedParameter<double>("dMinMetSig",0.0);

	BookHistograms();
	uFailNjetCut = 0; uFailMinHTCut = 0; uFailMinPFMHTCut = 0; uFailDphiCut = 0;
	uFailMetSigCut = 0;
}


MetSigForQcd::~MetSigForQcd()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MetSigForQcd::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
	processed->Fill(1);
	kRunLumiEvent.run  = iEvent.id().run();
	kRunLumiEvent.lumi = iEvent.id().luminosityBlock(); 
	kRunLumiEvent.evt  = iEvent.id().event();

	if (iVerbose) std::cout << __LINE__<< ":: Processing event: "
			<< kRunLumiEvent.run << ":" << kRunLumiEvent.lumi 
			<< ":" << kRunLumiEvent.evt << std::endl;

	//vertex cut is applied in the RA2Cleanning process.
	Handle<reco::VertexCollection> vertexHandle;
	//iEvent.getByLabel("goodVerticesRA2", vertexHandle);
	iEvent.getByLabel("offlinePrimaryVertices", vertexHandle);
   if (! vertexHandle.isValid())
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":vertexHandle handle not found!" << std::endl;
		assert(false);
	}

	Weight = 1;

	/* PU S3 reweighting for QCD MC sample
	*/
	double lumiWeight = 1;
/*	if ( doLumiWeighing )
	{
		lumiWeight = LumiWeights_.weight( iEvent );
		//std::cout << "lum wgt = " << lumiWeight << std::endl;
		sumLumiWeights += lumiWeight;
		Weight *= lumiWeight;
	}
*/
	//event weights for flat QCD samples
	double storedWeight = 1;
	if ( doEventWeighing )
	{
		edm::Handle<GenEventInfoProduct> genEvtInfoHandle;
		iEvent.getByLabel("generator", genEvtInfoHandle);
		storedWeight = genEvtInfoHandle->weight();
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
	}


	/*  Get all the handles needed and check their
	 *  validity.
	 */
	iEvent.getByLabel(htInputTag_, htHandle);
	if (! htHandle.isValid()) 
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":htHandle handle not found!" << std::endl;
		assert(false);
	}


	iEvent.getByLabel(patJetsPFPt50Eta25InputTag_, pfpt50eta25JetHandle);
	if (! pfpt50eta25JetHandle.isValid()) 
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":pfpt50eta25JetHandle handle not found!" << std::endl;
		assert(false);
	}

	iEvent.getByLabel(patJetsPFPt30Eta50InputTag_, pfpt30eta50JetHandle);
	if (! pfpt30eta50JetHandle.isValid()) 
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":pfpt30eta50JetHandle handle not found!" << std::endl;
		assert(false);
	}

	/////// PfMET information /////
	edm::Handle<reco::PFMETCollection> pfmet;
	iEvent.getByLabel("pfMet", pfmet);
	if (pfmet->size() == 0) {
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":pfMET handle not found!" << std::endl;
		assert(false);
	}


	iEvent.getByLabel(mhtInputTag_, mhtHandle);
	if (! mhtHandle.isValid()) 
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":MHT handle not found!" << std::endl;
		assert(false);
	}


	//recalc ht/mht
	float htTemp = 0;
	for (unsigned i=0; i < pfpt50eta25JetHandle->size(); ++i) 
	{ 
		htTemp += (*pfpt50eta25JetHandle)[i].pt(); 
	}

	const float myMht    = (*mhtHandle)[0].pt();
	const float myMhtPhi = (*mhtHandle)[0].phi();
	const float myHt     = htTemp;
	//const float myHt     = (*htHandle)[0];
	const unsigned myNjet= pfpt50eta25JetHandle->size(); 
	const float myMetSig = (*pfmet)[0].significance();

	//Check if RA2 dphi cuts are satisfied
	//dphi cuts for each of the jets, starting from highest pt jet
	vector<float> dPhiCut;
	dPhiCut.push_back(0.5); //1st jet
	dPhiCut.push_back(0.5); //2nd jet
	dPhiCut.push_back(0.3); //rest of the jets

	bool bPassDphiCut = true;

	for (unsigned i=0; i < pfpt30eta50JetHandle->size(); ++i) 
	{ 
		//this logic will work even in the case of exclusive njets
		//bcos nminjets == nmaxjets in such cases
		if ( i < (unsigned) dMinNjet_ )  
		{
			float dPhi = 0.3; //for rest of the jets.
			if (i < dPhiCut.size()) dPhi = dPhiCut.at(i);
			
			bPassDphiCut = bPassDphiCut && (std::abs(reco::deltaPhi((*pfpt30eta50JetHandle)[i].phi(), myMhtPhi)) > dPhi);

		} else break;
	}


	//APPLY cuts

	if ( bApplyNjetCut && (myNjet < (unsigned) dMinNjet_ || myNjet > (unsigned) dMaxNjet_ ) ) { ++uFailNjetCut;  return; }
	if ( bApplyHtCut   && myHt < dMinHT_                ) { ++uFailMinHTCut; return; }
	if ( bApplyMhtCut  && myMht < dMinMHT_              ) { ++uFailMinPFMHTCut; return; }
	if ( bApplyDphiCut && ! bPassDphiCut                ) { ++uFailDphiCut;   return; }
	if ( myMetSig < dMinMetSig_                         ) { ++uFailMetSigCut; return; }



	//loop over all ht bins
	for (unsigned i =0; i < htBins_.size()-1; ++i)
	{
		//stringstream dirName;
		//dirName << "HT" << htBins_.at(i) << "to" << htBins_.at(i+1);
	
		if (myHt >= htBins_.at(i) && myHt < htBins_.at(i+1)) 
		{
			FillHistograms(
				pfmet,
				vertexHandle,
				pfpt30eta50JetHandle,
				pfpt50eta25JetHandle,
				myHt, myMht, myMhtPhi,
				prescaleWeight,
				Weight,
				vHist.at(i)
				);
			break;
		}
	}

	passed->Fill(1);
}


// ------------ method called once each job just before starting event loop  ------------
void 
MetSigForQcd::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
void 
MetSigForQcd::endJob() 
{
	std::cout << "------------- " << __FUNCTION__ << ": summary " << std::endl;
	std::cout <<"MSQ: Events Processed  = " << processed->GetEntries() << std::endl;
	std::cout <<"MSQ: Events Passed     = " << passed->GetEntries() << std::endl;
	std::cout <<"MSQ: MetSig min        = " << dMinMetSig_ << std::endl;
	std::cout <<"MSQ: HT min            = " << dMinHT_ << std::endl;
	std::cout <<"MSQ: MHT min           = " << dMinMHT_ << std::endl;
	std::cout <<"MSQ: Njet min          = " << dMinNjet_ << std::endl;
	std::cout <<"MSQ: Njet max          = " << dMaxNjet_ << std::endl;
	std::cout <<"MSQ: ApplyNjetCut      = " << bApplyNjetCut << std::endl;
	std::cout <<"MSQ: ApplyHtCut        = " << bApplyHtCut << std::endl;
	std::cout <<"MSQ: ApplyMhtCut       = " << bApplyMhtCut << std::endl;
	std::cout <<"MSQ: ApplyDphiCut      = " << bApplyDphiCut << std::endl;
}

// ------------ method called when starting to processes a run  ------------
void 
MetSigForQcd::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
MetSigForQcd::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MetSigForQcd::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MetSigForQcd::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MetSigForQcd::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void MetSigForQcd::BookHistograms()
{
	edm::Service<TFileService> fs;
	assert (htBins_.size()>1 && "MetSigForQcd:: htBins_ size must be >1!");
	
	processed = fs->make<TH1F> ("processed"  ,"Events Processed;", 10, 0, 10);
	passed    = fs->make<TH1F> ("passed"  ,"Events Passed;", 10, 0, 10);

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
void MetSigForQcd::BookCommonHistograms(TFileDirectory& dir, const float htMin, 
								const float htMax, Hist_t& hist)
{

	//ht/mht/jet et,eta,phi,dphi-jet-mht, jet mass

	const double evt_mht_max = 1500, evt_mht_bins = 750;
	const double evt_ht_max = 4000, evt_ht_bins = 800;
	hist.evt.nvtx   = dir.make<TH1F> ("nvtx"  ,"RA2 Good Vertices; N Vtx;Events;", 40, 0, 40);
	hist.evt.mht    = dir.make<TH1F> ("mht"  ,"Calculated MHT;MHT [GeV];Events;", evt_mht_bins, 0, evt_mht_max);
	hist.evt.ht     = dir.make<TH1F> ("ht"   ,"HT from pfHT Handle ;pfHT [GeV];Events;", evt_ht_bins, 0, evt_ht_max);
	hist.evt.njet30 = dir.make<TH1F> ("njet30" ,"Njets (Et>30 GeV && | #Eta |<5.0;NJETS;Events;", 20, 0, 20);
	hist.evt.njet50 = dir.make<TH1F> ("njet50" ,"Njets (Et>50 GeV && | #Eta |<2.4;NJETS;Events;", 20, 0, 20);
	hist.evt.meff   = dir.make<TH1F> ("meff" ,"RA2:;MEff;Events;", 50, 0, 5000);
	hist.evt.prescaleWeights = dir.make<TH1F> ("prescaleWeights","Prescale Weights",2000,0,2000);
	hist.evt.eventWeights    = dir.make<TH1F> ("totalEventWeights","Total Event Weights",2000,0,2000);

	hist.evt.mht->Sumw2();
	hist.evt.ht->Sumw2();
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
	
	hist.evt.metsig = dir.make<TH1F> ("metsig","pf MET SIG", 100,0,500); 
	hist.evt.metprob = dir.make<TH1F> ("metsigprob","pf MET SIG Prob", 100,0,1); 
	hist.evt.metsigerr = dir.make<TH1F> ("metsigerr","pf MET SIG Err", 200,0,20); 
	hist.evt.metsigVsSumEt = dir.make<TH2F> ("metsigVsSumEt","pf MET SIG vs sqrt(SumEt)",500,0,1000, 1000,0,100); 
	hist.evt.metprobVsSumEt = dir.make<TH2F> ("metprobVsSumEt","pf MET Prob vs sqrt(SumEt)",500,0,1000, 100,0,1); 
	hist.evt.metsigVsSumEt_prof = dir.make<TProfile> ("metsigVsSumEtProf","pf MET Prob Vs sqrt(SumEt)",50,0,50,0,100); 
	hist.evt.metprobVsSumEt_prof = dir.make<TProfile> ("metprobVsSumEtProf","pf MET Prob Vs sqrt(SumEt)",50,0,50,0,1); 

	hist.evt.metsigVsNjet30_prof = dir.make<TProfile> ("metsigVsNjet30Prof","pf MET Sig Vs Njet (Et>30)",50,0,50,0,100); 
	hist.evt.metprobVsNjet30_prof = dir.make<TProfile> ("metprobVsNjet30Prof","pf MET Prob Vs Njet (Et>30)",50,0,50,0,1); 

	hist.evt.metsigVsNjet50_prof = dir.make<TProfile> ("metsigVsNjet50Prof","pf MET Sig Vs Njet (Et>50)",50,0,50,0,100); 
	hist.evt.metprobVsNjet50_prof = dir.make<TProfile> ("metprobVsNjet50Prof","pf MET Prob Vs Njet (Et>50)",50,0,50,0,1); 

	hist.evt.metsigVsHT = dir.make<TProfile> ("metsigVsHT","pf METSig Vs HT",300,0,3000,0,100); 
	hist.evt.metprobVsHT = dir.make<TProfile> ("metprobVsHT","pf METProb Vs HT",300,0,3000,0,1); 
	hist.evt.metsigVsMHT = dir.make<TProfile> ("metsigVsMHT","pf METSig Vs MHT",100,0,1000,0,100); 
	hist.evt.metprobVsMHT = dir.make<TProfile> ("metprobVsMHT","pf METProb Vs MHT",100,0,1000,0,1); 

}
void MetSigForQcd::FillHistograms(
	  	      edm::Handle<reco::PFMETCollection> pfmet,
				edm::Handle<reco::VertexCollection> vertexHandle,
				edm::Handle<std::vector<pat::Jet> > pf30eta50JetHandle,
				edm::Handle<std::vector<pat::Jet> > pf50eta25JetHandle,
				const float ht, const float& mht, const float mhtphi,
				const float prescaleWeight,
				const float Weight,
				Hist_t& hist
				)
{

	hist.evt.prescaleWeights->Fill(prescaleWeight);
	hist.evt.eventWeights->Fill(Weight);

	hist.evt.nvtx->Fill(vertexHandle->size(), Weight);

	hist.evt.njet50->Fill(pfpt50eta25JetHandle->size(), Weight);
	hist.evt.njet30->Fill(pfpt30eta50JetHandle->size(), Weight);

	hist.evt.meff->Fill(mht+ht, Weight);
	hist.evt.mht->Fill(mht, Weight);
	hist.evt.ht->Fill(ht, Weight);
	
	const int njets = (*pf30eta50JetHandle).size();
	const int maxJets =  (njets>10) ? 10 : njets;
	for (int i=0; i < maxJets; ++i)
	{
		const float pt  = (*pf30eta50JetHandle)[i].pt();
		const float phi = (*pf30eta50JetHandle)[i].phi();
		const float eta = (*pf30eta50JetHandle)[i].eta();
		const float dphi= fabs(TVector2::Phi_mpi_pi(phi - mhtphi));

		hist.jet[i].pt->Fill(pt, Weight);
		hist.jet[i].phi->Fill(phi, Weight);
		hist.jet[i].eta->Fill(eta, Weight);
		hist.jet[i].delphi->Fill(dphi, Weight);

	}

	float mpfMET_  = -1, mpfSumET_ = -1, mpfMETSign_ = -1, mpfMETSign_err = -1,
						ms_njets = -1, ms_neles = -1, ms_nmuons = -1;

	mpfMET_   = (*pfmet)[0].et();
	mpfSumET_ = (*pfmet)[0].sumEt();
	mpfMETSign_ = (*pfmet)[0].significance();
//	mpfMETSign_err = (*pfmet)[0].error();
//	ms_njets = (*pfmet)[0].getNumberOfJets();
//	ms_neles = (*pfmet)[0].getNumberOfElectrons();
//	ms_nmuons = (*pfmet)[0].getNumberOfMuons();

	const float prob = TMath::Prob(mpfMETSign_, met_ndof);

	hist.evt.metsig->Fill(mpfMETSign_, Weight);
	hist.evt.metprob->Fill(prob, Weight);
	//hist.evt.metsigerr->Fill(mpfMETSign_err, Weight);

	hist.evt.metsigVsSumEt->Fill(sqrt(mpfSumET_), mpfMETSign_, Weight);
	hist.evt.metprobVsSumEt->Fill(sqrt(mpfSumET_), prob, Weight);

	hist.evt.metprobVsSumEt_prof->Fill(sqrt(mpfSumET_), prob, Weight);
	hist.evt.metsigVsSumEt_prof->Fill(sqrt(mpfSumET_), mpfMETSign_, Weight);

	hist.evt.metprobVsNjet30_prof->Fill(pfpt30eta50JetHandle->size(), prob, Weight);
	hist.evt.metsigVsNjet30_prof->Fill(pfpt30eta50JetHandle->size(), mpfMETSign_, Weight);

	hist.evt.metprobVsNjet50_prof->Fill(pfpt50eta25JetHandle->size(), prob, Weight);
	hist.evt.metsigVsNjet50_prof->Fill(pfpt50eta25JetHandle->size(), mpfMETSign_, Weight);
	
	hist.evt.metsigVsHT->Fill(ht, mpfMETSign_, Weight);
	hist.evt.metprobVsHT->Fill(ht, prob, Weight);

	hist.evt.metsigVsMHT->Fill(mht, mpfMETSign_, Weight);
	hist.evt.metprobVsMHT->Fill(mht, prob, Weight);

}


void MetSigForQcd::PrintHeader()
{
	std::cout  << " >>>>>>>>>>>>>>>  Run:Lumi:Event: "
			<< kRunLumiEvent.run << ":" << kRunLumiEvent.lumi 
			<< ":" << kRunLumiEvent.evt 
			<< "<<<<<<<<<<<<<<" << std::endl;
}



//define this as a plug-in
DEFINE_FWK_MODULE(MetSigForQcd);

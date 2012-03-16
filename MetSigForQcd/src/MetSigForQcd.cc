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
// $Id$
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
			TH1F* processed;
			TH1F* passed;
			TH1F* nvtx;
			TH1F* mht;
			TH1F* pfmht;
			TH1F* ht;	//calculated by hand
			TH1F* pfht; //ht in PAT
			TH1F* njet30;
			TH1F* njet50;
			TH1F* meff;
			TH1F* metsig;
			TH1F* metprob;
			TH2F* metsigVsSumEt;
			TH2F* metprobVsSumEt;
			TProfile* metprobVsSumEt_prof;
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
				edm::Handle<edm::View<reco::MET> > mhtHandle,
				edm::Handle<double> htHandle, const float& myht, const float& mymht,
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

		unsigned int uProcessed; // number of processed events
		unsigned int uPassed; //number of events passed the filter
		RunLumiEvt_t kRunLumiEvent;
		int iVerbose; // control print levels
		double dMinHT;
		double dMinMHT;
		edm::InputTag patJetsPFPt50Eta25InputTag_, patJetsPFPt30Eta50InputTag_;
		edm::Handle<std::vector<pat::Jet> > pfpt50eta25JetHandle, pfpt30eta50JetHandle;
		unsigned uFailMinHTCut, uFailMinPFMHTCut;

		edm::LumiReWeighting LumiWeights_;
		std::vector<float> vMCNvtxDist, vDATANvtxDist;
		double sumLumiWeights; //sum of lumi weights for cross checks
		double Weight; //total weight for an event
		int doLumiWeighing, doEventWeighing;
		edm::Handle<double> prescaleWeightHandle;
		bool usePrescaleWeight;

		edm::Handle<std::vector<reco::PFMET> >pfMetHandle;
    	std::vector<double> htBins_;
		std::vector<Hist_t> vHist;
		float myHt, myMht;
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
	dMinHT          = iConfig.getUntrackedParameter<double>("dMinHT",0.0);
	dMinMHT         = iConfig.getUntrackedParameter<double>("dMinMHT",0.0);
	iVerbose        = iConfig.getUntrackedParameter<int>("verbose",0);
	doLumiWeighing  = iConfig.getUntrackedParameter<int>("ApplyLumiWeighing",0);
	doEventWeighing = iConfig.getUntrackedParameter<int>("ApplyEventWeighing",0);
	prescaleWeightInputTag = iConfig.getParameter<edm::InputTag>("prescaleWeight");
	usePrescaleWeight = iConfig.getUntrackedParameter<int>("usePrescaleWeight",0);
	sumLumiWeights   = 0;
	Weight           = 1;
	htBins_ = iConfig.getParameter<std::vector<double > >("htBins");

	BookHistograms();
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
	if ( doLumiWeighing )
	{
		lumiWeight = LumiWeights_.weight( iEvent );
		//std::cout << "lum wgt = " << lumiWeight << std::endl;
		sumLumiWeights += lumiWeight;
		Weight *= lumiWeight;
	}

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

	//APPLY RA2 cuts

	if ( pfpt50eta25JetHandle->size() <3 ) return;
	if ( (*htHandle) < dMinHT) { ++uFailMinHTCut; return; }
	if ( (*mhtHandle)[0].pt() < dMinMHT ) {++uFailMinPFMHTCut; return; }

	//recalc ht/mht
	float myHt = 0, myMht = 0;
	for (unsigned i=0; i<pfpt50eta25JetHandle->size(); ++i) myHt += (*pfpt50eta25JetHandle)[i].pt();

	TLorentzVector tlVec(0,0,0,0);
	for (unsigned i=0; i<pfpt30eta50JetHandle->size(); ++i)
	{
		const TLorentzVector iJetVec(
					(*pfpt30eta50JetHandle)[i].px(), (*pfpt30eta50JetHandle)[i].py(),
					(*pfpt30eta50JetHandle)[i].pz(), (*pfpt30eta50JetHandle)[i].energy()
					);
		tlVec += iJetVec;
	}
	myMht = tlVec.Pt();


	//Check if RA2 dphi cuts are satisfied
	bool bPassRA2dphiCut = true;
	bPassRA2dphiCut = bPassRA2dphiCut && (std::abs(reco::deltaPhi((*pfpt30eta50JetHandle)[0].phi(), (*mhtHandle)[0].phi())) > 0.5);
	bPassRA2dphiCut = bPassRA2dphiCut && (std::abs(reco::deltaPhi((*pfpt30eta50JetHandle)[1].phi(), (*mhtHandle)[0].phi())) > 0.5);
	bPassRA2dphiCut = bPassRA2dphiCut && (std::abs(reco::deltaPhi((*pfpt30eta50JetHandle)[2].phi(), (*mhtHandle)[0].phi())) > 0.3);


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
				mhtHandle,
				htHandle,
				myHt, 
				myMht,
				prescaleWeight,
				Weight,
				vHist.at(i)
				);
			break;
		}
	}

	++uPassed;



}


// ------------ method called once each job just before starting event loop  ------------
void 
MetSigForQcd::beginJob()
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
MetSigForQcd::endJob() 
{
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
	hist.evt.pfmht  = dir.make<TH1F> ("pfmht","MHT from PFmetHandle;MHT [GeV];Events;", evt_mht_bins, 0, evt_mht_max);
	hist.evt.mht    = dir.make<TH1F> ("mht"  ,"Calculated MHT;MHT [GeV];Events;", evt_mht_bins, 0, evt_mht_max);
	hist.evt.ht     = dir.make<TH1F> ("ht"   ,"HT from pfHT Handle ;pfHT [GeV];Events;", evt_ht_bins, 0, evt_ht_max);
	hist.evt.pfht   = dir.make<TH1F> ("pfht" ,"HT from Jets ET>50 GeV && |#Eta|<2.4;HT [GeV];Events;", evt_ht_bins, 0, evt_ht_max);
	hist.evt.njet30 = dir.make<TH1F> ("njet30" ,"Njets (Et>30 GeV && | #Eta |<5.0;NJETS;Events;", 20, 0, 20);
	hist.evt.njet50 = dir.make<TH1F> ("njet50" ,"Njets (Et>50 GeV && | #Eta |<2.4;NJETS;Events;", 20, 0, 20);
	hist.evt.meff   = dir.make<TH1F> ("meff" ,"RA2:;MEff;Events;", 50, 0, 5000);
	hist.evt.prescaleWeights = dir.make<TH1F> ("prescaleWeights","Prescale Weights",2000,0,2000);
	hist.evt.eventWeights    = dir.make<TH1F> ("totalEventWeights","Total Event Weights",2000,0,2000);

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
	
	hist.evt.metsig = dir.make<TH1F> ("metsig","pf MET SIG", 1000,0,100); 
	hist.evt.metprob = dir.make<TH1F> ("metsidprob","pf MET SIG Prob", 100,0,1); 
	hist.evt.metsigerr = dir.make<TH1F> ("metsigerr","pf MET SIG Err", 200,0,20); 
	hist.evt.metsigVsSumEt = dir.make<TH2F> ("metsigVsSumEt","pf MET SIG vs sqrt(SumEt)",500,0,1000, 1000,0,100); 
	hist.evt.metprobVsSumEt = dir.make<TH2F> ("metprobVsSumEt","pf MET Prob vs sqrt(SumEt)",500,0,1000, 100,0,1); 
	hist.evt.metprobVsSumEt_prof = dir.make<TProfile> ("metprobVsSumEtProf","pf MET Prob Vs sqrt(SumEt)",500,0,1000,0,1); 

}
void MetSigForQcd::FillHistograms(
	  	      edm::Handle<reco::PFMETCollection> pfmet,
				edm::Handle<reco::VertexCollection> vertexHandle,
				edm::Handle<std::vector<pat::Jet> > pf30eta50JetHandle,
				edm::Handle<std::vector<pat::Jet> > pf50eta25JetHandle,
				edm::Handle<edm::View<reco::MET> > mhtHandle,
				edm::Handle<double> htHandle, const float& myht, const float& mymht,
				const float prescaleWeight,
				const float Weight,
				Hist_t& hist
				)
{

	hist.evt.prescaleWeights->Fill(prescaleWeight);
	hist.evt.eventWeights->Fill(Weight);

	hist.evt.nvtx->Fill(vertexHandle->size(), Weight);

	const float mht = (*mhtHandle)[0].pt();
	const float ht  = (*htHandle);

	hist.evt.njet50->Fill(pfpt50eta25JetHandle->size(), Weight);
	hist.evt.njet30->Fill(pfpt30eta50JetHandle->size(), Weight);

	hist.evt.pfmht->Fill(mht, Weight);
	hist.evt.pfht->Fill(ht, Weight);
	hist.evt.meff->Fill(mht+ht, Weight);
	hist.evt.mht->Fill(mymht, Weight);
	hist.evt.ht->Fill(myht, Weight);
	
	const int njets = (*pf30eta50JetHandle).size();
	const int maxJets =  (njets>10) ? 10 : njets;
	for (int i=0; i < maxJets; ++i)
	{
		const float pt  = (*pf30eta50JetHandle)[0].pt();
		const float phi = (*pf30eta50JetHandle)[0].phi();
		const float eta = (*pf30eta50JetHandle)[0].eta();
		const float dphi= fabs(TVector2::Phi_mpi_pi(phi - (*mhtHandle)[0].phi()));

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

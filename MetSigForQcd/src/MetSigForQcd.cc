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
// $Id: MetSigForQcd.cc,v 1.4 2012/07/25 19:53:07 samantha Exp $
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
const static unsigned maxAllowedMSbins = 100;
const static unsigned maxAllowedJets = 12;

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
			TH1F* met;
			TH1F* ht;
			TH1F* njet30;
			TH1F* njet50;
			TH1F* metsig;
			TH1F* metsigprob;
			TH2F* metsigVsMet;
			TH2F* metsigVsMht;  //this can be compared to mhtsigVsMht plot check the deviation againts comman x-axis
			TH1F* mhtsig;
			TH1F* mhtsigprob;
			TH2F* mhtsigVsMht;
			TH1F* metAndmhtsigComp;   //1- mhtsig/metsig;
			TH2F* metsigVsmhtsig;    //just 2d plot
			TH1F* metSig_alt1;
			TH1F* mhtSig_alt1;
//			TH1F* metsigerr;
			TH1F* prescaleWeights;
			TH1F* eventWeights;
			TProfile* mhtsigVsMHT_prof;
			TProfile* metsigVsMET_prof;
			TProfile* mhtsigVsNjet50_prof;
			TProfile* mhtsigVsNjet30_prof;
			TProfile* metsigVsNjet50_prof;
			TProfile* metsigVsNjet30_prof;
		};

		struct JetHist_t{
			TH1F* pt;
			TH1F* eta;
			TH1F* phi;
			TH1F* delphi;  //delphi(jet,MHT)
		};

		struct MSigHist_t { //these are for each of the metsig/mhtsig bins
			TH1F* mhtormet;
			TH1F* significance;
			TH1F* njet30;
			TH1F* njet50;
			TH1F* dphimin;
		};

		struct Hist_t {  //hists for each inclusive jet category
			EventHist_t evt;
			JetHist_t jet[maxAllowedJets];
			MSigHist_t metsig[maxAllowedMSbins];
			MSigHist_t mhtsig[maxAllowedMSbins];
			MSigHist_t metsig_alt[maxAllowedMSbins];
			MSigHist_t mhtsig_alt[maxAllowedMSbins];
			MSigHist_t mhtcut[maxAllowedMSbins];
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
		void BookMetSigHistograms(TFileDirectory& dir, const float metsig, 
								MSigHist_t& hist);
		void BookMhtSigHistograms(TFileDirectory& dir, const float mhtsig, 
								MSigHist_t& hist);

		void FillHistograms(
				edm::Handle<reco::VertexCollection> vertexHandle,
				edm::Handle<std::vector<pat::Jet> > pf30eta50JetHandle,
				edm::Handle<std::vector<pat::Jet> > pf50eta25JetHandle,
	  	      const float& pfmet,
				const float& ht, const float& mht, const float& mhtphi,
				const float& metsig, const float& metsigprob, 
				const float& mhtsig, const float& mhtsigprob, 
				const float& metSig_alt1,
				const float& mhtSig_alt1,
				const float& prescaleWeight,
				const float& Weight,
				Hist_t& hist
				);
		void BookMSigHistograms(TFileDirectory& dir, const float& mxtsig, 
								MSigHist_t& hist, const int& metormht);

		void FillMSHistograms(
						const float& mhtormet,
						const float& significance,
						const unsigned& njet30,
						const unsigned& njet50,
						const float& dphimin,
						MSigHist_t& hist,
						const float Weight
						);

		float DelPhiMin(edm::Handle<std::vector<pat::Jet> > jetHandle
				, edm::Handle<edm::View<reco::MET> > mhtHandle);
				// ----------member data ---------------------------
		//jet collections
		edm::InputTag mhtInputTag_, htInputTag_;
		edm::Handle<edm::View<reco::MET> > mhtHandle;
		edm::Handle<double> htHandle;
		edm::InputTag prescaleWeightInputTag;
		edm::InputTag mhtSigProducerInputTag;

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
		edm::Handle<std::vector<reco::MET> > mhtSigProducerHandle;
		bool usePrescaleWeight;
		bool bApplyNjetCut, bApplyHtCut, bApplyMhtCut, bApplyDphiCut;

		edm::Handle<std::vector<reco::PFMET> >pfMetHandle;
    	std::vector<double> htBins_, mhtBins_, msigBins_, msigAltBins_;
		bool cutOnMHTorMET; //safety to make sure cut is made only one of them at a time.
		std::vector<Hist_t> vHist;
		float myHt, myMht;
		TH1F* processed;
		TH1F* passed;
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
	mhtBins_ = iConfig.getParameter<std::vector<double > >("mhtBins");
	msigBins_ = iConfig.getParameter<std::vector<double > >("msigBins");
	msigAltBins_ = iConfig.getParameter<std::vector<double > >("msigAltBins");
	bApplyNjetCut = iConfig.getUntrackedParameter<int>("applyNjetCut",1);
	bApplyHtCut   = iConfig.getUntrackedParameter<int>("applyHtCut",1);
	bApplyMhtCut  = iConfig.getUntrackedParameter<int>("applyMhtCut",1);
	bApplyDphiCut = iConfig.getUntrackedParameter<int>("applyDphiCut",1);
	dMinMetSig_   = iConfig.getUntrackedParameter<double>("dMinMetSig",0.0);
	mhtSigProducerInputTag = iConfig.getParameter<edm::InputTag>("mhtSigProducer");

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
	edm::Handle<reco::PFMETCollection> pfmetHandle;
	iEvent.getByLabel("pfMet", pfmetHandle);
	if (pfmetHandle->size() == 0) {
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":pfMET handle not found!" << std::endl;
		assert(false);
	}


	iEvent.getByLabel(mhtInputTag_, mhtHandle);
	if (! mhtHandle.isValid()) 
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":MHT handle not found!" << std::endl;
		assert(false);
	}


	iEvent.getByLabel(mhtSigProducerInputTag, mhtSigProducerHandle);
	if ( ! mhtSigProducerHandle.isValid()) 
	{ 
		std::cout << "mhtSigProducer not found!" << std::endl;	
		assert (false);
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
	const unsigned myNjet50= pfpt50eta25JetHandle->size(); 
	const unsigned myNjet30= pfpt30eta50JetHandle->size(); 
	const float myMetSig = (*pfmetHandle)[0].significance();

	//Check if RA2 dphi cuts are satisfied
	//dphi cuts for each of the jets, starting from highest pt jet
	bool bPassDphiCut = true;
	if (bApplyDphiCut)
	{
//		vector<float> dPhiCut;
//		dPhiCut.push_back(0.5); //1st jet
//		dPhiCut.push_back(0.5); //2nd jet
//		dPhiCut.push_back(0.3); //rest of the jets
//
//		for (unsigned i=0; i < pfpt30eta50JetHandle->size(); ++i) 
//		{ 
//			//this logic will work even in the case of exclusive njets
//			//bcos nminjets == nmaxjets in such cases
//			if ( i < (unsigned) dMinNjet_ )  
//			{
//				float dPhi = 0.3; //for rest of the jets.
//				if (i < dPhiCut.size()) dPhi = dPhiCut.at(i);
//
//				bPassDphiCut = bPassDphiCut && (std::abs(reco::deltaPhi((*pfpt30eta50JetHandle)[i].phi(), myMhtPhi)) > dPhi);
//
//			} else break;
//		}
		
		//apply dphi cut only to first 3 jets.
				
		if (pfpt30eta50JetHandle->size()>0) bPassDphiCut = bPassDphiCut && (std::abs(reco::deltaPhi((*pfpt30eta50JetHandle)[0].phi(), myMhtPhi)) > 0.5);
		if (pfpt30eta50JetHandle->size()>1) bPassDphiCut = bPassDphiCut && (std::abs(reco::deltaPhi((*pfpt30eta50JetHandle)[1].phi(), myMhtPhi)) > 0.5);
		if (pfpt30eta50JetHandle->size()>2) bPassDphiCut = bPassDphiCut && (std::abs(reco::deltaPhi((*pfpt30eta50JetHandle)[2].phi(), myMhtPhi)) > 0.3);

	}


	//APPLY cuts

	if ( bApplyNjetCut && (myNjet50 < (unsigned) dMinNjet_ || myNjet50 > (unsigned) dMaxNjet_ ) ) { ++uFailNjetCut;  return; }
	if ( bApplyHtCut   && myHt < dMinHT_                ) { ++uFailMinHTCut; return; }
	if ( bApplyDphiCut && ! bPassDphiCut                ) { ++uFailDphiCut;   return; }
	//do not need a MET cut. Just use MHT only. But use MET in plots.
	if ( bApplyMhtCut  && myMht < dMinMHT_              ) { ++uFailMinPFMHTCut; return; }
	
	//this is not necessary if I am using msig bins
	//if ( myMetSig < dMinMetSig_                         ) { ++uFailMetSigCut; return; }


	const float pfMet    = (*pfmetHandle)[0].et();
	const float pfMetHt  = (*pfmetHandle)[0].sumEt(); //don't use this. I do not understand how it is calculated.
	const float pfMetSig = (*pfmetHandle)[0].significance();
	const float pfMetSigProb = TMath::Prob(pfMet, met_ndof);
//	const float pfMetSigerr = (*pfmetHandle)[0].error();
//	const float ms_njets    = (*pfmetHandle)[0].getNumberOfJets();
//	const float ms_neles    = (*pfmetHandle)[0].getNumberOfElectrons();
//	const float ms_nmuons   = (*pfmetHandle)[0].getNumberOfMuons();
	const float dPhiMin = DelPhiMin(pfpt30eta50JetHandle,	mhtHandle);

	//calculate various METsig defnitions	
	const float metSig_alt1 = pfMet/sqrt(myHt);
	const float mhtSig_alt1 = myMht/sqrt(myHt);
	const float mhtSig      = (*mhtSigProducerHandle)[0].significance(); //lhx's metsig
	const float mhtSigProb  = TMath::Prob(mhtSig, met_ndof);

	//loop over all ht bins and fille histograms
	for (unsigned i =0; i < (htBins_.size()-1); ++i)
	{
		if (myHt >= htBins_.at(i) && myHt < htBins_.at(i+1)) 
		{
			FillHistograms(
					vertexHandle,
					pfpt30eta50JetHandle,
					pfpt50eta25JetHandle,
					pfMet,
					myHt, myMht, myMhtPhi,
					pfMetSig, pfMetSigProb,
					mhtSig, mhtSigProb,
					metSig_alt1,
					mhtSig_alt1,
					prescaleWeight,
					Weight,
					vHist.at(i)
					);

			//loop over metsig bins	
			for (unsigned j =0; j < (msigBins_.size()-1); ++j)
			{
				if (pfMetSig > msigBins_.at(j)) 
				{
					FillMSHistograms(pfMet, pfMetSig, myNjet30, myNjet50, dPhiMin, vHist.at(i).metsig[j], Weight);
				}
				if (mhtSig > msigBins_.at(j)) 
				{
					FillMSHistograms(myMht, mhtSig, myNjet30, myNjet50, dPhiMin, vHist.at(i).mhtsig[j], Weight);
				}
			}
			for (unsigned j =0; j < (msigAltBins_.size()-1); ++j)
			{
				if (metSig_alt1 > msigAltBins_.at(j)) 
				{
					FillMSHistograms(pfMet, metSig_alt1, myNjet30, myNjet50, dPhiMin, vHist.at(i).metsig_alt[j], Weight);
				}
				if (mhtSig_alt1 > msigAltBins_.at(j)) 
				{
					FillMSHistograms(myMht, mhtSig_alt1, myNjet30, myNjet50, dPhiMin, vHist.at(i).mhtsig_alt[j], Weight);
				}
			}

			for (unsigned j =0; j < (mhtBins_.size()-1); ++j)
			{
				if (myMht > mhtBins_.at(j)) 
				{
					//put 5.5 inplace for MetSig to avoid any misunderstanding when looking
					//at significane for MHT cut when there is none
					FillMSHistograms(myMht,5.5, myNjet30, myNjet50, dPhiMin, vHist.at(i).mhtcut[j], Weight);
				}
			}

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
	
	processed = fs->make<TH1F> ("processed"  ,"Events Processed;", 2, 0, 2);
	passed    = fs->make<TH1F> ("passed"  ,"Events Passed;", 2, 0, 2);

	if ((msigBins_.size()-1) > maxAllowedMSbins )
	{
		cout << __FUNCTION__ << ": Maximum MSig bins allowed is " << maxAllowedMSbins << "!!!" << endl;
		assert (false);
	}
	if ((msigAltBins_.size()-1) > maxAllowedMSbins )
	{
		cout << __FUNCTION__ << ": Maximum MSigAlt bins allowed is " << maxAllowedMSbins << "!!!" << endl;
		assert (false);
	}

	for (unsigned i =0; i < htBins_.size()-1; ++i)
	{
		stringstream dirName;
		dirName << "HT" << htBins_.at(i) << "to" << htBins_.at(i+1);
		TFileDirectory subDir = fs->mkdir(dirName.str());
		TFileDirectory subDir1 = subDir.mkdir("MetSig");
		TFileDirectory subDir2 = subDir.mkdir("MhtSig");
		TFileDirectory subDir3 = subDir.mkdir("MetSigAlt");
		TFileDirectory subDir4 = subDir.mkdir("MhtSigAlt");
		TFileDirectory subDir5 = subDir.mkdir("MhtCut");
		Hist_t hist;
		vHist.push_back(hist);
		BookCommonHistograms(subDir, htBins_.at(i), htBins_.at(i+1), vHist.at(i));

		for (unsigned j = 0; j < msigBins_.size()-1; ++j)
		{
			BookMSigHistograms(subDir1, msigBins_.at(j), vHist.at(i).metsig[j], 1);
			BookMSigHistograms(subDir2, msigBins_.at(j), vHist.at(i).mhtsig[j], 2);
		}
		for (unsigned j = 0; j < msigAltBins_.size()-1; ++j)
		{
			BookMSigHistograms(subDir3, msigAltBins_.at(j), vHist.at(i).metsig_alt[j], 3);
			BookMSigHistograms(subDir4, msigAltBins_.at(j), vHist.at(i).mhtsig_alt[j], 4);
		}
		for (unsigned j = 0; j < mhtBins_.size()-1; ++j)
		{
			BookMSigHistograms(subDir5, mhtBins_.at(j), vHist.at(i).mhtcut[j], 5);
		}

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
	hist.evt.met    = dir.make<TH1F> ("met"  ,"Calculated MET;MET [GeV];Events;", evt_mht_bins, 0, evt_mht_max);
	hist.evt.ht     = dir.make<TH1F> ("ht"   ,"HT from pfHT Handle ;pfHT [GeV];Events;", evt_ht_bins, 0, evt_ht_max);
	hist.evt.njet30 = dir.make<TH1F> ("njet30" ,"Njets (Et>30 GeV && | #Eta |<5.0;NJETS;Events;", 20, 0, 20);
	hist.evt.njet50 = dir.make<TH1F> ("njet50" ,"Njets (Et>50 GeV && | #Eta |<2.4;NJETS;Events;", 20, 0, 20);
	hist.evt.prescaleWeights = dir.make<TH1F> ("prescaleWeights","Prescale Weights",2000,0,2000);
	hist.evt.eventWeights    = dir.make<TH1F> ("totalEventWeights","Total Event Weights",2000,0,2000);

	hist.evt.mht->Sumw2();
	hist.evt.met->Sumw2();
	hist.evt.ht->Sumw2();
	hist.evt.njet30->Sumw2();
	hist.evt.njet50->Sumw2();
	hist.evt.prescaleWeights->Sumw2();
	hist.evt.eventWeights->Sumw2();

	const double pt_bins = 150, pt_max = 1500;
	for (int i =0; i < 7; ++i)
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
	
	hist.evt.metsig = dir.make<TH1F> ("metsig","pf MET SIG", 100,0,1000); 
	hist.evt.metsigprob = dir.make<TH1F> ("metsigprob","pf MET SIG Prob", 100,0,1); 
//	hist.evt.metsigerr = dir.make<TH1F> ("metsigerr","pf MET SIG Err", 200,0,20); 
	hist.evt.metsigVsMet = dir.make<TH2F> ("metsigVsMet","METSig Vs MET;MET;METSig;",100,0,1000,120,0,600); 
	
	hist.evt.mhtsig = dir.make<TH1F> ("mhtsig","MHT SIG", 100,0,1000); 
	hist.evt.mhtsigprob = dir.make<TH1F> ("mhtsigprob","MHT SIGprob", 100,0,1); 
	hist.evt.mhtsigVsMht = dir.make<TH2F> ("mhtsigVsMht","MHTSig Vs MHT;MHT;METSig;",100,0,1000,120,0,600); 
	hist.evt.metsigVsMht = dir.make<TH2F> ("metsigVsMht","METSig Vs MHT;MHT;METSig;",100,0,1000,120,0,600); 
	hist.evt.metsigVsmhtsig = dir.make<TH2F> ("metsigVsmhtsig","METSig Vs MHTSig;MHTSig;METSig;",100,0,1000,120,0,600); 
	hist.evt.metAndmhtsigComp = dir.make<TH1F> ("metAndmhtsigComp","1-MHTSig/METSig;;Events;",100,-5,5); 

	hist.evt.metSig_alt1 = dir.make<TH1F> ("metSig_alt1","METsig;METsig (MET/ #sqrt{HT});Events;", 100,0,100); 
	hist.evt.mhtSig_alt1 = dir.make<TH1F> ("mhtSig_alt1","MHTsig;MHTsig (MHT/ #sqrt{HT});Events;", 100,0,100); 

	hist.evt.mhtsigVsMHT_prof     = dir.make<TProfile> ("mhtsigVsMHT_prof",";MHT;MhtSig;", evt_mht_bins, 0, evt_mht_max, 0, 500);
	hist.evt.metsigVsMET_prof     = dir.make<TProfile> ("mhtsigVsMET_prof",";MET;MetSig;", evt_mht_bins, 0, evt_mht_max, 0, 500);
	hist.evt.mhtsigVsNjet50_prof  = dir.make<TProfile> ("mhtsig_njet50_prof",";Njets [p_{T}>50 GeV, | #eta |<2.5];MHTSig;", 30, 0, 30, 0, 500);
	hist.evt.mhtsigVsNjet30_prof  = dir.make<TProfile> ("mhtsig_njet30_prof",";Njets [p_{T}>30 GeV, | #eta |<5.0];MHTSig;", 30, 0, 30, 0, 500);
	hist.evt.metsigVsNjet50_prof  = dir.make<TProfile> ("metsig_njet50_prof",";Njets [p_{T}>50 GeV, | #eta |<2.5];METSig;", 30, 0, 30, 0, 500);
	hist.evt.metsigVsNjet30_prof  = dir.make<TProfile> ("metsig_njet30_prof",";Njets [p_{T}>30 GeV, | #eta |<5.0];METSig;", 30, 0, 30, 0, 500);
}

void MetSigForQcd::BookMSigHistograms(TFileDirectory& dir, const float& mxtsig, 
								MSigHist_t& hist, const int& metsigtype)
{
	//metormht, 0==met , 1==mht
	stringstream name, title, signame, sigtitle, dphimintitle;
	stringstream nj30name, nj30title, nj50name, nj50title, dphiminname;
	int ms_max = 800, ms_bins = 400; 

	name << "met_msigmin_" << mxtsig;
	nj30name << "njet30_msigmin_" << mxtsig;
	nj50name << "njet50_msigmin_" << mxtsig;
	signame << "metsig_msigmin_" << mxtsig;
	dphiminname << "dphimin_msigmin_" << mxtsig;

	if (metsigtype == 1) 
	{ 
		title << "METSig>" << mxtsig << ";MET [GeV];Events;" << mxtsig;
		sigtitle << "METSig>" << mxtsig << ";MET-Significance;Events;";
		nj30title << "METSig>" << mxtsig << ";Njets [E_{T}>30 GeV];Events;";
		nj50title << "METSig>" << mxtsig << ";Njets [E_{T}>50 GeV];Events;";
		dphimintitle << "METSig>" << mxtsig << ";#Delta #Phi _{min};Events;";

	} else if (metsigtype == 2)
	{
		title << "MHTSig>" << mxtsig << ";MHT [GeV];Events;";
		sigtitle << "MHTSig>" << mxtsig << ";MHT-Significance;Events;";
		nj30title << "MHTSig>" << mxtsig << ";Njets [E_{T}>30 GeV];Events;";
		nj50title << "MHTSig>" << mxtsig << ";Njets [E_{T}>50 GeV];Events;";
		dphimintitle << "MHTSig>" << mxtsig << ";#Delta #Phi _{min};Events;";
	} else if (metsigtype == 3)
	{
		ms_max = 100;
		ms_bins = 100; 
		const string def("MetSig= MET/#sqrt{pfHT}");
		title << def << " >" << mxtsig << ";MET [GeV];Events;";
		sigtitle << def << " >" << mxtsig << ";MET-Significance;Events;";
		nj30title << def << " >" << mxtsig << ";Njets [E_{T}>30 GeV];Events;";
		nj50title << def << " >" << mxtsig << ";Njets [E_{T}>50 GeV];Events;";
		dphimintitle << def << " >" << mxtsig << ";#Delta #Phi _{min};Events;";

	} else if (metsigtype == 4)
	{
		ms_max = 100;
		ms_bins = 100; 
		const string def("MhtSig= MHT/#sqrt{HT}");
		title << def << " >" << mxtsig << ";MHT [GeV];Events;";
		sigtitle << def << " >" << mxtsig << ";MHT-Significance;Events;";
		nj30title << def << " >" << mxtsig << ";Njets [E_{T}>30 GeV];Events;";
		nj50title << def << " >" << mxtsig << ";Njets [E_{T}>50 GeV];Events;";
		dphimintitle << def << " >" << mxtsig << ";#Delta #Phi _{min};Events;";
	} else if (metsigtype == 5)
	{
		ms_max = 1000;
		ms_bins = 200; 
		const string def("MHT Cut");
		title << def << " >" << mxtsig << ";MHT [GeV];Events;";
		sigtitle << "This had no meaning! Just to get event count! :" << def << " >" << mxtsig << ";MHT Cut;Events;";
		nj30title << def << " >" << mxtsig << ";Njets [E_{T}>30 GeV];Events;";
		nj50title << def << " >" << mxtsig << ";Njets [E_{T}>50 GeV];Events;";
		dphimintitle << def << " >" << mxtsig << ";#Delta #Phi _{min};Events;";
	}

	const double evt_mht_max = 1500, evt_mht_bins = 750;
	hist.mhtormet     = dir.make<TH1F> (name.str().c_str(), title.str().c_str(), evt_mht_bins, 0, evt_mht_max);
	hist.significance = dir.make<TH1F> (signame.str().c_str(), sigtitle.str().c_str(), ms_bins, 0, ms_max);
	hist.njet30       = dir.make<TH1F> (nj30name.str().c_str(), nj30title.str().c_str(), 20, 0, 20);
	hist.njet50       = dir.make<TH1F> (nj50name.str().c_str(), nj50title.str().c_str(), 20, 0, 20);
	hist.dphimin      = dir.make<TH1F> (dphiminname.str().c_str(), dphimintitle.str().c_str(), 125, 0, 2.5);

	hist.mhtormet->Sumw2();
	hist.significance->Sumw2();
	hist.njet30->Sumw2();
	hist.njet50->Sumw2();
	hist.dphimin->Sumw2();


}

void MetSigForQcd::FillHistograms(
				edm::Handle<reco::VertexCollection> vertexHandle,
				edm::Handle<std::vector<pat::Jet> > pf30eta50JetHandle,
				edm::Handle<std::vector<pat::Jet> > pf50eta25JetHandle,
	  	      const float& pfmet,
				const float& ht, const float& mht, const float& mhtphi,
				const float& metsig, const float& metsigprob, 
				const float& mhtsig, const float& mhtsigprob, 
				const float& metSig_alt1,
				const float& mhtSig_alt1,
				const float& prescaleWeight,
				const float& Weight,
				Hist_t& hist
				)
{

	hist.evt.prescaleWeights->Fill(prescaleWeight);
	hist.evt.eventWeights->Fill(Weight);

	hist.evt.nvtx->Fill(vertexHandle->size(), Weight);
	hist.evt.njet50->Fill(pfpt50eta25JetHandle->size(), Weight);
	hist.evt.njet30->Fill(pfpt30eta50JetHandle->size(), Weight);
	hist.evt.mht->Fill(mht, Weight);
	hist.evt.met->Fill(pfmet, Weight);
	hist.evt.ht->Fill(ht, Weight);
	
	const int njets = (*pf30eta50JetHandle).size();
	const int maxJets =  (njets>7) ? 7 : njets;
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

	hist.evt.metsig->Fill(metsig, Weight);
	hist.evt.metsigprob->Fill(metsigprob, Weight);
	//hist.evt.metsigerr->Fill(mpfMETSign_err, Weight);
	hist.evt.metsigVsMet->Fill(pfmet, metsig, Weight);

	hist.evt.mhtsig->Fill(mhtsig, Weight);
	hist.evt.mhtsigprob->Fill(mhtsigprob, Weight);
	hist.evt.mhtsigVsMht->Fill(mht, mhtsig, Weight);
	hist.evt.metsigVsMht->Fill(mht, metsig, Weight);
	hist.evt.metsigVsmhtsig->Fill(metsig, mhtsig, Weight);

	if (metsig>0) hist.evt.metAndmhtsigComp->Fill(1.0-(mhtsig/metsig), Weight);

	hist.evt.mhtSig_alt1->Fill(mhtSig_alt1, Weight);
	hist.evt.metSig_alt1->Fill(metSig_alt1, Weight);

	hist.evt.mhtsigVsMHT_prof->Fill(mht, mhtsig, Weight);
	hist.evt.metsigVsMET_prof->Fill(pfmet, metsig, Weight);
	hist.evt.mhtsigVsNjet50_prof->Fill(pfpt50eta25JetHandle->size(), mhtsig);
	hist.evt.mhtsigVsNjet30_prof->Fill(pfpt30eta50JetHandle->size(), mhtsig);
	hist.evt.metsigVsNjet50_prof->Fill(pfpt50eta25JetHandle->size(), metsig);
	hist.evt.metsigVsNjet30_prof->Fill(pfpt30eta50JetHandle->size(), metsig);

}


void MetSigForQcd::FillMSHistograms(
	  	      const float& mhtormet,
	  	      const float& significance,
				const unsigned& njet30,
				const unsigned& njet50,
				const float& dphimin,
				MSigHist_t& hist,
				const float Weight
				)
{
	hist.mhtormet->Fill(mhtormet, Weight);
	hist.significance->Fill(significance, Weight);
	hist.njet30->Fill(njet30, Weight);
	hist.njet50->Fill(njet50, Weight);
	hist.dphimin->Fill(dphimin, Weight);

}
void MetSigForQcd::PrintHeader()
{
	std::cout  << " >>>>>>>>>>>>>>>  Run:Lumi:Event: "
			<< kRunLumiEvent.run << ":" << kRunLumiEvent.lumi 
			<< ":" << kRunLumiEvent.evt 
			<< "<<<<<<<<<<<<<<" << std::endl;
}


float MetSigForQcd::DelPhiMin(edm::Handle<std::vector<pat::Jet> > jetHandle, 
					edm::Handle<edm::View<reco::MET> > mhtHandle
				)
{
	//PrintHeader();
	std::vector<float> vDelPhi_jetmht;
	const float mht = (*mhtHandle)[0].pt();

	for (unsigned i = 0 ; i < jetHandle->size() ; ++i)
	{
		if (i>2) break; //use only three leading jets
		const float delphi_jetmht = fabs(TVector2::Phi_mpi_pi((*jetHandle)[i].phi() - (*mhtHandle)[0].phi()));
		vDelPhi_jetmht.push_back(delphi_jetmht);
	}

	std::sort(vDelPhi_jetmht.begin(), vDelPhi_jetmht.end(), sort_using_less_than);	

	assert (vDelPhi_jetmht.size() == 3 && "ERROR: More than 3 dPhiMin calculations found!");
	//std::cout << "vDelPhi_jetmht size = " << vDelPhi_jetmht.size() << std::endl;
	
	const float dPhiMin = vDelPhi_jetmht.at(0);
	return dPhiMin;
}

//define this as a plug-in
DEFINE_FWK_MODULE(MetSigForQcd);

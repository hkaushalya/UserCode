// -*- C++ -*-
//
// Package:    RA2Example
// Class:      RA2Example
// 
/**\class RA2Example RA2Example.cc SusyAnalysis/RA2Example/src/RA2Example.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Samantha Hewamanage
// $Id: RA2Example.cc,v 1.4 2013/07/25 19:16:16 samantha Exp $
//
//


// system include files
#include <memory>
#include <algorithm>
#include <vector>
#include <string>

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "DataFormats/Math/interface/deltaR.h"

// TFile Service
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

using namespace std;

class RA2Example : public edm::EDAnalyzer {

	public:

		explicit RA2Example(const edm::ParameterSet & iConfig);
		~RA2Example();

	private:

		void beginJob() ;
		void endJob() ;
		void analyze(const edm::Event&, const edm::EventSetup&);
		void BookHistograms();
		void clearTreeVectors();  
  		vector<double> generateWeights(const TH1* data_npu_estimated) const;
		double getPUWeight(int npu) const; 

		bool debug_;
		edm::InputTag vtxSrc_, HBHENoiseFiltSrc_;
		edm::InputTag firedTrigNamesSrc_, firedTrigPrescaleSrc_;
		bool doPUReWeight_;
		int doEventWeighing_;
		edm::InputTag puWeigthSrc_, puWeigthABSrc_, puWeigthABCSrc_, puWeigthRA2Src_;
		edm::InputTag pfMetSrc_;
		edm::InputTag jetAllsrc_, genjetAllsrc_ ,genParticleSrc_;
		edm::InputTag mhtSrc_, htSrc_;
		bool mcFlag_;
		double minPFJetPt_, minGENJetPt_; 
	   std::vector<double> genParticleList_;

		edm::Service<TFileService> fs;
		TTree* tree;
		unsigned int         t_EvtRun, t_EvtLS, t_EvtEvent;
		int                  t_NVertices, t_tru_Npv, t_avg_Npv;
		int                  t_MCflag;
		double               t_PUWeight, t_PUWeightAB, t_PUWeightABC, t_PUWeightRA2;
		double               t_EvtWeight; //for flat sample weighing
		double               t_PFMetPx, t_PFMetPy;
		std::vector<double> *t_PFJetPt,  *t_PFJetEta,  *t_PFJetPhi,  *t_PFJetE;
  		std::vector<double> *t_PFJetBTag,*t_PFJetNHF,  *t_PFJetEMF, *t_PFbDiscriminator1, *t_PFbDiscriminator2, *t_PFbDiscriminator3, *t_PFbDiscriminator4;
		std::vector<double> *t_PFJetHOEne;
		double               t_PFht, t_PFmht;
		int                  t_NJetsPt30Eta2p5, t_NJetsPt30Eta5p0, t_NJetsPt50Eta2p5, t_NJetsPt50Eta5p0;

		std::vector<double > *t_genJetPt, *t_genJetEta, *t_genJetPhi, *t_genJetE;
		std::vector<double > *t_genParPt, *t_genParEta, *t_genParPhi, *t_genParE, *t_genParStatus, *t_genParID;
		double               t_minPFJetPt, t_minGenJetPt;
		int 					  t_allFilters;
		int                  t_beamHaloFilter, t_eeBadScFilter, t_eeNoiseFilter, t_greedyMuons;
		int 					  t_hcalLaserEventFilter, t_inconsistentMuons, t_ra2EcalBEFilter;
		int 					  t_ra2EcalTPFilter, t_trackingFailureFilter, t_HBHENoiseFilterRA2;
		int 						t_ecalLaserCorrFilter;
		std::vector<string > *t_firedTrigs;
		std::vector<double> *t_firedTrigsPrescale;

		bool no_beamHaloFilter_, no_eeBadScFilter_, no_eeNoiseFilter_, no_greedyMuonsFilter_;
		bool no_hcalLaserEventFilter_, no_inconsistentFilter_, no_ra2EcalBEFilter_, no_ra2EcalTPFilter_;
		bool no_trackingFailureFilter_, no_HBHENoiseFilter_, no_ecalLaserCorrFilter_;
  		std::vector<double> _puWeigths;
  		std::string   btagname_;
};


RA2Example::RA2Example(const edm::ParameterSet & iConfig) {

	debug_       = iConfig.getParameter<bool>("Debug");
	vtxSrc_      = iConfig.getParameter<edm::InputTag>("VertexSource");
	doPUReWeight_= iConfig.getParameter<bool>("DoPUReweight");
	doEventWeighing_ = iConfig.getUntrackedParameter<int>("ApplyEventWeighing",0);
	puWeigthSrc_    = iConfig.getParameter<edm::InputTag>("PUWeigthSource");
	puWeigthABSrc_  = iConfig.getParameter<edm::InputTag>("PUWeigthABSource");
	puWeigthABCSrc_ = iConfig.getParameter<edm::InputTag>("PUWeigthABCSource");
	//puWeigthRA2Src_ = iConfig.getParameter<edm::InputTag>("PUWeigthRA2Source");
	pfMetSrc_     = iConfig.getParameter<edm::InputTag>("PFMetSource");
	jetAllsrc_    = iConfig.getParameter<edm::InputTag>("JetAllSource");
  	btagname_       = iConfig.getParameter<std::string>  ("bTagName");
	genjetAllsrc_ = iConfig.getParameter<edm::InputTag>("genJetAllSource");
	genParticleSrc_ = iConfig.getParameter<edm::InputTag>("genParticleSource");
	genParticleList_ = iConfig.getParameter<std::vector<double> > ("genParticleList");
	mhtSrc_       = iConfig.getParameter<edm::InputTag>("MHTSource");
	mcFlag_        = iConfig.getParameter<bool>("MCflag");
	htSrc_        = iConfig.getParameter<edm::InputTag>("HTSource");
	minPFJetPt_   = iConfig.getUntrackedParameter<double>("MinPFJetPt",0.0);
	minGENJetPt_  = iConfig.getUntrackedParameter<double>("MinGENJetPt",0.0);
	HBHENoiseFiltSrc_ = iConfig.getParameter<edm::InputTag>("HBHENoiseFiltSrc");

	no_beamHaloFilter_        = iConfig.getParameter<bool>("no_beamHaloFilter");
	no_eeBadScFilter_         = iConfig.getParameter<bool>("no_eeBadScFilter");
	no_eeNoiseFilter_         = iConfig.getParameter<bool>("no_eeNoiseFilter");
	no_greedyMuonsFilter_     = iConfig.getParameter<bool>("no_greedyMuonsFilter");
	no_hcalLaserEventFilter_  = iConfig.getParameter<bool>("no_hcalLaserEventFilter");
	no_inconsistentFilter_    = iConfig.getParameter<bool>("no_inconsistentFilter");
	no_ra2EcalBEFilter_       = iConfig.getParameter<bool>("no_ra2EcalBEFilter");
	no_ra2EcalTPFilter_       = iConfig.getParameter<bool>("no_ra2EcalTPFilter");
	no_trackingFailureFilter_ = iConfig.getParameter<bool>("no_trackingFailureFilter");
	no_HBHENoiseFilter_       = iConfig.getParameter<bool>("no_HBHENoiseFilter");
	no_ecalLaserCorrFilter_   = iConfig.getParameter<bool>("no_ecalLaserCorrFilter");

	firedTrigNamesSrc_  = iConfig.getParameter<edm::InputTag>("firedTrigNamesSource");
	firedTrigPrescaleSrc_  = iConfig.getParameter<edm::InputTag>("firedTrigsPrescalesSource");


   std::string fileNamePU = iConfig.getParameter<std::string> ("FileNamePUDataDistribution");
	if (fileNamePU.length() != 0 && fileNamePU != "NONE") {
	//	edm::FileInPath filePUDataDistr(fileNamePU);
	//	std::cout << "  Reading PU scenario from '" << filePUDataDistr.fullPath() << "'" << std::endl;
		std::cout << __FUNCTION__ << ": Reading PU scenario from " << fileNamePU << std::endl;
		//TFile file(filePUDataDistr.fullPath().c_str(), "READ");
		TFile file(fileNamePU.c_str(), "READ");
		if (file.IsZombie())
		{
			cout << __FUNCTION__ << ": PU ROOT file, " << fileNamePU << " not found!" << endl;
			assert(false);
		}
		TH1 *h = 0;
		file.GetObject("pileup", h);
		if (h) {
			h->SetDirectory(0);
			_puWeigths = generateWeights(h);
		} else {
			//cout << __FUNCTION__ << ": Hist nameed 'pileup' not found in " << filePUDataDistr.fullPath() << endl;
			cout << __FUNCTION__ << ": Hist nameed 'pileup' not found in " << fileNamePU << endl;
			assert(false);
		}
	}


}

RA2Example::~RA2Example() {
}


void RA2Example::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

	using namespace edm;
	using namespace std;

	clearTreeVectors();
	bool saveEvent = true;

	// fill event ID 
	t_EvtRun   = iEvent.id().run();
	t_EvtLS    = iEvent.luminosityBlock();
	t_EvtEvent = iEvent.id().event();

//	cout << __FILE__ << endl;
//	cout << "========================== run:lumi:evt=" <<  iEvent.id().run() << ":" <<  iEvent.luminosityBlock() << ":" << iEvent.id().event() << endl;
	// save number of vertices and position of primary vertex
	edm::Handle< std::vector<reco::Vertex> > vertices;
	iEvent.getByLabel(vtxSrc_, vertices);
	if(vertices->size()<1) std::cout << "No vertices are reconstructed - check StdCleaning ?" << std::endl;

	t_NVertices = vertices->size();

	//not sure if I need this nor why others saved it. for now just 
	//keeping it.
	t_avg_Npv = 0;

	// save pileup weight
	double pu_event_wt = 1.0;
	double pu_event_wtAB = 1.0;
	double pu_event_wtABC = 1.0;
	double pu_event_wtRA2 = 1.0;
	edm::Handle<double> puweight;
	edm::Handle<double> puweightAB;
	edm::Handle<double> puweightABC;
	edm::Handle<double> puweightRA2;
	if( doPUReWeight_ && mcFlag_) {
//		iEvent.getByLabel(puWeigthSrc_, puweight);
//		pu_event_wt = *puweight;
//		iEvent.getByLabel(puWeigthABSrc_, puweightAB);
//		pu_event_wtAB = *puweightAB;
//		iEvent.getByLabel(puWeigthABCSrc_, puweightABC);
//		pu_event_wtABC = *puweightABC;
//iEvent.getByLabel(puWeigthRA2Src_, puweightRA2);
//pu_event_wtRA2 = *puweightRA2;

		edm::Handle<std::vector<PileupSummaryInfo> > puInfo;
		iEvent.getByLabel("addPileupInfo", puInfo);
		int npu = 0;
		if (puInfo.isValid()) {
			std::vector<PileupSummaryInfo>::const_iterator puIt;
			int n = 0;
			for (puIt = puInfo->begin(); puIt != puInfo->end(); ++puIt, ++n) {
				//std::cout << " Pileup Information: bunchXing, nvtx: " << puIt->getBunchCrossing() << " " << puIt->getPU_NumInteractions() << std::endl;
				if (puIt->getBunchCrossing() == 0) { // Select in-time bunch crossing
					npu = puIt->getTrueNumInteractions();
					t_tru_Npv = npu;
					break;
				}
			}
			pu_event_wt = getPUWeight(npu);
		}
	}
	t_PUWeight = pu_event_wt;
	t_PUWeightAB = pu_event_wtAB;
	t_PUWeightABC = pu_event_wtABC;
	//t_PUWeightRA2 = pu_event_wtRA2;

	if(debug_)std::cout << "t_NVertices " << t_NVertices << "  t_PUWeight " << t_PUWeight << std::endl;
	//std::cout << "t_NVertices " << t_NVertices  << " (" << t_tru_Npv << ")" << "  t_PUWeight " << t_PUWeight << std::endl;


	//event weights for flat QCD samples
	double storedWeight = 1;
	t_EvtWeight = storedWeight;
	if ( doEventWeighing_ && mcFlag_)
	{
		edm::Handle<GenEventInfoProduct> genEvtInfoHandle;
		iEvent.getByLabel("generator", genEvtInfoHandle);
		storedWeight = genEvtInfoHandle->weight();
		//std::cout << "storedWeight = " << storedWeight << std::endl;
	}
	t_EvtWeight = storedWeight;


	// save missing transverse energy
	edm::Handle < std::vector<pat::MET> > met;
	iEvent.getByLabel(pfMetSrc_, met);
	for (std::vector<pat::MET>::const_iterator it = met->begin(); it != met->end(); ++it) {
		t_PFMetPx = it->px();
		t_PFMetPy = it->py();
	}

	// save all jets along with btag information
	edm::Handle<edm::View<pat::Jet> > jetsAll;
	iEvent.getByLabel(jetAllsrc_, jetsAll);

	t_NJetsPt30Eta2p5 = 0;
	t_NJetsPt30Eta5p0 = 0; 
	t_NJetsPt50Eta2p5 = 0;
	t_NJetsPt50Eta5p0 = 0;

	reco::MET::LorentzVector recomht(0,0,0,0);
	//std::cout << "jetsAll->size() "<< jetsAll->size() << std::endl;
	for(unsigned int ijet=0; ijet<jetsAll->size(); ijet++) {
		//std::cout << (*jetsAll)[ijet].pt() << " ,  " << std::endl;
		if( (*jetsAll)[ijet].pt() > minPFJetPt_ ) {
			t_PFJetPt  ->push_back((*jetsAll)[ijet].pt() );
			t_PFJetEta ->push_back((*jetsAll)[ijet].eta());    
			t_PFJetPhi ->push_back((*jetsAll)[ijet].phi());
			t_PFJetE   ->push_back((*jetsAll)[ijet].energy()  );
    		t_PFJetBTag  ->push_back((*jetsAll)[ijet].bDiscriminator(btagname_.c_str()) );
			t_PFJetNHF   ->push_back((*jetsAll)[ijet].neutralHadronEnergyFraction());
			t_PFJetEMF   ->push_back((*jetsAll)[ijet].photonEnergyFraction());

			//calculate energy deposited in HO for each jet
			const std::vector<reco::PFCandidatePtr> & pfcands = (*jetsAll)[ijet].getPFConstituents();
			double hoEne   = 0.0;
			std::vector<reco::PFCandidatePtr>::const_iterator itr;
			for(itr = pfcands.begin(); itr != pfcands.end(); itr++) {
				hoEne   += (*itr)->hoEnergy();
			}
			t_PFJetHOEne ->push_back(hoEne);


			if( (*jetsAll)[ijet].pt()>30.0 && std::fabs((*jetsAll)[ijet].eta())<5.0 )
			{
				t_NJetsPt30Eta5p0++;
				recomht -= (*jetsAll)[ijet].p4();
			}
			if( (*jetsAll)[ijet].pt()>30.0 && std::abs((*jetsAll)[ijet].eta())<2.5 )
				t_NJetsPt30Eta2p5++;
			if( (*jetsAll)[ijet].pt()>50.0 && std::abs((*jetsAll)[ijet].eta())<5.0 )
				t_NJetsPt50Eta5p0++;
			if( (*jetsAll)[ijet].pt()>50.0 && std::abs((*jetsAll)[ijet].eta())<2.5 )
				t_NJetsPt50Eta2p5++;
		}
	}

	// save ra2 ht/mht mainly for debugging only
	edm::Handle<double> ht;
	iEvent.getByLabel(htSrc_, ht);
	t_PFht  = *ht;

	edm::Handle<edm::View<reco::MET> > mht;
	iEvent.getByLabel(mhtSrc_, mht);
	t_PFmht = (*mht)[0].pt() ;

	if(debug_)std::cout << "ht " << t_PFht << " mht " << t_PFmht << std::endl;
	//	cout << "recomht/jetsmht = " << recomht.pt() << "/" << t_PFmht << endl; 

	/*******************************************************
	 * save gen info
	 ******************************************************/
	t_MCflag = mcFlag_;
	if (mcFlag_)
	{
		//get gen-jets info
		edm::Handle<edm::View<reco::GenJet> > genjetsAll;
		iEvent.getByLabel(genjetAllsrc_, genjetsAll);

		//std::cout << "genjetsAll->size() "<< genjetsAll->size() << std::endl;
		for(unsigned int ijet=0; ijet<genjetsAll->size(); ijet++) {
			//std::cout << (*genjetsAll)[ijet].pt() << " ,  " << std::endl;
			if( (*genjetsAll)[ijet].pt() > minGENJetPt_ ) {
				t_genJetPt  ->push_back((*genjetsAll)[ijet].pt() );
				t_genJetEta ->push_back((*genjetsAll)[ijet].eta());    
				t_genJetPhi ->push_back((*genjetsAll)[ijet].phi());
				t_genJetE   ->push_back((*genjetsAll)[ijet].energy()  );
			}
		}

		//get gen particle info
     edm::Handle<std::vector<reco::GenParticle > > genParticles;
		iEvent.getByLabel(genParticleSrc_, genParticles);

		for (unsigned ig=0; ig<genParticles->size(); ig++) 
		{
			const reco::GenParticle& gen = genParticles->at(ig);
			const int pdgId = gen.pdgId();
			if (find(genParticleList_.begin(), genParticleList_.end(),(double) pdgId) == genParticleList_.end()) continue; 
			const int pdgStatus = gen.status();
	/*		cout << __FUNCTION__ << ":Gen par ig=id/stat/pt/eta/phi/e =\t"
					<< ig << "=\t" << pdgId << "\t" << pdgStatus 
					<< "\t" << gen.pt() 
					<< "\t" << gen.eta() 
					<< "\t" << gen.phi() 
					<< "\t" << gen.energy()
					<< endl;
	*/			
			t_genParPt ->push_back(gen.pt());
			t_genParEta->push_back(gen.eta());    
			t_genParPhi->push_back(gen.phi());
			t_genParE  ->push_back(gen.energy());
			t_genParStatus ->push_back((double)pdgStatus);
			t_genParID ->push_back((double)pdgId);
			//bool frombjet = find_mother( &gen, 5 );
			//bool frombbarjet = find_mother( &gen, -5 );

		}


	}

	/*******************************************************
	 * get cleaning filter status
	 ******************************************************/
	edm::Handle<bool> beamHaloVal, eeBadVal, eeNoiseVal, greedyMuonsVal;
	edm::Handle<bool> hcalLaserEventVal, inconsMuonsVal, ra2EcalBEVal, ecalLaserCorrVal;
	edm::Handle<bool> ra2EcalTPVal, trackFailureVal, HBHENoiseVal;


	//if filter is not found in the PAT set it to 2
	//so by default it will be true when cast to bool
	//during the readout
	
	const int FILTER_NOT_SET = 2;
	t_beamHaloFilter = t_eeBadScFilter = t_eeNoiseFilter = FILTER_NOT_SET;
	t_greedyMuons = t_hcalLaserEventFilter = t_inconsistentMuons = FILTER_NOT_SET;
	t_ra2EcalBEFilter = t_ra2EcalTPFilter = t_trackingFailureFilter  = FILTER_NOT_SET;
	t_HBHENoiseFilterRA2 = FILTER_NOT_SET;

	if (! no_beamHaloFilter_) {
		iEvent.getByLabel("beamHaloFilter", beamHaloVal);
		t_beamHaloFilter  = *beamHaloVal;
	}
	if (! no_eeBadScFilter_) {
		iEvent.getByLabel("eeBadScFilter", eeBadVal);
		t_eeBadScFilter   = *eeBadVal;
	}
	if (! no_eeNoiseFilter_) {
		iEvent.getByLabel("eeNoiseFilter", eeNoiseVal);
		t_eeNoiseFilter   = *eeNoiseVal;
	}
	if (! no_greedyMuonsFilter_) {
		iEvent.getByLabel("greedyMuons", greedyMuonsVal);
		t_greedyMuons     = *greedyMuonsVal;
	}
	if (! no_hcalLaserEventFilter_) {
		iEvent.getByLabel("hcalLaserEventFilter", hcalLaserEventVal);
		t_hcalLaserEventFilter = *hcalLaserEventVal;
	}
	if (! no_inconsistentFilter_) {
		iEvent.getByLabel("inconsistentMuons", inconsMuonsVal);
		t_inconsistentMuons = *inconsMuonsVal;
	}
	if (! no_ra2EcalBEFilter_) {
		iEvent.getByLabel("ra2EcalBEFilter", ra2EcalBEVal);
		t_ra2EcalBEFilter   = *ra2EcalBEVal;
	}
	if (! no_ra2EcalTPFilter_) {
		iEvent.getByLabel("ra2EcalTPFilter", ra2EcalTPVal);
		t_ra2EcalTPFilter   = *ra2EcalTPVal;
	}
	if (! no_trackingFailureFilter_) {
		iEvent.getByLabel("trackingFailureFilter", trackFailureVal);
		t_trackingFailureFilter = *trackFailureVal;
	}
	//this is not found in 52X skims
	if (! no_HBHENoiseFilter_) {
		iEvent.getByLabel(HBHENoiseFiltSrc_, HBHENoiseVal);
		t_HBHENoiseFilterRA2 = *HBHENoiseVal;
	}

	if (! no_ecalLaserCorrFilter_) {
		iEvent.getByLabel("ecalLaserCorrFilter", ecalLaserCorrVal);
		t_ecalLaserCorrFilter   = *ecalLaserCorrVal;
	}


	//get trigger info
	
	if (! mcFlag_)
	{
		edm::Handle<vector<std::string> > firedTrigNamesHandle;
		edm::Handle<vector<unsigned> > firedTrigPrescalesHandle;
		iEvent.getByLabel(firedTrigNamesSrc_, firedTrigNamesHandle);
		iEvent.getByLabel(firedTrigPrescaleSrc_, firedTrigPrescalesHandle);

		//cout << "trig size = " << (*firedTrigNamesHandle).size() << endl;
		for (unsigned i=0 ;i < (*firedTrigNamesHandle).size(); ++i)
		{
			//cout << (*firedTrigNamesHandle).at(i) << "-> prescale = " << (*firedTrigPrescalesHandle).at(i) << endl;
			t_firedTrigs->push_back((*firedTrigNamesHandle).at(i));
			t_firedTrigsPrescale->push_back((*firedTrigPrescalesHandle).at(i));
		}

		if (t_firedTrigs->size()<1) saveEvent = saveEvent && false;
	}
	//save only trigger selected data
	if (saveEvent) tree->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void RA2Example::beginJob() {

	BookHistograms();

}

// ------------ method called once each job just after ending the event loop  ------------
void RA2Example::endJob() {

	cout << "---------" << __FILE__ << ":" << __FUNCTION__ << "---------" << endl;
	cout << "\t Saved Events = " << tree->GetEntries() << endl;
	if (mcFlag_) 
	{
		cout <<"\t Gen Particles saved = ";
		for (unsigned i=0; i< genParticleList_.size();++i)
		{
			cout << genParticleList_.at(i) << ", ";
		}
		cout << endl;
	}


}

void RA2Example::BookHistograms() {

	// book tree here
	tree = fs->make<TTree>("tree", "tree");
	tree->SetAutoSave(10000);

	tree->Branch("t_EvtRun",    &t_EvtRun,   "t_EvtRun/i");
	tree->Branch("t_EvtLS",     &t_EvtLS,    "t_EvtLS/i");
	tree->Branch("t_EvtEvent",  &t_EvtEvent, "t_EvtEvent/i");
	tree->Branch("t_NVertices", &t_NVertices,"t_NVertices/I");
	tree->Branch("t_tru_Npv",   &t_tru_Npv,  "t_tru_Npv/I");
	tree->Branch("t_avg_Npv",   &t_avg_Npv,  "t_avg_Npv/I");
	tree->Branch("t_MCflag",    &t_MCflag,   "t_MCflag/I");
	tree->Branch("t_PUWeight",  &t_PUWeight, "t_PUWeight/D");
	tree->Branch("t_PUWeightAB",   &t_PUWeightAB,  "t_PUWeightAB/D");
	tree->Branch("t_PUWeightABC",  &t_PUWeightABC, "t_PUWeightABC/D");
	tree->Branch("t_PUWeightRA2",  &t_PUWeightRA2, "t_PUWeightRA2/D");
	tree->Branch("t_EvtWeight", &t_EvtWeight, "t_EvtWeight/D");
	tree->Branch("t_PFMetPx",   &t_PFMetPx,  "t_PFMetPx/D");
	tree->Branch("t_PFMetPy",   &t_PFMetPy,  "t_PFMetPy/D");

	t_PFJetPt  = new std::vector<double>();
	t_PFJetEta = new std::vector<double>();    
	t_PFJetPhi = new std::vector<double>();
	t_PFJetE   = new std::vector<double>();
	t_PFJetBTag  = new std::vector<double>();
	t_PFJetNHF   = new std::vector<double>();
	t_PFJetEMF   = new std::vector<double>();
	t_PFJetHOEne = new std::vector<double>();
	tree->Branch("t_PFJetPt",  "vector<double>", &t_PFJetPt );
	tree->Branch("t_PFJetEta", "vector<double>", &t_PFJetEta);
	tree->Branch("t_PFJetPhi", "vector<double>", &t_PFJetPhi);
	tree->Branch("t_PFJetE",   "vector<double>", &t_PFJetE);
  tree->Branch("t_PFJetBTag",  "vector<double>", &t_PFJetBTag);
  tree->Branch("t_PFJetNHF",   "vector<double>", &t_PFJetNHF);
  tree->Branch("t_PFJetEMF",   "vector<double>", &t_PFJetEMF);
  tree->Branch("t_PFJetHOEne", "vector<double>", &t_PFJetHOEne);
	tree->Branch("t_NJetsPt30Eta2p5", &t_NJetsPt30Eta2p5, "t_NJetsPt30Eta2p5/I");
	tree->Branch("t_NJetsPt30Eta5p0", &t_NJetsPt30Eta5p0, "t_NJetsPt30Eta5p0/I"); 
	tree->Branch("t_NJetsPt50Eta2p5", &t_NJetsPt50Eta2p5, "t_NJetsPt50Eta2p5/I"); 
	tree->Branch("t_NJetsPt50Eta5p0", &t_NJetsPt50Eta5p0, "t_NJetsPt50Eta5p0/I");

	tree->Branch("t_PFht",      &t_PFht,      "t_PFht/D");
	tree->Branch("t_PFmht",     &t_PFmht,     "t_PFmht/D");

	t_genJetPt  = new std::vector<double>();
	t_genJetEta = new std::vector<double>();    
	t_genJetPhi = new std::vector<double>();
	t_genJetE   = new std::vector<double>();
	tree->Branch("t_genJetPt",  "vector<double>", &t_genJetPt );
	tree->Branch("t_genJetEta", "vector<double>", &t_genJetEta);
	tree->Branch("t_genJetPhi", "vector<double>", &t_genJetPhi);
	tree->Branch("t_genJetE",   "vector<double>", &t_genJetE  );

	t_genParPt  = new std::vector<double>();
	t_genParEta = new std::vector<double>();    
	t_genParPhi = new std::vector<double>();
	t_genParE   = new std::vector<double>();
	t_genParStatus  = new std::vector<double>();
	t_genParID  = new std::vector<double>();
	tree->Branch("t_genParPt",  "vector<double>", &t_genParPt );
	tree->Branch("t_genParEta", "vector<double>", &t_genParEta);
	tree->Branch("t_genParPhi", "vector<double>", &t_genParPhi);
	tree->Branch("t_genParE",   "vector<double>", &t_genParE  );
	tree->Branch("t_genParStatus",  "vector<double>", &t_genParStatus );
	tree->Branch("t_genParID",  "vector<double>", &t_genParID );


	tree->Branch("t_beamHaloFilter", &t_beamHaloFilter, "t_beamHaloFilter/I");
	tree->Branch("t_eeBadScFilter", &t_eeBadScFilter, "t_eeBadScFilter/I");
	tree->Branch("t_eeNoiseFilter", &t_eeNoiseFilter, "t_eeNoiseFilter/I");
	tree->Branch("t_greedyMuons", &t_greedyMuons, "t_greedyMuons/I");
	tree->Branch("t_hcalLaserEventFilter", &t_hcalLaserEventFilter, "t_hcalLaserEventFilter/I");
	tree->Branch("t_inconsistentMuons", &t_inconsistentMuons, "t_inconsistentMuons/I");
	tree->Branch("t_ra2EcalBEFilter", &t_ra2EcalBEFilter, "t_ra2EcalBEFilter/I");
	tree->Branch("t_ra2EcalTPFilter", &t_ra2EcalTPFilter, "t_ra2EcalTPFilter/I");
	tree->Branch("t_trackingFailureFilter", &t_trackingFailureFilter, "t_trackingFailureFilter/I");
	tree->Branch("t_HBHENoiseFilterRA2", &t_HBHENoiseFilterRA2, "t_HBHENoiseFilterRA2/I");
	tree->Branch("t_ecalLaserCorrFilter", &t_ecalLaserCorrFilter, "t_ecalLaserCorrFilter/I");

	t_firedTrigs  = new std::vector<string>();
	tree->Branch("t_firedTrigs",  "vector<string>", &t_firedTrigs );
	t_firedTrigsPrescale   = new std::vector<double>();
	tree->Branch("t_firedTrigsPrescale",  "vector<double>", &t_firedTrigsPrescale );

}

void RA2Example::clearTreeVectors() {

	t_PFJetPt->clear();
	t_PFJetEta->clear();    
	t_PFJetPhi->clear();
	t_PFJetE->clear();
	t_PFJetBTag->clear();
	t_PFJetNHF->clear();
	t_PFJetEMF->clear();
	t_PFJetHOEne->clear();

	t_genJetPt->clear();
	t_genJetEta->clear();    
	t_genJetPhi->clear();
	t_genJetE->clear();

	t_genParPt->clear();
	t_genParEta->clear();    
	t_genParPhi->clear();
	t_genParE->clear();
	t_genParStatus->clear();
	t_genParID->clear();

	t_firedTrigs->clear();
	t_firedTrigsPrescale->clear();
}


// Get weight factor dependent on number of
// added PU interactions
// --------------------------------------------------
double RA2Example::getPUWeight(int npu) const {
   double w = 1.;
   if (npu < static_cast<int> (_puWeigths.size())) {
      w = _puWeigths.at(npu);
   } else {
      std::cerr << "WARNING in WeightProcessor::getPUWeight: Number of PU vertices = " << npu
            << " out of histogram binning." << std::endl;
   }

   return w;
}

// Generate weights for given data PU distribution
// Scenarios from: https://twiki.cern.ch/twiki/bin/view/CMS/Pileup_MC_Gen_Scenarios
// Code adapted from: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupReweighting
// --------------------------------------------------
std::vector<double> RA2Example::generateWeights(const TH1* data_npu_estimated) const
{

  unsigned int nMaxPU = 0;
  double *npuProbs = 0;

    nMaxPU = 60;
    double npuSummer12_S10[60] = {
      2.560E-06,
      5.239E-06,
      1.420E-05,
      5.005E-05,
      1.001E-04,
      2.705E-04,
      1.999E-03,
      6.097E-03,
      1.046E-02,
      1.383E-02,
      1.685E-02,
      2.055E-02,
      2.572E-02,
      3.262E-02,
      4.121E-02,
      4.977E-02,
      5.539E-02,
      5.725E-02,
      5.607E-02,
      5.312E-02,
      5.008E-02,
      4.763E-02,
      4.558E-02,
      4.363E-02,
      4.159E-02,
      3.933E-02,
      3.681E-02,
      3.406E-02,
      3.116E-02,
      2.818E-02,
      2.519E-02,
      2.226E-02,
      1.946E-02,
      1.682E-02,
      1.437E-02,
      1.215E-02,
      1.016E-02,
      8.400E-03,
      6.873E-03,
      5.564E-03,
      4.457E-03,
      3.533E-03,
      2.772E-03,
      2.154E-03,
      1.656E-03,
      1.261E-03,
      9.513E-04,
      7.107E-04,
      5.259E-04,
      3.856E-04,
      2.801E-04,
      2.017E-04,
      1.439E-04,
      1.017E-04,
      7.126E-05,
      4.948E-05,
      3.405E-05,
      2.322E-05,
      1.570E-05,
      5.005E-06};
    npuProbs = npuSummer12_S10;

  std::vector<double> result(nMaxPU);
  double s = 0.0;
  for(unsigned int npu = 0; npu < nMaxPU; ++npu) {
    double npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));
    result[npu] = npu_estimated / npuProbs[npu];
    s += npu_estimated;
  }
  // normalize weights such that the total sum of weights over thw whole sample is 1.0, i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
  for (unsigned int npu = 0; npu < nMaxPU; ++npu) {
    result[npu] /= s;
  }

  return result;
}



#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(RA2Example);

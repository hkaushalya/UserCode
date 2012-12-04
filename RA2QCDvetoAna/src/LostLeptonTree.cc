// -*- C++ -*-
//
// Package:    LostLeptonTree
// Class:      LostLeptonTree
// 
/**\class LostLeptonTree LostLeptonTree.cc SusyAnalysis/LostLeptonTree/src/LostLeptonTree.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Samantha Hewamanage
// $Id: LostLeptonTree.cc,v 1.1 2012/11/21 23:35:31 samantha Exp $
//
//


// system include files
#include <memory>
#include <algorithm>
#include <vector>

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


class LostLeptonTree : public edm::EDAnalyzer {

  public:

  explicit LostLeptonTree(const edm::ParameterSet & iConfig);
  ~LostLeptonTree();

  private:

  void beginJob() ;
  void endJob() ;
  void analyze(const edm::Event&, const edm::EventSetup&);
  void BookHistograms();
  void clearTreeVectors();  

  bool debug_;
  edm::InputTag vtxSrc_, HBHENoiseFiltSrc_;
  bool doPUReWeight_;
  int doEventWeighing_;
  edm::InputTag puWeigthSrc_, puWeigthABSrc_, puWeigthABCSrc_;
  edm::InputTag pfMetSrc_;
  edm::InputTag jetAllsrc_, genjetAllsrc_;
  edm::InputTag mhtSrc_, htSrc_;
  double minPFJetPt_, minGENJetPt_; 

  edm::Service<TFileService> fs;
  TTree* tree;
  unsigned int         t_EvtRun, t_EvtLS, t_EvtEvent;
  int                  t_NVertices;
  double               t_PUWeight, t_PUWeightAB, t_PUWeightABC;
  double               t_EvtWeight; //for flat sample weighing
  double               t_PFMetPx, t_PFMetPy;
  std::vector<double> *t_PFJetPt,  *t_PFJetEta,  *t_PFJetPhi,  *t_PFJetE,  *t_PFJetBTag;
  double               t_PFht, t_PFmht;
  int                  t_NJetsPt30Eta2p5, t_NJetsPt30Eta5p0, t_NJetsPt50Eta2p5, t_NJetsPt50Eta5p0;

  std::vector<double > *t_genJetPt, *t_genJetEta, *t_genJetPhi, *t_genJetE;
  double               t_minPFJetPt, t_minGenJetPt;
  int 					  t_allFilters;
  int                  t_beamHaloFilter, t_eeBadScFilter, t_eeNoiseFilter, t_greedyMuons;
  int 					  t_hcalLaserEventFilter, t_inconsistentMuons, t_ra2EcalBEFilter;
  int 					  t_ra2EcalTPFilter, t_trackingFailureFilter, t_HBHENoiseFilterRA2;

};


LostLeptonTree::LostLeptonTree(const edm::ParameterSet & iConfig) {
  debug_       = iConfig.getParameter<bool>("Debug");
  vtxSrc_      = iConfig.getParameter<edm::InputTag>("VertexSource");
  doPUReWeight_= iConfig.getParameter<bool>("DoPUReweight");
  doEventWeighing_ = iConfig.getUntrackedParameter<int>("ApplyEventWeighing",0);
  puWeigthSrc_ = iConfig.getParameter<edm::InputTag>("PUWeigthSource");
  puWeigthABSrc_ = iConfig.getParameter<edm::InputTag>("PUWeigthABSource");
  puWeigthABCSrc_ = iConfig.getParameter<edm::InputTag>("PUWeigthABCSource");
  pfMetSrc_    = iConfig.getParameter<edm::InputTag>("PFMetSource");
  jetAllsrc_   = iConfig.getParameter<edm::InputTag>("JetAllSource");
  genjetAllsrc_   = iConfig.getParameter<edm::InputTag>("genJetAllSource");
  mhtSrc_      = iConfig.getParameter<edm::InputTag>("MHTSource");
  htSrc_       = iConfig.getParameter<edm::InputTag>("HTSource");
  minPFJetPt_   = iConfig.getUntrackedParameter<double>("MinPFJetPt",0.0);
  minGENJetPt_   = iConfig.getUntrackedParameter<double>("MinGENJetPt",0.0);
  HBHENoiseFiltSrc_ = iConfig.getParameter<edm::InputTag>("HBHENoiseFiltSrc");
}

LostLeptonTree::~LostLeptonTree() {
}


void LostLeptonTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  using namespace edm;
  using namespace std;

  clearTreeVectors();

  // fill event ID 
  t_EvtRun   = iEvent.id().run();
  t_EvtLS    = iEvent.luminosityBlock();
  t_EvtEvent = iEvent.id().event();
  
  // save number of vertices and position of primary vertex
  edm::Handle< std::vector<reco::Vertex> > vertices;
  iEvent.getByLabel(vtxSrc_, vertices);
  if(vertices->size()<1) std::cout << "No vertices are reconstructed - check StdCleaning ?" << std::endl;

  t_NVertices = vertices->size();
  
  // save pileup weight
  double pu_event_wt = 1.0;
  double pu_event_wtAB = 1.0;
  double pu_event_wtABC = 1.0;
  edm::Handle<double> puweight;
  edm::Handle<double> puweightAB;
  edm::Handle<double> puweightABC;
  if( doPUReWeight_ ) {
    iEvent.getByLabel(puWeigthSrc_, puweight);
    pu_event_wt = *puweight;
    iEvent.getByLabel(puWeigthABSrc_, puweightAB);
    pu_event_wtAB = *puweightAB;
    iEvent.getByLabel(puWeigthABCSrc_, puweightABC);
    pu_event_wtABC = *puweightABC;
  }
  t_PUWeight = pu_event_wt;
  t_PUWeightAB = pu_event_wtAB;
  t_PUWeightABC = pu_event_wtABC;

  if(debug_)std::cout << "t_NVertices " << t_NVertices << "  t_PUWeight " << t_PUWeight << std::endl;


	//event weights for flat QCD samples
	double storedWeight = 1;
	t_EvtWeight = storedWeight;
	if ( doEventWeighing_ )
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
  
  //========= save gen info
  //edm::Handle<edm::View<pat::Jet> > genjetsAll;
	//vector<reco::GenJet>                  "ak5GenJets"  
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

	//get cleaning filter status
  edm::Handle<bool> beamHaloVal, eeBadVal, eeNoiseVal, greedyMuonsVal;
  edm::Handle<bool> hcalLaserEventVal, inconsMuonsVal, ra2EcalBEVal;
  edm::Handle<bool> ra2EcalTPVal, trackFailureVal, HBHENoiseVal;
  
  iEvent.getByLabel("beamHaloFilter", beamHaloVal);
  iEvent.getByLabel("eeBadScFilter", eeBadVal);
  iEvent.getByLabel("eeNoiseFilter", eeNoiseVal);
  iEvent.getByLabel("greedyMuons", greedyMuonsVal);
  iEvent.getByLabel("hcalLaserEventFilter", hcalLaserEventVal);
  iEvent.getByLabel("inconsistentMuons", inconsMuonsVal);
  iEvent.getByLabel("ra2EcalBEFilter", ra2EcalBEVal);
  iEvent.getByLabel("ra2EcalTPFilter", ra2EcalTPVal);
  iEvent.getByLabel("trackingFailureFilter", trackFailureVal);
  //this is not found in 52X skims
  //iEvent.getByLabel(HBHENoiseFiltSrc_, HBHENoiseVal);

  t_beamHaloFilter  = *beamHaloVal;
  t_eeBadScFilter   = *eeBadVal;
  t_eeNoiseFilter   = *eeNoiseVal;
  t_greedyMuons     = *greedyMuonsVal;
  t_hcalLaserEventFilter = *hcalLaserEventVal;
  t_inconsistentMuons = *inconsMuonsVal;
  t_ra2EcalBEFilter   = *ra2EcalBEVal;
  t_ra2EcalTPFilter   = *ra2EcalTPVal;
  t_trackingFailureFilter = *trackFailureVal; 
  //t_HBHENoiseFilterRA2 = *HBHENoiseVal;

  tree->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void LostLeptonTree::beginJob() {

  BookHistograms();

}

// ------------ method called once each job just after ending the event loop  ------------
void LostLeptonTree::endJob() {

}

void LostLeptonTree::BookHistograms() {

  // book tree here
  tree = fs->make<TTree>("tree", "tree");
  tree->SetAutoSave(10000);

  tree->Branch("t_EvtRun",    &t_EvtRun,   "t_EvtRun/i");
  tree->Branch("t_EvtLS",     &t_EvtLS,    "t_EvtLS/i");
  tree->Branch("t_EvtEvent",  &t_EvtEvent, "t_EvtEvent/i");
  tree->Branch("t_NVertices", &t_NVertices,"t_NVertices/I");
  tree->Branch("t_PUWeight",  &t_PUWeight, "t_PUWeight/D");
  tree->Branch("t_PUWeightAB",  &t_PUWeightAB, "t_PUWeightAB/D");
  tree->Branch("t_PUWeightABC",  &t_PUWeightABC, "t_PUWeightABC/D");
  tree->Branch("t_EvtWeight",  &t_EvtWeight, "t_EvtWeight/D");
  tree->Branch("t_PFMetPx",   &t_PFMetPx,  "t_PFMetPx/D");
  tree->Branch("t_PFMetPy",   &t_PFMetPy,  "t_PFMetPy/D");

  t_PFJetPt  = new std::vector<double>();
  t_PFJetEta = new std::vector<double>();    
  t_PFJetPhi = new std::vector<double>();
  t_PFJetE   = new std::vector<double>();
  tree->Branch("t_PFJetPt",  "vector<double>", &t_PFJetPt );
  tree->Branch("t_PFJetEta", "vector<double>", &t_PFJetEta);
  tree->Branch("t_PFJetPhi", "vector<double>", &t_PFJetPhi);
  tree->Branch("t_PFJetE",   "vector<double>", &t_PFJetE);
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

  tree->Branch("t_beamHaloFilter", &t_beamHaloFilter, "t_beamHaloFilter/I");
  tree->Branch("t_eeBadScFilter", &t_eeBadScFilter, "t_eeBadScFilter/I");
  tree->Branch("t_eeNoiseFilter", &t_eeNoiseFilter, "t_eeNoiseFilter/I");
  tree->Branch("t_greedyMuons", &t_greedyMuons, "t_greedyMuons/I");
  tree->Branch("t_hcalLaserEventFilter", &t_hcalLaserEventFilter, "t_hcalLaserEventFilter/I");
  tree->Branch("t_inconsistentMuons", &t_inconsistentMuons, "t_inconsistentMuons/I");
  tree->Branch("t_ra2EcalBEFilter", &t_ra2EcalBEFilter, "t_ra2EcalBEFilter/I");
  tree->Branch("t_ra2EcalTPFilter", &t_ra2EcalTPFilter, "t_ra2EcalTPFilter/I");
  tree->Branch("t_trackingFailureFilter", &t_trackingFailureFilter, "t_trackingFailureFilter/I");
//  tree->Branch("t_HBHENoiseFilterRA2", &t_HBHENoiseFilterRA2, "t_HBHENoiseFilterRA2/I");


}

void LostLeptonTree::clearTreeVectors() {

  t_PFJetPt  ->clear();
  t_PFJetEta ->clear();    
  t_PFJetPhi ->clear();
  t_PFJetE   ->clear();

  t_genJetPt  ->clear();
  t_genJetEta ->clear();    
  t_genJetPhi ->clear();
  t_genJetE   ->clear();

}



#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(LostLeptonTree);

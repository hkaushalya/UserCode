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
// $Id: RA2Example.cc,v 1.1 2013/07/26 19:04:37 samantha Exp $
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

		bool debug_;
		edm::InputTag vtxSrc_, HBHENoiseFiltSrc_;
		edm::InputTag firedTrigNamesSrc_, firedTrigPrescaleSrc_;
		bool doPUReWeight_;
		int doEventWeighing_;
		edm::InputTag puWeigthSrc_, puWeigthABSrc_, puWeigthABCSrc_, puWeigthRA2Src_;
		edm::InputTag pfMetSrc_;
		edm::InputTag jetAllsrc_, chsjetAllsrc_, genjetAllsrc_ ,genParticleSrc_;
		edm::InputTag mhtSrc_, htSrc_;
		bool mcFlag_;
		double minPFJetPt_, minGENJetPt_; 
	   std::vector<double> genParticleList_;
		TH1F *hist_nvtx, *hist_njet, *hist_jetpt;

		edm::Service<TFileService> fs;
};


RA2Example::RA2Example(const edm::ParameterSet & iConfig) {

	//this info comes from the pythos/ra2example_cfi
	//(or from the runRA2Example_condor_cfg.py)
	debug_           = iConfig.getParameter<bool>("Debug");
	vtxSrc_          = iConfig.getParameter<edm::InputTag>("VertexSource");
	pfMetSrc_        = iConfig.getParameter<edm::InputTag>("PFMetSource");
	jetAllsrc_       = iConfig.getParameter<edm::InputTag>("JetAllSource");
	chsjetAllsrc_    = iConfig.getParameter<edm::InputTag>("chsJetAllSource");
	genjetAllsrc_    = iConfig.getParameter<edm::InputTag>("genJetAllSource");
	genParticleSrc_  = iConfig.getParameter<edm::InputTag>("genParticleSource");
	genParticleList_ = iConfig.getParameter<std::vector<double> > ("genParticleList");
	mhtSrc_          = iConfig.getParameter<edm::InputTag>("MHTSource");
	mcFlag_          = iConfig.getParameter<bool>("MCflag");
	htSrc_           = iConfig.getParameter<edm::InputTag>("HTSource");
	minPFJetPt_      = iConfig.getUntrackedParameter<double>("MinPFJetPt",0.0);
	minGENJetPt_     = iConfig.getUntrackedParameter<double>("MinGENJetPt",0.0);
}

RA2Example::~RA2Example() {
}


void RA2Example::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

	using namespace edm;
	using namespace std;

	// event ID 
	const double t_EvtRun   = iEvent.id().run();
	const double t_EvtLS    = iEvent.luminosityBlock();
	const double t_EvtEvent = iEvent.id().event();

	cout << "========================== run:lumi:evt=" <<  t_EvtRun 
			<< ":" <<  t_EvtLS << ":" << t_EvtEvent << endl;

	//primary vertex
	edm::Handle< std::vector<reco::Vertex> > vertices;
	iEvent.getByLabel(vtxSrc_, vertices);
	const double t_NVertices = vertices->size();
	if(debug_)std::cout << "t_NVertices " << t_NVertices << std::endl;
	hist_nvtx->Fill(t_NVertices);


	// missing transverse energy
	edm::Handle < std::vector<pat::MET> > met;
	iEvent.getByLabel(pfMetSrc_, met);
	for (std::vector<pat::MET>::const_iterator it = met->begin(); it != met->end(); ++it) {
		const double t_PFMetPx = it->px();
		const double t_PFMetPy = it->py();
	}

	// jets information
	edm::Handle<edm::View<pat::Jet> > jetsAll;
	iEvent.getByLabel(jetAllsrc_, jetsAll);

	//std::cout << "jetsAll->size() "<< jetsAll->size() << std::endl;
	hist_njet->Fill(jetsAll->size());
	for (unsigned int ijet=0; ijet<jetsAll->size(); ijet++) 
	{
		if ( (*jetsAll)[ijet].pt() > minPFJetPt_ )
		{
			const double t_PFJetPt = (*jetsAll)[ijet].pt();
			hist_jetpt->Fill(t_PFJetPt);

			cout << __FUNCTION__ << ":" << __LINE__ << ":: reco jet [" << ijet << "]=" << t_PFJetPt  << endl;

			const double t_PFJetEta   = (*jetsAll)[ijet].eta();    
			const double t_PFJetPhi   = (*jetsAll)[ijet].phi();
			const double t_PFJetE     = (*jetsAll)[ijet].energy();
			const double t_PFJetNHF   = (*jetsAll)[ijet].neutralHadronEnergyFraction();
			const double t_PFJetEMF   = (*jetsAll)[ijet].photonEnergyFraction();

			/****************************************************
			 * This is how you can access the PFJet constituents
			 * As and eg. I am calculating total energy deposited 
			 * in HO for each jet
			 ****************************************************/
			const std::vector<reco::PFCandidatePtr> & pfcands = (*jetsAll)[ijet].getPFConstituents();
			double hoEne   = 0.0;
			std::vector<reco::PFCandidatePtr>::const_iterator itr;
			for (itr = pfcands.begin(); itr != pfcands.end(); itr++) 
			{
				hoEne   += (*itr)->hoEnergy();
			}
		}
	}

	// save ra2 ht/mht mainly for debugging only
	edm::Handle<double> ht;
	iEvent.getByLabel(htSrc_, ht);
	const double t_PFht  = *ht;

	edm::Handle<edm::View<reco::MET> > mht;
	iEvent.getByLabel(mhtSrc_, mht);
	const double t_PFmht = (*mht)[0].pt() ;

	if(debug_)std::cout << "ht " << t_PFht << " mht " << t_PFmht << std::endl;

	/*******************************************************
	 * Accessing GEN-JETS and GEN-PARTICLES
	 ******************************************************/
	if (mcFlag_)
	{
		edm::Handle<edm::View<reco::GenJet> > genjetsAll;
		iEvent.getByLabel(genjetAllsrc_, genjetsAll);

		//std::cout << "genjetsAll->size() "<< genjetsAll->size() << std::endl;
		for (unsigned int ijet=0; ijet<genjetsAll->size(); ijet++)
		{
			if ( (*genjetsAll)[ijet].pt() > minGENJetPt_ )
			{
				const double t_genJetPt = (*genjetsAll)[ijet].pt();
				const double t_genJetEta = (*genjetsAll)[ijet].eta();    
				const double t_genJetPhi = (*genjetsAll)[ijet].phi();
				const double t_genJetE   = (*genjetsAll)[ijet].energy();
				cout << __FUNCTION__ << ":" << __LINE__ << ":: gen jet [" << ijet << "]=" << t_genJetPt  << endl;
			}
		} //end gen-jet loop

		//get gen particle info
		edm::Handle<std::vector<reco::GenParticle > > genParticles;
		iEvent.getByLabel(genParticleSrc_, genParticles);

		for (unsigned ig=0; ig<genParticles->size(); ig++) 
		{
			const reco::GenParticle& gen = genParticles->at(ig);
			const int pdgId = gen.pdgId();
			if (find(genParticleList_.begin(), genParticleList_.end(),(double) pdgId) == genParticleList_.end()) continue; 
			const int pdgStatus = gen.status();
			/*	cout << __FUNCTION__ << ":Gen par ig=id/stat/pt/eta/phi/e =\t"
				<< ig << "=\t" << pdgId << "\t" << pdgStatus 
				<< "\t" << gen.pt() 
				<< "\t" << gen.eta() 
				<< "\t" << gen.phi() 
				<< "\t" << gen.energy()
				<< endl;
				*/		
			const double t_genParPt     = gen.pt();
			const double t_genParEta    = gen.eta();    
			const double t_genParPhi    = gen.phi();
			const double t_genParE      = gen.energy();
			const double t_genParStatus = (double)pdgStatus;
			const double t_genParID     = (double)pdgId;

		} //end genParticle loop

	} //end mcflag

}

// ------------ method called once each job just before starting event loop  ------------
void RA2Example::beginJob() {

	BookHistograms();

}

// ------------ method called once each job just after ending the event loop  ------------
void RA2Example::endJob() {


}

void RA2Example::BookHistograms() {

	hist_nvtx  = fs->make<TH1F>("nvtx", "nvertices",20,0,20);
	hist_njet  = fs->make<TH1F>("njet", "njet",20,0,20);
	hist_jetpt = fs->make<TH1F>("jetpt", "jet pt",100,0,1000);
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(RA2Example);

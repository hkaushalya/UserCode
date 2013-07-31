// -*- C++ -*-
//
// Package:    MuonEffRecoId
// Class:      MuonEffRecoId
// 
/**\class MuonEffRecoId MuonEffRecoId.cc UserCode/MuonEff/src/MuonEffRecoId.cc

 Description: Studies the DelPhi(jet,MET) for QCD veto.

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  samantha hewamanage
//         Created:  Tue Jul 26 22:15:44 CDT 2011
// $Id: MuonEffRecoId.cc,v 1.6 2012/06/07 18:31:25 samantha Exp $
//
//


// system include files
#include <memory>
#include <vector>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

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
#include "DataFormats/PatCandidates/interface/Muon.h"

//jet collections
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

//ROOT
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "assert.h"
#include "TVector3.h"
#include <iomanip>
#include "TPad.h"
#include <utility>
#include "TCanvas.h"
#include "TDirectory.h"

using namespace std;

const static unsigned nJetsToPlot = 6;
//
// class declaration
//

void Print(const TLorentzVector& v)
{
	cout << "pt/eta/phi/e=" << v.Pt() << "/" << v.Eta() << "/" << v.Phi() << "/" << v.E() << endl;
}

class MuonEffRecoId : public edm::EDAnalyzer {
   public:
      explicit MuonEffRecoId(const edm::ParameterSet&);
      ~MuonEffRecoId();

//      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      void beginJob() ;
      void analyze(const edm::Event&, const edm::EventSetup&);
      void endJob() ;
      
      bool beginRun(edm::Run&, edm::EventSetup const&);
      bool endRun(edm::Run&, edm::EventSetup const&);
      bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
		unsigned GetVectorIndex(const vector<pair<double,double> > & binEdges, const double& val);

      // ----------member data ---------------------------
		//jet collections
		edm::InputTag mhtInputTag_, htInputTag_, metInputTag_, recoMuonSrc_, recoisoMuonSrc_ ,genParticleSrc_;
		edm::Handle<edm::View<reco::MET> > mhtHandle;
		edm::Handle<double> htHandle;
		edm::Handle<std::vector<reco::PFMET> >pfMetHandle;

		struct RunLumiEvt_t
		{
			unsigned run;
			unsigned lumi;
			unsigned evt;
		};

		
		int inMinVtx;
		unsigned int iProcessed; // number of processed events
		unsigned int iPassed; //number of events passed the filter
		int inMinNdofVtx; //minimum number of dof for the vtx
		RunLumiEvt_t kRunLumiEvent;
		int iVerbose; // control print levels
		double dMaxPrimVtxZ;
		double dMaxPrimVtxRho;
		double dMinHT, dMinMHt;
		double dMinJetEt4MHt, dMaxJetEta4MHt;

		struct EventHist_t {
			TH1F* nvtx;
			TH1F* mht;
			TH1F* ht;
			TH1F* njet_et30eta5p0;
			TH1F* njet_et50eta2p5;
			TH1F* nRecoMuons;
			TH1F* nGenMuons;
			TH1F* matched_MuPtRatio;  //mathced reco/gen muon info
			TH1F* matched_MuERatio;
			TH1F* matched_MuEDiff;
			TH1F* matched_MuDelR;
			TH1F* matched_MuDelEta;
			TH2D* muEff_reco_pteta_num; //intermidiate plot
			TH2D* muEff_gen_pteta_num; //intermidiate plot
			TH2D* muEff_recoiso_pteta_num; //intermidiate plot
//			TH2D* muEff_pteta;		//end product// this cannot be just num/den if I am spliiting the jobs, but the average will be.
			TH2D* allgenmuon_pteta;
			TH2D* allrecomuon_pteta;
			TH1D* allmatchingmuon_delr;
		};

		struct JetHist_t{
			TH1F* pt;
			TH1F* eta;
			TH1F* phi;
			TH1F* delphi;  //delphi(jet,MHT)
		};

		struct Hist_t {
			EventHist_t evt;
			vector<JetHist_t> jet;
		};


		vector <pair<double, double> > vJetBins_, vHtBins_;
		vector<vector<Hist_t> > Hist;
		
		edm::InputTag recoJetAllSource_;
//		edm::InputTag patJetsPFPt50Eta25InputTag_;
		edm::Handle<std::vector<pat::Jet> > recoJetHandle;
		edm::Handle<std::vector<pat::Jet> > pfpt50Eta25JetHandle;

		void PrintHeader();
		TLorentzVector vMHtVec;

};


//
// constructors and destructor
//
MuonEffRecoId::MuonEffRecoId(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
	recoMuonSrc_          = iConfig.getParameter<edm::InputTag>("recoMuonSource");
	recoisoMuonSrc_          = iConfig.getParameter<edm::InputTag>("recoisoMuonSource");
	genParticleSrc_   = iConfig.getParameter<edm::InputTag>("genParticleSource");
	recoJetAllSource_ = iConfig.getParameter<edm::InputTag>("recoJetAllSource");
	htInputTag_       = iConfig.getParameter<edm::InputTag>("HTSource");
	dMinHT            = iConfig.getUntrackedParameter<double>("MinHT",0.0);
	iVerbose          = iConfig.getUntrackedParameter<int>("verbose",0);
	iProcessed = 0;
	iPassed    = 0;
	inMinVtx   = 1;
	vMHtVec.SetPxPyPzE(0,0,0,0);
	
	//for now
	vJetBins_.push_back(make_pair(3,1000));
	vJetBins_.push_back(make_pair(2,2));
	vJetBins_.push_back(make_pair(3,5));
	vJetBins_.push_back(make_pair(6,7));
	vJetBins_.push_back(make_pair(8,1000));
	vHtBins_.push_back(make_pair(500,8000));
}


MuonEffRecoId::~MuonEffRecoId()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
void
MuonEffRecoId::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
	++iProcessed;
	kRunLumiEvent.run  = iEvent.id().run();
	kRunLumiEvent.lumi = iEvent.id().luminosityBlock(); 
	kRunLumiEvent.evt  = iEvent.id().event();

	
	if (iVerbose) PrintHeader();

	Handle<reco::VertexCollection> vertexHandle;
	iEvent.getByLabel("goodVertices", vertexHandle);
   if (! vertexHandle.isValid())
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":vertexHandle handle not found!" << std::endl;
		assert(false);
	}

	iEvent.getByLabel(htInputTag_, htHandle);
	if (! htHandle.isValid()) 
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":htHandle handle not found!" << std::endl;
		assert(false);
	}
	
	iEvent.getByLabel(recoJetAllSource_, recoJetHandle);
	if (! recoJetHandle.isValid()) 
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":recoJetHandle handle not found!" << std::endl;
		assert(false);
	}

	/*************
	 *before any cuts fill the hists to check
	 * PAT contents
	 *************/

	const double ht     = (*htHandle);
	const double njet   = recoJetHandle->size();
	const double nvtx = vertexHandle->size();

	edm::Handle<edm::View<pat::Muon> > reco_muons; //muons that are reconstructed to measure reco eff
	iEvent.getByLabel(recoMuonSrc_, reco_muons);

	edm::Handle<edm::View<pat::Muon> > recoiso_muons;  //muons reconstructed and passed isolation: to measure iso eff
	iEvent.getByLabel(recoisoMuonSrc_, recoiso_muons);

   edm::Handle<std::vector<reco::GenParticle > > genParticles;
	iEvent.getByLabel(genParticleSrc_, genParticles);

	const int muonPDGcode = 13;
	const double minMuonPt = 10.0;
	const int wbosonPDGcode = 24;
  //if (nvtx > 0 && recoiso_muons->size()==1) 
  if (nvtx > 0) 
  {
	  bool foundGenMuon = false; 
	  bool foundRecoMatch = false;
	  TLorentzVector matched_genMuonVec(0,0,0,0), matched_recoMuonVec(0,0,0,0);
	  for (unsigned ig=0; ig<genParticles->size(); ig++) 
	  {
		  const reco::GenParticle& gen = genParticles->at(ig);
		  const int pdgId     = gen.pdgId();
		  const int pdgStatus = gen.status();
		  const double genPt  = gen.pt();
		  const double genE  = gen.energy();

		  bool hasWmother = false;
		  double mompt = 0;
		  double mome = 0;
		  for (unsigned idx=0; idx< gen.numberOfMothers(); ++idx)
		  {
			  if ( std::abs(gen.mother(idx)->pdgId()) == wbosonPDGcode )
			  {
				  hasWmother = true;
				  mompt = gen.mother(idx)->pt();
				  mome = gen.mother(idx)->energy();
				  break;
			  }
		  }
		  
		  if (fabs(pdgId) != muonPDGcode || genPt<minMuonPt || ! hasWmother) continue;
		  //cout << "gen[pdg/pt/e] = " << pdgId << "/" << genPt  << "/" << genE << endl;
		  //cout << "mom[pdg/pt/e] = " << wbosonPDGcode << "/" << mompt << "/" << mome << endl;

		  foundGenMuon = true;  //found a gen muon from a W decay

		  TLorentzVector genmuVec(0,0,0,0);
		  genmuVec.SetPxPyPzE(gen.px(), gen.py(), gen.pz(), gen.energy());
		  matched_genMuonVec = genmuVec; //this is fine as long as I am exiting the loops after the first match!!!

		  Hist.at(0).at(0).evt.allgenmuon_pteta->Fill(genmuVec.Pt(), genmuVec.Eta());

			/* try to find a matching RECO muon */
		  for (edm::View<pat::Muon>::const_iterator m = reco_muons->begin(); m != reco_muons->end(); ++m) 
		  {
			  TLorentzVector recomuVec(0,0,0,0);
			  recomuVec.SetPxPyPzE(m->px(), m->py(), m->pz(), m->energy());
			  const double dR = genmuVec.DeltaR(recomuVec);
			  Hist.at(0).at(0).evt.allmatchingmuon_delr->Fill(dR);
			  Hist.at(0).at(0).evt.allrecomuon_pteta->Fill(recomuVec.Pt(), recomuVec.Eta());

			  if (dR>0.3) continue;
			  foundRecoMatch = true;
			  matched_recoMuonVec = recomuVec;

			  /*  cout << __FUNCTION__ << ":Gen par ig=id/stat/pt/eta/phi/e =\t"
					<< ig << "=\t" << pdgId << "\t" << pdgStatus 
					<< "\t" << gen.pt() 
					<< "\t" << gen.eta() 
					<< "\t" << gen.phi() 
					<< "\t" << gen.energy()
					<< endl;
					*/
			  break; //one match
		  }

		  break; //looking for only one leading muon
	  } //end genparticle loop


		/* now see if the recoMuon can be matched to a iso muon */
		/* should I consider only event with a gen muon for this? */
		/* should I just loop over all */ 
		
		TLorentzVector matched_recoIsoMuonVec(0,0,0,0);
		bool foundRecoIsoMatch = false;
	  /* try to find a matching RECO muon to Iso muon*/
	  for (edm::View<pat::Muon>::const_iterator m = recoiso_muons->begin(); m != recoiso_muons->end(); ++m) 
	  {
		  TLorentzVector recoisomuVec(0,0,0,0);
		  recoisomuVec.SetPxPyPzE(m->px(), m->py(), m->pz(), m->energy());
		  const double dR = matched_recoMuonVec.DeltaR(recoisomuVec);

		  if (dR>0.3) continue;
		  foundRecoIsoMatch = true;
		  matched_recoIsoMuonVec = recoisomuVec;
		  break; //one match
	  }


	  //now we have a gen muon and/or matching reco muon
	  //calcualted event kinematics and fill eff hists

	  reco::MET::LorentzVector recomht(0,0,0,0);
	  unsigned reco_njet_pt30eta5p0 = 0, reco_njet_pt50eta2p5 = 0;
	  double reco_ht = 0;

	  //check if the the mathced muon is part of jet collection. I am sure it is.
	  //Just make sure.
	  for(unsigned int ijet=0; ijet<recoJetHandle->size(); ijet++) 
	  {
		  const double pt  = (*recoJetHandle)[ijet].pt();
		  const double eta = std::fabs((*recoJetHandle)[ijet].eta());

		  if( pt>30.0 && fabs(eta)<5.0 )
		  {
			  ++reco_njet_pt30eta5p0;
			  recomht -= (*recoJetHandle)[ijet].p4();
		  }
		  if( pt>50.0 && fabs(eta)<2.5 )
		  {
			  ++reco_njet_pt50eta2p5;
			  reco_ht += pt;
		  }
	  } //end of jets counting

	  //loose jet id stuff
	  /*		if ((*jetHandle)[i].neutralHadronEnergyFraction() >=0.99) continue;
				if ((*jetHandle)[i].neutralEmEnergyFraction() >=0.99) continue;
				if (((*jetHandle)[i].getPFConstituents()).size() <=1) continue;
				if ((*jetHandle)[i].chargedHadronEnergyFraction() <=0) continue;
				if ((*jetHandle)[i].chargedMultiplicity() <=0) continue;
				if ((*jetHandle)[i].chargedEmEnergyFraction() >=0.99) continue;
				*/


	  const unsigned iHtBin  = GetVectorIndex(vHtBins_, reco_ht);
	  const unsigned iJetBin = GetVectorIndex(vJetBins_, reco_njet_pt50eta2p5);
	  //cout << "njet/ht=" << reco_njet_pt50eta2p5 << "/" << reco_ht << endl;
	  //cout << "iHtBin/iJetBin = " << iHtBin << "/" << iJetBin << endl;

	  //discard event if it does not fall into any njet/ht bin
	  if (iHtBin > vHtBins_.size()) return;
	  if (iJetBin > vJetBins_.size()) return;

	  //Print(matched_genMuonVec);
	  //Print(matched_recoMuonVec);
	  Hist.at(iJetBin).at(iHtBin).evt.nvtx->Fill(nvtx);
	  Hist.at(iJetBin).at(iHtBin).evt.mht->Fill(recomht.Pt());
	  Hist.at(iJetBin).at(iHtBin).evt.ht->Fill(reco_ht);
	  Hist.at(iJetBin).at(iHtBin).evt.njet_et50eta2p5->Fill(reco_njet_pt50eta2p5);
	  Hist.at(iJetBin).at(iHtBin).evt.njet_et30eta5p0->Fill(reco_njet_pt30eta5p0);


	  if (foundGenMuon)
	  {
		  const double matched_genMuon_e  = matched_genMuonVec.E();
		  const double matched_genMuon_pt  = matched_genMuonVec.Pt();
		  const double matched_genMuon_eta = matched_genMuonVec.Eta();
		  Hist.at(iJetBin).at(iHtBin).evt.muEff_gen_pteta_num->Fill(matched_genMuon_pt, matched_genMuon_eta);
		  if (foundRecoMatch)
		  {
			  Hist.at(iJetBin).at(iHtBin).evt.muEff_reco_pteta_num->Fill(matched_genMuon_pt, matched_genMuon_eta);

			  const double matched_recoMuon_e  = matched_recoMuonVec.E();
			  const double matched_recoMuon_pt  = matched_recoMuonVec.Pt();
			  const double matched_recoMuon_eta = matched_recoMuonVec.Eta();
			  const double delR = matched_genMuonVec.DeltaR(matched_recoMuonVec);
			  const double pt_ratio = matched_recoMuon_pt/matched_genMuon_pt;
			  const double e_ratio = matched_recoMuon_e/matched_genMuon_e;
			  const double e_diff = matched_genMuon_e - matched_recoMuon_e;
			  Hist.at(iJetBin).at(iHtBin).evt.matched_MuPtRatio->Fill(pt_ratio);
			  Hist.at(iJetBin).at(iHtBin).evt.matched_MuERatio->Fill(e_ratio);
			  Hist.at(iJetBin).at(iHtBin).evt.matched_MuEDiff->Fill(e_diff);
			  Hist.at(iJetBin).at(iHtBin).evt.matched_MuDelR->Fill(delR);
			  //Hist.at(iJetBin).at(iHtBin).evt.matched_MuDelEta->Fill();
			  
			  if (foundRecoIsoMatch)
			  {
				  Hist.at(iJetBin).at(iHtBin).evt.muEff_recoiso_pteta_num->Fill(matched_genMuon_pt, matched_genMuon_eta);
			  }
		  }
	  }

  }//end of nvtx>=1


	  ++iPassed;
}

// ------------ method called once each job just before starting event loop  ------------
void 
MuonEffRecoId::beginJob()
{

	//generate hists
	edm::Service<TFileService> fs;
	TFileDirectory dir = fs->mkdir("Hist");

	for (unsigned jetbin=0; jetbin<vJetBins_.size(); ++jetbin)
	{
		unsigned minjet= vJetBins_.at(jetbin).first;
		unsigned maxjet= vJetBins_.at(jetbin).second;
		vector<Hist_t> vHist;
		for (unsigned htbin=0; htbin<vHtBins_.size(); ++htbin)
		{
			Hist_t hist;
			unsigned minht= vHtBins_.at(htbin).first;
			unsigned maxht= vHtBins_.at(htbin).second;
			stringstream folder, header;
			folder << "Njet" << minjet << "to" << maxjet << "Ht" << minht << "to" << maxht;
			TFileDirectory subDir = dir.mkdir(folder.str());
			//gDirectory->pwd();
			stringstream nvtx_name, mht_name, ht_name, njet_name, nRecoMuons_name, nGenMuons_name;
			stringstream nvtx_title, mht_title, ht_title, njet_title, nRecoMuons_title, nGenMuons_title;
			stringstream njet_et30eta5p0_title, njet_et50eta2p5_title;

			header << "Njets[" << minjet << "-" << maxjet << "], HT[" << minht << "-" << maxht<< "]";
			nvtx_title << header.str() << ";N vtx; Events;";
			ht_title << header.str() << ";HT[GeV]; Events;";
			mht_title << header.str() << ";MHT[GeV]; Events;";
			njet_title << header.str() << ";Njets (all reco jets); Events;";
			njet_et30eta5p0_title << header.str() << ";Njets [P_{T}>30 GeV, |#eta |<5.0]; Events;";
			njet_et50eta2p5_title << header.str() << ";Njets [P_{T}>50 GeV, |#eta |<2.5]; Events;";
			nRecoMuons_title << header.str() << ";N-Reco Muons; Events;";
			nGenMuons_title << header.str() << ";N-Gen Muons; Events;";


			hist.evt.nvtx = subDir.make<TH1F> ("nvtx" , nvtx_title.str().c_str(), 35, 0, 35);
			hist.evt.ht = subDir.make<TH1F> ("ht" , ht_title.str().c_str(), 320, 0, 8000);
			hist.evt.mht = subDir.make<TH1F> ("mht" , mht_title.str().c_str(), 60, 0, 1500);
			hist.evt.njet_et30eta5p0 = subDir.make<TH1F> ("njet_et30eta5p0" , njet_et30eta5p0_title.str().c_str(), 15, 0, 15);
			hist.evt.njet_et50eta2p5 = subDir.make<TH1F> ("njet_et50eta2p5" , njet_et50eta2p5_title.str().c_str(), 15, 0, 15);
//			hist.evt.nRecoMuons = subDir.make<TH1F> ("nrecoMuons" , nRecoMuons_title.str().c_str(), 15, 0, 15);
//			hist.evt.nGenMuons = subDir.make<TH1F> ("nGenMuons" , nGenMuons_title.str().c_str(), 15, 0, 15);

			//const int n_eff_ptbins  = 13; 
			//const double eff_ptbins[]  = {0,20,30,40,50,70,100,150,200,300,500,700,1000,1500};
			const int n_eff_ptbins  = 12; 
			const double eff_ptbins[]  = {0,5,10,15,20,30,40,50,75,100,300,500,1000};
			const int n_eff_etabins = 8;
			const double eff_etabins[] = {0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,5.0};

			stringstream allgenmuon_pteta_title, allrecomuon_pteta_title, allmatching_delr_title;
			allgenmuon_pteta_title << header.str() << ":all gen muons pt/eta before any cuts;Muon PT;Muon #eta;";
			allrecomuon_pteta_title << header.str() << ":all gen muons pt/eta before any cuts;Muon PT;Muon #eta;";
			allmatching_delr_title << header.str() << ":all matching delR before any cuts;DelR;Events;";
			hist.evt.allgenmuon_pteta= subDir.make<TH2D> ("allgenmuon_pteta" , allgenmuon_pteta_title.str().c_str(), 1000,0,1000,100,0,5);
			hist.evt.allrecomuon_pteta= subDir.make<TH2D> ("allrecomuon_pteta" , allrecomuon_pteta_title.str().c_str(), 1000,0,1000,100,0,5);
			hist.evt.allmatchingmuon_delr= subDir.make<TH1D> ("allmatching_delR" , allmatching_delr_title.str().c_str(), 100,0,5);

			stringstream muEff_genmuon_pteta_title, muEff_recomuon_pteta_title, muEff_recoisomuon_pteta_title;
			muEff_genmuon_pteta_title << header.str() << " (Events with Generator Muon from W decay) ;Gen-Muon P_{T}; Gen-Muon #eta ;";
			muEff_recomuon_pteta_title << header.str() << " (Event with a Reco muon matching to a Generator Muon from W decay) ;Gen-Muon P_{T}; Gen-Muon #eta ;";
			muEff_recoisomuon_pteta_title << header.str() << " (Events with Reco muon passing Iso cut;Gen-Muon P_{T}; Gen-Muon #eta;";
			hist.evt.muEff_gen_pteta_num = subDir.make<TH2D> ("muEff_genMuon_den" ,muEff_genmuon_pteta_title.str().c_str(), n_eff_ptbins, eff_ptbins, n_eff_etabins, eff_etabins);
			hist.evt.muEff_reco_pteta_num = subDir.make<TH2D> ("muEff_recoMuon_num" ,muEff_recomuon_pteta_title.str().c_str(), n_eff_ptbins, eff_ptbins, n_eff_etabins, eff_etabins);
			hist.evt.muEff_recoiso_pteta_num = subDir.make<TH2D> ("muEff_recoisoMuon_num" ,muEff_recoisomuon_pteta_title.str().c_str(), n_eff_ptbins, eff_ptbins, n_eff_etabins, eff_etabins);

			stringstream matched_MuPtRatio_title, matched_MuERatio_title, matched_MuEDiff_title, matched_MuDelR_title;
			matched_MuPtRatio_title << header.str() << ": Matched muons pt ratio; P_{T}^{Reco Muon}/P_{T}^{Gen Muon};Events;";
			matched_MuERatio_title << header.str() << ": Matched muons E ratio; E_{T}^{Reco Muon}/E_{T}^{Gen Muon};Events;";
			matched_MuEDiff_title << header.str() << ": Matched muons E difference; |E^{Gen Muon} - E^{Reco Muon}|;Events;";
			matched_MuDelR_title << header.str() << ": Matched muons #Delta R;#Delta R;Events;";

			hist.evt.matched_MuPtRatio = subDir.make<TH1F> ("matched_MuPtRatio" , matched_MuPtRatio_title.str().c_str(), 80, 0, 4);
			hist.evt.matched_MuERatio  = subDir.make<TH1F> ("matched_MuERatio" , matched_MuERatio_title.str().c_str(), 80, 0, 4);
			hist.evt.matched_MuEDiff   = subDir.make<TH1F> ("matched_MuEDiff" , matched_MuEDiff_title.str().c_str(), 50, 0, 100);
			hist.evt.matched_MuDelR    = subDir.make<TH1F> ("matched_MuDelR" , matched_MuDelR_title.str().c_str(), 30, 0, 3);
			//gDirectory->ls();

			vHist.push_back(hist);
		}
		Hist.push_back(vHist);
	}
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonEffRecoId::endJob() {

	std::cout << __FILE__ << "::" << __FUNCTION__  << "\n" << std::endl;
	std::cout << "[ATM:00] Job settings " << std::endl;
	std::cout << "[ATM:01] Events Processed --- = " << iProcessed << std::endl;
	std::cout << "[ATM:02] Events Passed ------ = " << iPassed << std::endl;
	std::cout << "[ATM:10] Primary Vtx Min ---- = " << inMinVtx << std::endl; 
	std::cout << "[ATM:31] Minimum HT --------- = " << dMinHT << std::endl;
	std::cout << "[ATM:32] Minimum MHT -------- = " << dMinMHt << std::endl;
	std::cout << "[ATM:33] Filter Summary------------- " << std::endl;
}

// ------------ method called when starting to processes a run  ------------
bool 
MuonEffRecoId::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
MuonEffRecoId::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

void MuonEffRecoId::PrintHeader()
{
	std::cout  << "======= Run/Event: " << kRunLumiEvent.run
			<< ":" << kRunLumiEvent.evt << std::endl;
}

unsigned MuonEffRecoId::GetVectorIndex(const vector<pair<double, double> >& binEdges, const double& val)
{
	for (unsigned bin=0; bin < binEdges.size(); ++bin)
	{
		const double min = binEdges.at(bin).first;
		const double max = binEdges.at(bin).second;
		//cout << "val->min-max=" << val << "->" << min << "-" << max << endl;
		if (val>=min && val<=max) return bin;
	}
	return 999;
}
//define this as a plug-in
DEFINE_FWK_MODULE(MuonEffRecoId);

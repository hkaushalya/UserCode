// -*- C++ -*-
//
// Package:    recoJERStudy
// Class:      recoJERStudy
// 
/**\class recoJERStudy recoJERStudy.cc UserCode/recoJERStudy/src/recoJERStudy.cc

 Description: [JET ENERGY RESOLUTION STUDY]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  samantha hewamanage
//         Created:  Mon Aug 15 11:32:01 CDT 2011
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
#include "DataFormats/Math/interface/deltaR.h"

//ROOT
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "assert.h"
#include "TVector3.h"
#include <iomanip>

//
// class declaration
//

class recoJERStudy : public edm::EDFilter {
   public:
      explicit recoJERStudy(const edm::ParameterSet&);
      ~recoJERStudy();

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
		struct JetHist_t
		{
			TH1F* pt;
			TH1F* eta;
			TH1F* phi;
			TH1F* match_pt_res;
			TH1F* match_pt_ratio;
			TH1F* match_delR;
			TH1F* match_delEta;
			TH1F* match_delPhi;
		};
		TH1F* h_minDelR;
		TH1F* h_minDelPhi;
		TH1F* h_minDelEta;


		float eta_bin_step, pt_bin_step;
		std::vector<float> eta_bins, pt_bins;
		std::vector<JetHist_t> JERHist; 
		void BookHistograms(JetHist_t& hist);
		void FillHistograms(JetHist_t& hist);

		edm::InputTag detJetsInputTag_;
		edm::InputTag genJetsInputTag_;
		edm::Handle<reco::PFJetCollection> detJetHandle;
		edm::Handle<std::vector<reco::GenJet> > genJetHandle;

		int iVerbose; // control print levels
		unsigned int iProcessed; // number of processed events
		unsigned int iPassed; //number of events passed the filter

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
static bool sort_using_less_than(double u, double v)
{
	   return u < v;
}
//
// constructors and destructor
//
recoJERStudy::recoJERStudy(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed

	detJetsInputTag_ = iConfig.getParameter<edm::InputTag>("detJetsInputTag");
	genJetsInputTag_ = iConfig.getParameter<edm::InputTag>("genJetsInputTag");
	//mhtInputTag_ = iConfig.getParameter<edm::InputTag>("mhtInputTag");
	//htInputTag_ = iConfig.getParameter<edm::InputTag>("htInputTag");
	//inMinVtx = iConfig.getUntrackedParameter<int>("nMinVtx",1);
	//inMinNdofVtx = iConfig.getUntrackedParameter<int>("ndofVtx",0);
	//dMaxPrimVtxZ = iConfig.getUntrackedParameter<double>("maxDelzVtx",24.0);
	iVerbose = iConfig.getUntrackedParameter<int>("verbose",0);
	iProcessed = 0;
	iPassed = 0;


	eta_bin_step = 0.1;
	pt_bin_step = 10;
	const float eta_max = 5.0;
	const float pt_max = 1000;
	for (float f=0.0; f<= eta_max; f+=eta_bin_step)
	{
		std::cout << "eta bin = " << f << std::endl; 
		eta_bins.push_back(f);
		std::cout << "pt bins = " ; 
		for (float pt=0.0; pt<= pt_max; pt+=pt_bin_step)
		{
			std::cout <<  pt << " | "; 
			pt_bins.push_back(f);
		}
		std::cout << std::endl; 
	}


	//generate hists
	edm::Service<TFileService> fs;
//	BookHistograms(fs, Hist[0], 60, 80); 

	h_minDelR = fs->make<TH1F> ("minDelR" ,"min DelR", 100, 0, 10);
	h_minDelEta = fs->make<TH1F> ("minDelEta" ,"min Del Eta", 100, 0, 10);
	h_minDelPhi = fs->make<TH1F> ("minDelPhi" ,"min Del Phi", 50, 0, 5);

}


recoJERStudy::~recoJERStudy()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
recoJERStudy::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
	++iProcessed;

	iEvent.getByLabel(detJetsInputTag_, detJetHandle);
	if (! detJetHandle.isValid()) 
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":detJetHandle handle not found!" << std::endl;
		assert(false);
	}
	iEvent.getByLabel(genJetsInputTag_, genJetHandle);
	if (! genJetHandle.isValid()) 
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":genJetHandle handle not found!" << std::endl;
		assert(false);
	}


	std::vector<float> vDelR, vDelPhi, vDelEta;
	unsigned numJets = detJetHandle->size();
	unsigned numGenJets = genJetHandle->size();
	unsigned nMatched = 0, nDetJets = 0, nGenJets=0;
	const float minPt = 100.0;
	const float maxEta = 5.0;

	for (unsigned j = 0 ; j < genJetHandle->size() ; ++j)
	{
			if ((*genJetHandle)[j].pt()<minPt || fabs((*genJetHandle)[j].eta())> maxEta) continue;
			++nGenJets;
	}

	for (unsigned i = 0 ; i < numJets ; ++i)
	{
		if ((*detJetHandle)[i].pt()<minPt || fabs((*detJetHandle)[i].eta())> maxEta) continue;
		++nDetJets;
		const TLorentzVector detJetVec((*detJetHandle)[i].px(),(*detJetHandle)[i].py(),(*detJetHandle)[i].pz(),(*detJetHandle)[i].energy());

		for (unsigned j = 0 ; j < genJetHandle->size() ; ++j)
		{
			const TLorentzVector genJetVec((*genJetHandle)[j].px(),(*genJetHandle)[j].py(),(*genJetHandle)[j].pz(),(*genJetHandle)[j].energy());
			const float delphi = fabs(TVector2::Phi_mpi_pi(detJetVec.Phi() - genJetVec.Phi()));
			const float deleta = fabs(detJetVec.Eta() - genJetVec.Eta());
			const float delR = detJetVec.DeltaR(genJetVec);
			vDelR.push_back(delR);
			vDelPhi.push_back(delphi);
			vDelEta.push_back(deleta);

			h_minDelPhi->Fill(delphi);
			h_minDelEta->Fill(deleta);
			h_minDelR->Fill(delR);
			if (delR<0.4) 
			{
				++nMatched;
				std::cout << "match jets [" << i << "," << j << "][" << (*detJetHandle)[i].pt() << "/" << (*genJetHandle)[j].pt() << "/" << std::endl;  
				(*genJetHandle)[j].print();
				break;
			}
		}
	}

	std::sort(vDelPhi.begin(), vDelPhi.end(), sort_using_less_than);	
	std::sort(vDelEta.begin(), vDelEta.end(), sort_using_less_than);	
	std::sort(vDelR.begin(), vDelR.end(), sort_using_less_than);	

	if (nMatched != numJets || nMatched != numGenJets)
	{
		std::cout << "njets det/gen/matched = " << nDetJets << "/ " << nGenJets << " / " << nMatched << std::endl;
	}
/*	if (vDelPhi.size()) 
	{
		std::cout << "min DelR/phi/eta = " << vDelR.at(0) << "/" << vDelPhi.at(0) << "/ " << vDelEta.at(0) << std::endl;  
		h_minDelPhi->Fill(vDelPhi.at(0));
		h_minDelEta->Fill(vDelEta.at(0));
		h_minDelR->Fill(vDelR.at(0));
	}
*/
	++iPassed;
   return true;
}

// ------------ method called once each job just before starting event loop  ------------
void 
recoJERStudy::beginJob()
{


}

// ------------ method called once each job just after ending the event loop  ------------
void 
recoJERStudy::endJob() {
	std::cout << __FILE__ << "::" << __FUNCTION__  << "\n" << std::endl;
	std::cout << "[JER:00] Job settings " << std::endl;
	std::cout << "[JER:01] Events Processed --- = " << iProcessed << std::endl;
	std::cout << "[JER:02] Events Passed ------ = " << iPassed << std::endl;
	//std::cout << "[JER:03] Primary Vtx Min ---- = " << inMinVtx << std::endl; 
	//std::cout << "[JER:04] Primary Vtx Min ndof = " << inMinNdofVtx << std::endl;
}

// ------------ method called when starting to processes a run  ------------
bool 
recoJERStudy::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
recoJERStudy::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
recoJERStudy::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
recoJERStudy::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
recoJERStudy::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(recoJERStudy);

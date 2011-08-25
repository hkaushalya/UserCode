// -*- C++ -*-
//
// Package:    JERStudy
// Class:      JERStudy
// 
/**\class JERStudy JERStudy.cc UserCode/JERStudy/src/JERStudy.cc

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

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
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

class JERStudy : public edm::EDFilter {
   public:
      explicit JERStudy(const edm::ParameterSet&);
      ~JERStudy();

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
		struct MatchJetHist_t
		{
			TH1F* gen_pt;
			TH1F* gen_eta;
			TH1F* gen_phi;
			TH1F* det_pt;
			TH1F* det_eta;
			TH1F* det_phi;
			TH1F* match_E_res;
			TH1F* match_pt_res;
			TH1F* match_pt_ratio;
			TH1F* match_E_ratio;
			TH1F* match_delR;
			TH1F* match_delEta;
			TH1F* match_delPhi;
		};

		TH1F* h_minDelR;
		TH1F* h_minDelPhi;
		TH1F* h_minDelEta;
		TH1F* h_evtWeights;

		struct EtaPtBin_t
		{
			float eta_min;
			float eta_max;
			float pt_min;
			float pt_max;
			std::string str_eta_range;
			std::string str_pt_range;
		};

		float eta_bin_step, pt_bin_step;
//		std::vector<float> eta_bins, pt_bins;
		std::vector<EtaPtBin_t> vEtaPtBins;
		std::vector<MatchJetHist_t> vJERHists; 
		void BookHistograms(edm::Service<TFileService>& fs
					, std::vector<MatchJetHist_t>::iterator it
					, const EtaPtBin_t& etpt);

		unsigned GetJetPtEtaBin(const TLorentzVector& genVec);

		void FillMatchJetHists(const TLorentzVector& genJetVec, 
						const TLorentzVector& detJetVec); 

		edm::InputTag patJetsInputTag_;
		edm::InputTag genJetsInputTag_;
		edm::Handle<std::vector<pat::Jet> > detJetHandle;
		edm::Handle<std::vector<reco::GenJet> > genJetHandle;

		double storedWeight;
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
static void Dump(const TLorentzVector& tl)
{
	std::cout << "pt/ eta/ phi " << std::setw(10) << tl.Pt() 
		<< std::setw(10) << tl.Eta() << std::setw(10) << tl.Phi()
		<< std::endl;
}
static void Dump(pat::Jet  patJet)
{
	std::cout << "pt/ eta/ phi " << std::setw(10) << patJet.pt() 
		<< std::setw(10) << patJet.eta() << std::setw(10) << patJet.phi()
		<< std::endl;
}
static void Dump(reco::GenJet  genJet)
{
	std::cout << "pt/ eta/ phi " << std::setw(10) << genJet.pt() 
		<< std::setw(10) << genJet.eta() << std::setw(10) << genJet.phi()
		<< std::endl;
}

//
// constructors and destructor
//
JERStudy::JERStudy(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed

	patJetsInputTag_ = iConfig.getParameter<edm::InputTag>("patJetsInputTag");
	genJetsInputTag_ = iConfig.getParameter<edm::InputTag>("genJetsInputTag");
	//mhtInputTag_ = iConfig.getParameter<edm::InputTag>("mhtInputTag");
	//htInputTag_ = iConfig.getParameter<edm::InputTag>("htInputTag");
	//inMinVtx = iConfig.getUntrackedParameter<int>("nMinVtx",1);
	//inMinNdofVtx = iConfig.getUntrackedParameter<int>("ndofVtx",0);
	//dMaxPrimVtxZ = iConfig.getUntrackedParameter<double>("maxDelzVtx",24.0);
	iVerbose = iConfig.getUntrackedParameter<int>("verbose",0);
	iProcessed = 0;
	iPassed = 0;


	edm::Service<TFileService> fs;
	eta_bin_step = 0.1;
	const float eta_max = 5.0;
	const float pt_max = 1500;
	
	//JET ENERGY reolution Et@50 GeV ~0.2%, Et@150 GeV ~0.1% and slowly falttens

	const float pt_lessthan100_resolution = 0.2;
	const float pt_morethan100_resolution = 0.1;

	for (float f=0.0; f<= eta_max; f+=eta_bin_step)
	{
		if (iVerbose) std::cout << "eta bin = " << f << std::endl; 
		pt_bin_step = 20;
		for (float pt=0.0; pt<= pt_max; pt+=pt_bin_step)
		{
			if (iVerbose)std::cout <<  pt << " | "; 
			if (pt>=960) pt_bin_step = 100.0;
			else if (pt>=800) pt_bin_step = 80.0;
			else if (pt>=620) pt_bin_step = 60.0;
			else if (pt>=420) pt_bin_step = 40.0;
			else if (pt>=300) pt_bin_step = 30.0;

			EtaPtBin_t etpt;
			etpt.eta_min = f; 
			etpt.eta_max = f + eta_bin_step;
			etpt.pt_min = pt; 
			etpt.pt_max = pt+pt_bin_step;

			std::stringstream eta_range, pt_range;
			eta_range << etpt.eta_min << "<|#Eta|<" << etpt.eta_max;
			pt_range << etpt.pt_min << "<pt<" << etpt.pt_max;
			etpt.str_eta_range = eta_range.str();
			etpt.str_pt_range = pt_range.str();

			vEtaPtBins.push_back(etpt);

			MatchJetHist_t hist;
			vJERHists.push_back(hist);
			BookHistograms(fs, vJERHists.end()-1, etpt);
		}
		if (iVerbose) std::cout << std::endl; 
	}


	//generate hists
//	BookHistograms(fs, Hist[0], 60, 80); 

	h_minDelR = fs->make<TH1F> ("minDelR" ,"min DelR", 100, 0, 10);
	h_minDelEta = fs->make<TH1F> ("minDelEta" ,"min Del Eta", 100, 0, 10);
	h_minDelPhi = fs->make<TH1F> ("minDelPhi" ,"min Del Phi", 50, 0, 5);
	h_evtWeights = fs->make<TH1F> ("eventWeights" ,"Event Weights", 10000, 0, 0.0001);

}


JERStudy::~JERStudy()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
JERStudy::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
	++iProcessed;

	// Get pthat event weight
	edm::Handle<GenEventInfoProduct> genEvtInfoHandle;
	iEvent.getByLabel("generator", genEvtInfoHandle);
	if (genEvtInfoHandle.isValid()) {
		storedWeight = genEvtInfoHandle->weight();
	} else 
	{
		std::cout << __FUNCTION__ << ":" << __LINE__
				<<": generator evt handle not found!" << std::endl;
		assert(false);
	}

//	std::cout << "weight = " << storedWeight << std::endl;
	h_evtWeights->Fill(storedWeight);


	iEvent.getByLabel(patJetsInputTag_, detJetHandle);
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
	unsigned nMatched = 0, nDetJets = 0, nGenJets = 0;
	const float minPt = 30.0;
	const float maxEta = 5.0;

	if (iVerbose)
	{
		std::cout << "=======================" << std::endl;
		std::cout << "********** all det jet minpt >" << minPt << std::endl;
	}
	for (unsigned i = 0 ; i < numJets ; ++i)
	{
		if ((*detJetHandle)[i].pt()<minPt || fabs((*detJetHandle)[i].eta())> maxEta) continue;
		++nDetJets;
		if(iVerbose) Dump((*detJetHandle)[i]);
	}

	if (iVerbose) std::cout << "********** all gen jet minpt >" << minPt << std::endl;

	for (unsigned j = 0 ; j < genJetHandle->size() ; ++j)
	{
		if ((*genJetHandle)[j].pt()<minPt || fabs((*genJetHandle)[j].eta())> maxEta) continue;
		++nGenJets;
		if (iVerbose) Dump((*genJetHandle)[j]);
	}


	if (iVerbose)	std::cout << "**********" << std::endl;

	for (unsigned i = 0 ; i < numJets ; ++i)
	{
		if ((*detJetHandle)[i].pt()<minPt || fabs((*detJetHandle)[i].eta())> maxEta) continue;
		const TLorentzVector iJetVec((*detJetHandle)[i].px(),(*detJetHandle)[i].py(),(*detJetHandle)[i].pz(),(*detJetHandle)[i].energy());
		for (unsigned j = 0 ; j < numJets ; ++j)
		{
			const TLorentzVector jJetVec((*detJetHandle)[j].px(),(*detJetHandle)[j].py(),(*detJetHandle)[j].pz(),(*detJetHandle)[j].energy());
			if (iJetVec.DeltaR(jJetVec)<0.3)
			{
				//std::cout << "self match jets [" << i << "," << j << "::pt ratio=" << iJetVec.Pt()/jJetVec.Pt() << std::endl;
				//Dump(iJetVec);
				//Dump(jJetVec);
			}
		}
	}

	if (iVerbose) std::cout << "**********" << std::endl;

	for (unsigned j = 0 ; j < genJetHandle->size() ; ++j)
	{
		if ((*genJetHandle)[j].pt()<minPt || fabs((*genJetHandle)[j].eta())> maxEta) continue;
		const TLorentzVector genJetVec((*genJetHandle)[j].px(),(*genJetHandle)[j].py(),(*genJetHandle)[j].pz(),(*genJetHandle)[j].energy());

		for (unsigned i = 0 ; i < detJetHandle->size() ; ++i)
		{
			//if ((*detJetHandle)[i].pt()<minPt || fabs((*detJetHandle)[i].eta())> maxEta) continue;
			const TLorentzVector detJetVec((*detJetHandle)[i].px(),(*detJetHandle)[i].py(),(*detJetHandle)[i].pz(),(*detJetHandle)[i].energy());

			const float delphi = fabs(TVector2::Phi_mpi_pi(detJetVec.Phi() - genJetVec.Phi()));
			const float deleta = fabs(detJetVec.Eta() - genJetVec.Eta());
			const float delR = detJetVec.DeltaR(genJetVec);
			vDelR.push_back(delR);
			vDelPhi.push_back(delphi);
			vDelEta.push_back(deleta);

			h_minDelPhi->Fill(delphi);
			h_minDelEta->Fill(deleta);
			h_minDelR->Fill(delR);
			if (delR<0.3) 
			{
				++nMatched;
				//std::cout << "match jets [" << i << "," << j << "::pt ratio=" << detJetVec.Pt()/genJetVec.Pt() << std::endl;
				//Dump(genJetVec);
				//Dump(detJetVec);
				FillMatchJetHists(genJetVec, detJetVec);
				break;
			}
		}
	}

	std::sort(vDelPhi.begin(), vDelPhi.end(), sort_using_less_than);	
	std::sort(vDelEta.begin(), vDelEta.end(), sort_using_less_than);	
	std::sort(vDelR.begin(), vDelR.end(), sort_using_less_than);	

	if (nMatched != nDetJets || nMatched != nGenJets)
	{
		if (iVerbose) std::cout << "################### njets det/gen/matched = " << nDetJets << "/ " << nGenJets << " / " << nMatched << std::endl;
		//if ((nGenJets-nMatched)>2) assert(false);
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
JERStudy::beginJob()
{


}

// ------------ method called once each job just after ending the event loop  ------------
void 
JERStudy::endJob() {
	std::cout << __FILE__ << "::" << __FUNCTION__  << "\n" << std::endl;
	std::cout << "[JER:00] Job settings " << std::endl;
	std::cout << "[JER:01] Events Processed --- = " << iProcessed << std::endl;
	std::cout << "[JER:02] Events Passed ------ = " << iPassed << std::endl;
	//std::cout << "[JER:03] Primary Vtx Min ---- = " << inMinVtx << std::endl; 
	//std::cout << "[JER:04] Primary Vtx Min ndof = " << inMinNdofVtx << std::endl;
}

// ------------ method called when starting to processes a run  ------------
bool 
JERStudy::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
JERStudy::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
JERStudy::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
JERStudy::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
JERStudy::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


unsigned JERStudy::GetJetPtEtaBin(const TLorentzVector& genVec)
{

	for (unsigned bin=0; bin < vEtaPtBins.size(); ++bin)
	{
		if (genVec.Pt()>vEtaPtBins.at(bin).pt_min 
			&& genVec.Pt() <= vEtaPtBins.at(bin).pt_max
			&& fabs(genVec.Eta())> vEtaPtBins.at(bin).eta_min
			&& fabs(genVec.Eta())<= vEtaPtBins.at(bin).eta_max)
		{
			if (iVerbose)
			{	
				std::cout << __FUNCTION__ << ": pt bin found "
					<< bin << " pt[" << vEtaPtBins.at(bin).pt_min << ","
					<< vEtaPtBins.at(bin).pt_max << "] eta["
					<< vEtaPtBins.at(bin).eta_min << ", "
					<< vEtaPtBins.at(bin).eta_max << "]" << std::endl;
			}
			return bin;
		}
	}
	
	std::cout << __FUNCTION__ << ": NO pt/eta(" << genVec.Pt() << ", " << genVec.Eta() << ") bin found!" << std::endl;
//	assert (false);
	return 999999;

}

void JERStudy::FillMatchJetHists(const TLorentzVector& genJetVec, 
						const TLorentzVector& detJetVec)
{
	const unsigned bin = GetJetPtEtaBin(genJetVec);
	if (bin>99999) return; //no bin is found

	vJERHists.at(bin).gen_pt->Fill(genJetVec.Pt(), storedWeight);
	vJERHists.at(bin).gen_eta->Fill(genJetVec.Eta(), storedWeight);
	vJERHists.at(bin).gen_phi->Fill(genJetVec.Phi(), storedWeight);

	vJERHists.at(bin).match_delR->Fill(genJetVec.DeltaR(detJetVec), storedWeight);
	vJERHists.at(bin).match_pt_res->Fill(((detJetVec.Pt()-genJetVec.Pt())/genJetVec.Pt()), storedWeight);
	vJERHists.at(bin).match_pt_ratio->Fill((detJetVec.Pt()/genJetVec.Pt()), storedWeight);
	vJERHists.at(bin).match_E_res->Fill(((detJetVec.E()-genJetVec.E())/genJetVec.E()), storedWeight);
	vJERHists.at(bin).match_E_ratio->Fill((detJetVec.E()/genJetVec.E()), storedWeight);

}

void JERStudy::BookHistograms(edm::Service<TFileService>& fs
					, std::vector<MatchJetHist_t>::iterator it
					, const EtaPtBin_t& etpt)
{
	const float pt_bins = 200, pt_max = 2000;
	const float eta_bins = 100, eta_min = -5, eta_max = 5;
	const float phi_bins = 160, phi_min = -8, phi_max = 8;

	//need to finish this
		std::stringstream etpt_range_name, etpt_range_title;	
		etpt_range_name << "eta_" << etpt.eta_min << "_" << etpt.eta_max << "_pt_" << etpt.pt_min 
					<< "_" << etpt.pt_max;
		etpt_range_title << etpt.eta_min << "<|#eta|<" << etpt.eta_max << ", " << etpt.pt_min 
					<< "<pt<" << etpt.pt_max;

		std::stringstream name_pt, name_eta, name_phi, title_pt, title_eta, title_phi;
		std::stringstream name_delR, name_res_pt, name_res_e, title_delR, title_res_pt, title_res_e;
		std::stringstream name_ratio_pt, title_ratio_pt, name_ratio_e, title_ratio_e;

		name_pt << etpt_range_name.str() << "_genJet_pt";
		name_eta << etpt_range_name.str() << "_genJet_eta";
		name_phi << etpt_range_name.str() << "_genJet_phi";
		name_delR << etpt_range_name.str() << "_delR";
		name_res_pt << etpt_range_name.str() << "_res_pt";
		name_ratio_pt << etpt_range_name.str() << "_ratio_pt";
		name_ratio_e << etpt_range_name.str() << "_ratio_e";
		name_res_e << etpt_range_name.str() << "_res_e";
		title_pt << etpt_range_title.str()  << ";pt [GeV];Arbitrary;";
		title_eta << etpt_range_title.str() << ";#eta;Arbitrary;";
		title_phi << etpt_range_title.str() << ";#phi;Arbitrary;";
		title_delR << etpt_range_title.str() << ";#delta R;Arbitrary;";
		title_res_pt << etpt_range_title.str() << ";pt Resolution: (det-gen)/gen;Arbitrary;";
		title_ratio_pt << etpt_range_title.str() << ";pt Ratio: det/gen;Arbitrary;";
		title_res_e << etpt_range_title.str() << ";E Resolution: (det-gen)/gen;Arbitrary;";
		title_ratio_e << etpt_range_title.str() << ";E Ratio: det/gen;Arbitrary;";

		it->gen_pt = fs->make<TH1F> (name_pt.str().c_str(),title_pt.str().c_str(), pt_bins, 0, pt_max);
		it->gen_eta = fs->make<TH1F> (name_eta.str().c_str(),title_eta.str().c_str(), eta_bins, eta_min, eta_max);
		it->gen_phi = fs->make<TH1F> (name_phi.str().c_str(),title_phi.str().c_str(), phi_bins, phi_min, phi_max);
		it->match_delR = fs->make<TH1F> (name_delR.str().c_str(),title_delR.str().c_str(), 100, 0, 5);
		it->match_pt_ratio = fs->make<TH1F> (name_ratio_pt.str().c_str(),title_ratio_pt.str().c_str(), 120, -1, 5);
		it->match_E_ratio = fs->make<TH1F> (name_ratio_e.str().c_str(),title_ratio_e.str().c_str(), 120, -1, 5);
		it->match_pt_res = fs->make<TH1F> (name_res_pt.str().c_str(),title_res_pt.str().c_str(), 120, -1, 5);
		it->match_E_res = fs->make<TH1F> (name_res_e.str().c_str(),title_res_e.str().c_str(), 120, -1, 5);

}
//define this as a plug-in
DEFINE_FWK_MODULE(JERStudy);

// -*- C++ -*-
//
// Package:    Factorization
// Class:      Factorization
// 
/**\class Factorization Factorization.cc UserCode/Factorization/src/Factorization.cc

 Description: Repeat of 2010 Factorization Study for 2011.

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  samantha hewamanage
//         Created:  Tue Jul 26 22:15:44 CDT 2011
// $Id: Factorization.cc,v 1.1 2011/11/15 21:55:50 samantha Exp $
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
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
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

class Factorization : public edm::EDFilter {
   public:
      explicit Factorization(const edm::ParameterSet&);
      ~Factorization();

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
		unsigned int uProcessed; // number of processed events
		unsigned int uPassed; //number of events passed the filter
		RunLumiEvt_t kRunLumiEvent;
		int iVerbose; // control print levels
		double dMinHT;
		double dMinMHT;

		struct EventHist_t {
			TH1F* mht;
			TH1F* ht;
			TH1F* njet;
			TH1F* meff;
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
			TH1F* delphiMin_jetmht[3];
		};

		TH1F* hDelPhiMin_mht[8];
		TH1F* hPass[16];
		TH1F* hFail[41];
		TH1F* hSignalRegion[27];

		EventHist_t evtHist;
		JetHist_t pf30_jet1Hist, pf30_jet2Hist, pf30_jet3Hist; //for upto 5 jets may be
		JetHist_t pf50eta25_jet1Hist, pf50eta25_jet2Hist, pf50eta25_jet3Hist; //for upto 5 jets may be

		edm::InputTag patJetsPFPt50Eta25InputTag_;
		edm::Handle<std::vector<pat::Jet> > pfpt50eta25JetHandle;

		TH1F* MHT_by_phislice[6];
		void DoDelMinStudy(edm::Handle<std::vector<pat::Jet> > pfpt50eta25JetHandle
				, edm::Handle<edm::View<reco::MET> > mhtHandle);
		void PrintHeader();
		TLorentzVector vMetVec;
		edm::Handle<std::vector<reco::PFMET> >pfMetHandle;
		unsigned uFailMinHTCut, uFailMinPFMHTCut;

		edm::LumiReWeighting LumiWeights_;
		std::vector<float> vMCNvtxDist, vDATANvtxDist;
		double sumLumiWeights; //sum of lumi weights for cross checks
		double Weight; //total weight for an event
		int doLumiWeighing, doEventWeighing;


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
Factorization::Factorization(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
	patJetsPFPt50Eta25InputTag_ = iConfig.getParameter<edm::InputTag>("patJetsPFPt50Eta25InputTag");
	mhtInputTag_    = iConfig.getParameter<edm::InputTag>("mhtInputTag");
	htInputTag_     = iConfig.getParameter<edm::InputTag>("htInputTag");
	//outputFile      = iConfig.getUntrackedParameter<double>("outputFile","Default.root");
	dMinHT          = iConfig.getUntrackedParameter<double>("dMinHT",0.0);
	dMinMHT         = iConfig.getUntrackedParameter<double>("dMinMHT",0.0);
	iVerbose        = iConfig.getUntrackedParameter<int>("verbose",0);
	doLumiWeighing  = iConfig.getUntrackedParameter<int>("ApplyLumiWeighing",0);
	doEventWeighing = iConfig.getUntrackedParameter<int>("ApplyEventWeighing",0);
	uProcessed       = 0;
	uPassed          = 0;
	uFailMinHTCut    = 0;
	uFailMinPFMHTCut = 0;
	sumLumiWeights   = 0;
	Weight           = 1;


	//generate hists
	edm::Service<TFileService> fs;

	const float met_min = 0, met_max=500, met_bins=100;
	MHT_by_phislice[0] = fs->make<TH1F> ("mht_phislice_lt0.1" ,"MHT (#Delta#Phi_{min}<0.1);MHT [GeV];Events;", met_bins, 0, met_max);
	MHT_by_phislice[1] = fs->make<TH1F> ("mht_phislice_lt0.2" ,"MHT (0.1<#Delta#Phi_{min}<0.2);MHT [GeV];Events;", met_bins, 0, met_max);
	MHT_by_phislice[2] = fs->make<TH1F> ("mht_phislice_lt0.3" ,"MHT (0.2<#Delta#Phi_{min}<0.3);MHT [GeV];Events;", met_bins, 0, met_max);
	MHT_by_phislice[3] = fs->make<TH1F> ("mht_phislice_lt0.5" ,"MHT (0.3<#Delta#Phi_{min}<0.5);MHT [GeV];Events;", met_bins, 0, met_max);
	MHT_by_phislice[4] = fs->make<TH1F> ("mht_phislice_lt0.8" ,"MHT (0.5<#Delta#Phi_{min}<0.8);MHT [GeV];Events;", met_bins, 0, met_max);
	MHT_by_phislice[5] = fs->make<TH1F> ("mht_phislice_gt0.8" ,"MHT (#Delta#Phi_{min}>0.8);MHT [GeV];Events;", met_bins, 0, met_max);

	//these are general event hist to check the PAT tuples cuts
	const double evt_met_max = 800, evt_met_bins = 400;
	const double evt_ht_max = 4000, evt_ht_bins = 80;

	evtHist.mht  = fs->make<TH1F> ("mht" ,"RA2: (MHT from PFmetHandle);MHT [GeV];Events;", evt_met_bins, 0, evt_met_max);
	evtHist.ht   = fs->make<TH1F> ("ht" ,"RA2: HT from Jets ET>50 && |#Eta|<2.4;HT [GeV];Events;", evt_ht_bins, 0, evt_ht_max);
	evtHist.njet = fs->make<TH1F> ("njet_et50eta24" ,"RA2: Njets (Et>50 && |#Eta|<2.4;NJETS;Events;", 10, 0, 10);
	evtHist.meff = fs->make<TH1F> ("meteff" ,"RA2:;MEff;Events;", 50, 0, 5000);


	const double pt_bins = 200, pt_max = 2000;
	pf30_jet1Hist.pt  = fs->make<TH1F> ("pf30_jet1_pt"  ,"RA2: PF30-Jet1 pt;pt [GeV];Events;", pt_bins, 0, pt_max);
	pf30_jet1Hist.eta = fs->make<TH1F> ("pf30_jet1_eta" ,"RA2: PF30-Jet1 eta;eta;Events;", 100, -5, 5);
	pf30_jet1Hist.phi = fs->make<TH1F> ("pf30_jet1_phi" ,"RA2: PF30-Jet1 phi;phi;Events;", 160, -8, 8);
//	pf30_jet1Hist.delphi = fs->make<TH1F> ("pf30_jet1_delphi" ,"RA2: PF30-Jet1: delphi;delphi;Events;", 160, -8, 8);

	pf30_jet2Hist.pt  = fs->make<TH1F> ("pf30_jet2_pt"  ,"RA2: PF30-Jet2 pt;pt [GeV];Events;", pt_bins, 0, pt_max);
	pf30_jet2Hist.eta = fs->make<TH1F> ("pf30_jet2_eta" ,"RA2: PF30-Jet2 eta;eta;Events;", 100, -5, 5);
	pf30_jet2Hist.phi = fs->make<TH1F> ("pf30_jet2_phi" ,"RA2: PF30-Jet2 phi;phi;Events;", 160, -8, 8);
//	pf30_jet2Hist.delphi = fs->make<TH1F> ("pf30_jet2_delphi" ,"RA2: PF30-Jet2: delphi;delphi;Events;", 160, -8, 8);

	pf30_jet3Hist.pt = fs->make<TH1F> ("pf30_jet3_pt" ,"RA2: PF30-Jet3 pt;pt [GeV];Events;", pt_bins, 0, pt_max);
	pf30_jet3Hist.eta = fs->make<TH1F> ("pf30_jet3_eta" ,"RA2: PF30-Jet3 eta;eta;Events;", 100, -5, 5);
	pf30_jet3Hist.phi = fs->make<TH1F> ("pf30_jet3_phi" ,"RA2: PF30-Jet3 phi;phi;Events;", 160, -8, 8);
//	pf30_jet3Hist.delphi = fs->make<TH1F> ("pf30_jet3_delphi" ,"RA2: PF30-Jet3: delphi;delphi;Events;", 160, -8, 8);

	pf30_jet1Hist.pt->Sumw2();
	pf30_jet2Hist.pt->Sumw2();
	pf30_jet3Hist.pt->Sumw2();
	pf30_jet1Hist.eta->Sumw2();
	pf30_jet2Hist.eta->Sumw2();
	pf30_jet3Hist.eta->Sumw2();
	pf30_jet1Hist.phi->Sumw2();
	pf30_jet2Hist.phi->Sumw2();
	pf30_jet3Hist.phi->Sumw2();

	//const float npassFailHistBins = 16;
	//const float passFailHistBins[] = {0,20,40,60,80,100,120,140,160,180,200,250,300,350,400,500,600};
	//const float npassFailHistBins = 23;
	//const float passFailHistBins[] = {50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,130,140,150,175,200,250,350,600,1000};
	//const float npassFailHistBins = 14;
	//const float passFailHistBins[] = {50,60,70,80,90,100,110,120,140,160,200,250,350,600,1000};
	const float npassFailHistBins = 13;
	const float passFailHistBins[] = {50,60,70,80,90,100,110,120,140,160,200,300,500,1000};

	hPass[0]  = fs->make<TH1F> ("Pass_0","PASS from #Delta#Phi_{min} cut >0.15", npassFailHistBins, passFailHistBins);
	hPass[1]  = fs->make<TH1F> ("Pass_1","PASS from #Delta#Phi_{min} cui >0.2", npassFailHistBins, passFailHistBins);
	hPass[2]  = fs->make<TH1F> ("Pass_2","PASS from #Delta#Phi_{min} cut >0.25", npassFailHistBins, passFailHistBins);
	hPass[3]  = fs->make<TH1F> ("Pass_3","PASS from #Delta#Phi_{min} cut >0.3", npassFailHistBins, passFailHistBins);
	hPass[4]  = fs->make<TH1F> ("Pass_4","PASS from #Delta#Phi_{min} cut >0.35", npassFailHistBins, passFailHistBins);
	hPass[5]  = fs->make<TH1F> ("Pass_5","PASS from #Delta#Phi_{min} cut >0.4", npassFailHistBins, passFailHistBins);
	hPass[6]  = fs->make<TH1F> ("Pass_RA2dphi","PASS from RA2 dPhi Selection: dphi (j1 & j1 , mht) >0.5 and dphi (j3, mht)>0.3)", npassFailHistBins, passFailHistBins);
	hPass[7]  = fs->make<TH1F> ("Pass_RA2dphi_HT500" ,"PASS from RA2 dPhi Selection: dphi (j1 & j1 , mht) >0.5 and dphi (j3, mht)>0.3) && HT>500 GeV", npassFailHistBins, passFailHistBins);
	hPass[8]  = fs->make<TH1F> ("Pass_RA2dphi_HT800" ,"PASS from RA2 dPhi Selection: dphi (j1 & j1 , mht) >0.5 and dphi (j3, mht)>0.3) && HT>800 GeV", npassFailHistBins, passFailHistBins);
	hPass[9]  = fs->make<TH1F> ("Pass_RA2dphi_HT1000","PASS from RA2 dPhi Selection: dphi (j1 & j1 , mht) >0.5 and dphi (j3, mht)>0.3) && HT>1000 GeV", npassFailHistBins, passFailHistBins);
	hPass[10] = fs->make<TH1F> ("Pass_RA2dphi_HT1200","PASS from RA2 dPhi Selection: dphi (j1 & j1 , mht) >0.5 and dphi (j3, mht)>0.3) && HT>1200 GeV", npassFailHistBins, passFailHistBins);
	hPass[11] = fs->make<TH1F> ("Pass_RA2dphi_HT1400","PASS from RA2 dPhi Selection: dphi (j1 & j1 , mht) >0.5 and dphi (j3, mht)>0.3) && HT>1400 GeV", npassFailHistBins, passFailHistBins);
	hPass[12] = fs->make<TH1F> ("Pass_RA2dphi_500HT800"  ,"PASS from RA2 dPhi Selection: dphi (j1 & j1 , mht) >0.5 and dphi (j3, mht)>0.3) && 500<HT<800 GeV", npassFailHistBins, passFailHistBins);
	hPass[13] = fs->make<TH1F> ("Pass_RA2dphi_800HT1000" ,"PASS from RA2 dPhi Selection: dphi (j1 & j1 , mht) >0.5 and dphi (j3, mht)>0.3) && 800<HT<1000 GeV", npassFailHistBins, passFailHistBins);
	hPass[14] = fs->make<TH1F> ("Pass_RA2dphi_1000HT1200","PASS from RA2 dPhi Selection: dphi (j1 & j1 , mht) >0.5 and dphi (j3, mht)>0.3) && 1000<HT<1200 GeV", npassFailHistBins, passFailHistBins);
	hPass[15] = fs->make<TH1F> ("Pass_RA2dphi_1200HT1400","PASS from RA2 dPhi Selection: dphi (j1 & j1 , mht) >0.5 and dphi (j3, mht)>0.3) && 1200<HT<1400 GeV", npassFailHistBins, passFailHistBins);
	for (int i = 0; i < 16; ++i ) { hPass[i]->Sumw2();}

	hFail[0]  = fs->make<TH1F> ("Fail_0","FAIL from #Delta#Phi_{min} cut <0.15", npassFailHistBins, passFailHistBins);
	hFail[1]  = fs->make<TH1F> ("Fail_1","FAIL from #Delta#Phi_{min} cut <0.2", npassFailHistBins, passFailHistBins);
	hFail[2]  = fs->make<TH1F> ("Fail_2","FAIL from #Delta#Phi_{min} cut <0.25", npassFailHistBins, passFailHistBins);
	hFail[3]  = fs->make<TH1F> ("Fail_3","FAIL from #Delta#Phi_{min} cut <0.3", npassFailHistBins, passFailHistBins);
	hFail[4]  = fs->make<TH1F> ("Fail_4","FAIL from #Delta#Phi_{min} cut <0.35", npassFailHistBins, passFailHistBins);
	hFail[5]  = fs->make<TH1F> ("Fail_5","FAIL from #Delta#Phi_{min} cut <0.4", npassFailHistBins, passFailHistBins);
	hFail[6]  = fs->make<TH1F> ("Fail_lt_point2","Fail from #Delta#Phi_{min} cut <0.2", 1500, 0, 1500);
	hFail[7]  = fs->make<TH1F> ("Fail_lt_point3","Fail from #Delta#Phi_{min} cut <0.3", 1500, 0, 1500);
	hFail[8]  = fs->make<TH1F> ("Fail_lt_point2_HT500","Fail from #Delta#Phi_{min} cut <0.2 && HT>500 GeV", 1500, 0, 1500);
	hFail[9]  = fs->make<TH1F> ("Fail_lt_point2_HT800","Fail from #Delta#Phi_{min} cut <0.2 && HT>800 GeV", 1500, 0, 1500);
	hFail[10] = fs->make<TH1F> ("Fail_lt_point2_HT1000","Fail from #Delta#Phi_{min} cut <0.2 && HT>1000 GeV", 1500, 0, 1500);
	hFail[11] = fs->make<TH1F> ("Fail_lt_point2_HT1200","Fail from #Delta#Phi_{min} cut <0.2 && HT>1200 GeV", 1500, 0, 1500);
	hFail[12] = fs->make<TH1F> ("Fail_lt_point2_HT1400","Fail from #Delta#Phi_{min} cut <0.2 && HT>1400 GeV", 1500, 0, 1500);
	hFail[13] = fs->make<TH1F> ("Fail_lt_point2_500HT800"  ,"Fail from #Delta#Phi_{min} cut <0.2 && 500<HT<800 GeV", 1500, 0, 1500);
	hFail[14] = fs->make<TH1F> ("Fail_lt_point2_800HT1000" ,"Fail from #Delta#Phi_{min} cut <0.2 && 800<HT<1000 GeV", 1500, 0, 1500);
	hFail[15] = fs->make<TH1F> ("Fail_lt_point2_1000HT1200","Fail from #Delta#Phi_{min} cut <0.2 && 1000<HT<1200 GeV", 1500, 0, 1500);
	hFail[16] = fs->make<TH1F> ("Fail_lt_point2_1200HT1400","Fail from #Delta#Phi_{min} cut <0.2 && 1200<HT<1400 GeV", 1500, 0, 1500);
	hFail[17] = fs->make<TH1F> ("Syst1_Fail_ht500","HT>500: FAIL from dphi (j1 & j1 , mht) >0.5 and dphi (j3, mht)>0.1 ", npassFailHistBins, passFailHistBins);
	hFail[18] = fs->make<TH1F> ("Syst1_Fail_ht500_fineBin","FineBin:HT>500: FAIL from dphi (j1 & j1 , mht) >0.5 and dphi (j3, mht)>0.1 ", 1500, 0, 1500);
	hFail[19] = fs->make<TH1F> ("Syst1_Fail_ht600","HT>600: FAIL from dphi (j1 & j1 , mht) >0.5 and dphi (j3, mht)>0.1 ", npassFailHistBins, passFailHistBins);
	hFail[20] = fs->make<TH1F> ("Syst1_Fail_ht600_fineBin","FineBin: HT>600: FAIL from dphi (j1 & j1 , mht) >0.5 and dphi (j3, mht)>0.1 ", 1500, 0, 1500);
	hFail[21] = fs->make<TH1F> ("Syst2_Fail_ht500","HT>500: FAIL from dphi (j1 & j1 , mht) >0.4 and dphi (j3, mht)>0.2", npassFailHistBins, passFailHistBins);
	hFail[22] = fs->make<TH1F> ("Syst2_Fail_ht500_fineBin","FineBin: HT>500: FAIL from dphi (j1 & j1 , mht) >0.4 and dphi (j3, mht)>0.2 ", 1500, 0, 1500);
	hFail[23] = fs->make<TH1F> ("Syst2_Fail_ht600","HT>600: FAIL from dphi (j1 & j1 , mht) >0.4 and dphi (j3, mht)>0.2", npassFailHistBins, passFailHistBins);
	hFail[24] = fs->make<TH1F> ("Syst2_Fail_ht600_fineBin","FineBin: HT>600: FAIL from dphi (j1 & j1 , mht) >0.4 and dphi (j3, mht)>0.2 ", 1500, 0, 1500);

	hFail[25] = fs->make<TH1F> ("Syst1_Fail_HT800_fineBin"     ,"FineBin:HT>800: FAIL from dphi (j1 & j1 , mht) >0.5 and dphi (j3, mht)>0.1 " , 1500, 0, 1500);
	hFail[26] = fs->make<TH1F> ("Syst1_Fail_HT1000_fineBin"    ,"FineBin:HT>1000: FAIL from dphi (j1 & j1 , mht) >0.5 and dphi (j3, mht)>0.1 ", 1500, 0, 1500);
	hFail[27] = fs->make<TH1F> ("Syst1_Fail_HT1200_fineBin"    ,"FineBin:HT>1200: FAIL from dphi (j1 & j1 , mht) >0.5 and dphi (j3, mht)>0.1 ", 1500, 0, 1500);
	hFail[28] = fs->make<TH1F> ("Syst1_Fail_HT1400_fineBin"    ,"FineBin:HT>1400: FAIL from dphi (j1 & j1 , mht) >0.5 and dphi (j3, mht)>0.1 ", 1500, 0, 1500);
	hFail[29] = fs->make<TH1F> ("Syst1_Fail_500HT800_fineBin"  ,"FineBin:500<HT<800: FAIL from dphi (j1 & j1 , mht) >0.5 and dphi (j3, mht)>0.1 "  , 1500, 0, 1500);
	hFail[30] = fs->make<TH1F> ("Syst1_Fail_800HT1000_fineBin" ,"FineBin:800<HT<1000: FAIL from dphi (j1 & j1 , mht) >0.5 and dphi (j3, mht)>0.1 " , 1500, 0, 1500);
	hFail[31] = fs->make<TH1F> ("Syst1_Fail_1000HT1200_fineBin","FineBin:1000<HT<1200: FAIL from dphi (j1 & j1 , mht) >0.5 and dphi (j3, mht)>0.1 ", 1500, 0, 1500);
	hFail[32] = fs->make<TH1F> ("Syst1_Fail_1200HT1400_fineBin","FineBin:1200<HT<1400: FAIL from dphi (j1 & j1 , mht) >0.5 and dphi (j3, mht)>0.1 ", 1500, 0, 1500);

	hFail[33] = fs->make<TH1F> ("Syst2_Fail_HT800_fineBin"     ,"FineBin:HT>800: FAIL from dphi (j1 & j1 , mht) >0.4 and dphi (j3, mht)>0.2 " , 1500, 0, 1500);
	hFail[34] = fs->make<TH1F> ("Syst2_Fail_HT1000_fineBin"    ,"FineBin:HT>1000: FAIL from dphi (j1 & j1 , mht) >0.4 and dphi (j3, mht)>0.2 ", 1500, 0, 1500);
	hFail[35] = fs->make<TH1F> ("Syst2_Fail_HT1200_fineBin"    ,"FineBin:HT>1200: FAIL from dphi (j1 & j1 , mht) >0.4 and dphi (j3, mht)>0.2 ", 1500, 0, 1500);
	hFail[36] = fs->make<TH1F> ("Syst2_Fail_HT1400_fineBin"    ,"FineBin:HT>1400: FAIL from dphi (j1 & j1 , mht) >0.4 and dphi (j3, mht)>0.2 ", 1500, 0, 1500);
	hFail[37] = fs->make<TH1F> ("Syst2_Fail_500HT800_fineBin"  ,"FineBin:500<HT<800: FAIL from dphi (j1 & j1 , mht) >0.4 and dphi (j3, mht)>0.2 "  , 1500, 0, 1500);
	hFail[38] = fs->make<TH1F> ("Syst2_Fail_800HT1000_fineBin" ,"FineBin:800<HT<1000: FAIL from dphi (j1 & j1 , mht) >0.4 and dphi (j3, mht)>0.2 " , 1500, 0, 1500);
	hFail[39] = fs->make<TH1F> ("Syst2_Fail_1000HT1200_fineBin","FineBin:1000<HT<1200: FAIL from dphi (j1 & j1 , mht) >0.4 and dphi (j3, mht)>0.2 ", 1500, 0, 1500);
	hFail[40] = fs->make<TH1F> ("Syst2_Fail_1200HT1400_fineBin","FineBin:1200<HT<1400: FAIL from dphi (j1 & j1 , mht) >0.4 and dphi (j3, mht)>0.2 ", 1500, 0, 1500);

	for (int i = 0; i <41; ++i ) { hFail[i]->Sumw2();}

	hSignalRegion[0]  = fs->make<TH1F> ("Signal_HT500MHT200" ,"Signal Region: HT>500 GeV & MHT>200 GeV && pass RA2 #Delta#Phi_{min}", 1500, 0, 1500);
	hSignalRegion[1]  = fs->make<TH1F> ("Signal_HT500MHT350" ,"Signal Region: HT>500 GeV & MHT>350 GeV && pass RA2 #Delta#Phi_{min}", 1500, 0, 1500);
	hSignalRegion[2]  = fs->make<TH1F> ("Signal_HT500MHT500" ,"Signal Region: HT>500 GeV & MHT>500 GeV && pass RA2 #Delta#Phi_{min}", 1500, 0, 1500);
	hSignalRegion[3]  = fs->make<TH1F> ("Signal_HT800MHT200" ,"Signal Region: HT>800 GeV & MHT>200 GeV && pass RA2 #Delta#Phi_{min}", 1500, 0, 1500);
	hSignalRegion[4]  = fs->make<TH1F> ("Signal_HT800MHT350" ,"Signal Region: HT>800 GeV & MHT>350 GeV && pass RA2 #Delta#Phi_{min}", 1500, 0, 1500);
	hSignalRegion[5]  = fs->make<TH1F> ("Signal_HT800MHT500" ,"Signal Region: HT>800 GeV & MHT>500 GeV && pass RA2 #Delta#Phi_{min}", 1500, 0, 1500);
	hSignalRegion[6]  = fs->make<TH1F> ("Signal_HT1000MHT200","Signal Region: HT>1000 GeV & MHT>200 GeV && pass RA2 #Delta#Phi_{min}", 1500, 0, 1500);
	hSignalRegion[7]  = fs->make<TH1F> ("Signal_HT1000MHT350","Signal Region: HT>1000 GeV & MHT>350 GeV && pass RA2 #Delta#Phi_{min}", 1500, 0, 1500);
	hSignalRegion[8]  = fs->make<TH1F> ("Signal_HT1000MHT500","Signal Region: HT>1000 GeV & MHT>500 GeV && pass RA2 #Delta#Phi_{min}", 1500, 0, 1500);
	hSignalRegion[9]  = fs->make<TH1F> ("Signal_HT1200MHT200","Signal Region: HT>1200 GeV & MHT>200 GeV && pass RA2 #Delta#Phi_{min}", 1500, 0, 1500);
	hSignalRegion[10] = fs->make<TH1F> ("Signal_HT1200MHT350","Signal Region: HT>1200 GeV & MHT>350 GeV && pass RA2 #Delta#Phi_{min}", 1500, 0, 1500);
	hSignalRegion[11] = fs->make<TH1F> ("Signal_HT1200MHT500","Signal Region: HT>1200 GeV & MHT>500 GeV && pass RA2 #Delta#Phi_{min}", 1500, 0, 1500);
	hSignalRegion[12] = fs->make<TH1F> ("Signal_HT1400MHT200","Signal Region: HT>1400 GeV & MHT>200 GeV && pass RA2 #Delta#Phi_{min}", 1500, 0, 1500);
	hSignalRegion[13] = fs->make<TH1F> ("Signal_HT1400MHT350","Signal Region: HT>1400 GeV & MHT>350 GeV && pass RA2 #Delta#Phi_{min}", 1500, 0, 1500);
	hSignalRegion[14] = fs->make<TH1F> ("Signal_HT1400MHT500","Signal Region: HT>1400 GeV & MHT>500 GeV && pass RA2 #Delta#Phi_{min}", 1500, 0, 1500);
	

	hSignalRegion[15] = fs->make<TH1F> ("Signal_HT500to800MHT200"  ,"Signal Region: 500<HT<800 GeV & MHT>200 GeV && pass RA2 #Delta#Phi_{min}", 1500, 0, 1500);
	hSignalRegion[16] = fs->make<TH1F> ("Signal_HT500to800MHT350"  ,"Signal Region: 500<HT<800 GeV & MHT>350 GeV && pass RA2 #Delta#Phi_{min}", 1500, 0, 1500);
	hSignalRegion[17] = fs->make<TH1F> ("Signal_HT500to800MHT500"  ,"Signal Region: 500<HT<800 GeV & MHT>500 GeV && pass RA2 #Delta#Phi_{min}", 1500, 0, 1500);
	hSignalRegion[18] = fs->make<TH1F> ("Signal_HT800to1000MHT200" ,"Signal Region: 800<HT<1000 GeV & MHT>200 GeV && pass RA2 #Delta#Phi_{min}", 1500, 0, 1500);
	hSignalRegion[19] = fs->make<TH1F> ("Signal_HT800to1000MHT350" ,"Signal Region: 800<HT<1000 GeV & MHT>350 GeV && pass RA2 #Delta#Phi_{min}", 1500, 0, 1500);
	hSignalRegion[20] = fs->make<TH1F> ("Signal_HT800to1000MHT500" ,"Signal Region: 800<HT<1000 GeV & MHT>500 GeV && pass RA2 #Delta#Phi_{min}", 1500, 0, 1500);
	hSignalRegion[21] = fs->make<TH1F> ("Signal_HT1000to1200MHT200","Signal Region: 1000<HT<1200 GeV & MHT>200 GeV && pass RA2 #Delta#Phi_{min}", 1500, 0, 1500);
	hSignalRegion[22] = fs->make<TH1F> ("Signal_HT1000to1200MHT350","Signal Region: 1000<HT<1200 GeV & MHT>350 GeV && pass RA2 #Delta#Phi_{min}", 1500, 0, 1500);
	hSignalRegion[23] = fs->make<TH1F> ("Signal_HT1000to1200MHT500","Signal Region: 1000<HT<1200 GeV & MHT>500 GeV && pass RA2 #Delta#Phi_{min}", 1500, 0, 1500);
	hSignalRegion[24] = fs->make<TH1F> ("Signal_HT1200to1400MHT200","Signal Region: 1200<HT<1400 GeV & MHT>200 GeV && pass RA2 #Delta#Phi_{min}", 1500, 0, 1500);
	hSignalRegion[25] = fs->make<TH1F> ("Signal_HT1200to1400MHT350","Signal Region: 1200<HT<1400 GeV & MHT>350 GeV && pass RA2 #Delta#Phi_{min}", 1500, 0, 1500);
	hSignalRegion[26] = fs->make<TH1F> ("Signal_HT1200to1400MHT500","Signal Region: 1200<HT<1400 GeV & MHT>500 GeV && pass RA2 #Delta#Phi_{min}", 1500, 0, 1500);
	for (int i = 0; i <27; ++i ) { hSignalRegion[i]->Sumw2();}

	hDelPhiMin_mht[0] = fs->make<TH1F> ("delPhiMin_MHT_60TO80"  ,"RA2: #Delta#Phi_{min} for 60 <#slash{H}_{T}<80 GeV" , 150, 0, 3);
	hDelPhiMin_mht[1] = fs->make<TH1F> ("delPhiMin_MHT_80to100" ,"RA2: #Delta#Phi_{min} for 80 <#slash{H}_{T}<100 GeV", 150, 0, 3);
	hDelPhiMin_mht[2] = fs->make<TH1F> ("delPhiMin_MHT_100to120","RA2: #Delta#Phi_{min} for 100<#slash{H}_{T}<120 GeV", 150, 0, 3);
	hDelPhiMin_mht[3] = fs->make<TH1F> ("delPhiMin_MHT_120to140","RA2: #Delta#Phi_{min} for 120<#slash{H}_{T}<140 GeV", 150, 0, 3);
	hDelPhiMin_mht[4] = fs->make<TH1F> ("delPhiMin_MHT_140to170","RA2: #Delta#Phi_{min} for 140<#slash{H}_{T}<170 GeV", 150, 0, 3);
	hDelPhiMin_mht[5] = fs->make<TH1F> ("delPhiMin_MHT_170to200","RA2: #Delta#Phi_{min} for 170<#slash{H}_{T}<200 GeV", 150, 0, 3);
	hDelPhiMin_mht[6] = fs->make<TH1F> ("delPhiMin_MHT_200to250","RA2: #Delta#Phi_{min} for 200<#slash{H}_{T}<250 GeV", 150, 0, 3);
	hDelPhiMin_mht[7] = fs->make<TH1F> ("delPhiMin_MHT_250up"   ,"RA2: #Delta#Phi_{min} for #slash{H}_{T}>250 GeV"    , 150, 0, 3);

	hDelPhiMin_mht[1]->Sumw2();
	hDelPhiMin_mht[2]->Sumw2();
	hDelPhiMin_mht[3]->Sumw2();
	hDelPhiMin_mht[4]->Sumw2();
	hDelPhiMin_mht[5]->Sumw2();
	hDelPhiMin_mht[6]->Sumw2();
	hDelPhiMin_mht[7]->Sumw2();
}


Factorization::~Factorization()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
Factorization::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
	++uProcessed;
	kRunLumiEvent.run  = iEvent.id().run();
	kRunLumiEvent.lumi = iEvent.id().luminosityBlock(); 
	kRunLumiEvent.evt  = iEvent.id().event();

//    std::cout << __LINE__<< ":: Processing event: "
//		 	<<  kRunLumiEvent.run << ":" << kRunLumiEvent.lumi 
//			<< ":" << kRunLumiEvent.evt << std::endl;
//			return 0;

    if (iVerbose) std::cout << __LINE__<< ":: Processing event: "
		 	<<  kRunLumiEvent.run << ":" << kRunLumiEvent.lumi 
			<< ":" << kRunLumiEvent.evt << std::endl;


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
		//std::cout << "storedWeight = " << storedWeight << std::endl;
		Weight *= storedWeight;
	}

	//std::cout << "Weight = " << Weight << std::endl;


	

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

	iEvent.getByLabel(mhtInputTag_, mhtHandle);
	if (! mhtHandle.isValid()) 
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":MHT handle not found!" << std::endl;
		assert(false);
	}

	//APPLY RA2 cuts

	if ( pfpt50eta25JetHandle->size() <3 ) return 0;
	if ( (*htHandle) < dMinHT) { ++uFailMinHTCut; return 0; }
	if ( (*mhtHandle)[0].pt() < dMinMHT ) {++uFailMinPFMHTCut; return 0; }

	/*
	 * Fill hists after RA2b base selection
	 */
	evtHist.njet->Fill(pfpt50eta25JetHandle->size(), Weight);
	evtHist.mht->Fill( (*mhtHandle)[0].pt(), Weight);
	evtHist.ht->Fill((*htHandle), Weight);
	evtHist.meff->Fill( (*mhtHandle)[0].pt() + (*htHandle), Weight);

	pf30_jet1Hist.pt->Fill((*pfpt50eta25JetHandle)[0].pt(), Weight);
	pf30_jet2Hist.pt->Fill((*pfpt50eta25JetHandle)[1].pt(), Weight);
	pf30_jet3Hist.pt->Fill((*pfpt50eta25JetHandle)[2].pt(), Weight);
	pf30_jet1Hist.eta->Fill((*pfpt50eta25JetHandle)[0].eta(), Weight);
	pf30_jet2Hist.eta->Fill((*pfpt50eta25JetHandle)[1].eta(), Weight);
	pf30_jet3Hist.eta->Fill((*pfpt50eta25JetHandle)[2].eta(), Weight);

	DoDelMinStudy(pfpt50eta25JetHandle,	mhtHandle);

	++uPassed;
   return true;
}

// ------------ method called once each job just before starting event loop  ------------
void 
Factorization::beginJob()
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
Factorization::endJob() {

	std::cout << __FILE__ << "::" << __FUNCTION__  << "\n" << std::endl;
	std::cout << "[ATM:00] Job settings " << std::endl;
	std::cout << "[ATM:01] Events Processed --- = " << uProcessed << std::endl;
	std::cout << "[ATM:02] Events Passed ------ = " << uPassed << std::endl;
	std::cout << "[ATM:03] Lumi Weighing? ----- = " << doLumiWeighing << std::endl;
	std::cout << "[ATM:04] Event Weighing? ---- = " << doEventWeighing << std::endl;
	std::cout << "[ATM:31] Minimum HT --------- = " << dMinHT << std::endl;
	std::cout << "[ATM:32] Minimum MHT -------- = " << dMinMHT << std::endl;
	std::cout << "[ATM:40] PASS Summary --------- " << std::endl;
	std::cout << "[ATM:43] Pass ht ------------ = " << (uProcessed - uFailMinHTCut) 
																	<< " (" << uFailMinHTCut << ")" << std::endl;
	std::cout << "[ATM:44] Pass mht ----------- = " << (uProcessed - uFailMinHTCut - uFailMinPFMHTCut) 
																	<< " (" << uFailMinPFMHTCut << ")" << std::endl;
	std::cout << "[ATM:50] LumiWeights Avg ---- = ";
	if ( uPassed>0) std::cout << sumLumiWeights/(double)uPassed << std::endl;
	else std::cout << "0" << std::endl;

}

// ------------ method called when starting to processes a run  ------------
bool 
Factorization::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
Factorization::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
Factorization::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
Factorization::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Factorization::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void Factorization::DoDelMinStudy(edm::Handle<std::vector<pat::Jet> > pfpt50eta25JetHandle
				, edm::Handle<edm::View<reco::MET> > mhtHandle
				)
{
	//PrintHeader();
	std::vector<float> vDelPhi_jetmht;
	const float mht = (*mhtHandle)[0].pt();

	for (unsigned i = 0 ; i < pfpt50eta25JetHandle->size() ; ++i)
	{
		if (i>2) break; //use only three leading jets
		//std::cout << "i = " <<  i << std::endl;
		const float delphi_jetmht = fabs(TVector2::Phi_mpi_pi((*pfpt50eta25JetHandle)[i].phi() - (*mhtHandle)[0].phi()));
		vDelPhi_jetmht.push_back(delphi_jetmht);
	}

	std::sort(vDelPhi_jetmht.begin(), vDelPhi_jetmht.end(), sort_using_less_than);	

	assert (vDelPhi_jetmht.size() == 3 && "ERROR: More than 3 dPhiMin calculations found!");
	//std::cout << "vDelPhi_jetmht size = " << vDelPhi_jetmht.size() << std::endl;
	const float dPhiMin = vDelPhi_jetmht.at(0);

	// dphimin in slices of MHT
	if ( mht >= 60 && mht< 80 ) hDelPhiMin_mht[0]->Fill(dPhiMin, Weight);
	else if ( mht >= 80 && mht< 100 ) hDelPhiMin_mht[1]->Fill(dPhiMin, Weight);
	else if ( mht >= 100 && mht< 120 ) hDelPhiMin_mht[2]->Fill(dPhiMin, Weight);
	else if ( mht >= 120 && mht< 140 ) hDelPhiMin_mht[3]->Fill(dPhiMin, Weight);
	else if ( mht >= 140 && mht< 170 ) hDelPhiMin_mht[4]->Fill(dPhiMin, Weight);
	else if ( mht >= 170 && mht< 200 ) hDelPhiMin_mht[5]->Fill(dPhiMin, Weight);
	else if ( mht >= 200 && mht< 250 ) hDelPhiMin_mht[6]->Fill(dPhiMin, Weight);
	else if ( mht >= 250 ) hDelPhiMin_mht[7]->Fill(dPhiMin, Weight);

	// mht in slices of dphimin
	if ( dPhiMin < 0.1 ) MHT_by_phislice[0]->Fill(mht);
	else if ( dPhiMin >= 0.1 && dPhiMin < 0.2 ) MHT_by_phislice[1]->Fill(mht, Weight);
	else if ( dPhiMin >= 0.2 && dPhiMin < 0.3 ) MHT_by_phislice[2]->Fill(mht, Weight);
	else if ( dPhiMin >= 0.3 && dPhiMin < 0.5 ) MHT_by_phislice[3]->Fill(mht, Weight);
	else if ( dPhiMin >= 0.5 && dPhiMin < 0.8 ) MHT_by_phislice[4]->Fill(mht, Weight);
	else if ( dPhiMin > 0.8 ) MHT_by_phislice[5]->Fill(mht, Weight);


	//make PASS/FAIL plots with RA2 cut

	if (dPhiMin>0.15) hPass[0]->Fill(mht, Weight);
	else 	hFail[0]->Fill(mht, Weight);

	if (dPhiMin>0.2)
	{
		hPass[1]->Fill(mht, Weight);
	} else
	{
		hFail[1]->Fill(mht, Weight);
		hFail[6]->Fill(mht, Weight);
	}

	if (dPhiMin>0.25) hPass[2]->Fill(mht, Weight);
	else 	hFail[2]->Fill(mht, Weight);

	if (dPhiMin>=0.3)
	{
		hPass[3]->Fill(mht, Weight);
	} else
	{
		hFail[3]->Fill(mht, Weight);
		hFail[7]->Fill(mht, Weight);
	}

	if (dPhiMin>=0.35) hPass[4]->Fill(mht, Weight);
	else 	hFail[4]->Fill(mht, Weight);

	if (dPhiMin>=0.4) hPass[5]->Fill(mht, Weight);
	else 	hFail[5]->Fill(mht, Weight);


	if (dPhiMin<0.2)
	{
		if ( (*htHandle) > 500) hFail[8]->Fill(mht, Weight);
		if ( (*htHandle) > 800) hFail[9]->Fill(mht, Weight);
		if ( (*htHandle) > 1000) hFail[10]->Fill(mht, Weight);
		if ( (*htHandle) > 1200) hFail[11]->Fill(mht, Weight);
		if ( (*htHandle) > 1400) hFail[12]->Fill(mht, Weight);

		if ( (*htHandle) > 500  && (*htHandle) < 800) hFail[13]->Fill(mht, Weight);
		if ( (*htHandle) > 800  && (*htHandle) < 1000) hFail[14]->Fill(mht, Weight);
		if ( (*htHandle) > 1000 && (*htHandle) < 1200) hFail[15]->Fill(mht, Weight);
		if ( (*htHandle) > 1200 && (*htHandle) < 1400) hFail[16]->Fill(mht, Weight);
	}



	//Pass selection with RA2 dphi cuts
	bool passed = true;
	if (pfpt50eta25JetHandle->size() >= 1) passed = passed && (std::abs(reco::deltaPhi((*pfpt50eta25JetHandle)[0].phi(), (*mhtHandle)[0].phi())) > 0.5);
	if (pfpt50eta25JetHandle->size() >= 2) passed = passed && (std::abs(reco::deltaPhi((*pfpt50eta25JetHandle)[1].phi(), (*mhtHandle)[0].phi())) > 0.5);
	if (pfpt50eta25JetHandle->size() >= 3) passed = passed && (std::abs(reco::deltaPhi((*pfpt50eta25JetHandle)[2].phi(), (*mhtHandle)[0].phi())) > 0.3);

	if (passed) 
	{
		hPass[6]->Fill(mht, Weight);
		if ( (*htHandle) > 500)
		{
			hPass[7]->Fill(mht, Weight);
			if (mht>200) hSignalRegion[0]->Fill(mht,Weight);
			if (mht>350) hSignalRegion[1]->Fill(mht,Weight);
			if (mht>500) hSignalRegion[2]->Fill(mht,Weight);
		}
		if ( (*htHandle) > 800) 
		{
			hPass[8]->Fill(mht, Weight);
			if (mht>200) hSignalRegion[3]->Fill(mht,Weight);
			if (mht>350) hSignalRegion[4]->Fill(mht,Weight);
			if (mht>500) hSignalRegion[5]->Fill(mht,Weight);
		}
		if ( (*htHandle) > 1000) 
		{
			hPass[9]->Fill(mht, Weight);
			if (mht>200) hSignalRegion[6]->Fill(mht,Weight);
			if (mht>350) hSignalRegion[7]->Fill(mht,Weight);
			if (mht>500) hSignalRegion[8]->Fill(mht,Weight);
		}
		if ( (*htHandle) > 1200) 
		{
			hPass[10]->Fill(mht, Weight);
			if (mht>200) hSignalRegion[9]->Fill(mht,Weight);
			if (mht>350) hSignalRegion[10]->Fill(mht,Weight);
			if (mht>500) hSignalRegion[11]->Fill(mht,Weight);
		}
		if ( (*htHandle) > 1400) 
		{
			hPass[11]->Fill(mht, Weight);
			if (mht>200) hSignalRegion[12]->Fill(mht,Weight);
			if (mht>350) hSignalRegion[13]->Fill(mht,Weight);
			if (mht>500) hSignalRegion[14]->Fill(mht,Weight);
		}

		if ( (*htHandle) > 500  && (*htHandle) < 800 )
		{
			hPass[12]->Fill(mht, Weight);
			if (mht>200) hSignalRegion[15]->Fill(mht,Weight);
			if (mht>350) hSignalRegion[16]->Fill(mht,Weight);
			if (mht>500) hSignalRegion[17]->Fill(mht,Weight);
		}
		if ( (*htHandle) > 800  && (*htHandle) < 1000)
		{
			hPass[13]->Fill(mht, Weight);
			if (mht>200) hSignalRegion[18]->Fill(mht,Weight);
			if (mht>350) hSignalRegion[19]->Fill(mht,Weight);
			if (mht>500) hSignalRegion[20]->Fill(mht,Weight);
		}
		if ( (*htHandle) > 1000 && (*htHandle) < 1200)
		{
			hPass[14]->Fill(mht, Weight);
			if (mht>200) hSignalRegion[21]->Fill(mht,Weight);
			if (mht>350) hSignalRegion[22]->Fill(mht,Weight);
			if (mht>500) hSignalRegion[23]->Fill(mht,Weight);
		}
		if ( (*htHandle) > 1200 && (*htHandle) < 1400)
		{
			hPass[15]->Fill(mht, Weight);
			if (mht>200) hSignalRegion[24]->Fill(mht,Weight);
			if (mht>350) hSignalRegion[25]->Fill(mht,Weight);
			if (mht>500) hSignalRegion[26]->Fill(mht,Weight);
		}
	}


	/**************************************************/
	//                 for systematics
	/**************************************************/
	passed = true;
	if (pfpt50eta25JetHandle->size() >= 1) passed = passed && (std::abs(reco::deltaPhi((*pfpt50eta25JetHandle)[0].phi(), (*mhtHandle)[0].phi())) > 0.5);
	if (pfpt50eta25JetHandle->size() >= 2) passed = passed && (std::abs(reco::deltaPhi((*pfpt50eta25JetHandle)[1].phi(), (*mhtHandle)[0].phi())) > 0.5);
	if (pfpt50eta25JetHandle->size() >= 3) passed = passed && (std::abs(reco::deltaPhi((*pfpt50eta25JetHandle)[2].phi(), (*mhtHandle)[0].phi())) > 0.1);
	
	if (! passed)
	{
		if ( (*htHandle) > 500) 
		{
			hFail[17]->Fill(mht, Weight);
			hFail[18]->Fill(mht, Weight);
		}
		if ( (*htHandle) > 600)
		{
			hFail[19]->Fill(mht, Weight);
			hFail[20]->Fill(mht, Weight);
		}

		if ( (*htHandle) > 800) hFail[25]->Fill(mht, Weight);
		if ( (*htHandle) > 1000) hFail[26]->Fill(mht, Weight);
		if ( (*htHandle) > 1200) hFail[27]->Fill(mht, Weight);
		if ( (*htHandle) > 1400) hFail[28]->Fill(mht, Weight);

		if ( (*htHandle) > 500  && (*htHandle) < 800) hFail[29]->Fill(mht, Weight);
		if ( (*htHandle) > 800  && (*htHandle) < 1000) hFail[30]->Fill(mht, Weight);
		if ( (*htHandle) > 1000 && (*htHandle) < 1200) hFail[31]->Fill(mht, Weight);
		if ( (*htHandle) > 1200 && (*htHandle) < 1400) hFail[32]->Fill(mht, Weight);

	}


	passed = true;
	if (pfpt50eta25JetHandle->size() >= 1) passed = passed && (std::abs(reco::deltaPhi((*pfpt50eta25JetHandle)[0].phi(), (*mhtHandle)[0].phi())) > 0.4);
	if (pfpt50eta25JetHandle->size() >= 2) passed = passed && (std::abs(reco::deltaPhi((*pfpt50eta25JetHandle)[1].phi(), (*mhtHandle)[0].phi())) > 0.4);
	if (pfpt50eta25JetHandle->size() >= 3) passed = passed && (std::abs(reco::deltaPhi((*pfpt50eta25JetHandle)[2].phi(), (*mhtHandle)[0].phi())) > 0.2);
	
	if (! passed)
	{
		if ( (*htHandle) > 500) 
		{
			hFail[21]->Fill(mht, Weight);
			hFail[22]->Fill(mht, Weight);
		}
		if ( (*htHandle) > 600)
		{
			hFail[23]->Fill(mht, Weight);
			hFail[24]->Fill(mht, Weight);
		}

		if ( (*htHandle) > 800) hFail[33]->Fill(mht, Weight);
		if ( (*htHandle) > 1000) hFail[34]->Fill(mht, Weight);
		if ( (*htHandle) > 1200) hFail[35]->Fill(mht, Weight);
		if ( (*htHandle) > 1400) hFail[36]->Fill(mht, Weight);

		if ( (*htHandle) > 500  && (*htHandle) < 800) hFail[37]->Fill(mht, Weight);
		if ( (*htHandle) > 800  && (*htHandle) < 1000) hFail[38]->Fill(mht, Weight);
		if ( (*htHandle) > 1000 && (*htHandle) < 1200) hFail[39]->Fill(mht, Weight);
		if ( (*htHandle) > 1200 && (*htHandle) < 1400) hFail[40]->Fill(mht, Weight);
	}


	if (iVerbose)
	{
		std::cout << "delphi ordered (jetmet)= ";
		for (std::vector<float>::const_iterator it = vDelPhi_jetmht.begin(); it != vDelPhi_jetmht.end(); ++it)
		{
			std::cout << "  " << (*it);
		}
		std::cout << std::endl;
	}
}

void Factorization::PrintHeader()
{
	std::cout  << "======= Run/Event: " << kRunLumiEvent.run
			<< ":" << kRunLumiEvent.evt << std::endl;
}


//define this as a plug-in
DEFINE_FWK_MODULE(Factorization);

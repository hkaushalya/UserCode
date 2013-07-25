// -*- C++ -*-
//
// Package:    Factorization_SmearedJets
// Class:      Factorization_SmearedJets
// 
/**\class Factorization_SmearedJets Factorization_SmearedJets.cc UserCode/Factorization_SmearedJets/src/Factorization_SmearedJets.cc

 Description: Repeat of 2010 Factorization_SmearedJets Study for 2011.

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  samantha hewamanage
//         Created:  Tue Jul 26 22:15:44 CDT 2011
// $Id: Factorization_SmearedJets.cc,v 1.2 2012/05/03 14:51:52 samantha Exp $
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
#include "TRandom3.h"
#include "TF1.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "PhysicsTools/KinFitter/interface/TKinFitter.h"
#include "PhysicsTools/KinFitter/interface/TFitParticleEtEtaPhi.h"
#include "PhysicsTools/KinFitter/interface/TFitConstraintEp.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "CommonTools/Utils/interface/PtComparator.h"
//
// class declaration
//

using namespace std;

class Factorization_SmearedJets : public edm::EDFilter {
   public:
      explicit Factorization_SmearedJets(const edm::ParameterSet&);
      ~Factorization_SmearedJets();

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
		edm::InputTag prescaleWeightInputTag;

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

		struct CutsPassedHists_t
		{
			TH1F *processed;
			TH1F *njet; 
			TH1F *ht;
			TH1F *mht;
		};

		CutsPassedHists_t cutsHist;

		TH1F* hDelPhiMin_mht[8];
		TH1F* hPass[16];
		TH1F* hFail[41];
		TH1F* hSignalRegion[27];

		EventHist_t evtHist;
		JetHist_t pf30_jet1Hist, pf30_jet2Hist, pf30_jet3Hist; //for upto 5 jets may be
		JetHist_t pf50eta25_jet1Hist, pf50eta25_jet2Hist, pf50eta25_jet3Hist; //for upto 5 jets may be
		TH1F* hDelPhiMin;

		edm::InputTag smearedGenJetInputTag_;
		//edm::Handle<std::vector<reco::GenJet> > smearedJetsHandle;
		edm::Handle<std::vector<pat::Jet> > smearedJetsHandle;

		TH1F* MHT_by_phislice[6];
		void DoDelMinStudy(edm::Handle<std::vector<pat::Jet> > smearedJetsHandle, 
									const TLorentzVector &mhtVec);
		void CalcSmearedHTMHT(edm::Handle<std::vector<pat::Jet> > smearedJetsHandle, 
									int &smearedNJet, float &smearedHT, 
									float &smearedMHT, float &smearedMEFF, TLorentzVector &smearedMHTVEC);
		void GenerateJets(edm::Event& iEvent, const edm::EventSetup& iSetup);

		void PrintHeader();
		TLorentzVector vMetVec;
		edm::Handle<std::vector<reco::PFMET> >pfMetHandle;
		unsigned uFailNjetCut, uFailMinHTCut, uFailMinPFMHTCut;

		edm::LumiReWeighting LumiWeights_;
		std::vector<float> vMCNvtxDist, vDATANvtxDist;
		double sumLumiWeights; //sum of lumi weights for cross checks
		double Weight; //total weight for an event
		int doLumiWeighing, doEventWeighing;
		edm::Handle<double> prescaleWeightHandle;
		bool usePrescaleWeight;
		TH1F* hist_prescaleWeights;
		TH1F* hist_eventWeights;
		float smearedHT, smearedMHT, smearedMEFF; //generated using smeared jets
		int   smearedNJetEt50Eta25; 
		TLorentzVector smearedMHTVEC;
		int binningOption_;

		/**********************************************
		 *  R+S smearing stuff 
		 **********************************************/
      
      typedef math::XYZTLorentzVector LorentzVector;
      typedef std::vector<std::string>::const_iterator StrIter;

      double rebalancedJetPt_;
      std::string rebalanceMode_; // "MHTall", "MHThigh" or "MET" only for smearCollection = "Reco"

      int nSmearedJets_;
      double smearedJetPt_;
      std::vector<double> PtBinEdges_scaling_;
      std::vector<double> EtaBinEdges_scaling_;
      std::vector<double> AdditionalSmearing_;
      std::vector<double> LowerTailScaling_;
      std::vector<double> UpperTailScaling_;
      double AdditionalSmearing_variation_;
      double LowerTailScaling_variation_;
      double UpperTailScaling_variation_;
		
		std::string smearCollection_; // "Gen" or "Reco"
      edm::InputTag genjets_;
      edm::InputTag jets_;
      edm::InputTag weightName_;
      std::string jets_reb_;
      std::string met_reb_;
      std::string jets_smeared_;
      std::string genjets_smeared_;

      //// vector of response function
      std::vector<std::vector<std::vector<TH1F*> > > smearFunc;
      std::vector<std::vector<std::vector<TH1F*> > > smearFunc_Core;
      std::vector<std::vector<std::vector<TH1F*> > > smearFunc_LowerTail;
      std::vector<std::vector<std::vector<TH1F*> > > smearFunc_UpperTail;
      std::vector<std::vector<std::vector<TH1F*> > > smearFunc_scaled;
      std::vector<std::vector<TH1F*> > SigmaPtHist;
      std::vector<std::vector<TF1*> > SigmaPt;
      std::vector<std::vector<TH1F*> > SigmaPtHist_scaled;
      std::vector<std::vector<TF1*> > SigmaPt_scaled;
      std::vector<double> PtBinEdges_;
      std::vector<double> EtaBinEdges_;
      std::string inputhist1_;
      std::string inputhist2_;
      std::string inputhist3p_;
      std::string smearingfile_;
      std::string outputfile_;
      int NRebin_;
      bool controlPlots_;
      bool absoluteTailScaling_;
      bool cleverPrescaleTreating_;
      double A0RMS_;
      double A1RMS_;
      double probExtreme_;
      double MHTmin_;
      double MHTmax_;
      double HTmin_;
      double HTmax_;
      int NbinsMHT_;
      int NbinsHT_;
      int Ntries_;
      int NJets_;
      double JetsHTPt_;
      double JetsHTEta_;
      double JetsMHTPt_;
      double JetsMHTEta_;
      std::vector<double> JetDeltaMin_;
      double MHTcut_low_;
      double MHTcut_medium_;
      double MHTcut_high_;
      double HTcut_low_;
      double HTcut_medium_;
      double HTcut_high_;
      double HTcut_veryhigh_;
      double HTcut_extremehigh_;

      std::string uncertaintyName_;

      int plotindex_;
      TRandom3 *rand_;

      double JetResolution_Pt2(const double&, const double&, const int&);
      double JetResolution_Ptrel(const double&, const double&, const int&);
      double JetResolution_Eta2(const double&, const double&);
      double JetResolution_Phi2(const double&, const double&);
      double JetResolutionHist_Pt_Smear(const double&, const double&, const int&);
      double GetAdditionalSmearing(const double&, const double&);
      double GetLowerTailScaling(const double&, const double&);
      double GetUpperTailScaling(const double&, const double&);
      int GetIndex(const double&, const std::vector<double>*);
      void FoldWithGaussian(const TH1&, TH1&, const double&);
      void StretchHisto(const TH1&, TH1&, const double&);
      void FillPredictionHistos(const std::vector<pat::Jet>&, const int&, const double&);
      double calcHT(const std::vector<pat::Jet>&);
      double calcHT_gen(const std::vector<reco::GenJet>&);
      math::PtEtaPhiMLorentzVector calcMHT(const std::vector<pat::Jet>&);
      math::PtEtaPhiMLorentzVector calcMHT_gen(const std::vector<reco::GenJet>&);
      int calcNJets(const std::vector<pat::Jet>&);
      int calcNJets_gen(const std::vector<reco::GenJet>&);
      bool calcMinDeltaPhi(const std::vector<pat::Jet>&, math::PtEtaPhiMLorentzVector&);
      bool calcMinDeltaPhi_gen(const std::vector<reco::GenJet>&, math::PtEtaPhiMLorentzVector&);

      bool RebalanceJets_KinFitter(edm::View<pat::Jet>*, std::vector<pat::Jet> &);
      void SmearingJets(const std::vector<pat::Jet>&, std::vector<pat::Jet> &);
      void SmearingGenJets(edm::View<reco::GenJet>*, std::vector<reco::GenJet> &);
      void FillPredictionHistos_gen(const std::vector<reco::GenJet>&, const int&, const double&);

      std::string GetName(const std::string plot, const std::string uncert = "", const std::string ptbin = "") const;

      TH2F* h_RecJetRes_Pt;
      TH2F* h_RecJetRes_Eta;
      TH2F* h_RebJetRes_Pt;
      TH2F* h_RebJetRes_Eta;
      TH2F* h_SmearedJetRes_Pt;
      TH2F* h_SmearedJetRes_Eta;

      TH1F *h_HTall_gen, *h_HTall_rec, *h_HTall_smeared, *h_HTall_reb;
      TH1F *h_HThigh_gen, *h_HThigh_rec, *h_HThigh_smeared, *h_HThigh_reb;
      TH1F *h_MHTall_gen, *h_MHTall_rec, *h_MHTall_smeared, *h_MHTall_reb;
      TH1F *h_MHThigh_gen, *h_MHThigh_rec, *h_MHThigh_smeared, *h_MHThigh_reb;

      TH1F *h_fitProb;
      TH1F *h_weight;
      TH1F *h_weightedWeight;

      TH2F *h_presel_MHT_prediction;
      TH2F *h_presel_HT_prediction;
      TH2F *h_presel_Meff_prediction;
      TH2F *h_baselineNoDeltaPhi_MHT_prediction;
      TH2F *h_baselineNoDeltaPhi_HT_prediction;
      TH2F *h_baselineNoDeltaPhi_Meff_prediction;
      TH2F *h_lowHT_MHT_prediction;
      TH2F *h_lowHT_Meff_prediction;
      TH2F *h_lowMHT_HT_prediction;
      TH2F *h_mediumHT_MHT_prediction;
      TH2F *h_mediumHT_Meff_prediction;
      TH2F *h_mediumMHT_HT_prediction;
      TH2F *h_highHT_MHT_prediction;
      TH2F *h_highHT_Meff_prediction;
      TH2F *h_highMHT_HT_prediction;
      TH2F *h_veryhighHT_MHT_prediction;
      TH2F *h_veryhighHT_Meff_prediction;
      TH2F *h_extremehighHT_MHT_prediction;
      TH2F *h_extremehighHT_Meff_prediction;

      double weight_;

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
Factorization_SmearedJets::Factorization_SmearedJets(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
	smearedGenJetInputTag_ = iConfig.getParameter<edm::InputTag>("genJetSrc");
	//outputFile    = iConfig.getUntrackedParameter<double>("outputFile","Default.root");
	dMinHT          = iConfig.getUntrackedParameter<double>("dMinHT",0.0);
	dMinMHT         = iConfig.getUntrackedParameter<double>("dMinMHT",0.0);
	iVerbose        = iConfig.getUntrackedParameter<int>("verbose",0);
	doLumiWeighing  = iConfig.getUntrackedParameter<int>("ApplyLumiWeighing",0);
	doEventWeighing = iConfig.getUntrackedParameter<int>("ApplyEventWeighing",0);
	prescaleWeightInputTag = iConfig.getParameter<edm::InputTag>("prescaleWeight");
	usePrescaleWeight = iConfig.getUntrackedParameter<int>("usePrescaleWeight",0);
	binningOption_    = iConfig.getUntrackedParameter<int>("binningOption",0);
	uProcessed       = 0;
	uPassed          = 0;
	uFailNjetCut     = 0;
	uFailMinHTCut    = 0;
	uFailMinPFMHTCut = 0;
	sumLumiWeights   = 0;
	Weight           = 1;

	//inputs for R+S smearing
   rebalancedJetPt_ = iConfig.getParameter<double> ("RebalanceJetPt");
   rebalanceMode_ = iConfig.getParameter<std::string> ("RebalanceMode");
   //nSmearedJets_ = iConfig.getParameter<int> ("NSmearedJets");
   smearedJetPt_ = iConfig.getParameter<double> ("SmearedJetPt");
   PtBinEdges_scaling_ = iConfig.getParameter<std::vector<double> > ("PtBinEdges_scaling");
   EtaBinEdges_scaling_ = iConfig.getParameter<std::vector<double> > ("EtaBinEdges_scaling");
   AdditionalSmearing_ = iConfig.getParameter<std::vector<double> > ("AdditionalSmearing");
   LowerTailScaling_ = iConfig.getParameter<std::vector<double> > ("LowerTailScaling");
   UpperTailScaling_ = iConfig.getParameter<std::vector<double> > ("UpperTailScaling");
   AdditionalSmearing_variation_ = iConfig.getParameter<double> ("AdditionalSmearing_variation");
   LowerTailScaling_variation_ = iConfig.getParameter<double> ("LowerTailScaling_variation");
   UpperTailScaling_variation_ = iConfig.getParameter<double> ("UpperTailScaling_variation");
   smearCollection_ = iConfig.getParameter<std::string> ("SmearCollection");

   genjets_ = iConfig.getParameter<edm::InputTag> ("genjetCollection");
   jets_ = iConfig.getParameter<edm::InputTag> ("jetCollection");
   jets_reb_ = iConfig.getParameter<std::string> ("jetCollection_reb");
   jets_smeared_ = iConfig.getParameter<std::string> ("jetCollection_smeared");
   genjets_smeared_ = iConfig.getParameter<std::string> ("genjetCollection_smeared");
   uncertaintyName_ = iConfig.getParameter<std::string> ("uncertaintyName");
   inputhist1_ = iConfig.getParameter<std::string> ("InputHisto1");
   inputhist2_ = iConfig.getParameter<std::string> ("InputHisto2");
   inputhist3p_ = iConfig.getParameter<std::string> ("InputHisto3p");
   controlPlots_ = iConfig.getParameter<bool> ("ControlPlots");
   smearingfile_ = iConfig.getParameter<std::string> ("SmearingFile");
   outputfile_ = iConfig.getParameter<std::string> ("OutputFile");
   NRebin_ = iConfig.getParameter<int> ("NRebin");
   weightName_ = iConfig.getParameter<edm::InputTag> ("weightName");
   absoluteTailScaling_ = iConfig.getParameter<bool> ("absoluteTailScaling");
   cleverPrescaleTreating_ = iConfig.getParameter<bool> ("cleverPrescaleTreating");
   A0RMS_ = iConfig.getParameter<double> ("A0RMS");
   A1RMS_ = iConfig.getParameter<double> ("A1RMS");
   probExtreme_ = iConfig.getParameter<double> ("probExtreme");
   MHTmin_ = iConfig.getParameter<double> ("MHTmin");
   MHTmax_ = iConfig.getParameter<double> ("MHTmax");
   HTmin_ = iConfig.getParameter<double> ("HTmin");
   HTmax_ = iConfig.getParameter<double> ("HTmax");
   NbinsMHT_ = iConfig.getParameter<int> ("NbinsMHT");
   NbinsHT_ = iConfig.getParameter<int> ("NbinsHT");
   Ntries_ = iConfig.getParameter<int> ("Ntries");
   NJets_ = iConfig.getParameter<int> ("NJets");
   JetsHTPt_ = iConfig.getParameter<double> ("JetsHTPt");
   JetsHTEta_ = iConfig.getParameter<double> ("JetsHTEta");
   JetsMHTPt_ = iConfig.getParameter<double> ("JetsMHTPt");
   JetsMHTEta_ = iConfig.getParameter<double> ("JetsMHTEta");
   JetDeltaMin_ = iConfig.getParameter<std::vector<double> > ("JetDeltaMin");
   MHTcut_low_ = iConfig.getParameter<double> ("MHTcut_low");
   MHTcut_medium_ = iConfig.getParameter<double> ("MHTcut_medium");
   MHTcut_high_ = iConfig.getParameter<double> ("MHTcut_high");
   HTcut_low_ = iConfig.getParameter<double> ("HTcut_low");
   HTcut_medium_ = iConfig.getParameter<double> ("HTcut_medium");
   HTcut_high_ = iConfig.getParameter<double> ("HTcut_high");
   HTcut_veryhigh_ = iConfig.getParameter<double> ("HTcut_veryhigh");
   HTcut_extremehigh_ = iConfig.getParameter<double> ("HTcut_extremehigh");

   PtBinEdges_ = iConfig.getParameter<std::vector<double> > ("PtBinEdges");
   EtaBinEdges_ = iConfig.getParameter<std::vector<double> > ("EtaBinEdges");

   unsigned int needed_dim = (PtBinEdges_scaling_.size() - 1) * (EtaBinEdges_scaling_.size() - 1);
   if (AdditionalSmearing_.size() != needed_dim) {
      throw edm::Exception(edm::errors::Configuration, "AdditionalSmearing has not correct dimension");
   }
   if (LowerTailScaling_.size() != needed_dim) {
      throw edm::Exception(edm::errors::Configuration, "LowerTailScaling has not correct dimension");
   }
   if (UpperTailScaling_.size() != needed_dim) {
      throw edm::Exception(edm::errors::Configuration, "UpperTailScaling has not correct dimension");
   }

   // Different seed per initialization
   gRandom->SetSeed(0);
   rand_ = new TRandom3(0);

}


Factorization_SmearedJets::~Factorization_SmearedJets()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

	//R+S Smearing stuff
   for (std::vector<std::vector<std::vector<TH1F*> > >::iterator it = smearFunc.begin(); it != smearFunc.end(); ++it) {
      for (std::vector<std::vector<TH1F*> >::iterator jt = it->begin(); jt != it->end(); ++jt) {
         for (std::vector<TH1F*>::iterator kt = jt->begin(); kt != jt->end(); ++kt) {
            delete *kt;
         }
      }
   }
   smearFunc.clear();

   for (std::vector<std::vector<std::vector<TH1F*> > >::iterator it = smearFunc_Core.begin(); it
         != smearFunc_Core.end(); ++it) {
      for (std::vector<std::vector<TH1F*> >::iterator jt = it->begin(); jt != it->end(); ++jt) {
         for (std::vector<TH1F*>::iterator kt = jt->begin(); kt != jt->end(); ++kt) {
            delete *kt;
         }
      }
   }
   smearFunc_Core.clear();

   for (std::vector<std::vector<std::vector<TH1F*> > >::iterator it = smearFunc_LowerTail.begin(); it
         != smearFunc_LowerTail.end(); ++it) {
      for (std::vector<std::vector<TH1F*> >::iterator jt = it->begin(); jt != it->end(); ++jt) {
         for (std::vector<TH1F*>::iterator kt = jt->begin(); kt != jt->end(); ++kt) {
            delete *kt;
         }
      }
   }
   smearFunc_LowerTail.clear();

   for (std::vector<std::vector<std::vector<TH1F*> > >::iterator it = smearFunc_UpperTail.begin(); it
         != smearFunc_UpperTail.end(); ++it) {
      for (std::vector<std::vector<TH1F*> >::iterator jt = it->begin(); jt != it->end(); ++jt) {
         for (std::vector<TH1F*>::iterator kt = jt->begin(); kt != jt->end(); ++kt) {
            delete *kt;
         }
      }
   }
   smearFunc_UpperTail.clear();

   for (std::vector<std::vector<std::vector<TH1F*> > >::iterator it = smearFunc_scaled.begin(); it
         != smearFunc_scaled.end(); ++it) {
      for (std::vector<std::vector<TH1F*> >::iterator jt = it->begin(); jt != it->end(); ++jt) {
         for (std::vector<TH1F*>::iterator kt = jt->begin(); kt != jt->end(); ++kt) {
            delete *kt;
         }
      }
   }
   smearFunc_scaled.clear();

   PtBinEdges_.clear();
   EtaBinEdges_.clear();

   if (rand_)
      delete rand_;

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
Factorization_SmearedJets::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	std::cout << __LINE__ << std::endl;
   using namespace edm;
	++uProcessed;
	cutsHist.processed->Fill(1);
	kRunLumiEvent.run  = iEvent.id().run();
	kRunLumiEvent.lumi = iEvent.id().luminosityBlock(); 
	kRunLumiEvent.evt  = iEvent.id().event();

//    std::cout << __LINE__<< ":: Processing event: "
//		 	<<  kRunLumiEvent.run << ":" << kRunLumiEvent.lumi 
//			<< ":" << kRunLumiEvent.evt << std::endl;
//			return 0;

	smearedHT   = 0;
	smearedMHT  = 0;
	smearedMEFF = 0 ;
	smearedMHTVEC.SetPxPyPzE(0,0,0,0);

    if (iVerbose) std::cout << __LINE__<< ":: Processing event: "
		 	<<  kRunLumiEvent.run << ":" << kRunLumiEvent.lumi 
			<< ":" << kRunLumiEvent.evt << std::endl;

	std::cout << __LINE__ << std::endl;

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

	std::cout << __LINE__ << std::endl;
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

	std::cout << __LINE__ << std::endl;
	/* 
	 * prescaled trigger weights
	 */
	double prescaleWeight = 1;
	if ( iEvent.isRealData() && usePrescaleWeight)
	{
		//iEvent.getByLabel("prescaleweightProducer:weight", prescaleWeightHandle);
		iEvent.getByLabel(prescaleWeightInputTag, prescaleWeightHandle);
	   //iEvent.getByLabel("prescaleweightProducer", prescaleWeightHandle);
      //std::cout << "prescaleWeight = " << (*prescaleWeightHandle) << std::endl;
		if ( ! prescaleWeightHandle.isValid()) 
		{ 
			std::cout << "prescaleWeightHandle found!" << std::endl;	
			assert (false);
		}
		prescaleWeight = (*prescaleWeightHandle);
		hist_prescaleWeights->Fill(prescaleWeight);
		Weight *= prescaleWeight;
		//std::cout << "Weight = " << Weight << std::endl;
	}
	std::cout << __LINE__ << std::endl;

	hist_eventWeights->Fill(Weight);
	//std::cout << "Weight = " << Weight << std::endl;

	/*  Get all the handles needed and check their
	 *  validity.
	 */

	std::cout << __LINE__ << std::endl;
	const bool bUseRSMod = true;
	if (bUseRSMod)
	{
		iEvent.getByLabel(smearedGenJetInputTag_, smearedJetsHandle);
		//iEvent.getByLabel("QCDfromSmearing:smearedJets", smearedJetsHandle);
		//std::cout << __FILE__ << ":" << __FUNCTION__ << ":: smearedGenJetInputTag_ = " << smearedGenJetInputTag_ << std::endl;
		smearedNJetEt50Eta25 = 0;
		smearedHT = 0;
		smearedMHT = 0;
		smearedMEFF = 0;
		smearedMHTVEC.SetPxPyPzE(0,0,0,0);

		if (! smearedJetsHandle.isValid()) 
		{
			std::cout << __FUNCTION__ << ":" << __LINE__ << ":smearedJetsHandle handle not found!" << std::endl;
			assert(false);
		} else {
			CalcSmearedHTMHT(smearedJetsHandle, smearedNJetEt50Eta25,
					smearedHT, smearedMHT, smearedMEFF, smearedMHTVEC);
		}
	} else 
	{ //this is to do the smearing mutiple times per event
	  // all i need to do is use rebalance jets and do mutiple
	  // experiments per event

		
		GenerateJets(iEvent, iSetup);
		return 0;
		iEvent.getByLabel(smearedGenJetInputTag_, smearedJetsHandle);
		smearedNJetEt50Eta25 = 0;
		smearedHT = 0;
		smearedMHT = 0;
		smearedMEFF = 0;
		smearedMHTVEC.SetPxPyPzE(0,0,0,0);
		CalcSmearedHTMHT(smearedJetsHandle, smearedNJetEt50Eta25,
				smearedHT, smearedMHT, smearedMEFF, smearedMHTVEC);
	}

	//APPLY RA2 cuts
	//std::cout << "smeared jet size = " << smearedJetsHandle->size() << std::endl;
	// need to count number of jets after smearing. Ideally, I want to remove RA2 preselection and use a
								  // skim without njet selection so the jets flutuating upward can be included. For now I do not have 
								  // skim as I am using Seema's dataset.
	
	if ( smearedNJetEt50Eta25 < 3 ) return 0;
	cutsHist.njet->Fill(1);
	if ( smearedHT < dMinHT) { ++uFailMinHTCut; return 0; }
	cutsHist.ht->Fill(1);
	if ( smearedMHT < dMinMHT ) {++uFailMinPFMHTCut; return 0; }
	cutsHist.mht->Fill(1);

	evtHist.njet->Fill(smearedNJetEt50Eta25, Weight);
	evtHist.mht->Fill(smearedMHT, Weight);
	evtHist.ht->Fill(smearedHT, Weight);
	evtHist.meff->Fill(smearedMEFF, Weight);

	DoDelMinStudy(smearedJetsHandle,	smearedMHTVEC);

	++uPassed;
   return true;
}

// ------------ method called once each job just before starting event loop  ------------
void 
Factorization_SmearedJets::beginJob()
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


	//generate hists
	edm::Service<TFileService> fs;

	cutsHist.processed  = fs->make<TH1F> ("evts_processed" ,"Number of events pass the processed cut", 5, 0, 5);
	cutsHist.njet  = fs->make<TH1F> ("evtsPassedNjetCut" ,"Number of events pass the njet cut", 5, 0, 5);
	cutsHist.ht    = fs->make<TH1F> ("evtsPassedHtCut" ,"Number of events pass the ht cut", 5, 0, 5);
	cutsHist.mht   = fs->make<TH1F> ("evtsPassedMhtCut" ,"Number of events pass the mht cut", 5, 0, 5);

	const float met_min = 0, met_max=500, met_bins=100;
	MHT_by_phislice[0] = fs->make<TH1F> ("mht_phislice_lt0.1" ,"MHT (#Delta#Phi_{min}<0.1);MHT [GeV];Events;", met_bins, met_min, met_max);
	MHT_by_phislice[1] = fs->make<TH1F> ("mht_phislice_lt0.2" ,"MHT (0.1<#Delta#Phi_{min}<0.2);MHT [GeV];Events;", met_bins, met_min, met_max);
	MHT_by_phislice[2] = fs->make<TH1F> ("mht_phislice_lt0.3" ,"MHT (0.2<#Delta#Phi_{min}<0.3);MHT [GeV];Events;", met_bins, met_min, met_max);
	MHT_by_phislice[3] = fs->make<TH1F> ("mht_phislice_lt0.5" ,"MHT (0.3<#Delta#Phi_{min}<0.5);MHT [GeV];Events;", met_bins, met_min, met_max);
	MHT_by_phislice[4] = fs->make<TH1F> ("mht_phislice_lt0.8" ,"MHT (0.5<#Delta#Phi_{min}<0.8);MHT [GeV];Events;", met_bins, met_min, met_max);
	MHT_by_phislice[5] = fs->make<TH1F> ("mht_phislice_gt0.8" ,"MHT (#Delta#Phi_{min}>0.8);MHT [GeV];Events;", met_bins, met_min, met_max);

	//these are general event hist to check the PAT tuples cuts
	const double evt_met_max = 800, evt_met_bins = 400;
	const double evt_ht_max = 4000, evt_ht_bins = 80;

	evtHist.mht  = fs->make<TH1F> ("mht" ,"MHT from smeared jets;MHT [GeV];Events;", evt_met_bins, 0, evt_met_max);
	evtHist.ht   = fs->make<TH1F> ("ht" ,"HT from Jets ET>50 && |#eta|<2.5;HT [GeV];Events;", evt_ht_bins, 0, evt_ht_max);
	evtHist.njet = fs->make<TH1F> ("njet_et50eta24" ,"Njets (Et>50 && |#eta|<2.5;NJETS;Events;", 10, 0, 10);
	evtHist.meff = fs->make<TH1F> ("meteff" ,"MEFF;MEff;Events;", 50, 0, 5000);


	const double pt_bins = 200, pt_max = 2000;
	pf30_jet1Hist.pt  = fs->make<TH1F> ("pf30_jet1_pt"  ,"RA2: PF30-Jet1 pt;pt [GeV];Events;", pt_bins, 0, pt_max);
	pf30_jet1Hist.eta = fs->make<TH1F> ("pf30_jet1_eta" ,"RA2: PF30-Jet1 eta;eta;Events;", 100, -5, 5);
	pf30_jet1Hist.phi = fs->make<TH1F> ("pf30_jet1_phi" ,"RA2: PF30-Jet1 phi;phi;Events;", 160, -8, 8);
	pf30_jet1Hist.delphi = fs->make<TH1F> ("pf30_jet1_delphi" ,"RA2: PF30-Jet1: #Delta#Phi; #Delta#Phi (jet, MHT) ;Events;", 40, 0, 4);

	pf30_jet2Hist.pt  = fs->make<TH1F> ("pf30_jet2_pt"  ,"RA2: PF30-Jet2 pt;pt [GeV];Events;", pt_bins, 0, pt_max);
	pf30_jet2Hist.eta = fs->make<TH1F> ("pf30_jet2_eta" ,"RA2: PF30-Jet2 eta;eta;Events;", 100, -5, 5);
	pf30_jet2Hist.phi = fs->make<TH1F> ("pf30_jet2_phi" ,"RA2: PF30-Jet2 phi;phi;Events;", 160, -8, 8);
	pf30_jet2Hist.delphi = fs->make<TH1F> ("pf30_jet2_delphi" ,"RA2: PF30-Jet2: #Delta#Phi; #Delta#Phi (jet, MHT);Events;", 40, 0, 4);

	pf30_jet3Hist.pt = fs->make<TH1F> ("pf30_jet3_pt" ,"RA2: PF30-Jet3 pt;pt [GeV];Events;", pt_bins, 0, pt_max);
	pf30_jet3Hist.eta = fs->make<TH1F> ("pf30_jet3_eta" ,"RA2: PF30-Jet3 eta;eta;Events;", 100, -5, 5);
	pf30_jet3Hist.phi = fs->make<TH1F> ("pf30_jet3_phi" ,"RA2: PF30-Jet3 phi;phi;Events;", 160, -8, 8);
	pf30_jet3Hist.delphi = fs->make<TH1F> ("pf30_jet3_delphi" ,"RA2: PF30-Jet3: #Delta#Phi; #Delta#Phi (jet, MHT);Events;", 40, 0, 4);

	pf30_jet1Hist.pt->Sumw2();
	pf30_jet2Hist.pt->Sumw2();
	pf30_jet3Hist.pt->Sumw2();
	pf30_jet1Hist.eta->Sumw2();
	pf30_jet2Hist.eta->Sumw2();
	pf30_jet3Hist.eta->Sumw2();
	pf30_jet1Hist.phi->Sumw2();
	pf30_jet2Hist.phi->Sumw2();
	pf30_jet3Hist.phi->Sumw2();
	pf30_jet1Hist.delphi->Sumw2();
	pf30_jet2Hist.delphi->Sumw2();
	pf30_jet3Hist.delphi->Sumw2();

	/*std::vector<float> vBins;

	if (binningOption_ == 1) //to check event counts in each of the signal bins
	{
		const float nHistBins = 7;
		const float HistBins[] = {50,100,150,200,350,500,600,1000};
		for (int i=0; i < nHistBins; ++i) vBins.push_back(i);
	} else 
	{
		const float nHistBins = 13;
		const float HistBins[] = {50,60,70,80,90,100,110,120,140,160,200,300,500,1000};
		for (int i=0; i < nHistBins; ++i) vBins.push_back(i);
	}

	const float npassFailHistBins = (float) vBins.size();
	const int iBins = (int) vBins.size();
	float passFailHistBins[iBins];
	for (unsigned i=0; i < vBins.size(); ++i) passFailHistBins[i] = vBins.at(i);
	*/

	const float npassFailHistBins = 13;
	const float passFailHistBins[] = {50,60,70,80,90,100,110,120,140,160,200,300,500,1000};
	//to check event counts in each of the signal bins
	//const float npassFailHistBins = 7;
	//const float passFailHistBins[] = {50,100,150,200,350,500,600,1000};

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

	hDelPhiMin = fs->make<TH1F> ("delPhiMin"  ,"RA2: #Delta#Phi_{min}" , 150, 0, 3);
	hDelPhiMin->Sumw2();
	
	hist_prescaleWeights = fs->make<TH1F> ("prescaleWeights","Prescale Weights",2000,0,2000);
	hist_prescaleWeights->Sumw2();
	hist_eventWeights    = fs->make<TH1F> ("totalEventWeights","Total Event Weights",2000,0,2000);
	hist_eventWeights->Sumw2();
	
	/******************************
	* R+S PLOTS
	******************************/

   if (controlPlots_) {
      h_RecJetRes_Pt = fs->make<TH2F> ("RecJetRes_Pt", "RecJetRes_Pt", 100, 0., 1000., 100, 0., 3.);
      h_RecJetRes_Pt->Sumw2();
      h_RecJetRes_Eta = fs->make<TH2F> ("RecJetRes_Eta", "RecJetRes_Eta", 100, -5., 5., 100, 0., 3.);
      h_RecJetRes_Eta->Sumw2();
      h_RebJetRes_Pt = fs->make<TH2F> ("RebJetRes_Pt", "RebJetRes_Pt", 100, 0., 1000., 100, 0., 3.);
      h_RebJetRes_Pt->Sumw2();
      h_RebJetRes_Eta = fs->make<TH2F> ("RebJetRes_Eta", "RebJetRes_Eta", 100, -5., 5., 100, 0., 3.);
      h_RebJetRes_Eta->Sumw2();
      h_SmearedJetRes_Pt = fs->make<TH2F> ("SmearedJetRes_Pt", "SmearedJetRes_Pt", 100, 0., 1000., 100, 0., 3.);
      h_SmearedJetRes_Pt->Sumw2();
      h_SmearedJetRes_Eta = fs->make<TH2F> ("SmearedJetRes_Eta", "SmearedJetRes_Eta", 100, -5., 5., 100, 0., 3.);
      h_SmearedJetRes_Eta->Sumw2();
      h_HTall_gen = fs->make<TH1F> ("HTall_gen", "HTall_gen", NbinsHT_, HTmin_, HTmax_);
      h_HTall_gen->Sumw2();
      h_HTall_rec = fs->make<TH1F> ("HTall_rec", "HTall_rec", NbinsHT_, HTmin_, HTmax_);
      h_HTall_rec->Sumw2();
      h_HTall_smeared = fs->make<TH1F> ("HTall_smeared", "HTall_smeared", NbinsHT_, HTmin_, HTmax_);
      h_HTall_smeared->Sumw2();
      h_HTall_reb = fs->make<TH1F> ("HTall_reb", "HTall_reb", NbinsHT_, HTmin_, HTmax_);
      h_HTall_reb->Sumw2();
      h_HThigh_gen = fs->make<TH1F> ("HThigh_gen", "HThigh_gen", NbinsHT_, HTmin_, HTmax_);
      h_HThigh_gen->Sumw2();
      h_HThigh_rec = fs->make<TH1F> ("HThigh_rec", "HThigh_rec", NbinsHT_, HTmin_, HTmax_);
      h_HThigh_rec->Sumw2();
      h_HThigh_smeared = fs->make<TH1F> ("HThigh_smeared", "HThigh_smeared", NbinsHT_, HTmin_, HTmax_);
      h_HThigh_smeared->Sumw2();
      h_HThigh_reb = fs->make<TH1F> ("HThigh_reb", "HThigh_reb", NbinsHT_, HTmin_, HTmax_);
      h_HThigh_reb->Sumw2();
      h_MHTall_gen = fs->make<TH1F> ("MHTall_gen", "MHTall_gen", NbinsMHT_, MHTmin_, MHTmax_);
      h_MHTall_gen->Sumw2();
      h_MHTall_rec = fs->make<TH1F> ("MHTall_rec", "MHTall_rec", NbinsMHT_, MHTmin_, MHTmax_);
      h_MHTall_rec->Sumw2();
      h_MHTall_smeared = fs->make<TH1F> ("MHTall_smeared", "MHTall_smeared", NbinsMHT_, MHTmin_, MHTmax_);
      h_MHTall_smeared->Sumw2();
      h_MHTall_reb = fs->make<TH1F> ("MHTall_reb", "MHTall_reb", NbinsMHT_, MHTmin_, MHTmax_);
      h_MHTall_reb->Sumw2();
      h_MHThigh_gen = fs->make<TH1F> ("MHThigh_gen", "MHThigh_gen", NbinsMHT_, MHTmin_, MHTmax_);
      h_MHThigh_gen->Sumw2();
      h_MHThigh_rec = fs->make<TH1F> ("MHThigh_rec", "MHThigh_rec", NbinsMHT_, MHTmin_, MHTmax_);
      h_MHThigh_rec->Sumw2();
      h_MHThigh_smeared = fs->make<TH1F> ("MHThigh_smeared", "MHThigh_smeared", NbinsMHT_, MHTmin_, MHTmax_);
      h_MHThigh_smeared->Sumw2();
      h_MHThigh_reb = fs->make<TH1F> ("MHThigh_reb", "MHThigh_reb", NbinsMHT_, MHTmin_, MHTmax_);
      h_MHThigh_reb->Sumw2();
      h_fitProb = fs->make<TH1F> ("h_fitProb", "h_fitProb", 100, 0., 1.);
      h_fitProb->Sumw2();
      h_weight = fs->make<TH1F> ("h_weight", "h_weight", 70, -1., 6.);
      h_weight->Sumw2();
      h_weightedWeight = fs->make<TH1F> ("h_weightedWeight", "h_weightedWeight", 70, -1., 6.);
      h_weightedWeight->Sumw2();
   }

   //// Bootstrap plots

   //// Preselection
   h_presel_MHT_prediction = fs->make<TH2F> ("presel_MHT_prediction", "presel_MHT_prediction", NbinsMHT_, MHTmin_,
         MHTmax_, Ntries_, 0.5, Ntries_ + 0.5);
   h_presel_MHT_prediction->Sumw2();
   h_presel_HT_prediction = fs->make<TH2F> ("presel_HT_prediction", "presel_HT_prediction", NbinsHT_, HTmin_, HTmax_,
         Ntries_, 0.5, Ntries_ + 0.5);
   h_presel_HT_prediction->Sumw2();
   h_presel_Meff_prediction = fs->make<TH2F> ("presel_Meff_prediction", "presel_Meff_prediction", NbinsHT_, HTmin_,
         HTmax_, Ntries_, 0.5, Ntries_ + 0.5);
   h_presel_Meff_prediction->Sumw2();

   //// Baseline, no DeltaPhi
   h_baselineNoDeltaPhi_MHT_prediction = fs->make<TH2F> ("baselineNoDeltaPhi_MHT_prediction",
         "baselineNoDeltaPhi_MHT_prediction", NbinsMHT_, MHTmin_, MHTmax_, Ntries_, 0.5, Ntries_ + 0.5);
   h_baselineNoDeltaPhi_MHT_prediction->Sumw2();
   h_baselineNoDeltaPhi_HT_prediction = fs->make<TH2F> ("baselineNoDeltaPhi_HT_prediction",
         "baselineNoDeltaPhi_HT_prediction", NbinsHT_, HTmin_, HTmax_, Ntries_, 0.5, Ntries_ + 0.5);
   h_baselineNoDeltaPhi_HT_prediction->Sumw2();
   h_baselineNoDeltaPhi_Meff_prediction = fs->make<TH2F> ("baselineNoDeltaPhi_Meff_prediction",
         "baselineNoDeltaPhi_Meff_prediction", NbinsHT_, HTmin_, HTmax_, Ntries_, 0.5, Ntries_ + 0.5);
   h_baselineNoDeltaPhi_Meff_prediction->Sumw2();

   //// low HT
   h_lowHT_MHT_prediction = fs->make<TH2F> ("lowHT_MHT_prediction", "lowHT_MHT_prediction", NbinsMHT_, MHTmin_,
         MHTmax_, Ntries_, 0.5, Ntries_ + 0.5);
   h_lowHT_MHT_prediction->Sumw2();
   h_lowHT_Meff_prediction = fs->make<TH2F> ("lowHT_Meff_prediction", "lowHT_Meff_prediction", NbinsHT_, HTmin_,
         HTmax_, Ntries_, 0.5, Ntries_ + 0.5);
   h_lowHT_Meff_prediction->Sumw2();

   //// low MHT
   h_lowMHT_HT_prediction = fs->make<TH2F> ("lowMHT_HT_prediction", "lowMHT_HT_prediction", NbinsHT_, HTmin_, HTmax_,
         Ntries_, 0.5, Ntries_ + 0.5);
   h_lowMHT_HT_prediction->Sumw2();

   //// medium HT
   h_mediumHT_MHT_prediction = fs->make<TH2F> ("mediumHT_MHT_prediction", "mediumHT_MHT_prediction", NbinsMHT_,
         MHTmin_, MHTmax_, Ntries_, 0.5, Ntries_ + 0.5);
   h_mediumHT_MHT_prediction->Sumw2();
   h_mediumHT_Meff_prediction = fs->make<TH2F> ("mediumHT_Meff_prediction", "mediumHT_Meff_prediction", NbinsHT_,
         HTmin_, HTmax_, Ntries_, 0.5, Ntries_ + 0.5);
   h_mediumHT_Meff_prediction->Sumw2();

   //// medium MHT
   h_mediumMHT_HT_prediction = fs->make<TH2F> ("mediumMHT_HT_prediction", "mediumMHT_HT_prediction", NbinsHT_, HTmin_,
         HTmax_, Ntries_, 0.5, Ntries_ + 0.5);
   h_mediumMHT_HT_prediction->Sumw2();

   //// high HT
   h_highHT_MHT_prediction = fs->make<TH2F> ("highHT_MHT_prediction", "highHT_MHT_prediction", NbinsMHT_, MHTmin_,
         MHTmax_, Ntries_, 0.5, Ntries_ + 0.5);
   h_highHT_MHT_prediction->Sumw2();
   h_highHT_Meff_prediction = fs->make<TH2F> ("highHT_Meff_prediction", "highHT_Meff_prediction", NbinsHT_, HTmin_,
         HTmax_, Ntries_, 0.5, Ntries_ + 0.5);
   h_highHT_Meff_prediction->Sumw2();

   //// high MHT
   h_highMHT_HT_prediction = fs->make<TH2F> ("highMHT_HT_prediction", "highMHT_HT_prediction", NbinsHT_, HTmin_,
         HTmax_, Ntries_, 0.5, Ntries_ + 0.5);
   h_highMHT_HT_prediction->Sumw2();

   //// very high HT
   h_veryhighHT_MHT_prediction = fs->make<TH2F> ("veryhighHT_MHT_prediction", "veryhighHT_MHT_prediction", NbinsMHT_,
         MHTmin_, MHTmax_, Ntries_, 0.5, Ntries_ + 0.5);
   h_veryhighHT_MHT_prediction->Sumw2();
   h_veryhighHT_Meff_prediction = fs->make<TH2F> ("veryhighHT_Meff_prediction", "veryhighHT_Meff_prediction", NbinsHT_,
         HTmin_, HTmax_, Ntries_, 0.5, Ntries_ + 0.5);
   h_veryhighHT_Meff_prediction->Sumw2();

   //// extreme high HT
   h_extremehighHT_MHT_prediction = fs->make<TH2F> ("extremehighHT_MHT_prediction", "extremehighHT_MHT_prediction",
         NbinsMHT_, MHTmin_, MHTmax_, Ntries_, 0.5, Ntries_ + 0.5);
   h_extremehighHT_MHT_prediction->Sumw2();
   h_extremehighHT_Meff_prediction = fs->make<TH2F> ("extremehighHT_Meff_prediction", "extremehighHT_Meff_prediction",
         NbinsHT_, HTmin_, HTmax_, Ntries_, 0.5, Ntries_ + 0.5);
   h_extremehighHT_Meff_prediction->Sumw2();



	/*************************
	 *  R+S Plots
	 ************************/
   //// load response histos

   //// open root file/tree and create SmearingFunction histo
   TFile *f1 = new TFile(smearingfile_.c_str(), "READ", "", 0);

   smearFunc.resize(3); //// three bins for jet rank
   for (std::vector<std::vector<std::vector<TH1F*> > >::iterator it = smearFunc.begin(); it != smearFunc.end(); ++it) {
      it->resize(EtaBinEdges_.size() - 1);
      for (std::vector<std::vector<TH1F*> >::iterator jt = it->begin(); jt != it->end(); ++jt) {
         jt->resize(PtBinEdges_.size() - 1);
      }
   }

   smearFunc_Core.resize(3); //// three bins for jet rank
   for (std::vector<std::vector<std::vector<TH1F*> > >::iterator it = smearFunc_Core.begin(); it
         != smearFunc_Core.end(); ++it) {
      it->resize(EtaBinEdges_.size() - 1);
      for (std::vector<std::vector<TH1F*> >::iterator jt = it->begin(); jt != it->end(); ++jt) {
         jt->resize(PtBinEdges_.size() - 1);
      }
   }

   smearFunc_LowerTail.resize(3); //// three bins for jet rank
   for (std::vector<std::vector<std::vector<TH1F*> > >::iterator it = smearFunc_LowerTail.begin(); it
         != smearFunc_LowerTail.end(); ++it) {
      it->resize(EtaBinEdges_.size() - 1);
      for (std::vector<std::vector<TH1F*> >::iterator jt = it->begin(); jt != it->end(); ++jt) {
         jt->resize(PtBinEdges_.size() - 1);
      }
   }

   smearFunc_UpperTail.resize(3); //// three bins for jet rank
   for (std::vector<std::vector<std::vector<TH1F*> > >::iterator it = smearFunc_UpperTail.begin(); it
         != smearFunc_UpperTail.end(); ++it) {
      it->resize(EtaBinEdges_.size() - 1);
      for (std::vector<std::vector<TH1F*> >::iterator jt = it->begin(); jt != it->end(); ++jt) {
         jt->resize(PtBinEdges_.size() - 1);
      }
   }

   smearFunc_scaled.resize(3); //// three bins for jet rank
   for (std::vector<std::vector<std::vector<TH1F*> > >::iterator it = smearFunc_scaled.begin(); it
         != smearFunc_scaled.end(); ++it) {
      it->resize(EtaBinEdges_.size() - 1);
      for (std::vector<std::vector<TH1F*> >::iterator jt = it->begin(); jt != it->end(); ++jt) {
         jt->resize(PtBinEdges_.size() - 1);
      }
   }

   SigmaPtHist_scaled.resize(3);
   for (int i_jet = 0; i_jet < 3; ++i_jet) {
      SigmaPtHist_scaled.at(i_jet).resize(EtaBinEdges_.size() - 1);
      for (unsigned int i_eta = 0; i_eta < EtaBinEdges_.size() - 1; ++i_eta) {
         char hname[100];
         sprintf(hname, "SigmaPtHist_scaled_Eta%i_Jet%i", i_eta, i_jet + 1);
         SigmaPtHist_scaled.at(i_jet).at(i_eta) = fs->make<TH1F> (hname, hname, PtBinEdges_.size() - 1,
               &(PtBinEdges_.at(0)));
         SigmaPtHist_scaled.at(i_jet).at(i_eta)->Sumw2();
      }
   }

   SigmaPtHist.resize(3);
   for (int i_jet = 0; i_jet < 3; ++i_jet) {
      SigmaPtHist.at(i_jet).resize(EtaBinEdges_.size() - 1);
      for (unsigned int i_eta = 0; i_eta < EtaBinEdges_.size() - 1; ++i_eta) {
         char hname[100];
         sprintf(hname, "SigmaPtHist_Eta%i_Jet%i", i_eta, i_jet + 1);
         SigmaPtHist.at(i_jet).at(i_eta) = fs->make<TH1F> (hname, hname, PtBinEdges_.size() - 1, &(PtBinEdges_.at(0)));
         SigmaPtHist.at(i_jet).at(i_eta)->Sumw2();
      }
   }

   SigmaPt_scaled.resize(3);
   for (unsigned int i_jet = 0; i_jet < 3; ++i_jet) {
      SigmaPt_scaled.at(i_jet).resize(EtaBinEdges_.size() - 1);
   }

   SigmaPt.resize(3);
   for (unsigned int i_jet = 0; i_jet < 3; ++i_jet) {
      SigmaPt.at(i_jet).resize(EtaBinEdges_.size() - 1);
   }

   //// Fetch histos and fit gaussian core
   for (unsigned int i_Pt = 0; i_Pt < PtBinEdges_.size() - 1; ++i_Pt) {
      for (unsigned int i_eta = 0; i_eta < EtaBinEdges_.size() - 1; ++i_eta) {
         //// Get the histos
         char hname[100];
         sprintf(hname, "%s_Pt%i_Eta%i", inputhist1_.c_str(), i_Pt, i_eta);
         smearFunc.at(0).at(i_eta).at(i_Pt) = (TH1F*) f1->FindObjectAny(hname);
         sprintf(hname, "%s_Pt%i_Eta%i", inputhist2_.c_str(), i_Pt, i_eta);
         smearFunc.at(1).at(i_eta).at(i_Pt) = (TH1F*) f1->FindObjectAny(hname);
         sprintf(hname, "%s_Pt%i_Eta%i", inputhist3p_.c_str(), i_Pt, i_eta);
         smearFunc.at(2).at(i_eta).at(i_Pt) = (TH1F*) f1->FindObjectAny(hname);
         for (unsigned int i_jet = 0; i_jet < 3; ++i_jet) {
            smearFunc.at(i_jet).at(i_eta).at(i_Pt)->Rebin(NRebin_);
            if (probExtreme_ > 0) {
               double p = probExtreme_ * smearFunc.at(i_jet).at(i_eta).at(i_Pt)->Integral();
               smearFunc.at(i_jet).at(i_eta).at(i_Pt)->SetBinContent(1, p);
            }
            //// Get width of gaussian core
            if (smearFunc.at(i_jet).at(i_eta).at(i_Pt)->GetEntries() > 50) {
               double RMS = smearFunc.at(i_jet).at(i_eta).at(i_Pt)->GetRMS();
               double MEAN = smearFunc.at(i_jet).at(i_eta).at(i_Pt)->GetMean();
               TF1* fitfunction = new TF1("f", "gaus(0)", MEAN - 1 * RMS, MEAN + 1 * RMS);
               fitfunction->SetParameters(smearFunc.at(i_jet).at(i_eta).at(i_Pt)->GetMaximum(), MEAN, RMS);
               smearFunc.at(i_jet).at(i_eta).at(i_Pt)->Fit(fitfunction, "LLRQN");
               double Pt = SigmaPtHist_scaled.at(i_jet).at(i_eta)->GetBinCenter(i_Pt);
               double eta = (EtaBinEdges_.at(i_eta) + EtaBinEdges_.at(i_eta + 1)) / 2;
               double f = GetAdditionalSmearing(Pt, eta);
               SigmaPtHist.at(i_jet).at(i_eta)->SetBinContent(i_Pt + 1, std::abs(fitfunction->GetParameter(2)));
               SigmaPtHist.at(i_jet).at(i_eta)->SetBinError(i_Pt + 1, fitfunction->GetParError(2));
               SigmaPtHist_scaled.at(i_jet).at(i_eta)->SetBinContent(i_Pt + 1, std::abs(fitfunction->GetParameter(2))
                     * f);
               SigmaPtHist_scaled.at(i_jet).at(i_eta)->SetBinError(i_Pt + 1, fitfunction->GetParError(2) * f);

               //// Split smearFunc in Core and Tail
               TH1F* hResponseFit = new TH1F(*smearFunc.at(i_jet).at(i_eta).at(i_Pt));
               hResponseFit->Reset();
               for (int i = 0; i < hResponseFit->GetNbinsX(); ++i) {
                  hResponseFit->SetBinContent(i, fitfunction->Eval(hResponseFit->GetBinCenter(i)));
               }

               //// Split lower tail
               smearFunc_LowerTail.at(i_jet).at(i_eta).at(i_Pt) = new TH1F(*smearFunc.at(i_jet).at(i_eta).at(i_Pt));
               smearFunc_LowerTail.at(i_jet).at(i_eta).at(i_Pt)->Add(hResponseFit, -1.);
               for (int i = 0; i <= smearFunc_LowerTail.at(i_jet).at(i_eta).at(i_Pt)->GetNbinsX(); ++i) {
                  double tmp_v = smearFunc_LowerTail.at(i_jet).at(i_eta).at(i_Pt)->GetBinContent(i);
                  double tmp_e = smearFunc_LowerTail.at(i_jet).at(i_eta).at(i_Pt)->GetBinError(i);
                  //// by definition a tail has positive entries
                  if (tmp_v < 0) {
                     tmp_v = 0;
                     tmp_e = 0;
                  } else {
                     //// suppress everything except for low response tail
                     double scale = 1;
                     double x = smearFunc_LowerTail.at(i_jet).at(i_eta).at(i_Pt)->GetBinCenter(i);
                     if (x > MEAN - 1 * RMS)
                        scale = 0.;
                     tmp_v = scale * smearFunc_LowerTail.at(i_jet).at(i_eta).at(i_Pt)->GetBinContent(i);
                     tmp_e = scale * smearFunc_LowerTail.at(i_jet).at(i_eta).at(i_Pt)->GetBinError(i);
                  }
                  smearFunc_LowerTail.at(i_jet).at(i_eta).at(i_Pt)->SetBinContent(i, tmp_v);
                  smearFunc_LowerTail.at(i_jet).at(i_eta).at(i_Pt)->SetBinError(i, tmp_e);
               }

               //// Split upper tail
               smearFunc_UpperTail.at(i_jet).at(i_eta).at(i_Pt) = new TH1F(*smearFunc.at(i_jet).at(i_eta).at(i_Pt));
               smearFunc_UpperTail.at(i_jet).at(i_eta).at(i_Pt)->Add(hResponseFit, -1.);
               for (int i = 0; i <= smearFunc_UpperTail.at(i_jet).at(i_eta).at(i_Pt)->GetNbinsX(); ++i) {
                  double tmp_v = smearFunc_UpperTail.at(i_jet).at(i_eta).at(i_Pt)->GetBinContent(i);
                  double tmp_e = smearFunc_UpperTail.at(i_jet).at(i_eta).at(i_Pt)->GetBinError(i);
                  //// by definition a tail has positive entries
                  if (tmp_v < 0) {
                     tmp_v = 0;
                     tmp_e = 0;
                  } else {
                     //// suppress everything except for low response tail
                     double scale = 1;
                     double x = smearFunc_UpperTail.at(i_jet).at(i_eta).at(i_Pt)->GetBinCenter(i);
                     if (x < MEAN + 1 * RMS)
                        scale = 0.;
                     tmp_v = scale * smearFunc_UpperTail.at(i_jet).at(i_eta).at(i_Pt)->GetBinContent(i);
                     tmp_e = scale * smearFunc_UpperTail.at(i_jet).at(i_eta).at(i_Pt)->GetBinError(i);
                  }
                  smearFunc_UpperTail.at(i_jet).at(i_eta).at(i_Pt)->SetBinContent(i, tmp_v);
                  smearFunc_UpperTail.at(i_jet).at(i_eta).at(i_Pt)->SetBinError(i, tmp_e);
               }

               smearFunc_Core.at(i_jet).at(i_eta).at(i_Pt) = new TH1F(*smearFunc.at(i_jet).at(i_eta).at(i_Pt));
               smearFunc_Core.at(i_jet).at(i_eta).at(i_Pt)->Add(smearFunc_LowerTail.at(i_jet).at(i_eta).at(i_Pt), -1.);
               smearFunc_Core.at(i_jet).at(i_eta).at(i_Pt)->Add(smearFunc_UpperTail.at(i_jet).at(i_eta).at(i_Pt), -1.);

            } else {
               //// Set core and tail if needed
               smearFunc_Core.at(i_jet).at(i_eta).at(i_Pt) = new TH1F(*smearFunc.at(i_jet).at(i_eta).at(i_Pt));
               smearFunc_LowerTail.at(i_jet).at(i_eta).at(i_Pt) = new TH1F(*smearFunc.at(i_jet).at(i_eta).at(i_Pt));
               smearFunc_LowerTail.at(i_jet).at(i_eta).at(i_Pt)->Reset();
               smearFunc_UpperTail.at(i_jet).at(i_eta).at(i_Pt) = new TH1F(*smearFunc.at(i_jet).at(i_eta).at(i_Pt));
               smearFunc_UpperTail.at(i_jet).at(i_eta).at(i_Pt)->Reset();
            }
         }
      }
   }

   //// Fit scaled gaussian sigma as function of pt
   for (unsigned int i_jet = 0; i_jet < 3; ++i_jet) {
      for (unsigned int i_eta = 0; i_eta < EtaBinEdges_.size() - 1; ++i_eta) {
         char fname[100];
         sprintf(fname, "SigmaPtScaled_Eta%i_Jet%i", i_eta, i_jet + 1);
         bool first = false;
         int FirstBin = 1;
         int LastBin = SigmaPtHist_scaled.at(i_jet).at(i_eta)->GetNbinsX();
         for (int j = 1; j <= SigmaPtHist_scaled.at(i_jet).at(i_eta)->GetNbinsX(); ++j) {
            if (!first && SigmaPtHist_scaled.at(i_jet).at(i_eta)->GetBinContent(j) > 0) {
               first = true;
               FirstBin = j;
            }
            if (first && SigmaPtHist_scaled.at(i_jet).at(i_eta)->GetBinContent(j) < 0.001) {
               LastBin = j - 1;
               break;
            }
         }
         SigmaPt_scaled.at(i_jet).at(i_eta) = new TF1(fname,
               "TMath::Sqrt(sign([0])*pow([0]/x,2)+pow([1],2)*pow(x,[2]-1.)+pow([3],2))",
               SigmaPtHist_scaled.at(i_jet).at(i_eta)->GetBinCenter(FirstBin),
               SigmaPtHist_scaled.at(i_jet).at(i_eta)->GetBinCenter(LastBin));
         SigmaPt_scaled.at(i_jet).at(i_eta)->SetParameters(1.2, 0., 0.03);
         SigmaPtHist_scaled.at(i_jet).at(i_eta)->Fit(SigmaPt_scaled.at(i_jet).at(i_eta), "LLRQ");
      }
   }

   //// Fit gaussian sigma as function of pt
   for (unsigned int i_jet = 0; i_jet < 3; ++i_jet) {
      for (unsigned int i_eta = 0; i_eta < EtaBinEdges_.size() - 1; ++i_eta) {
         char fname[100];
         sprintf(fname, "SigmaPt_Eta%i_Jet%i", i_eta, i_jet + 1);
         bool first = false;
         int FirstBin = 1;
         int LastBin = SigmaPtHist.at(i_jet).at(i_eta)->GetNbinsX();
         for (int j = 1; j <= SigmaPtHist.at(i_jet).at(i_eta)->GetNbinsX(); ++j) {
            if (!first && SigmaPtHist.at(i_jet).at(i_eta)->GetBinContent(j) > 0) {
               first = true;
               FirstBin = j;
            }
            if (first && SigmaPtHist.at(i_jet).at(i_eta)->GetBinContent(j) < 0.001) {
               LastBin = j - 1;
               break;
            }
         }
         SigmaPt.at(i_jet).at(i_eta) = new TF1(fname,
               "TMath::Sqrt(sign([0])*pow([0]/x,2)+pow([1],2)*pow(x,[2]-1.)+pow([3],2))", SigmaPtHist.at(i_jet).at(
                     i_eta)->GetBinCenter(FirstBin), SigmaPtHist.at(i_jet).at(i_eta)->GetBinCenter(LastBin));
         SigmaPt.at(i_jet).at(i_eta)->SetParameters(1.2, 0., 0.03);
         SigmaPtHist.at(i_jet).at(i_eta)->Fit(SigmaPt.at(i_jet).at(i_eta), "LLRQ");
      }
   }

   //// Book and fill histograms for smeared and scaled response functions
   for (unsigned int i_jet = 0; i_jet < 3; ++i_jet) {
      for (unsigned int i_Pt = 0; i_Pt < PtBinEdges_.size() - 1; ++i_Pt) {
         for (unsigned int i_eta = 0; i_eta < EtaBinEdges_.size() - 1; ++i_eta) {
            char hname[100];
            sprintf(hname, "SmearedAndScaledResolution_Pt%i_Eta%i_Jet%i", i_Pt, i_eta, i_jet + 1);
            smearFunc_scaled.at(i_jet).at(i_eta).at(i_Pt) = fs->make<TH1F> (hname, hname,
                  smearFunc.at(i_jet).at(i_eta).at(i_Pt)->GetNbinsX(),
                  smearFunc.at(i_jet).at(i_eta).at(i_Pt)->GetXaxis()->GetXmin(),
                  smearFunc.at(i_jet).at(i_eta).at(i_Pt)->GetXaxis()->GetXmax());

            if (smearFunc.at(i_jet).at(i_eta).at(i_Pt)->GetEntries() > 50) {
               //// fold core and tail with additional gaussian
               TH1F smearFunc_Core_tmp(*smearFunc_Core.at(i_jet).at(i_eta).at(i_Pt));
               TH1F smearFunc_LowerTail_tmp(*smearFunc_LowerTail.at(i_jet).at(i_eta).at(i_Pt));
               TH1F smearFunc_UpperTail_tmp(*smearFunc_UpperTail.at(i_jet).at(i_eta).at(i_Pt));
               double AddSmear = GetAdditionalSmearing(PtBinEdges_.at(i_Pt), EtaBinEdges_.at(i_eta));
               if (AddSmear > 1) {
                  //// additional sigma from (1+x)*sigma = sqrt(sigma^2+add_sigma^2)
                  //// or from sigma' = sqrt((sigma'/(1+x))^2+add_sigma^2)
                  double sigma = SigmaPtHist_scaled.at(i_jet).at(i_eta)->GetBinContent(i_Pt + 1);
                  //// if no sigma was fitted use the extrapolation
                  if (sigma == 0)
                     sigma = SigmaPt_scaled.at(i_jet).at(i_eta)->Eval(
                           SigmaPtHist_scaled.at(i_jet).at(i_eta)->GetBinCenter(i_Pt + 1));
                  double AdditionalSigma = TMath::Sqrt(1 - 1 / pow(AddSmear, 2)) * sigma;
                  smearFunc_Core_tmp.Reset();
                  smearFunc_LowerTail_tmp.Reset();
                  smearFunc_UpperTail_tmp.Reset();
                  FoldWithGaussian(*smearFunc_Core.at(i_jet).at(i_eta).at(i_Pt), smearFunc_Core_tmp, AdditionalSigma);
                  FoldWithGaussian(*smearFunc_LowerTail.at(i_jet).at(i_eta).at(i_Pt), smearFunc_LowerTail_tmp,
                        AdditionalSigma);
                  FoldWithGaussian(*smearFunc_UpperTail.at(i_jet).at(i_eta).at(i_Pt), smearFunc_UpperTail_tmp,
                        AdditionalSigma);
               } else if (AddSmear < 1) {
                  smearFunc_Core_tmp.Reset();
                  StretchHisto(*smearFunc_Core.at(i_jet).at(i_eta).at(i_Pt), smearFunc_Core_tmp, AddSmear);
               }

               //// Scale tails
               double LowerTailScale = GetLowerTailScaling(PtBinEdges_.at(i_Pt), EtaBinEdges_.at(i_eta));
               double UpperTailScale = GetUpperTailScaling(PtBinEdges_.at(i_Pt), EtaBinEdges_.at(i_eta));
               //cout << "absolute scaling factor: " << TailScale << endl;
               if (absoluteTailScaling_) {
                  double RMS = smearFunc.at(i_jet).at(i_eta).at(i_Pt)->GetRMS();
                  //cout << "Integral from " << 1-A1RMS_*RMS << " to " << 1-A0RMS_*RMS << endl;

                  //// get integral of tails
                  int i_min = smearFunc.at(i_jet).at(i_eta).at(i_Pt)->FindBin(1 - A1RMS_ * RMS);
                  int i_max = smearFunc.at(i_jet).at(i_eta).at(i_Pt)->FindBin(1 - A0RMS_ * RMS);
                  double RLowerTail = smearFunc_LowerTail.at(i_jet).at(i_eta).at(i_Pt)->Integral(i_min, i_max);
                  double Rcore = smearFunc_Core.at(i_jet).at(i_eta).at(i_Pt)->Integral(i_min, i_max);
                  if (RLowerTail > 0)
                     LowerTailScale = (LowerTailScale * (RLowerTail + Rcore) - Rcore) / RLowerTail;

                  i_min = smearFunc.at(i_jet).at(i_eta).at(i_Pt)->FindBin(1 + A0RMS_ * RMS);
                  i_max = smearFunc.at(i_jet).at(i_eta).at(i_Pt)->FindBin(1 + A1RMS_ * RMS);
                  double RUpperTail = smearFunc_UpperTail.at(i_jet).at(i_eta).at(i_Pt)->Integral(i_min, i_max);
                  Rcore = smearFunc_Core.at(i_jet).at(i_eta).at(i_Pt)->Integral(i_min, i_max);
                  if (RUpperTail > 0)
                     UpperTailScale = (UpperTailScale * (RUpperTail + Rcore) - Rcore) / RUpperTail;

               }
               //cout << "tail scaling factor: " << TailScale << endl;
               smearFunc_scaled.at(i_jet).at(i_eta).at(i_Pt)->Add(&smearFunc_Core_tmp);
               smearFunc_scaled.at(i_jet).at(i_eta).at(i_Pt)->Add(&smearFunc_LowerTail_tmp, LowerTailScale);
               smearFunc_scaled.at(i_jet).at(i_eta).at(i_Pt)->Add(&smearFunc_UpperTail_tmp, UpperTailScale);
            } else {
               //// Replace Histograms with only few entries by gaussians
               double N = 1;
               if (smearFunc.at(i_jet).at(i_eta).at(i_Pt)->GetEntries() > 0) {
                  N = smearFunc.at(i_jet).at(i_eta).at(i_Pt)->Integral();
               }
               cout << "Too few entries for (i_Pt, i_eta, i_jet): " << i_Pt << ", " << i_eta << ", " << i_jet
                     << ", entries = " << smearFunc.at(i_jet).at(i_eta).at(i_Pt)->GetEntries() << endl;
               for (int j = 1; j <= smearFunc_scaled.at(i_jet).at(i_eta).at(i_Pt)->GetNbinsX(); ++j) {
                  double pt = (PtBinEdges_.at(i_Pt) + PtBinEdges_.at(i_Pt + 1)) / 2;
                  double g = N * TMath::Gaus(smearFunc_scaled.at(i_jet).at(i_eta).at(i_Pt)->GetBinCenter(j), 1.,
                        SigmaPt_scaled.at(i_jet).at(i_eta)->Eval(pt));
                  smearFunc_scaled.at(i_jet).at(i_eta).at(i_Pt)->SetBinContent(j, g);
                  smearFunc_Core.at(i_jet).at(i_eta).at(i_Pt)->SetBinContent(j, g);
                  smearFunc_LowerTail.at(i_jet).at(i_eta).at(i_Pt)->SetBinContent(j, 0.);
                  smearFunc_UpperTail.at(i_jet).at(i_eta).at(i_Pt)->SetBinContent(j, 0.);
               }
            }
         }
      }
   }


}

// ------------ method called once each job just after ending the event loop  ------------
void 
Factorization_SmearedJets::endJob() {

	std::cout << __FILE__ << "::" << __FUNCTION__  << "\n" << std::endl;
	std::cout << "[ATM:00] Job settings " << std::endl;
	std::cout << "[ATM:01] Events Processed --- = " << uProcessed << std::endl;
	std::cout << "[ATM:02] Events Passed ------ = " << uPassed << std::endl;
	std::cout << "[ATM:03] Lumi Weighing? ----- = " << doLumiWeighing << std::endl;
	std::cout << "[ATM:04] Event Weighing? ---- = " << doEventWeighing << std::endl;
	std::cout << "[ATM:05] Prescale Weighing? - = " << usePrescaleWeight << std::endl;
	std::cout << "[ATM:31] Minimum HT --------- = " << dMinHT << std::endl;
	std::cout << "[ATM:32] Minimum MHT -------- = " << dMinMHT << std::endl;
	std::cout << "[ATM:40] PASS Summary --------- " << std::endl;
	std::cout << "[ATM:41] Pass NJet ---------- = " << (uProcessed - uFailNjetCut) << std::endl;
	std::cout << "[ATM:42] Pass ht ------------ = " << (uProcessed - uFailNjetCut - uFailMinHTCut) << std::endl;
	std::cout << "[ATM:43] Pass mht ----------- = " << (uProcessed - uFailNjetCut - uFailMinHTCut - uFailMinPFMHTCut)  << std::endl;
	std::cout << "[ATM:50] LumiWeights Avg ---- = ";
	if ( uPassed>0) std::cout << sumLumiWeights/(double)uPassed << std::endl;
	else std::cout << "0" << std::endl;

}

// ------------ method called when starting to processes a run  ------------
bool 
Factorization_SmearedJets::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
Factorization_SmearedJets::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
Factorization_SmearedJets::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
Factorization_SmearedJets::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Factorization_SmearedJets::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void Factorization_SmearedJets::DoDelMinStudy(edm::Handle<std::vector<pat::Jet> > smearedJetsHandle
				, const TLorentzVector &mhtVec
				)
{
	//PrintHeader();
	std::vector<float> vDelPhi_jetmht;
	const float mht = mhtVec.Pt();
	std::vector<unsigned> vLeadJetIndex;

	for (unsigned i = 0 ; i < smearedJetsHandle->size() ; ++i)
	{
		const float pt  = (*smearedJetsHandle)[i].pt();
		const float eta = (*smearedJetsHandle)[i].eta(); 
		if (pt < 50.0 || fabs(eta) > 2.5) continue;
		vLeadJetIndex.push_back(i);
		//std::cout << "i = " <<  i << std::endl;
		const float delphi_jetmht = fabs(TVector2::Phi_mpi_pi((*smearedJetsHandle)[i].phi() - mhtVec.Phi()));
		vDelPhi_jetmht.push_back(delphi_jetmht);
		if (vLeadJetIndex.size() >=3 ) break; //use only three leading jets
	}

	std::sort(vDelPhi_jetmht.begin(), vDelPhi_jetmht.end(), sort_using_less_than);	
	assert (vDelPhi_jetmht.size() == 3 && "ERROR: More/less than 3 dPhiMin calculations found!");
	//std::cout << "vDelPhi_jetmht size = " << vDelPhi_jetmht.size() << std::endl;
	
	const float dPhiMin = vDelPhi_jetmht.at(0);
	const float jet1_pt  = (*smearedJetsHandle)[vLeadJetIndex.at(0)].pt();
	const float jet2_pt  = (*smearedJetsHandle)[vLeadJetIndex.at(1)].pt();
	const float jet3_pt  = (*smearedJetsHandle)[vLeadJetIndex.at(2)].pt();
	const float jet1_phi = (*smearedJetsHandle)[vLeadJetIndex.at(0)].phi();
	const float jet2_phi = (*smearedJetsHandle)[vLeadJetIndex.at(1)].phi();
	const float jet3_phi = (*smearedJetsHandle)[vLeadJetIndex.at(2)].phi();
	const float jet1_eta = (*smearedJetsHandle)[vLeadJetIndex.at(0)].eta();
	const float jet2_eta = (*smearedJetsHandle)[vLeadJetIndex.at(1)].eta();
	const float jet3_eta = (*smearedJetsHandle)[vLeadJetIndex.at(2)].eta();
	const float j1mht_dphi = fabs(TVector2::Phi_mpi_pi(jet1_phi - mhtVec.Phi()));
	const float j2mht_dphi = fabs(TVector2::Phi_mpi_pi(jet2_phi - mhtVec.Phi()));
	const float j3mht_dphi = fabs(TVector2::Phi_mpi_pi(jet3_phi - mhtVec.Phi()));

	pf30_jet1Hist.pt->Fill(jet1_pt, Weight);
	pf30_jet2Hist.pt->Fill(jet2_pt, Weight);
	pf30_jet3Hist.pt->Fill(jet3_pt, Weight);
	pf30_jet1Hist.phi->Fill(jet1_phi, Weight);
	pf30_jet2Hist.phi->Fill(jet2_phi, Weight);
	pf30_jet3Hist.phi->Fill(jet3_phi, Weight);
	pf30_jet1Hist.eta->Fill(jet1_eta, Weight);
	pf30_jet2Hist.eta->Fill(jet2_eta, Weight);
	pf30_jet3Hist.eta->Fill(jet3_eta, Weight);
	pf30_jet1Hist.delphi->Fill(j1mht_dphi, Weight);
	pf30_jet2Hist.delphi->Fill(j2mht_dphi, Weight);
	pf30_jet3Hist.delphi->Fill(j3mht_dphi, Weight);

	hDelPhiMin->Fill(dPhiMin);


	// dphimin in slices of MHT
	if      ( mht >= 60  && mht< 80  ) hDelPhiMin_mht[0]->Fill(dPhiMin, Weight);
	else if ( mht >= 80  && mht< 100 ) hDelPhiMin_mht[1]->Fill(dPhiMin, Weight);
	else if ( mht >= 100 && mht< 120 ) hDelPhiMin_mht[2]->Fill(dPhiMin, Weight);
	else if ( mht >= 120 && mht< 140 ) hDelPhiMin_mht[3]->Fill(dPhiMin, Weight);
	else if ( mht >= 140 && mht< 170 ) hDelPhiMin_mht[4]->Fill(dPhiMin, Weight);
	else if ( mht >= 170 && mht< 200 ) hDelPhiMin_mht[5]->Fill(dPhiMin, Weight);
	else if ( mht >= 200 && mht< 250 ) hDelPhiMin_mht[6]->Fill(dPhiMin, Weight);
	else if ( mht >= 250             ) hDelPhiMin_mht[7]->Fill(dPhiMin, Weight);

	// mht in slices of dphimin
	if      ( dPhiMin <  0.1                  ) MHT_by_phislice[0]->Fill(mht, Weight);
	else if ( dPhiMin >= 0.1 && dPhiMin < 0.2 ) MHT_by_phislice[1]->Fill(mht, Weight);
	else if ( dPhiMin >= 0.2 && dPhiMin < 0.3 ) MHT_by_phislice[2]->Fill(mht, Weight);
	else if ( dPhiMin >= 0.3 && dPhiMin < 0.5 ) MHT_by_phislice[3]->Fill(mht, Weight);
	else if ( dPhiMin >= 0.5 && dPhiMin < 0.8 ) MHT_by_phislice[4]->Fill(mht, Weight);
	else if ( dPhiMin >  0.8                  ) MHT_by_phislice[5]->Fill(mht, Weight);


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
		if (smearedHT >  500.0) hFail[8]->Fill(mht, Weight);
		if (smearedHT >  800.0) hFail[9]->Fill(mht, Weight);
		if (smearedHT > 1000.0) hFail[10]->Fill(mht, Weight);
		if (smearedHT > 1200.0) hFail[11]->Fill(mht, Weight);
		if (smearedHT > 1400.0) hFail[12]->Fill(mht, Weight);

		if (smearedHT >  500.0 && smearedHT <  800.0) hFail[13]->Fill(mht, Weight);
		if (smearedHT >  800.0 && smearedHT < 1000.0) hFail[14]->Fill(mht, Weight);
		if (smearedHT > 1000.0 && smearedHT < 1200.0) hFail[15]->Fill(mht, Weight);
		if (smearedHT > 1200.0 && smearedHT < 1400.0) hFail[16]->Fill(mht, Weight);
	}



	//Pass selection with RA2 dphi cuts
	bool passed = true;
	passed = passed && (std::abs(reco::deltaPhi(jet1_phi, mhtVec.Phi())) > 0.5);
	passed = passed && (std::abs(reco::deltaPhi(jet2_phi, mhtVec.Phi())) > 0.5);
	passed = passed && (std::abs(reco::deltaPhi(jet3_phi, mhtVec.Phi())) > 0.3);

	if (passed) 
	{
		hPass[6]->Fill(mht, Weight);
		if (smearedHT > 500)
		{
			hPass[7]->Fill(mht, Weight);
			if (mht>200) hSignalRegion[0]->Fill(mht,Weight);
			if (mht>350) hSignalRegion[1]->Fill(mht,Weight);
			if (mht>500) hSignalRegion[2]->Fill(mht,Weight);
		}
		if (smearedHT > 800) 
		{
			hPass[8]->Fill(mht, Weight);
			if (mht>200) hSignalRegion[3]->Fill(mht,Weight);
			if (mht>350) hSignalRegion[4]->Fill(mht,Weight);
			if (mht>500) hSignalRegion[5]->Fill(mht,Weight);
		}
		if (smearedHT > 1000) 
		{
			hPass[9]->Fill(mht, Weight);
			if (mht>200) hSignalRegion[6]->Fill(mht,Weight);
			if (mht>350) hSignalRegion[7]->Fill(mht,Weight);
			if (mht>500) hSignalRegion[8]->Fill(mht,Weight);
		}
		if (smearedHT > 1200) 
		{
			hPass[10]->Fill(mht, Weight);
			if (mht>200) hSignalRegion[9]->Fill(mht,Weight);
			if (mht>350) hSignalRegion[10]->Fill(mht,Weight);
			if (mht>500) hSignalRegion[11]->Fill(mht,Weight);
		}
		if (smearedHT > 1400) 
		{
			hPass[11]->Fill(mht, Weight);
			if (mht>200) hSignalRegion[12]->Fill(mht,Weight);
			if (mht>350) hSignalRegion[13]->Fill(mht,Weight);
			if (mht>500) hSignalRegion[14]->Fill(mht,Weight);
		}

		if (smearedHT > 500  && smearedHT < 800 )
		{
			hPass[12]->Fill(mht, Weight);
			if (mht>200) hSignalRegion[15]->Fill(mht,Weight);
			if (mht>350) hSignalRegion[16]->Fill(mht,Weight);
			if (mht>500) hSignalRegion[17]->Fill(mht,Weight);
		}
		if (smearedHT > 800  && smearedHT < 1000)
		{
			hPass[13]->Fill(mht, Weight);
			if (mht>200) hSignalRegion[18]->Fill(mht,Weight);
			if (mht>350) hSignalRegion[19]->Fill(mht,Weight);
			if (mht>500) hSignalRegion[20]->Fill(mht,Weight);
		}
		if (smearedHT > 1000 && smearedHT < 1200)
		{
			hPass[14]->Fill(mht, Weight);
			if (mht>200) hSignalRegion[21]->Fill(mht,Weight);
			if (mht>350) hSignalRegion[22]->Fill(mht,Weight);
			if (mht>500) hSignalRegion[23]->Fill(mht,Weight);
		}
		if (smearedHT > 1200 && smearedHT < 1400)
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
	passed = passed && (std::abs(reco::deltaPhi(jet1_phi, mhtVec.Phi())) > 0.5);
	passed = passed && (std::abs(reco::deltaPhi(jet2_phi, mhtVec.Phi())) > 0.5);
	passed = passed && (std::abs(reco::deltaPhi(jet3_phi, mhtVec.Phi())) > 0.1);
	
	if (! passed)
	{
		if (smearedHT > 500) 
		{
			hFail[17]->Fill(mht, Weight);
			hFail[18]->Fill(mht, Weight);
		}
		if (smearedHT > 600)
		{
			hFail[19]->Fill(mht, Weight);
			hFail[20]->Fill(mht, Weight);
		}

		if (smearedHT > 800) hFail[25]->Fill(mht, Weight);
		if (smearedHT > 1000) hFail[26]->Fill(mht, Weight);
		if (smearedHT > 1200) hFail[27]->Fill(mht, Weight);
		if (smearedHT > 1400) hFail[28]->Fill(mht, Weight);

		if (smearedHT > 500  && smearedHT < 800) hFail[29]->Fill(mht, Weight);
		if (smearedHT > 800  && smearedHT < 1000) hFail[30]->Fill(mht, Weight);
		if (smearedHT > 1000 && smearedHT < 1200) hFail[31]->Fill(mht, Weight);
		if (smearedHT > 1200 && smearedHT < 1400) hFail[32]->Fill(mht, Weight);

	}


	passed = true;
	passed = passed && (std::abs(reco::deltaPhi(jet1_phi, mhtVec.Phi())) > 0.4);
	passed = passed && (std::abs(reco::deltaPhi(jet2_phi, mhtVec.Phi())) > 0.4);
	passed = passed && (std::abs(reco::deltaPhi(jet3_phi, mhtVec.Phi())) > 0.2);

	if (! passed)
	{
		if (smearedHT > 500) 
		{
			hFail[21]->Fill(mht, Weight);
			hFail[22]->Fill(mht, Weight);
		}
		if (smearedHT > 600)
		{
			hFail[23]->Fill(mht, Weight);
			hFail[24]->Fill(mht, Weight);
		}

		if (smearedHT >  800) hFail[33]->Fill(mht, Weight);
		if (smearedHT > 1000) hFail[34]->Fill(mht, Weight);
		if (smearedHT > 1200) hFail[35]->Fill(mht, Weight);
		if (smearedHT > 1400) hFail[36]->Fill(mht, Weight);

		if (smearedHT >  500 && smearedHT <  800) hFail[37]->Fill(mht, Weight);
		if (smearedHT >  800 && smearedHT < 1000) hFail[38]->Fill(mht, Weight);
		if (smearedHT > 1000 && smearedHT < 1200) hFail[39]->Fill(mht, Weight);
		if (smearedHT > 1200 && smearedHT < 1400) hFail[40]->Fill(mht, Weight);
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

void Factorization_SmearedJets::PrintHeader()
{
	std::cout  << "======= Run/Event: " << kRunLumiEvent.run
			<< ":" << kRunLumiEvent.evt << std::endl;
}
void Factorization_SmearedJets::CalcSmearedHTMHT(
				edm::Handle<std::vector<pat::Jet> > smearedJetsHandle, 
				int &NJET, float &HT, float &MHT, float &MEFF, TLorentzVector &MHTVEC)
{
	int njet = 0;
	float ht = 0;
	TLorentzVector tlSumJets(0,0,0,0);
	PrintHeader();

	for (unsigned i = 0 ; i < smearedJetsHandle->size() ; ++i)
	{
		const float pt  = (*smearedJetsHandle)[i].pt();
		const float eta = (*smearedJetsHandle)[i].eta(); 
		const float px  = (*smearedJetsHandle)[i].px();
		const float py  = (*smearedJetsHandle)[i].py();
		const float pz  = (*smearedJetsHandle)[i].pz();
		const float en  = (*smearedJetsHandle)[i].energy();
		const TLorentzVector jetvec(px, py, pz, en);
		if (i<5) std::cout << "i = " <<  i << ":" << pt << "\t" << eta << std::endl;
		if (pt > 50.0 && fabs(eta) < 2.5)
		{
			++njet;
			ht += pt;
		}
		if (pt > 30.0 && fabs(eta) < 5.0) { tlSumJets += jetvec; }
	}

	NJET = njet;
	HT   = ht;
	MHT  = tlSumJets.Pt();
	MEFF = HT + MHT;
	MHTVEC = tlSumJets;

}

void Factorization_SmearedJets::GenerateJets(edm::Event& iEvent, 
				const edm::EventSetup& iSetup)
{

   //Weight
   edm::Handle<double> event_weight;
   iEvent.getByLabel(weightName_, event_weight);
   weight_ = (event_weight.isValid() ? (*event_weight) : 1.0);
   //if (!event_weight.isValid()) cout << "weight not found" << endl;

   if (controlPlots_) {
      h_weight->Fill(log10(weight_));
      h_weightedWeight->Fill(log10(weight_), weight_);
   }

   //GenJets
   edm::Handle<edm::View<reco::GenJet> > gj;
   edm::View<reco::GenJet> Jets_gen;
   bool gj_present = iEvent.getByLabel(genjets_, gj);
   if (gj_present) {
      Jets_gen = *gj;
   }

   //PATJets
   edm::Handle<edm::View<pat::Jet> > Jets;
   iEvent.getByLabel(jets_, Jets);
   edm::View<pat::Jet> Jets_rec = *Jets;

   // collection of rebalanced jets
   std::auto_ptr<vector<pat::Jet> > Jets_reb(new vector<pat::Jet> );

   // collection of smeared jets
   std::auto_ptr<vector<pat::Jet> > Jets_smeared(new vector<pat::Jet> );

   // collection of smeared gen jets
   std::auto_ptr<vector<reco::GenJet> > GenJets_smeared(new vector<reco::GenJet> );

   double HTall_rec = 0;
   double HTlow_rec = 0;
   double HThigh_rec = 0;
   math::PtEtaPhiMLorentzVector vMHTall_rec(0., 0., 0., 0.);
   math::PtEtaPhiMLorentzVector vMHTlow_rec(0., 0., 0., 0.);
   math::PtEtaPhiMLorentzVector vMHThigh_rec(0., 0., 0., 0.);
   //// Fill measured particles to vector
   for (edm::View<pat::Jet>::const_iterator it = Jets_rec.begin(); it != Jets_rec.end(); ++it) {
      if (it->pt() < JetsHTPt_) {
         vMHTall_rec -= it->p4();
         vMHTlow_rec -= it->p4();
         HTall_rec += it->pt();
         HTlow_rec += it->pt();
      } else {
         vMHTall_rec -= it->p4();
         vMHThigh_rec -= it->p4();
         HTall_rec += it->pt();
         HThigh_rec += it->pt();
      }
   }
   if (controlPlots_) {
      h_HTall_rec->Fill(HTall_rec, weight_);
      h_HThigh_rec->Fill(HThigh_rec, weight_);
      h_MHTall_rec->Fill(vMHTall_rec.pt(), weight_);
      h_MHThigh_rec->Fill(vMHThigh_rec.pt(), weight_);
   }

   double HTall_gen = 0;
   double HTlow_gen = 0;
   double HThigh_gen = 0;
   math::PtEtaPhiMLorentzVector vMHTall_gen(0., 0., 0., 0.);
   math::PtEtaPhiMLorentzVector vMHTlow_gen(0., 0., 0., 0.);
   math::PtEtaPhiMLorentzVector vMHThigh_gen(0., 0., 0., 0.);
   //// Fill measured particles to vector
   for (edm::View<reco::GenJet>::const_iterator it = Jets_gen.begin(); it != Jets_gen.end(); ++it) {
      if (it->pt() < JetsHTPt_) {
         vMHTall_gen -= it->p4();
         vMHTlow_gen -= it->p4();
         HTall_gen += it->pt();
         HTlow_gen += it->pt();
      } else {
         vMHTall_gen -= it->p4();
         vMHThigh_gen -= it->p4();
         HTall_gen += it->pt();
         HThigh_gen += it->pt();
      }
   }
   if (controlPlots_) {
      h_HTall_gen->Fill(HTall_gen, weight_);
      h_HThigh_gen->Fill(HThigh_gen, weight_);
      h_MHTall_gen->Fill(vMHTall_gen.pt(), weight_);
      h_MHThigh_gen->Fill(vMHThigh_gen.pt(), weight_);
   }
   //cout << "HT gen  = " << HThigh_gen << endl;

   //// Pt resolution of reconstructed CaloJets
   if (gj_present) {
      for (edm::View<reco::GenJet>::const_iterator it = Jets_gen.begin(); it != Jets_gen.end(); ++it) {
         double dRmin = 100;
         const reco::Jet* matchedJet = 0;
         for (edm::View<pat::Jet>::const_iterator jt = Jets_rec.begin(); jt != Jets_rec.end(); ++jt) {
            double dR = deltaR(*jt, *it);
            if (dR < dRmin) {
               dRmin = dR;
               matchedJet = &(*jt);
            }
         }
         if (controlPlots_) {
            if (dRmin < 0.15) {
               if (fabs(it->eta()) < 1.5)
                  h_RecJetRes_Pt->Fill(it->pt(), matchedJet->pt() / it->pt(), weight_);
               if (it->pt() > 100.)
                  h_RecJetRes_Eta->Fill(it->eta(), matchedJet->pt() / it->pt(), weight_);
            }
         }
      }
   }

   Jets_reb->reserve(Jets_rec.size());
   Jets_smeared->reserve(Jets_rec.size());

   //
   // Rebalance multi jet system
   //
   if (smearCollection_ == "Reco") {
      bool isRebalanced = RebalanceJets_KinFitter(&Jets_rec, *(Jets_reb.get()));
      if (!isRebalanced) {
         cout << "Bad event: Not possible to rebalance!" << endl;
         weight_ = 0;
      }
      double HThigh_reb = 0;
      math::PtEtaPhiMLorentzVector vMHThigh_reb(0., 0., 0., 0.);
      double HTlow_reb = 0;
      math::PtEtaPhiMLorentzVector vMHTlow_reb(0., 0., 0., 0.);
      double HTall_reb = 0;
      math::PtEtaPhiMLorentzVector vMHTall_reb(0., 0., 0., 0.);
      for (vector<pat::Jet>::const_iterator it = Jets_reb-> begin(); it != Jets_reb->end(); ++it) {
         if (it->pt() > JetsHTPt_) {
            vMHThigh_reb -= it->p4();
            HThigh_reb += it->pt();
            vMHTall_reb -= it->p4();
            HTall_reb += it->pt();
         } else {
            vMHTlow_reb -= it->p4();
            HTlow_reb += it->pt();
            vMHTall_reb -= it->p4();
            HTall_reb += it->pt();
         }
      }
      if (!isRebalanced) {
         cout << "Bad event: Can't be rebalanced!!!" << endl;
         cout << "Reconstructed: HT, MHT = " << HThigh_rec << ", " << vMHThigh_rec.pt() << endl;
         cout << "Rebalanced: HT, MHT = " << HTall_reb << ", " << vMHTall_reb.pt() << endl;
      }
      if (controlPlots_) {
         h_HTall_reb->Fill(HTall_reb, weight_);
         h_HThigh_reb->Fill(HThigh_reb, weight_);
         //         if (abs(HThigh_reb - HThigh_rec) > 100) {
         //            cout << "WARNING!!!!!!!!!!!!!!!!!!!" << endl;
         //            cout << uncertaintyName_ << ": HTall reb = " << HTall_reb << " HThigh reb = " << HThigh_reb
         //                  << " HTall rec = " << HTall_rec << " HThigh rec = " << HThigh_rec << " weight = " << weight_ << endl;
         //         }
         h_MHTall_reb->Fill(vMHTall_reb.pt(), weight_);
         h_MHThigh_reb->Fill(vMHThigh_reb.pt(), weight_);
      }

      //// Pt resolution of rebalanced CaloJets
      if (gj_present) {
         for (edm::View<reco::GenJet>::const_iterator it = Jets_gen.begin(); it != Jets_gen.end(); ++it) {
            double dRmin = 100;
            const pat::Jet* matchedJet = 0;
            for (vector<pat::Jet>::const_iterator jt = Jets_reb-> begin(); jt != Jets_reb->end(); ++jt) {
               double dR = deltaR(*jt, *it);
               if (dR < dRmin) {
                  dRmin = dR;
                  matchedJet = &(*jt);
               }
            }
            if (controlPlots_) {
               if (dRmin < 0.15) {
                  if (fabs(it->eta()) < 1.5)
                     h_RebJetRes_Pt->Fill(it->pt(), matchedJet->pt() / it->pt(), weight_);
                  if (it->pt() > 100.)
                     h_RebJetRes_Eta->Fill(it->eta(), matchedJet->pt() / it->pt(), weight_);
               }
            }
         }
      }
   }

   //
   // Smear rebalanced multi jet system
   //
   if (smearCollection_ == "Reco") {
      SmearingJets(*(Jets_reb.get()), *(Jets_smeared.get()));
   } else if (smearCollection_ == "Gen") {
      SmearingGenJets(&Jets_gen, *(GenJets_smeared.get()));
   }

   double HThigh_smeared = 0;
   math::PtEtaPhiMLorentzVector vMHThigh_smeared(0., 0., 0., 0.);
   if (smearCollection_ == "Reco") {
      for (vector<pat::Jet>::const_iterator it = Jets_smeared-> begin(); it != Jets_smeared->end(); ++it) {
         vMHThigh_smeared -= it->p4();
         HThigh_smeared += it->pt();
      }
   } else if (smearCollection_ == "Gen") {
      for (vector<reco::GenJet>::const_iterator it = GenJets_smeared-> begin(); it != GenJets_smeared->end(); ++it) {
         vMHThigh_smeared -= it->p4();
         HThigh_smeared += it->pt();
      }
   }
   if (controlPlots_) {
      h_HTall_smeared->Fill(HTlow_rec + HThigh_smeared, weight_);
      h_HThigh_smeared->Fill(HThigh_smeared, weight_);
      h_MHTall_smeared->Fill((vMHTlow_rec + vMHThigh_smeared).pt(), weight_);
      h_MHThigh_smeared->Fill(vMHThigh_smeared.pt(), weight_);
   }

   //// Pt resolution of smeared Jets
   if (smearCollection_ == "Gen") {
      if (gj_present) {
         for (edm::View<reco::GenJet>::const_iterator it = Jets_gen.begin(); it != Jets_gen.end(); ++it) {
            double dRmin = 100;
            const reco::GenJet* matchedJet = 0;
            for (vector<reco::GenJet>::const_iterator jt = GenJets_smeared-> begin(); jt != GenJets_smeared->end(); ++jt) {
               double dR = deltaR(*jt, *it);
               if (dR < dRmin) {
                  dRmin = dR;
                  matchedJet = &(*jt);
               }
            }
            if (controlPlots_) {
               if (dRmin < 0.15) {
                  if (fabs(it->eta()) < 1.5)
                     h_SmearedJetRes_Pt->Fill(it->pt(), matchedJet->pt() / it->pt(), weight_);
                  if (it->pt() > 100.)
                     h_SmearedJetRes_Eta->Fill(it->eta(), matchedJet->pt() / it->pt(), weight_);
               }
            }
         }
      }
   } else if (smearCollection_ == "Reco") {
      if (gj_present) {
         for (edm::View<reco::GenJet>::const_iterator it = Jets_gen.begin(); it != Jets_gen.end(); ++it) {
            double dRmin = 100;
            const pat::Jet* matchedJet = 0;
            for (vector<pat::Jet>::const_iterator jt = Jets_smeared-> begin(); jt != Jets_smeared->end(); ++jt) {
               double dR = deltaR(*jt, *it);
               if (dR < dRmin) {
                  dRmin = dR;
                  matchedJet = &(*jt);
               }
            }
            if (controlPlots_) {
               if (dRmin < 0.15) {
                  if (fabs(it->eta()) < 1.5)
                     h_SmearedJetRes_Pt->Fill(it->pt(), matchedJet->pt() / it->pt(), weight_);
                  if (it->pt() > 100.)
                     h_SmearedJetRes_Eta->Fill(it->eta(), matchedJet->pt() / it->pt(), weight_);
               }
            }
         }
      }
   }

   // put products into event

   if (smearCollection_ == "Reco") {
      GreaterByPt<pat::Jet> ptComparator_;
      std::sort(Jets_reb->begin(), Jets_reb->end(), ptComparator_);
      //iEvent.put(Jets_reb, jets_reb_ + uncertaintyName_);
      iEvent.put(Jets_reb, jets_reb_);
      std::sort(Jets_smeared->begin(), Jets_smeared->end(), ptComparator_);
      //iEvent.put(Jets_smeared, jets_smeared_ + uncertaintyName_);
      iEvent.put(Jets_smeared, jets_smeared_);
   }

   if (smearCollection_ == "Gen") {

      std::vector<pat::Jet> * sJets = new std::vector<pat::Jet>();
      std::auto_ptr<std::vector<pat::Jet> > smearedJets(sJets);

      edm::View<reco::GenJet>::const_iterator it = Jets_gen.begin();
      for (vector<reco::GenJet>::const_iterator jt = GenJets_smeared->begin(); jt != GenJets_smeared->end() && it
            != Jets_gen.end(); ++jt, ++it) {

         reco::GenJet Jet = *jt;
         reco::GenJet oJet = *it;

         //pat::Jet myNewJet(pat::Jet(reco::CaloJet(Jet.p4(), math::XYZPoint(0., 0., 0.), reco::CaloJet::Specific())));
         pat::Jet myNewJet(pat::Jet(reco::Jet(Jet.p4(), math::XYZPoint(0., 0., 0.))));
         myNewJet.setJetCharge(oJet.pt());

         sJets->push_back(myNewJet);
      }

      //iEvent.put(GenJets_smeared, genjets_smeared_);
      //iEvent.put(smearedJets, jets_smeared_);

   }
}


//--------------------------------------------------------------------------
// rebalance the events using a kinematic fit and transverse momentum balance
bool Factorization_SmearedJets::RebalanceJets_KinFitter(edm::View<pat::Jet>* Jets_rec, std::vector<pat::Jet> &Jets_reb) {

   bool result = true;

   //// Interface to KinFitter
   TKinFitter* myFit = new TKinFitter();

   std::vector<TLorentzVector*> lvec_m;

   std::vector<TMatrixD*> covMat_m;

   std::vector<TFitParticleEtEtaPhi*> fitted;
   std::vector<TFitParticleEtEtaPhi*> measured;
   std::map<int, const pat::Jet*> JetMap;
   double dPx = 0;
   double dPy = 0;
   double HTreco = 0;
   double HTreb = 0;
   double MHTx_low = 0;
   double MHTy_low = 0;
   double MHTx_high = 0;
   double MHTy_high = 0;

   //// Fill measured particles to vector
   int i = 0;
   for (edm::View<pat::Jet>::const_iterator it = Jets_rec-> begin(); it != Jets_rec->end(); ++it) {

      //if (it->pt() < rebalancedJetPt_ || abs(it->pt()) > 3.0) {
      //if (it->pt() < rebalancedJetPt_ || it->chargedEmEnergyFraction()>0.9 || it->muonEnergyFraction()>0.9) {
      if (it->pt() < rebalancedJetPt_) {
         if (rebalanceMode_ == "MHTall") {
            MHTx_low -= it->px();
            MHTy_low -= it->py();
            pat::Jet rebalancedJet(*((pat::Jet*) &(*it)));
            Jets_reb.push_back(rebalancedJet);
         }
      } else {
         MHTx_high -= it->px();
         MHTy_high -= it->py();
         JetMap[i] = &(*it);

         // The particles before fitting
         double tmppx, tmppy, tmppz, tmpe;
         tmppx = it->px();
         tmppy = it->py();
         tmppz = it->pz();
         tmpe = it->energy();

         TLorentzVector* lv = new TLorentzVector(tmppx, tmppy, tmppz, tmpe);
         lvec_m.push_back(lv);
         TMatrixD* cM = new TMatrixD(3, 3);
         (*cM)(0, 0) = JetResolution_Pt2(it->pt(), it->eta(), i);
         (*cM)(1, 1) = JetResolution_Eta2(it->energy(), it->eta());
         (*cM)(2, 2) = JetResolution_Phi2(it->energy(), it->eta());
         covMat_m.push_back(cM);
         char name[10];
         sprintf(name, "jet%i", i);
         TFitParticleEtEtaPhi* jet1 = new TFitParticleEtEtaPhi(name, name, lvec_m.back(), covMat_m.back());
         measured.push_back(jet1);
         TFitParticleEtEtaPhi* jet2 = new TFitParticleEtEtaPhi(name, name, lvec_m.back(), covMat_m.back());
         fitted.push_back(jet2);
         myFit->addMeasParticle(fitted.back());
         ++i;
      }
   }

   //// Add momentum constraints
   double MET_constraint_x = 0.;
   double MET_constraint_y = 0.;
   if (rebalanceMode_ == "MHTall") {
      //// rebalance MHT of all jets
      MET_constraint_x = MHTx_low;
      MET_constraint_y = MHTy_low;
   } else if (rebalanceMode_ == "MHThigh") {
      //// rebalance MHT of fitted jets
      MET_constraint_x = 0.;
      MET_constraint_y = 0.;
   } else {
      //// default: rebalance MHT of fitted jets
      MET_constraint_x = 0.;
      MET_constraint_y = 0.;
   }
   TFitConstraintEp* momentumConstr1 = new TFitConstraintEp("px", "px", 0, TFitConstraintEp::pX, MET_constraint_x);
   TFitConstraintEp* momentumConstr2 = new TFitConstraintEp("py", "py", 0, TFitConstraintEp::pY, MET_constraint_y);
   for (unsigned int i = 0; i < fitted.size(); ++i) {
      momentumConstr1->addParticle(fitted.at(i));
      momentumConstr2->addParticle(fitted.at(i));
   }
   myFit->addConstraint(momentumConstr1);
   myFit->addConstraint(momentumConstr2);

   //// Set fit parameters
   myFit->setVerbosity(0);
   myFit->setMaxNbIter(100);
   myFit->setMaxF(0.01 * 2);
   myFit->setMaxDeltaS(1.e-3);
   myFit->fit();
   //cout << "KinFitter: " << myFit->getStatus() << endl;
   int status = myFit->getStatus();

   double chi2 = 0;
   double F = 0;
   double prob = 0;
   //if (status == 0 || status == 1) {
   if (status == 0) {
      chi2 = myFit->getS();
      F = myFit->getF();
      int dof = myFit->getNDF();
      prob = TMath::Prob(chi2, dof);
      //if (prob < 1.e-8) result = false;
      //cout << "chi2, prop, F = " << chi2 << " " << prob << " " << F << endl;
   } else {
      chi2 = 99999;
      prob = 0;
      F = 99999;
      result = false;
      //cout << "chi2, prop, F = " << chi2 << " " << prob << " " << F << endl;
   }
   if (controlPlots_)
      h_fitProb->Fill(prob);

   //// Get the output of KinFitter
   for (unsigned int i = 0; i < measured.size(); ++i) {
      // create new rebalanced Jet
      pat::Jet::LorentzVector newP4(fitted.at(i)->getCurr4Vec()->Px(), fitted.at(i)->getCurr4Vec()->Py(),
            fitted.at(i)->getCurr4Vec()->Pz(), fitted.at(i)->getCurr4Vec()->E());
      pat::Jet rebalancedJet(*((pat::Jet*) JetMap[i]));
      HTreco += rebalancedJet.pt();
      //cout << "RECO: " << i << "th: pt = " << rebalancedJet.pt() << " phi = " << rebalancedJet.phi() << endl;
      rebalancedJet.setP4(newP4);
      //cout << "REB: " << i << "th: pt = " << rebalancedJet.pt() << " phi = " << rebalancedJet.phi() << endl;
      Jets_reb.push_back(rebalancedJet);
      dPx -= newP4.Px() - measured.at(i)->getCurr4Vec()->Px();
      dPy -= newP4.Py() - measured.at(i)->getCurr4Vec()->Py();
      HTreb += rebalancedJet.pt();
   }
   //cout << "HT reco = " << HTreco << endl;
   //cout << "HT reb  = " << HTreb << endl;

   delete myFit;
   for (unsigned int i = 0; i < measured.size(); ++i) {
      delete lvec_m.at(i);
      delete covMat_m.at(i);
      delete measured.at(i);
      delete fitted.at(i);
   }
   delete momentumConstr1;
   delete momentumConstr2;

   return result;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
void Factorization_SmearedJets::SmearingJets(const std::vector<pat::Jet> &Jets_reb, std::vector<pat::Jet> &Jets_smeared) {

   double dPx = 0;
   double dPy = 0;

   for (int i = 1; i <= Ntries_; ++i) {
      int Ntries2 = 1;
      double w = weight_;
      if (cleverPrescaleTreating_ == true && weight_ > 1) {
         Ntries2 = (int) weight_;
         w = weight_ / Ntries2;
      }
      for (int j = 1; j <= Ntries2; ++j) {
         Jets_smeared.clear();
         int i_jet = 0;
         for (std::vector<pat::Jet>::const_iterator it = Jets_reb.begin(); it != Jets_reb.end(); ++it) {
            if (it->pt() > smearedJetPt_) {
               double newPt = 0;
               newPt = it->pt() * JetResolutionHist_Pt_Smear(it->pt(), it->eta(), i_jet);
               //double newEta = rand_->Gaus(it->eta(), TMath::Sqrt(JetResolution_Eta2(it->energy(), it->eta())));
               //double newPhi = rand_->Gaus(it->phi(), TMath::Sqrt(JetResolution_Phi2(it->energy(), it->eta())));
               double newEta = it->eta();
               double newPhi = it->phi();
               pat::Jet::PolarLorentzVector newP4(newPt, newEta, newPhi, it->mass());
               pat::Jet smearedJet(*it);
               smearedJet.setP4(newP4);
               Jets_smeared.push_back(smearedJet);
               dPx -= newP4.Px() - it->px();
               dPy -= newP4.Py() - it->py();
               ++i_jet;
            } else {
               pat::Jet smearedJet(*it);
               Jets_smeared.push_back(smearedJet);
            }
         }
         GreaterByPt<reco::Candidate> ptComparator_;
         std::sort(Jets_smeared.begin(), Jets_smeared.end(), ptComparator_);
         //Fill HT and MHT prediction histos for i-th iteration of smearing
         FillPredictionHistos(Jets_smeared, i, w);
      }
   }

   return;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
void Factorization_SmearedJets::SmearingGenJets(edm::View<reco::GenJet>* Jets_gen, std::vector<reco::GenJet> &GenJets_smeared) {

   double dPx = 0;
   double dPy = 0;

   for (int i = 1; i <= Ntries_; ++i) {
      GenJets_smeared.clear();
      int i_jet = 0;

      for (edm::View<reco::GenJet>::const_iterator it = Jets_gen->begin(); it != Jets_gen->end(); ++it) {

         if (it->pt() > smearedJetPt_) {
            double newPt = 0;
            newPt = it->pt() * JetResolutionHist_Pt_Smear(it->pt(), it->eta(), i_jet);
            //double newEta = rand_->Gaus(it->eta(), TMath::Sqrt(JetResolution_Eta2(it->energy(), it->eta())));
            //double newPhi = rand_->Gaus(it->phi(), TMath::Sqrt(JetResolution_Phi2(it->energy(), it->eta())));
            double newEta = it->eta();
            double newPhi = it->phi();
            reco::GenJet::PolarLorentzVector newP4(newPt, newEta, newPhi, it->mass());
            reco::GenJet smearedJet(*it);
            smearedJet.setP4(newP4);
            GenJets_smeared.push_back(smearedJet);
            dPx -= newP4.Px() - it->px();
            dPy -= newP4.Py() - it->py();
            ++i_jet;
         } else {
            reco::GenJet smearedJet(*it);
            GenJets_smeared.push_back(smearedJet);
         }
      }
      GreaterByPt<reco::Candidate> ptComparator_;
      std::sort(GenJets_smeared.begin(), GenJets_smeared.end(), ptComparator_);
      //Fill HT and MHT prediction histos for i-th iteration of smearing
      FillPredictionHistos_gen(GenJets_smeared, i, weight_);
   }

   return;
}
//--------------------------------------------------------------------------
// pt resolution for KinFitter
double Factorization_SmearedJets::JetResolution_Pt2(const double& pt, const double& eta, const int& i) {
   int i_jet;
   i < 2 ? i_jet = i : i_jet = 2;
   int i_eta = GetIndex(eta, &EtaBinEdges_);
   //return pow(pt, 2) * pow(SigmaPt_scaled.at(i_jet).at(i_eta)->Eval(pt), 2);
   return pow(pt, 2) * pow(SigmaPt.at(i_jet).at(i_eta)->Eval(pt), 2);
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// relative pt resolution for KinFitter
double Factorization_SmearedJets::JetResolution_Ptrel(const double& pt, const double& eta, const int& i) {
   int i_jet;
   i < 2 ? i_jet = i : i_jet = 2;
   int i_eta = GetIndex(eta, &EtaBinEdges_);
   return SigmaPt_scaled.at(i_jet).at(i_eta)->Eval(pt);
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// eta resolution for KinFitter
double Factorization_SmearedJets::JetResolution_Eta2(const double& e, const double& eta) {
   //may be artifically reduced (no angular fit)
   return (pow(0.05 / TMath::Sqrt(e), 2) + pow(0.005, 2)) / 1.e6;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// phi resolution for KinFitter
double Factorization_SmearedJets::JetResolution_Phi2(const double& e, const double& eta) {
   //may be artifically reduced (no angular fit)
   return (pow(0.05 / TMath::Sqrt(e), 2) + pow(0.005, 2)) / 1.e6;
}

//--------------------------------------------------------------------------
void Factorization_SmearedJets::FillPredictionHistos_gen(const std::vector<reco::GenJet>& Jets_smeared, const int& i,
      const double& w) {

   int NJets = calcNJets_gen(Jets_smeared);
   double HT = calcHT_gen(Jets_smeared);
   math::PtEtaPhiMLorentzVector vMHT = calcMHT_gen(Jets_smeared);
   double MHT = vMHT.pt();
   bool minDeltaPhi = calcMinDeltaPhi_gen(Jets_smeared, vMHT);
   double Meff = HT + MHT;

   if (NJets >= NJets_) {

      h_presel_MHT_prediction->Fill(MHT, double(i), w);
      h_presel_HT_prediction->Fill(HT, double(i), w);
      h_presel_Meff_prediction->Fill(Meff, double(i), w);

      if (HT > HTcut_low_)
         h_baselineNoDeltaPhi_MHT_prediction->Fill(MHT, double(i), w);
      if (HT > HTcut_low_ && MHT > MHTcut_low_)
         h_baselineNoDeltaPhi_Meff_prediction->Fill(Meff, double(i), w);
      if (MHT > MHTcut_low_)
         h_baselineNoDeltaPhi_HT_prediction->Fill(HT, double(i), w);

      if (minDeltaPhi) {

         if (HT > HTcut_low_)
            h_lowHT_MHT_prediction->Fill(MHT, double(i), w);
         if (HT > HTcut_low_ && MHT > MHTcut_low_)
            h_lowHT_Meff_prediction->Fill(Meff, double(i), w);
         if (MHT > MHTcut_low_)
            h_lowMHT_HT_prediction->Fill(HT, double(i), w);

         if (HT > HTcut_medium_)
            h_mediumHT_MHT_prediction->Fill(MHT, double(i), w);
         if (HT > HTcut_medium_ && MHT > MHTcut_low_)
            h_mediumHT_Meff_prediction->Fill(Meff, double(i), w);
         if (MHT > MHTcut_medium_)
            h_mediumMHT_HT_prediction->Fill(HT, double(i), w);

         if (HT > HTcut_high_)
            h_highHT_MHT_prediction->Fill(MHT, double(i), w);
         if (HT > HTcut_high_ && MHT > MHTcut_low_)
            h_highHT_Meff_prediction->Fill(Meff, double(i), w);
         if (MHT > MHTcut_high_)
            h_highMHT_HT_prediction->Fill(HT, double(i), w);

         if (HT > HTcut_veryhigh_)
            h_veryhighHT_MHT_prediction->Fill(MHT, double(i), w);
         if (HT > HTcut_veryhigh_ && MHT > MHTcut_low_)
            h_veryhighHT_Meff_prediction->Fill(Meff, double(i), w);

         if (HT > HTcut_extremehigh_)
            h_extremehighHT_MHT_prediction->Fill(MHT, double(i), w);
         if (HT > HTcut_extremehigh_ && MHT > MHTcut_low_)
            h_extremehighHT_Meff_prediction->Fill(Meff, double(i), w);

      }

   }

   return;
}

//--------------------------------------------------------------------------
// pt resolution for smearing
double Factorization_SmearedJets::JetResolutionHist_Pt_Smear(const double& pt, const double& eta, const int& i) {
   int i_jet;
   i < 2 ? i_jet = i : i_jet = 2;
   int i_Pt = GetIndex(pt, &PtBinEdges_);
   int i_eta = GetIndex(eta, &EtaBinEdges_);

   double res = smearFunc_scaled.at(i_jet).at(i_eta).at(i_Pt)->GetRandom();

   return res;
}

//--------------------------------------------------------------------------
void Factorization_SmearedJets::FillPredictionHistos(const std::vector<pat::Jet>& Jets_smeared, const int& i, const double& w) {

   int NJets = calcNJets(Jets_smeared);
   double HT = calcHT(Jets_smeared);
   math::PtEtaPhiMLorentzVector vMHT = calcMHT(Jets_smeared);
   double MHT = vMHT.pt();
   bool minDeltaPhi = calcMinDeltaPhi(Jets_smeared, vMHT);
   double Meff = HT + MHT;

   if (NJets >= NJets_) {

      h_presel_MHT_prediction->Fill(MHT, double(i), w);
      h_presel_HT_prediction->Fill(HT, double(i), w);
      h_presel_Meff_prediction->Fill(Meff, double(i), w);

      if (HT > HTcut_low_)
         h_baselineNoDeltaPhi_MHT_prediction->Fill(MHT, double(i), w);
      if (HT > HTcut_low_ && MHT > MHTcut_low_)
         h_baselineNoDeltaPhi_Meff_prediction->Fill(Meff, double(i), w);
      if (MHT > MHTcut_low_)
         h_baselineNoDeltaPhi_HT_prediction->Fill(HT, double(i), w);

      if (minDeltaPhi) {

         if (HT > HTcut_low_)
            h_lowHT_MHT_prediction->Fill(MHT, double(i), w);
         if (HT > HTcut_low_ && MHT > MHTcut_low_)
            h_lowHT_Meff_prediction->Fill(Meff, double(i), w);
         if (MHT > MHTcut_low_)
            h_lowMHT_HT_prediction->Fill(HT, double(i), w);

         if (HT > HTcut_medium_)
            h_mediumHT_MHT_prediction->Fill(MHT, double(i), w);
         if (HT > HTcut_medium_ && MHT > MHTcut_low_)
            h_mediumHT_Meff_prediction->Fill(Meff, double(i), w);
         if (MHT > MHTcut_medium_)
            h_mediumMHT_HT_prediction->Fill(HT, double(i), w);

         if (HT > HTcut_high_)
            h_highHT_MHT_prediction->Fill(MHT, double(i), w);
         if (HT > HTcut_high_ && MHT > MHTcut_low_)
            h_highHT_Meff_prediction->Fill(Meff, double(i), w);
         if (MHT > MHTcut_high_)
            h_highMHT_HT_prediction->Fill(HT, double(i), w);

         if (HT > HTcut_veryhigh_)
            h_veryhighHT_MHT_prediction->Fill(MHT, double(i), w);
         if (HT > HTcut_veryhigh_ && MHT > MHTcut_low_)
            h_veryhighHT_Meff_prediction->Fill(Meff, double(i), w);

         if (HT > HTcut_extremehigh_)
            h_extremehighHT_MHT_prediction->Fill(MHT, double(i), w);
         if (HT > HTcut_extremehigh_ && MHT > MHTcut_low_)
            h_extremehighHT_Meff_prediction->Fill(Meff, double(i), w);

      }

   }

   return;
}
//--------------------------------------------------------------------------
int Factorization_SmearedJets::GetIndex(const double& x, const std::vector<double>* vec) {
   int index = -1;
   // this is a check
   //int index = 0;
   for (std::vector<double>::const_iterator it = vec->begin(); it != vec->end(); ++it) {
      if ((*it) > fabs(x))
         break;
      ++index;
   }
   if (index < 0)
      index = 0;
   if (index > (int) vec->size() - 2)
      index = vec->size() - 2;

   return index;
}
//--------------------------------------------------------------------------
int Factorization_SmearedJets::calcNJets(const std::vector<pat::Jet>& Jets_smeared) {
   int NJets = 0;
   for (vector<pat::Jet>::const_iterator it = Jets_smeared.begin(); it != Jets_smeared.end(); ++it) {
      if (it->pt() > JetsHTPt_ && std::abs(it->eta()) < JetsHTEta_) {
         ++NJets;
      }
   }
   return NJets;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
double Factorization_SmearedJets::calcHT(const std::vector<pat::Jet>& Jets_smeared) {
   double HT = 0;
   for (vector<pat::Jet>::const_iterator it = Jets_smeared.begin(); it != Jets_smeared.end(); ++it) {
      if (it->pt() > JetsHTPt_ && std::abs(it->eta()) < JetsHTEta_) {
         HT += it->pt();
      }
   }
   return HT;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
bool Factorization_SmearedJets::calcMinDeltaPhi(const std::vector<pat::Jet>& Jets_smeared, math::PtEtaPhiMLorentzVector& MHT) {
   bool result = true;
   unsigned int i = 0;
   for (vector<pat::Jet>::const_iterator it = Jets_smeared.begin(); it != Jets_smeared.end(); ++it) {
      if (it->pt() > JetsMHTPt_ && std::abs(it->eta()) < JetsMHTEta_) {
         if (i < JetDeltaMin_.size()) {
            if (std::abs(deltaPhi(MHT, *it)) < JetDeltaMin_.at(i))
               result = false;
            ++i;
         } else {
            break;
         }
      }
   }
   return result;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
math::PtEtaPhiMLorentzVector Factorization_SmearedJets::calcMHT(const std::vector<pat::Jet>& Jets_smeared) {
   math::PtEtaPhiMLorentzVector MHT(0, 0, 0, 0);
   for (vector<pat::Jet>::const_iterator it = Jets_smeared.begin(); it != Jets_smeared.end(); ++it) {
      if (it->pt() > JetsMHTPt_ && std::abs(it->eta()) < JetsMHTEta_) {
         MHT -= it->p4();
      }
   }
   return MHT;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
int Factorization_SmearedJets::calcNJets_gen(const std::vector<reco::GenJet>& Jets_smeared) {
   int NJets = 0;
   for (vector<reco::GenJet>::const_iterator it = Jets_smeared.begin(); it != Jets_smeared.end(); ++it) {
      if (it->pt() > JetsHTPt_ && std::abs(it->eta()) < JetsHTEta_) {
         ++NJets;
      }
   }
   return NJets;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
double Factorization_SmearedJets::calcHT_gen(const std::vector<reco::GenJet>& Jets_smeared) {
   double HT = 0;
   for (vector<reco::GenJet>::const_iterator it = Jets_smeared.begin(); it != Jets_smeared.end(); ++it) {
      if (it->pt() > JetsHTPt_ && std::abs(it->eta()) < JetsHTEta_) {
         HT += it->pt();
      }
   }
   return HT;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
bool Factorization_SmearedJets::calcMinDeltaPhi_gen(const std::vector<reco::GenJet>& Jets_smeared, math::PtEtaPhiMLorentzVector& MHT) {
   bool result = true;
   unsigned int i = 0;
   for (vector<reco::GenJet>::const_iterator it = Jets_smeared.begin(); it != Jets_smeared.end(); ++it) {
      if (it->pt() > JetsMHTPt_ && std::abs(it->eta()) < JetsMHTEta_) {
         if (i < JetDeltaMin_.size()) {
            if (std::abs(deltaPhi(MHT, *it)) < JetDeltaMin_.at(i))
               result = false;
            ++i;
         } else {
            break;
         }
      }
   }
   return result;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
math::PtEtaPhiMLorentzVector Factorization_SmearedJets::calcMHT_gen(const std::vector<reco::GenJet>& Jets_smeared) {
   math::PtEtaPhiMLorentzVector MHT(0, 0, 0, 0);
   for (vector<reco::GenJet>::const_iterator it = Jets_smeared.begin(); it != Jets_smeared.end(); ++it) {
      if (it->pt() > JetsMHTPt_ && std::abs(it->eta()) < JetsMHTEta_) {
         MHT -= it->p4();
      }
   }
   return MHT;
}
//--------------------------------------------------------------------------
// Fold histogram input with gaussian resolution
void Factorization_SmearedJets::StretchHisto(const TH1& input, TH1& output, const double& f) {

   if (input.Integral() > 0) {
      double mean = input.GetMean();
      for (int i = 0; i < 1000000; ++i) {
         double r = input.GetRandom();
         double rprime = mean + (r - mean) * f;
         output.Fill(rprime);
      }
      output.Scale(input.Integral() / output.Integral());
   }

}

//--------------------------------------------------------------------------
double Factorization_SmearedJets::GetLowerTailScaling(const double& pt, const double& eta) {
   int i_Pt = GetIndex(pt, &PtBinEdges_scaling_);
   int i_eta = GetIndex(eta, &EtaBinEdges_scaling_);
   double result = LowerTailScaling_.at(i_eta * (PtBinEdges_scaling_.size() - 1) + i_Pt) * LowerTailScaling_variation_;
   return result;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
double Factorization_SmearedJets::GetUpperTailScaling(const double& pt, const double& eta) {
   int i_Pt = GetIndex(pt, &PtBinEdges_scaling_);
   int i_eta = GetIndex(eta, &EtaBinEdges_scaling_);
   double result = UpperTailScaling_.at(i_eta * (PtBinEdges_scaling_.size() - 1) + i_Pt) * UpperTailScaling_variation_;
   return result;
}

//--------------------------------------------------------------------------
double Factorization_SmearedJets::GetAdditionalSmearing(const double& pt, const double& eta) {
   int i_Pt = GetIndex(pt, &PtBinEdges_scaling_);
   int i_eta = GetIndex(eta, &EtaBinEdges_scaling_);
   double result = AdditionalSmearing_.at(i_eta * (PtBinEdges_scaling_.size() - 1) + i_Pt)
         * AdditionalSmearing_variation_;
   //if (result < 1) result = 1;
   return result;
}

//--------------------------------------------------------------------------
// Fold histogram input with gaussian resolution
void Factorization_SmearedJets::FoldWithGaussian(const TH1& input, TH1& output, const double& sigma) {

   double min = input.GetXaxis()->GetXmin();
   double max = input.GetXaxis()->GetXmax();
   for (int i = 0; i < input.GetNbinsX(); ++i) {
      double weight = input.GetBinContent(i);
      double mean = input.GetBinCenter(i);
      TF1 gauss("gauss", "gaus(0)", min, max);
      gauss.SetParameters(weight * 1 / sigma / sqrt(2 * TMath::Pi()), mean, sigma);
      for (int j = 0; j < output.GetNbinsX(); ++j) {
         double xmin = output.GetBinLowEdge(j);
         double xmax = output.GetBinLowEdge(j) + output.GetBinWidth(j);
         output.AddBinContent(j, gauss.Integral(xmin, xmax));
      }
   }

}
//define this as a plug-in
DEFINE_FWK_MODULE(Factorization_SmearedJets);

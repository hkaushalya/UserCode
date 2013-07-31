// -*- C++ -*-
//
// Package:    EmptySelectors
// Class:      EmptySelectors
#include <memory>
#include <string>
#include <vector>
#include <cassert>
#include <cmath>
#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "recipeAUX/OxbridgeMT2/interface/ChengHanBisect_Mt2_332_Calculator.h"
#include "PhysicsTools/KinFitter/interface/TKinFitter.h"
#include "PhysicsTools/KinFitter/interface/TFitParticleEtEtaPhi.h"
#include "PhysicsTools/KinFitter/interface/TFitConstraintEp.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "CommonTools/Utils/interface/PtComparator.h"

#include "TRandom3.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TKey.h"
#include "TFile.h"
#include "TMath.h"
#include "TArray.h"

//#include "Utilities/Parang/interface/Paramatrix.h"
//#include "RA2Classic/QCDBkgRS/interface/SmearFunction.h"
#include "../interface/SmearFunctionChris.h"
#include <TROOT.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPostScript.h>


//includes needed for this file

using namespace edm;
using namespace reco;
using namespace std;

class RebalanceAndSmear : public edm::EDProducer {
   public:
      explicit RebalanceAndSmear(const edm::ParameterSet&);
      ~RebalanceAndSmear();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
   virtual void beginJob() ;
   virtual void produce(edm::Event&, const edm::EventSetup&);
   virtual void endJob() ;
   
   virtual void beginRun(edm::Run&, edm::EventSetup const&);
   virtual void endRun(edm::Run&, edm::EventSetup const&);
   virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
   virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
   double JetResolution_Pt2(const double&, const double&, const int&);
   double JetResolution_Ptrel(const double&, const double&, const int&);
   double JetResolution_Eta2(const double&, const double&);
   double JetResolution_Phi2(const double&, const double&);
   double JetResolutionHist_Pt_Smear(const double&, const double&, const int&,const double&, const int&);
   double GetHFProb(const int&, const double&, const int&);
   int GetIndex(const double&, const std::vector<double>*);
   void FillPredictions(const std::vector<pat::Jet>&, const int&, const double&);
   void FillPredictions_gen(const std::vector<reco::GenJet>&, const int&, const double&);
   double calcHT(const std::vector<pat::Jet>&);
   double calcHT_gen(const std::vector<reco::GenJet>&);
   math::PtEtaPhiMLorentzVector calcMHT(const std::vector<pat::Jet>&);
   math::PtEtaPhiMLorentzVector calcMHT_gen(const std::vector<reco::GenJet>&);
   int calcNJets(const std::vector<pat::Jet>&);
   int calcNJets_gen(const std::vector<reco::GenJet>&);
   bool calcMinDeltaPhi(const std::vector<pat::Jet>&, math::PtEtaPhiMLorentzVector&);
   bool calcMinDeltaPhi_gen(const std::vector<reco::GenJet>&, math::PtEtaPhiMLorentzVector&);
   void FillLeadingJetPredictions(const std::vector<pat::Jet>&); 
   void FillDeltaPhiPredictions(const std::vector<pat::Jet>&, math::PtEtaPhiMLorentzVector&); 
   void FillLeadingJetPredictions_gen(const std::vector<reco::GenJet>&); 
   void FillDeltaPhiPredictions_gen(const std::vector<reco::GenJet>&, math::PtEtaPhiMLorentzVector&); 
   double GetRebalanceCorrection(double jet_pt);

   
   bool RebalanceJets_KinFitter(edm::View<pat::Jet>*, std::vector<pat::Jet> &);
//   void SmearingJets(const std::vector<pat::Jet> &Jets_reb, std::vector<pat::Jet> &Jets_smeared, std::vector<std::vector<pat::Jet> >& Jets_keepAndSmeared, math::XYZTLorentzVector & unclusteredMET_reco, vector<reco::MET> & smearedMET, double & smearedWeight);


   void SmearingJets(const std::vector<pat::Jet> &Jets_reb, 
                                     math::XYZTLorentzVector & unclusteredMET_reco, 
                                     vector<vector<math::XYZTLorentzVector> > & Jets_smeared, 
                                     vector<math::XYZTLorentzVector> & METs_smeared,
                                     vector<math::XYZTLorentzVector> & MHTs_smeared,
                                     vector<vector<math::XYZTLorentzVector> > & triplet_smeared,
                                     vector<vector<math::XYZTLorentzVector> > & rSystem_smeared,
                                     vector<vector<math::XYZTLorentzVector> > & bJetsInR_smeared,
                                     vector<double> & dPhiJet1AndMET_smeared,
                                     vector<double> & dPhiJet2AndMET_smeared,
                                     vector<double> & dPhiJet3AndMET_smeared,
                                     vector<double> & m23OverM123_smeared,
                                     vector<double> & m123_smeared,
                                     vector<double> & MT2_smeared,
                                     vector<double> & MTt_smeared,
                                     vector<double> & MTb_smeared) ;

   void SmearingGenJets(edm::View<reco::GenJet>*, std::vector<reco::GenJet> &);
   
   std::string GetName(const std::string plot, const std::string uncert = "", const std::string ptbin = "") const;
   void TripletSelector(std::vector<pat::Jet> & , std::vector<pat::Jet> & , std::vector<pat::Jet> & , std::vector<pat::Jet> & , double & , double & );
   void MTMT2(reco::MET & MET, std::vector<pat::Jet> & triplet,std::vector<pat::Jet> & rSystem, std::vector<pat::Jet> & bJetsInR, double & MT2, double & MTt, double& MTb);

   typedef math::XYZTLorentzVector LorentzVector;
   typedef std::vector<std::string>::const_iterator StrIter;


   SmearFunctionChris *smearFunc_;        // Object of class SmearFunctionChris
   
   double rebalancedJetPt_;
   std::string rebalanceMode_; // "MHTall", "MHThigh" or "MET" only for smearCollection = "Reco"
   
   int nSmearedJets_;
   double smearedJetPt_;
   
   edm::InputTag vertices_;
   edm::InputTag genjets_;
   edm::InputTag genParticleSrc_;
   edm::InputTag pfParticleSrc_;
   edm::InputTag jets_;
   edm::InputTag METs_;
   edm::InputTag weightName_;
   std::string jets_reb_;
   std::string met_reb_;
   std::string jets_smeared_;
   
   //// vector of response function
   std::vector<std::vector<std::vector<TH1F*> > > smearFunc;
   std::vector<std::vector<std::vector<TH1F*> > > smearFunc_Core;
   std::vector<std::vector<std::vector<TH1F*> > > smearFunc_LowerTail;
   std::vector<std::vector<std::vector<TH1F*> > > smearFunc_UpperTail;
   std::vector<std::vector<std::vector<TH1F*> > > smearFunc_scaled;
   std::vector<std::vector<std::vector<TH1F*> > > smearFunc_HF;
   std::vector<std::vector<std::vector<TH1F*> > > smearFunc_Core_HF;
   std::vector<std::vector<std::vector<TH1F*> > > smearFunc_LowerTail_HF;
   std::vector<std::vector<std::vector<TH1F*> > > smearFunc_UpperTail_HF;
   std::vector<std::vector<std::vector<TH1F*> > > smearFunc_scaled_HF;
   std::vector<std::vector<std::vector<TH1F*> > > smearFunc_NoHF;
   std::vector<std::vector<std::vector<TH1F*> > > smearFunc_Core_NoHF;
   std::vector<std::vector<std::vector<TH1F*> > > smearFunc_LowerTail_NoHF;
   std::vector<std::vector<std::vector<TH1F*> > > smearFunc_UpperTail_NoHF;
   std::vector<std::vector<std::vector<TH1F*> > > smearFunc_scaled_NoHF;
   std::vector<std::vector<TH1F*> > SigmaPtHist;
   std::vector<std::vector<TF1*> > SigmaPt;
   std::vector<std::vector<TH1F*> > SigmaPtHist_scaled;
   std::vector<std::vector<TF1*> > SigmaPt_scaled;
   std::string RebalanceCorrectionFile_;
   bool isData_;
   bool cleverPrescaleTreating_;
   bool useRebalanceCorrectionFactors_;
   double MHTmin_;
   double MHTmax_;
   double HTmin_;
   double HTmax_;
   int Ntries_;
   int NJets_;
   double JetsHTPt_;
   double JetsHTEta_;
   double JetsMHTPt_;
   double JetsMHTEta_;
   vector<double> JetDeltaMin_;
   double MET_requirement_;
   unsigned randomNumberSeed_;

   int plotindex_;
   TRandom3 *rand_;
   double weight_;
   double nPassMHT_;

   TH1F * h_RebCorrectionFactor;

   UShort_t vtxN;
   UShort_t Njets_pred;
   UShort_t Ntries_pred;
   Float_t HT_seed;
   Float_t HT_pred;
   Float_t MHT_pred;
   Float_t weight;
   Float_t Jet1Pt_pred;
   Float_t Jet2Pt_pred;
   Float_t Jet3Pt_pred;
   Float_t Jet1Eta_pred;
   Float_t Jet2Eta_pred;
   Float_t Jet3Eta_pred;
   Float_t DeltaPhi1_pred;
   Float_t DeltaPhi2_pred;
   Float_t DeltaPhi3_pred;



   std::vector<double> PtBinEdges_scaling_;
   std::vector<double> EtaBinEdges_scaling_;
   std::vector<double> AdditionalSmearing_;
   std::vector<double> LowerTailScaling_;
   std::vector<double> UpperTailScaling_;
   double AdditionalSmearing_variation_;
   double LowerTailScaling_variation_;
   double UpperTailScaling_variation_;

   std::string inputhist1HF_;
   std::string inputhist2HF_;
   std::string inputhist3pHF_;
   std::string inputhist1NoHF_;
   std::string inputhist2NoHF_;
   std::string inputhist3pNoHF_;
   std::string smearingfile_;
   std::string bprobabilityfile_;
   std::string outputfile_;

   int NRebin_;
   bool absoluteTailScaling_;
   double A0RMS_;
   double A1RMS_;
   double probExtreme_;

   std::string uncertaintyName_;

   std::vector<double> PtBinEdges_;
   std::vector<double> EtaBinEdges_;

   reco::Candidate::Point METVertex_;


   double Rmin_;
   double Rmax_;
   double arctanmin_;
   double arctanmax_;
   double m23OverM123Cut_;
   double topMass_;
   double bJetPtCut_;
   double bJetEtaCut_;
   string bJetDisc_;
   double bJetDiscCut_;

   double mTop_;
   double mWMin_;
   double mWMax_;

   double tripletJetPtCut_;
   double tripletDRCut_;

};

RebalanceAndSmear::RebalanceAndSmear(const edm::ParameterSet& iConfig)
{
   rebalancedJetPt_ = iConfig.getParameter<double> ("RebalanceJetPt");
   rebalanceMode_ = iConfig.getParameter<std::string> ("RebalanceMode");
   smearedJetPt_ = iConfig.getParameter<double> ("SmearedJetPt");
   vertices_ = iConfig.getParameter<edm::InputTag>("VertexCollection");
   jets_ = iConfig.getParameter<edm::InputTag> ("jetCollection");
   genParticleSrc_ = iConfig.getParameter<edm::InputTag>("genParticleSrc");
   pfParticleSrc_ = iConfig.getParameter<edm::InputTag>("pfParticleSrc");
   METs_ = iConfig.getParameter<edm::InputTag> ("METCollection");
   jets_reb_ = iConfig.getParameter<std::string> ("jetCollection_reb");
   jets_smeared_ = iConfig.getParameter<std::string> ("jetCollection_smeared");
   RebalanceCorrectionFile_ = iConfig.getParameter<std::string> ("RebalanceCorrectionFile");
   isData_ = iConfig.getParameter<bool> ("IsData");
   weightName_ = iConfig.getParameter<edm::InputTag> ("weightName");
   cleverPrescaleTreating_ = iConfig.getParameter<bool> ("cleverPrescaleTreating");
   useRebalanceCorrectionFactors_ = iConfig.getParameter<bool> ("useRebalanceCorrectionFactors");
   Ntries_ = iConfig.getParameter<int> ("Ntries");
   NJets_ = iConfig.getParameter<int> ("NJets");
   JetsHTPt_ = iConfig.getParameter<double> ("JetsHTPt");
   JetsHTEta_ = iConfig.getParameter<double> ("JetsHTEta");
   JetsMHTPt_ = iConfig.getParameter<double> ("JetsMHTPt");
   JetsMHTEta_ = iConfig.getParameter<double> ("JetsMHTEta");
   JetDeltaMin_ = iConfig.getParameter<std::vector<double> > ("JetDeltaMin");
   PtBinEdges_ = iConfig.getParameter<std::vector<double> > ("PtBinEdges");
   EtaBinEdges_ = iConfig.getParameter<std::vector<double> > ("EtaBinEdges");

   MET_requirement_ = iConfig.getParameter<double> ("MET_requirement");
   randomNumberSeed_ = iConfig.getParameter<unsigned>("randomNumberSeed");



   PtBinEdges_scaling_ = iConfig.getParameter<std::vector<double> > ("PtBinEdges_scaling");
   EtaBinEdges_scaling_ = iConfig.getParameter<std::vector<double> > ("EtaBinEdges_scaling");
   LowerTailScaling_ = iConfig.getParameter<std::vector<double> > ("LowerTailScaling");
   UpperTailScaling_ = iConfig.getParameter<std::vector<double> > ("UpperTailScaling");
   AdditionalSmearing_ = iConfig.getParameter<std::vector<double> > ("AdditionalSmearing");
   LowerTailScaling_variation_ = iConfig.getParameter<double> ("LowerTailScaling_variation");
   UpperTailScaling_variation_ = iConfig.getParameter<double> ("UpperTailScaling_variation");
   AdditionalSmearing_variation_ = iConfig.getParameter<double> ("AdditionalSmearing_variation");
   NRebin_ = iConfig.getParameter<int> ("NRebin");
   uncertaintyName_ = iConfig.getParameter<std::string> ("uncertaintyName");
   inputhist1HF_ = iConfig.getParameter<std::string> ("InputHisto1_HF");
   inputhist2HF_ = iConfig.getParameter<std::string> ("InputHisto2_HF");
   inputhist3pHF_ = iConfig.getParameter<std::string> ("InputHisto3p_HF");
   inputhist1NoHF_ = iConfig.getParameter<std::string> ("InputHisto1_NoHF");
   inputhist2NoHF_ = iConfig.getParameter<std::string> ("InputHisto2_NoHF");
   inputhist3pNoHF_ = iConfig.getParameter<std::string> ("InputHisto3p_NoHF");
   smearingfile_ = iConfig.getParameter<std::string> ("SmearingFile");
   bprobabilityfile_ = iConfig.getParameter<std::string> ("BProbabilityFile");
   outputfile_ = iConfig.getParameter<std::string> ("OutputFile");
   absoluteTailScaling_ = iConfig.getParameter<bool> ("absoluteTailScaling");
   A0RMS_ = iConfig.getParameter<double> ("A0RMS");
   A1RMS_ = iConfig.getParameter<double> ("A1RMS");
   probExtreme_ = iConfig.getParameter<double> ("probExtreme");  
   PtBinEdges_ = iConfig.getParameter<std::vector<double> > ("PtBinEdges");
   EtaBinEdges_ = iConfig.getParameter<std::vector<double> > ("EtaBinEdges");


   Rmin_ = iConfig.getParameter<double> ("Rmin");
   Rmax_ = iConfig.getParameter<double> ("Rmax");
   arctanmin_ = iConfig.getParameter<double> ("arctanmin");
   arctanmax_ = iConfig.getParameter<double> ("arctanmax");
   m23OverM123Cut_ = iConfig.getParameter<double>("m23OverM123Cut");
   topMass_ = iConfig.getParameter<double> ("topMass");

   mTop_ = iConfig.getParameter<double>("mTop");
   mWMin_ = iConfig.getParameter<double>("mWMin");
   mWMax_ = iConfig.getParameter<double>("mWMax");

   tripletJetPtCut_ = iConfig.getParameter<double>("tripletJetPtCut");
   tripletDRCut_ = iConfig.getParameter<double>("tripletDRCut");
   bJetPtCut_ = iConfig.getParameter<double>("bJetPtCut");
   bJetEtaCut_ = iConfig.getParameter<double>("bJetEtaCut");
   bJetDisc_ = iConfig.getParameter<string>("bJetDisc");
   bJetDiscCut_ = iConfig.getParameter<double>("bJetDiscCut");

   gRandom->SetSeed(randomNumberSeed_);
   rand_ = new TRandom3(randomNumberSeed_);

   smearFunc_ = new SmearFunctionChris(iConfig); 

   produces<std::vector<LorentzVector> >("rebalancedJets");
   produces<LorentzVector >("rebalancedMET");
   produces<LorentzVector>("unclusteredMET");
   produces<std::vector<std::vector<LorentzVector> > >("smearedJets");
   produces<std::vector<LorentzVector> >("smearedMET");
   produces<std::vector<LorentzVector> >("smearedMHT");
   produces<std::vector<std::vector<LorentzVector> > >("smearedTriplet");
   produces<std::vector<std::vector<LorentzVector> > >("smearedRSystem");
   produces<std::vector<std::vector<LorentzVector> > >("smearedBJet");
   produces<std::vector<double> >("dPhiJet1AndMET");
   produces<std::vector<double> >("dPhiJet2AndMET");
   produces<std::vector<double> >("dPhiJet3AndMET");
   produces<std::vector<double> >("m23OverM123");
   produces<std::vector<double> >("m123");
   produces<std::vector<double> >("MTt");
   produces<std::vector<double> >("MTb");
   produces<std::vector<double> >("MT2");
  



}


RebalanceAndSmear::~RebalanceAndSmear()
{
 
   PtBinEdges_.clear();
   EtaBinEdges_.clear();

   delete smearFunc_;

   if (rand_)
      delete rand_;


   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

std::string RebalanceAndSmear::GetName(const std::string plot, const std::string uncert, const std::string ptbin) const {
   std::string result(plot);
   if (uncert != "" && uncert != " ")
      result += "_" + uncert;
   if (ptbin != "" && ptbin != " ")
      result += "_" + ptbin;
   return result;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
int RebalanceAndSmear::GetIndex(const double& x, const std::vector<double>* vec) {
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

//--------------------------------------------------------------------------
// pt resolution for KinFitter
double RebalanceAndSmear::JetResolution_Pt2(const double& pt, const double& eta, const int& i) {
   int i_jet;
   i < 2 ? i_jet = i : i_jet = 2;
   int i_eta = GetIndex(eta, &EtaBinEdges_);
   //return pow(pt, 2) * pow(smearFunc_->getSigmaPtScaledForRebalancing(i_jet, i_eta)->Eval(pt), 2);
   return pow(pt, 2) * pow(smearFunc_->getSigmaPtForRebalancing(i_jet, i_eta)->Eval(pt), 2);
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// relative pt resolution for KinFitter
double RebalanceAndSmear::JetResolution_Ptrel(const double& pt, const double& eta, const int& i) {
   int i_jet;
   i < 2 ? i_jet = i : i_jet = 2;
   int i_eta = GetIndex(eta, &EtaBinEdges_);
   return smearFunc_->getSigmaPtScaledForRebalancing(i_jet, i_eta)->Eval(pt);
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// eta resolution for KinFitter
double RebalanceAndSmear::JetResolution_Eta2(const double& e, const double& eta) {
   //may be artifically reduced (no angular fit)
   return (pow(0.05 / TMath::Sqrt(e), 2) + pow(0.005, 2)) / 1.e6;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// phi resolution for KinFitter
double RebalanceAndSmear::JetResolution_Phi2(const double& e, const double& eta) {
   //may be artifically reduced (no angular fit)
   return (pow(0.05 / TMath::Sqrt(e), 2) + pow(0.005, 2)) / 1.e6;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// pt resolution for smearing
double RebalanceAndSmear::JetResolutionHist_Pt_Smear(const double& pt, const double& eta, const int& i, const double& HT, const int& NJets) {
   int i_jet;
   i < 2 ? i_jet = i : i_jet = 2;
   int i_Pt = GetIndex(pt, &PtBinEdges_);
   int i_eta = GetIndex(eta, &EtaBinEdges_);

   int jet_rank = i+1;
   double x_rand = rand_->Rndm();
   double hf_prob = GetHFProb(NJets, HT, jet_rank);

   double res = 1.0;
   if( x_rand <= hf_prob ){
      // get heavy flavor smear function
      res = smearFunc_->getSmearFunc(1, i_jet, i_eta, i_Pt)->GetRandom();
   }
   else {
      // get no heavy flavor smear function
      res = smearFunc_->getSmearFunc(0, i_jet, i_eta, i_Pt)->GetRandom();
   }

   return res;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// HF probability for smearing
double RebalanceAndSmear::GetHFProb(const int& NJets, const double& HT, const int& jet_rank) {
  
   int i_bin;

   //Added by chris

/*   if( NJets == 2){
      i_bin = h_bProb_NJets2->FindBin(HT, jet_rank);
      return h_bProb_NJets2->GetBinContent(i_bin);
   }
   if( NJets == 3){
      i_bin = h_bProb_NJets3->FindBin(HT, jet_rank);
      return h_bProb_NJets3->GetBinContent(i_bin);
   }
   if( NJets == 4){
      i_bin = h_bProb_NJets4->FindBin(HT, jet_rank);
      return h_bProb_NJets4->GetBinContent(i_bin);
   }
   if( NJets == 5 || NJets == 6){
      i_bin = h_bProb_NJets5p6->FindBin(HT, jet_rank);
      return h_bProb_NJets5p6->GetBinContent(i_bin);
   }
   if( NJets >= 7){
      i_bin = h_bProb_NJets7p->FindBin(HT, jet_rank);
      return h_bProb_NJets7p->GetBinContent(i_bin);
   }
   else {
      return -1;
   }
*/
   return -1;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
double RebalanceAndSmear::GetRebalanceCorrection(double jet_pt)
{
   if( jet_pt > 1000. ) jet_pt = 999.;

   

   if ( jet_pt > rebalancedJetPt_ ) {
      int i_bin = h_RebCorrectionFactor->FindBin(jet_pt);

      //  cout << "Reco jet pt: " << jet_pt << endl;
      //cout << "Rebalance correction: " << h_RebCorrectionFactor->GetBinContent(i_bin) << endl;
  
      return h_RebCorrectionFactor->GetBinContent(i_bin);
   }

   else return 1;

   return 1;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// rebalance the events using a kinematic fit and transverse momentum balance
bool RebalanceAndSmear::RebalanceJets_KinFitter(edm::View<pat::Jet>* Jets_rec, std::vector<pat::Jet> &Jets_reb) {

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
        
         JetMap[i] = &(*it); 

         // The particles before fitting
         double tmppx, tmppy, tmppz, tmpe;

         if( useRebalanceCorrectionFactors_ ) {
            tmppx = it->px()/GetRebalanceCorrection( it->pt() );
            tmppy = it->py()/GetRebalanceCorrection( it->pt() );
            tmppz = it->pz()/GetRebalanceCorrection( it->pt() );
            tmpe = it->energy()/GetRebalanceCorrection( it->pt() );
         }
         else {
            tmppx = it->px();
            tmppy = it->py();
            tmppz = it->pz();
            tmpe = it->energy();
         }

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
   try
   {
      myFit->fit();
   }
   catch (...)
   {
      cout<<"Rebalance failed, could not fit KinFitter object"<<endl;
      return false;
   }
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
      //    if (prob < 0.01) result = false;
      //cout << "chi2, prop, F = " << chi2 << " " << prob << " " << F << endl;
   } else {
      chi2 = 99999;
      prob = 0;
      F = 99999;
      result = false;
      //cout << "chi2, prop, F = " << chi2 << " " << prob << " " << F << endl;
   }

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

void RebalanceAndSmear::SmearingJets(const std::vector<pat::Jet> &Jets_reb, 
                                     LorentzVector & unclusteredMET_reco, 
                                     vector<vector<LorentzVector> > & Jets_smeared, 
                                     vector<LorentzVector> & METs_smeared,
                                     vector<LorentzVector> & MHTs_smeared,
                                     vector<vector<LorentzVector> > & triplet_smeared,
                                     vector<vector<LorentzVector> > & rSystem_smeared,
                                     vector<vector<LorentzVector> > & bJetsInR_smeared,
                                     vector<double> & dPhiJet1AndMET_smeared,
                                     vector<double> & dPhiJet2AndMET_smeared,
                                     vector<double> & dPhiJet3AndMET_smeared,
                                     vector<double> & m23OverM123_smeared,
                                     vector<double> & m123_smeared,
                                     vector<double> & MT2_smeared,
                                     vector<double> & MTt_smeared,
                                     vector<double> & MTb_smeared) 
{

   double dPx = 0;
   double dPy = 0;

   double HT = calcHT(Jets_reb);
   int NJets_reb = calcNJets(Jets_reb);

   weight_ = 1.0;
   std::vector<pat::Jet> Jets_thisTry;
   std::vector<pat::Jet> Jets_thisTry_PtCut;
   std::vector<pat::Jet> triplet;
   std::vector<pat::Jet> rSystem;
   std::vector<pat::Jet> bJet;

   std::vector<LorentzVector> Jets_thisTry_P4s;   
   std::vector<LorentzVector> triplet_P4s;   
   std::vector<LorentzVector> rSystem_P4s;   
   std::vector<LorentzVector> bJet_P4s;   

   double dPhiJet1AndMET;
   double dPhiJet2AndMET;
   double dPhiJet3AndMET;
   double M23OverM123;
   double M123;
   double MT2;
   double MTt;
   double MTb;
   double MTt2MTb;
   


   for (int i = 1; i <= Ntries_; ++i) 
   {
      int Ntries2 = 1;
      double w = weight_;
      if (cleverPrescaleTreating_ == true && weight_ > 1) 
      {
         Ntries2 = (int) weight_;
         w = weight_ / Ntries2;
      }
      for (int j = 1; j <= Ntries2; ++j) 
      {
         Jets_thisTry.clear();
         Jets_thisTry_P4s.clear();
         Jets_thisTry_PtCut.clear();
         triplet.clear();
         rSystem.clear();
         bJet.clear();
         dPhiJet1AndMET = -1;
         dPhiJet2AndMET = -1;
         dPhiJet3AndMET = -1;
         int i_jet = 0;
         LorentzVector vMET_pred(0.,0.,0.,0.);
         LorentzVector vMHT_pred(0.,0.,0.,0.);
         

         for (std::vector<pat::Jet>::const_iterator it = Jets_reb.begin(); it != Jets_reb.end(); ++it) 
         {
            if (it->pt() > smearedJetPt_) 
            {
               double newPt = 0;
               newPt = it->pt() * JetResolutionHist_Pt_Smear(it->pt(), it->eta(), i_jet, HT, NJets_reb);
               //double newEta = rand_->Gaus(it->eta(), TMath::Sqrt(JetResolution_Eta2(it->energy(), it->eta())));
               //double newPhi = rand_->Gaus(it->phi(), TMath::Sqrt(JetResolution_Phi2(it->energy(), it->eta())));
               double newEta = it->eta();
               double newPhi = it->phi();
               pat::Jet::PolarLorentzVector newP4(newPt, newEta, newPhi, it->mass());
               pat::Jet smearedJet(*it);
               smearedJet.setP4(newP4);
               Jets_thisTry.push_back(smearedJet);
               dPx -= newP4.Px() - it->px();
               dPy -= newP4.Py() - it->py();
               ++i_jet;
               LorentzVector tempP4(newP4.px(), newP4.py(), 0 , sqrt(newP4.px()*newP4.px() + newP4.py()*newP4.py()));
               vMET_pred -= tempP4;
            } 
            else 
            {
               pat::Jet smearedJet(*it);
               Jets_thisTry.push_back(smearedJet);
            }

            

         }
         GreaterByPt<reco::Candidate> ptComparator_;
         std::sort(Jets_thisTry.begin(), Jets_thisTry.end(), ptComparator_);

         //Fill HT and MHT prediction histos for i-th iteration of smearing
         int NJets = calcNJets(Jets_thisTry);
            double px = vMET_pred.px();
            double py = vMET_pred.py();
            vMHT_pred.SetPxPyPzE(px, py, 0, sqrt(px*px + py*py));
            px = vMET_pred.px() + unclusteredMET_reco.px();
            py = vMET_pred.py() + unclusteredMET_reco.py();
            vMET_pred.SetPxPyPzE(px, py, 0, sqrt(px*px + py*py));
            reco::MET vMET(vMET_pred, METVertex_);
            
            
            
            for(unsigned jet_i = 0; jet_i < Jets_thisTry.size(); jet_i++)
            {
               Jets_thisTry_P4s.push_back(Jets_thisTry[jet_i].p4());
               if(jet_i == 0)
                  dPhiJet1AndMET = abs(deltaPhi(Jets_thisTry[jet_i].phi(), vMET.phi()));
               if(jet_i == 1)
                  dPhiJet2AndMET = abs(deltaPhi(Jets_thisTry[jet_i].phi(), vMET.phi()));
               if(jet_i == 2)
                  dPhiJet3AndMET = abs(deltaPhi(Jets_thisTry[jet_i].phi(), vMET.phi()));
               if(Jets_thisTry[jet_i].pt() > tripletJetPtCut_)
               {
                  Jets_thisTry_PtCut.push_back(Jets_thisTry[jet_i]);
               }
            }
            
            
            
            TripletSelector(Jets_thisTry_PtCut, triplet, rSystem, bJet, M23OverM123, M123);
            MTMT2(vMET,triplet,rSystem,bJet,MT2,MTt,MTb);
            
            
/*         cout<<"DPHI: "<<dPhiJet1AndMET<<" "<<dPhiJet2AndMET<<" "<<dPhiJet3AndMET<<endl;
           cout<<"TRIPLET SIZE: "<<triplet.size()<<endl;
           cout<<"RSYSTEM SIZE: "<<rSystem.size()<<endl;
           cout<<"BJET SIZE: "<<bJet.size()<<endl;
           cout<<"M23/M123: "<<M23OverM123<<endl;
           cout<<"M123: "<<M123<<endl;
           cout<<"MT2: "<<MT2<<endl;
           cout<<"MTt: "<<MTt<<endl;
           cout<<"MTb: "<<MTb<<endl;
*/
            for(unsigned jet_i = 0; jet_i < triplet.size(); jet_i++)
            {
               triplet_P4s.push_back(triplet[jet_i].p4());
            }
            
            for(unsigned jet_i = 0; jet_i < rSystem.size(); jet_i++)
            {
               rSystem_P4s.push_back(rSystem[jet_i].p4());
            }
            for(unsigned jet_i = 0; jet_i < bJet.size(); jet_i++)
            {
               bJet_P4s.push_back(bJet[jet_i].p4());
            }
            
            Jets_smeared.push_back(Jets_thisTry_P4s);
            METs_smeared.push_back(vMET_pred);
            MHTs_smeared.push_back(vMHT_pred);
            dPhiJet1AndMET_smeared.push_back(dPhiJet1AndMET);
            dPhiJet2AndMET_smeared.push_back(dPhiJet2AndMET);
            dPhiJet3AndMET_smeared.push_back(dPhiJet3AndMET);
            triplet_smeared.push_back(triplet_P4s);
            rSystem_smeared.push_back(rSystem_P4s);
            bJetsInR_smeared.push_back(bJet_P4s);
            m23OverM123_smeared.push_back(M23OverM123);
            m123_smeared.push_back(M123);
            MT2_smeared.push_back(MT2);
            MTt_smeared.push_back(MTt);
            MTb_smeared.push_back(MTb);
            nPassMHT_++;
            

 
         // save only events with MHT > 100. GeV for data
         
         // clean variables in tree
         weight = 0.;
         Ntries_pred = 0.;
         Njets_pred = 0;
         HT_pred = 0.;
         MHT_pred = 0.;
         Jet1Pt_pred = 0.;
         Jet2Pt_pred = 0.;
         Jet3Pt_pred = 0.;
         Jet1Eta_pred = 0.;
         Jet2Eta_pred = 0.;
         Jet3Eta_pred = 0.;
         DeltaPhi1_pred = 0.;
         DeltaPhi2_pred = 0.;
         DeltaPhi3_pred = 0.;
        
      }
   }
   return;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
void RebalanceAndSmear::SmearingGenJets(edm::View<reco::GenJet>* Jets_gen, std::vector<reco::GenJet> &GenJets_smeared) {

   double dPx = 0;
   double dPy = 0;

   // calculate quantities needed for smearing
   double HT; 
   int NJets_gen = 0;
   for (edm::View<reco::GenJet>::const_iterator it = Jets_gen->begin(); it != Jets_gen->end(); ++it) {
      if (it->pt() > JetsHTPt_ && std::abs(it->eta()) < JetsHTEta_) {
         ++NJets_gen;
         HT += it->pt();
      }
   }

   for (int i = 1; i <= Ntries_; ++i) {
      GenJets_smeared.clear();
      int i_jet = 0;

      for (edm::View<reco::GenJet>::const_iterator it = Jets_gen->begin(); it != Jets_gen->end(); ++it) {

         if (it->pt() > smearedJetPt_) {
            double newPt = 0;
            newPt = it->pt() * JetResolutionHist_Pt_Smear(it->pt(), it->eta(), i_jet, HT, NJets_gen);
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
      int NJets = calcNJets_gen(GenJets_smeared);
      if (NJets >= NJets_) {
         FillPredictions_gen(GenJets_smeared, i, weight_);

         // clean variables in tree
         weight = 0.;
         Ntries_pred = 0.;
         Njets_pred = 0;
         HT_pred = 0.;
         MHT_pred = 0.;
         Jet1Pt_pred = 0.;
         Jet2Pt_pred = 0.;
         Jet3Pt_pred = 0.;
         Jet1Eta_pred = 0.;
         Jet2Eta_pred = 0.;
         Jet3Eta_pred = 0.;
         DeltaPhi1_pred = 0.;
         DeltaPhi2_pred = 0.;
         DeltaPhi3_pred = 0.;
      }
   }

   return;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
int RebalanceAndSmear::calcNJets(const std::vector<pat::Jet>& Jets_smeared) {
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
double RebalanceAndSmear::calcHT(const std::vector<pat::Jet>& Jets_smeared) {
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
bool RebalanceAndSmear::calcMinDeltaPhi(const std::vector<pat::Jet>& Jets_smeared, math::PtEtaPhiMLorentzVector& MHT) {
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
math::PtEtaPhiMLorentzVector RebalanceAndSmear::calcMHT(const std::vector<pat::Jet>& Jets_smeared) {
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
void RebalanceAndSmear::FillLeadingJetPredictions(const std::vector<pat::Jet>& Jets_smeared) {
   int NJets = 0;
   for (vector<pat::Jet>::const_iterator it = Jets_smeared.begin(); it != Jets_smeared.end(); ++it) {
      if (it->pt() > JetsHTPt_ && std::abs(it->eta()) < JetsHTEta_) {
         ++NJets;
      
         if( NJets == 1 ) {
            Jet1Pt_pred = it->pt();
            Jet1Eta_pred = it->eta(); 
         }
         if( NJets == 2 ) {            
            Jet2Pt_pred = it->pt();
            Jet2Eta_pred = it->eta();  
         }
         if( NJets == 3 ) {
            Jet3Pt_pred = it->pt();
            Jet3Eta_pred = it->eta();
            break;    
         }
      }
   }
   if( NJets == 2 ) {
      Jet3Pt_pred = -1.;
      Jet3Eta_pred = 9999.;
   }

   return;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
void RebalanceAndSmear::FillDeltaPhiPredictions(const std::vector<pat::Jet>& Jets_smeared, math::PtEtaPhiMLorentzVector& vMHT) {
   int NJets = 0;
   for (vector<pat::Jet>::const_iterator it = Jets_smeared.begin(); it != Jets_smeared.end(); ++it) {
      if (it->pt() > JetsMHTPt_ && std::abs(it->eta()) < JetsMHTEta_) {
         ++NJets;
      
         if( NJets == 1 ) {
            DeltaPhi1_pred = std::abs(deltaPhi(vMHT, *it));
         }
         if( NJets == 2 ) {
            DeltaPhi2_pred = std::abs(deltaPhi(vMHT, *it));
         }
         if( NJets == 3 ) {
            DeltaPhi3_pred = std::abs(deltaPhi(vMHT, *it));
            break;
         }
      }
   }
   if( NJets == 2 ) {
      DeltaPhi3_pred = 9999.;
   }

   return;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
int RebalanceAndSmear::calcNJets_gen(const std::vector<reco::GenJet>& Jets_smeared) {
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
double RebalanceAndSmear::calcHT_gen(const std::vector<reco::GenJet>& Jets_smeared) {
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
bool RebalanceAndSmear::calcMinDeltaPhi_gen(const std::vector<reco::GenJet>& Jets_smeared, math::PtEtaPhiMLorentzVector& MHT) {
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
math::PtEtaPhiMLorentzVector RebalanceAndSmear::calcMHT_gen(const std::vector<reco::GenJet>& Jets_smeared) {
   math::PtEtaPhiMLorentzVector MHT(0, 0, 0, 0);
   for (vector<reco::GenJet>::const_iterator it = Jets_smeared.begin(); it != Jets_smeared.end(); ++it) {
      if (it->pt() > JetsMHTPt_ && std::abs(it->eta()) < JetsMHTEta_) {
         MHT -= it->p4();
      }
   }
   return MHT;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
void RebalanceAndSmear::FillLeadingJetPredictions_gen(const std::vector<reco::GenJet>& Jets_smeared) {
   int NJets = 0;
   for (vector<reco::GenJet>::const_iterator it = Jets_smeared.begin(); it != Jets_smeared.end(); ++it) {
      if (it->pt() > JetsHTPt_ && std::abs(it->eta()) < JetsHTEta_) {
         ++NJets;
      
         if( NJets == 1 ) {            
            Jet1Pt_pred = it->pt();
            Jet1Eta_pred = it->eta(); 
              
         }
         if( NJets == 2 ) {
            Jet2Pt_pred = it->pt();
            Jet2Eta_pred = it->eta(); 
         }
         if( NJets == 3 ) {
            Jet3Pt_pred = it->pt();
            Jet3Eta_pred = it->eta(); 
            break;  
         }
      }
   }
   if( NJets == 2 ) {
      Jet3Pt_pred = -1.;
      Jet3Eta_pred = 9999.; 
   }

   return;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
void RebalanceAndSmear::FillDeltaPhiPredictions_gen(const std::vector<reco::GenJet>& Jets_smeared, math::PtEtaPhiMLorentzVector& vMHT) {
   int NJets = 0;
   for (vector<reco::GenJet>::const_iterator it = Jets_smeared.begin(); it != Jets_smeared.end(); ++it) {
      if (it->pt() > JetsMHTPt_ && std::abs(it->eta()) < JetsMHTEta_) {
         ++NJets;
      
         if( NJets == 1 ) {
            DeltaPhi1_pred = std::abs(deltaPhi(vMHT, *it));
         }
         if( NJets == 2 ) {
            DeltaPhi2_pred = std::abs(deltaPhi(vMHT, *it));
         }
         if( NJets == 3 ) {
            DeltaPhi3_pred = std::abs(deltaPhi(vMHT, *it));
            break;
         }
      }
   }
   if( NJets == 2 ) {
      DeltaPhi3_pred = 9999.;
   }


   return;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
void RebalanceAndSmear::FillPredictions(const std::vector<pat::Jet>& Jets_smeared, const int& i, const double& w) {

   int NJets = calcNJets(Jets_smeared);
   double HT = calcHT(Jets_smeared);
   math::PtEtaPhiMLorentzVector vMHT = calcMHT(Jets_smeared);
   double MHT = vMHT.pt();
   //   bool minDeltaPhi = calcMinDeltaPhi(Jets_smeared, vMHT);
 
   weight = w;
   Ntries_pred = i;
   Njets_pred = NJets;
   HT_pred = HT;
   MHT_pred = MHT;
   FillDeltaPhiPredictions(Jets_smeared, vMHT);
   FillLeadingJetPredictions(Jets_smeared); 

   return;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
void RebalanceAndSmear::FillPredictions_gen(const std::vector<reco::GenJet>& Jets_smeared, const int& i,
      const double& w) {

   int NJets = calcNJets_gen(Jets_smeared);
   double HT = calcHT_gen(Jets_smeared);
   math::PtEtaPhiMLorentzVector vMHT = calcMHT_gen(Jets_smeared);
   double MHT = vMHT.pt();
   //  bool minDeltaPhi = calcMinDeltaPhi_gen(Jets_smeared, vMHT);

   weight = w;
   Ntries_pred = i;
   Njets_pred = NJets;
   HT_pred = HT;
   MHT_pred = MHT;
   FillDeltaPhiPredictions_gen(Jets_smeared, vMHT);
   FillLeadingJetPredictions_gen(Jets_smeared); 
    
   return;
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void
RebalanceAndSmear::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   //Weight
//   edm::Handle<double> event_weight;
//   iEvent.getByLabel(weightName_, event_weight);
//   weight_ = (event_weight.isValid() ? (*event_weight) : 1.0);
   //if (!event_weight.isValid()) cout << "weight not found" << endl;
   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByLabel(vertices_,vertices);
   if( vertices.isValid() ) 
   {
      vtxN = vertices->size();
   }  

   edm::Handle<edm::View<GenParticle> > genParticles;
   bool foundGen = iEvent.getByLabel(genParticleSrc_, genParticles);
   edm::Handle<edm::View<PFCandidate> > pfParticles;
   bool foundpf = iEvent.getByLabel(pfParticleSrc_, pfParticles);
   //PATJets
   edm::Handle<edm::View<pat::Jet> > Jets;
   iEvent.getByLabel(jets_, Jets);
   edm::View<pat::Jet> Jets_rec = *Jets;


   //PATMET
   edm::Handle<edm::View<pat::MET> > METs;
   iEvent.getByLabel(METs_, METs);

   reco::Candidate::Point METVertex_ = (*METs)[0].vertex();


   std::auto_ptr<vector<pat::Jet> > Jets_reb(new vector<pat::Jet> ); 
   std::auto_ptr<vector<reco::MET> > METs_reb(new vector<reco::MET> );
   std::auto_ptr<vector<LorentzVector> > Jets_reb_P4(new vector<LorentzVector> ); 
   std::auto_ptr<LorentzVector> MET_reb_P4(new LorentzVector );
   std::auto_ptr<LorentzVector> unclusteredMET_reb(new LorentzVector);
   std::auto_ptr<vector<vector<LorentzVector> > > Jets_smeared(new vector<vector<LorentzVector> >);
   std::auto_ptr<vector<LorentzVector> > METs_smeared(new vector<LorentzVector> );
   std::auto_ptr<vector<LorentzVector> > MHTs_smeared(new vector<LorentzVector> );
   std::auto_ptr<vector<vector<LorentzVector> > > triplet_smeared(new vector<vector<LorentzVector> > );
   std::auto_ptr<vector< vector<LorentzVector> > > rSystem_smeared(new vector<vector<LorentzVector> >);
   std::auto_ptr<vector< vector<LorentzVector> > >bJetsInR_smeared(new vector<vector<LorentzVector> > );
   std::auto_ptr<vector<double> > dPhiJet1AndMET(new vector<double>);
   std::auto_ptr<vector<double> > dPhiJet2AndMET(new vector<double>);
   std::auto_ptr<vector<double> > dPhiJet3AndMET(new vector<double>);
   std::auto_ptr<vector<double> > m23OverM123(new vector<double>);
   std::auto_ptr<vector<double> > m123(new vector<double>);
   std::auto_ptr<vector<double> > MTt(new vector<double>);
   std::auto_ptr<vector<double> > MTb(new vector<double>);
   std::auto_ptr<vector<double> > MT2(new vector<double>);



   LorentzVector unclusteredMET(0. ,0. ,0. ,0.);
   LorentzVector vMHT(0., 0., 0., 0.);
   
   nPassMHT_ = 0;

   for (edm::View<pat::Jet>::const_iterator it = Jets_rec.begin(); it != Jets_rec.end(); ++it)
   {
      if (it->pt() > JetsMHTPt_ && abs(it->eta()) < JetsMHTEta_) 
      {
         LorentzVector tempP4(it->px(), it->py(), 0, sqrt(it->px()*it->px() + it->py()*it->py()));
         vMHT -= tempP4;
      }
   }


   double px =(*METs)[0].px() - vMHT.px();
   double py =(*METs)[0].py() - vMHT.py();

   unclusteredMET.SetPxPyPzE(px, py, 0, sqrt(px*px+py*py));
   (*unclusteredMET_reb).SetPxPyPzE(px, py, 0, sqrt(px*px+py*py));

   if(unclusteredMET.pt() > 500)
   {
      if(foundGen)
      {
         cout<<"STATUS 3 PARTICLES: "<<endl;
         for(unsigned i = 0; i < genParticles->size(); i++)
         {
            if((*genParticles)[i].status() != 3) continue;
            cout<<"PDGID: "<<(*genParticles)[i].pdgId()<<" ETA: "<<(*genParticles)[i].eta()<<" PHI: "<<(*genParticles)[i].phi()<<" PT: "<<(*genParticles)[i].pt()<<endl;
         }
      }
      for(unsigned i = 0; i < pfParticles->size(); i++)
      {
         cout<<"PF - ETA: "<<(*pfParticles)[i].eta()<<" PHI: "<<(*pfParticles)[i].phi()<<" PT: "<<(*pfParticles)[i].pt()<<endl;
      }

      for (edm::View<pat::Jet>::const_iterator it = Jets_rec.begin(); it != Jets_rec.end(); ++it)
      {
         cout<<"JET - ETA: "<<it->eta()<<" PHI: "<<it->phi()<<" PT: "<<it->pt()<<" CSV: "<<it->bDiscriminator(bJetDisc_.c_str())<<endl;
      }
      cout<<"MET - ETA: "<<(*METs)[0].eta()<<" PHI: "<<(*METs)[0].phi()<<" PT: "<<(*METs)[0].pt()<<endl;
      cout<<"UCLUSTERED - ETA: "<<unclusteredMET.eta()<<" PHI: "<<unclusteredMET.phi()<<" PT: "<<unclusteredMET.pt()<<endl;
   }

   Jets_reb->reserve(Jets_rec.size());
   Jets_smeared->reserve(Jets_rec.size());

   bool isRebalanced = RebalanceJets_KinFitter(&Jets_rec, *(Jets_reb.get()));
   
   if (!isRebalanced) 
   {
      cout << "Bad event: Not possible to rebalance!" << endl;
   }
   
   GreaterByPt<pat::Jet> ptComparator_;
   std::sort(Jets_reb->begin(), Jets_reb->end(), ptComparator_);

   for(unsigned jet_i = 0; jet_i < Jets_reb->size(); jet_i++)
   {
      Jets_reb_P4->push_back((*Jets_reb)[jet_i].p4());
   }

   vMHT = LorentzVector(0., 0., 0., 0.);
   for (std::vector<pat::Jet>::const_iterator it = Jets_reb->begin(); it != Jets_reb->end(); ++it)
   {
      if (it->pt() > JetsMHTPt_ && abs(it->eta()) < JetsMHTEta_) 
      {
         LorentzVector tempP4(it->px(), it->py(), 0, sqrt(it->px()*it->px() + it->py()*it->py()));
         
         vMHT -= tempP4;
      }
   }
   

   px =unclusteredMET.px() + vMHT.px();
   py =unclusteredMET.py() + vMHT.py();

   
   LorentzVector rebalancedMET(px, py, 0, sqrt(px*px+py*py));

   METs_reb->push_back(reco::MET(rebalancedMET, (*METs)[0].vertex()));
   *MET_reb_P4 = (*METs_reb)[0].p4();
   vector<reco::MET> smearedMET;

   if(isRebalanced)
   {
      SmearingJets(*(Jets_reb.get()), 
                   unclusteredMET,
                   *(Jets_smeared.get()),
                   *(METs_smeared.get()), 
                   *(MHTs_smeared.get()), 
                   *(triplet_smeared.get()), 
                   *(rSystem_smeared.get()), 
                   *(bJetsInR_smeared.get()), 
                   *(dPhiJet1AndMET.get()),
                   *(dPhiJet2AndMET.get()),
                   *(dPhiJet3AndMET.get()),
                   *(m23OverM123.get()),
                   *(m123.get()),
                   *(MT2.get()),
                   *(MTt.get()),
                   *(MTb.get())                  
                  
         );
   }
   iEvent.put(Jets_reb_P4, "rebalancedJets");
   iEvent.put(MET_reb_P4, "rebalancedMET");
   iEvent.put(unclusteredMET_reb, "unclusteredMET");
   iEvent.put(Jets_smeared, "smearedJets");
   iEvent.put(METs_smeared, "smearedMET");
   iEvent.put(MHTs_smeared, "smearedMHT");
   iEvent.put(triplet_smeared,"smearedTriplet");
   iEvent.put(rSystem_smeared, "smearedRSystem");
   iEvent.put(bJetsInR_smeared,"smearedBJet");
   iEvent.put(dPhiJet1AndMET,"dPhiJet1AndMET");
   iEvent.put(dPhiJet2AndMET, "dPhiJet2AndMET");
   iEvent.put(dPhiJet3AndMET, "dPhiJet3AndMET");
   iEvent.put(m23OverM123, "m23OverM123");
   iEvent.put(m123, "m123");
   iEvent.put(MTt,"MTt");
   iEvent.put(MTb,"MTb");
   iEvent.put(MT2, "MT2");


}
//--------------------------------------------------------------------------


// ------------ method called once each job just before starting event loop  ------------
void RebalanceAndSmear::beginJob()
{

   TFile *f_rebCorr = new TFile(RebalanceCorrectionFile_.c_str(), "READ", "", 0);
   h_RebCorrectionFactor = (TH1F*) f_rebCorr->FindObjectAny("RebCorrection_vsReco_px"); 


}

// ------------ method called once each job just before starting event loop  ------------
void 
RebalanceAndSmear::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
RebalanceAndSmear::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
RebalanceAndSmear::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
RebalanceAndSmear::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
RebalanceAndSmear::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
RebalanceAndSmear::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
void RebalanceAndSmear::TripletSelector(std::vector<pat::Jet> & jets, std::vector<pat::Jet> & triplet, std::vector<pat::Jet> & rSystem, std::vector<pat::Jet> & bJet, double & M23OverM123, double & M123)
{
   

   double highPt = 0;
   double highCSV = -1;
   double pt;
   double csv;
   double disc;
   unsigned highPtIndex;
   unsigned highCSVIndex;
   
   vector<pat::Jet> bJets;


   for( unsigned i = 0; i<jets.size(); i++)
   {
      disc = jets[i].bDiscriminator(bJetDisc_.c_str());
      if(jets[i].pt() < bJetPtCut_ || abs(jets[i].eta()) > bJetEtaCut_ || 
         disc < bJetDiscCut_) continue;
      bJets.push_back(jets[i]);
   }
   
   if( bJets.size() == 0)
   {
      for( unsigned i = 0; i < jets.size(); i++)
      {
/*         csv = jets[i].bDiscriminator(bJetDisc_.c_str());
         cout<<"ETA: "<<jets[i].eta()<<" PHI: "<<jets[i].phi()<<" PT: "<<jets[i].pt()<<" CSV: "<<csv<<endl;
*/
         if(abs(jets[i].eta()) > bJetEtaCut_) continue;
         
         pt = jets[i].pt();
         csv = jets[i].bDiscriminator(bJetDisc_.c_str());
         
         if(pt > highPt)
         {
            highPt = pt;
            highPtIndex = i;
         }
         if(csv > highCSV)
         {
            highCSV = csv;
            highCSVIndex = i;
         }
      }
   }
   if(highCSV != -1 && highCSV != -10)
      bJets.push_back(jets[highCSVIndex]);

         
   
   LorentzVector p41;
   LorentzVector p42;
   LorentzVector p43;
   
   bool bJetPass;
   bool bJetInTrip;
   bool bJetOutTrip;
   bool passDijet;
   bool passM23OverM123;
   bool atLeastOnePassed = false;
   bool skipFailedM23OverM123 = false;

   unsigned selectedIndex;
   unsigned rank = 5;
   double bJetEta;
   double bJetPhi;
   double m123;
   double m12;
   double m13;
   double m23;
   double dTopMin = 999999;


   vector<double> m123s;
   vector<double> m12s;
   vector<double> m13s;
   vector<double> m23s;
   vector<vector<int> > bJetInRIndices;

   vector<bool> bJetsBehave;
   vector<bool> passesDijetCuts;
   vector<bool> passesM23OverM123Cut;

   auto_ptr<vector<bool> > selectedPassDijetCuts(new vector<bool>());
   auto_ptr<vector<BasicJet> > selectedTriplet(new vector<BasicJet>());
   auto_ptr<vector<BasicJet> > selectedRSystem(new vector<BasicJet>());
   auto_ptr<vector<pat::Jet> > selectedBJets(new vector<pat::Jet>());
   auto_ptr<vector<int> > selectedTripletIndex(new vector<int>());
   auto_ptr<vector<int> > selectedBJetsIndex(new vector<int>());
   auto_ptr<vector<double> > selectedM123(new vector<double>());
   auto_ptr<vector<double> > selectedM12(new vector<double>());
   auto_ptr<vector<double> > selectedM13(new vector<double>());
   auto_ptr<vector<double> > selectedM23(new vector<double>());
   auto_ptr<double>  selectedM23OverM123(new double());
   auto_ptr<vector<bool> > hasBJet(new vector<bool>());
   auto_ptr<vector<bool> > passDijetCuts(new vector<bool>());
   auto_ptr<vector<bool> > passM23OverM123Cut(new vector<bool>());


   vector<vector<int> > indices;


   for(int i = 0; i < int(jets.size()); i++)
   {
      for(int j = i+1; j < int(jets.size()); j++)
      {
         for(int k = j+1; k < int(jets.size()); k++)
         {

            LorentzVector tempLor = jets[i].p4() + jets[j].p4() + jets[k].p4();
            double dR1 = deltaR( jets[i].eta(), jets[i].phi(),
                                 tempLor.eta(), tempLor.phi());
            double dR2 = deltaR( jets[j].eta(), jets[j].phi(),
                                 tempLor.eta(), tempLor.phi());
            double dR3 = deltaR( jets[k].eta(), jets[k].phi(),
                                 tempLor.eta(), tempLor.phi());

            if(dR1 > tripletDRCut_) continue;
            if(dR2 > tripletDRCut_) continue;
            if(dR3 > tripletDRCut_) continue;
            
            vector<int> tempVec;     
            tempVec.push_back(i);
            tempVec.push_back(j);
            tempVec.push_back(k);
            indices.push_back(tempVec);
         }
      }
   }


   for(unsigned i = 0; i < indices.size(); i++)
   {
//      cout<<indices[i][0]<<" "<<indices[i][1]<<" "<<indices[i][2]<<" "<<endl;
      bJetPass = true;
      passDijet = false;
      passM23OverM123 = false;

      vector<int> tempBJetIndices;
      if(indices[i].size() != 3) 
         cout<<indices[i].size()<<endl;

      p41 = jets[indices[i][0]].p4();
      p42 = jets[indices[i][1]].p4();
      p43 = jets[indices[i][2]].p4();
      
      m123 = (p41 + p42 + p43).mass();
      m12 = (p41 + p42).mass();
      m13 = (p41 + p43).mass();
      m23 = (p42 + p43).mass();

      //Check if the bjet(s) is in the right place

      if( bJets.size() == 0)
      {
         bJetPass = false;
      }
      else if( bJets.size() == 1)
      {         
         bJetEta = bJets[0].eta();
         bJetPhi = bJets[0].phi();
         if( (deltaR(p41.eta(), p41.phi(), bJetEta, bJetPhi) < 0.01)   || 
             (deltaR(p42.eta(), p42.phi(), bJetEta, bJetPhi) < 0.01)   || 
             (deltaR(p43.eta(), p43.phi(), bJetEta, bJetPhi) < 0.01) )
         {
            bJetPass = false;
         }

         else
         {
            tempBJetIndices.push_back(0);
         }
      }
      else if( bJets.size() >= 2)
      {
         bJetInTrip = false;
         bJetOutTrip = false;
         for (unsigned j = 0; j < bJets.size(); j++)
         {
            bJetEta = bJets[j].eta();
            bJetPhi = bJets[j].phi();
            if( (deltaR(p41.eta(), p41.phi(), bJetEta, bJetPhi) < 0.01)  || 
                (deltaR(p42.eta(), p42.phi(), bJetEta, bJetPhi) < 0.01)  || 
                (deltaR(p43.eta(), p43.phi(), bJetEta, bJetPhi) < 0.01) )
            {
               bJetInTrip = true;
            }
            else
            {
               tempBJetIndices.push_back(j);
               bJetOutTrip = true;
            }
         }

         if( !bJetInTrip || !bJetOutTrip)
            bJetPass = false;

      }

      
      //check if it passes the dijetmass cuts
      
      if(arctanmin_ < atan(m13/m12) && atan(m13/m12) < arctanmax_ && 
         Rmin_ < m23/m123 && m23/m123 < Rmax_) 
         passDijet = true;

      if(Rmin_*Rmin_ * (1+m13*m13/(m12*m12))  < 1 - m23*m23/(m123*m123) &&
         Rmax_*Rmax_ * (1+m13*m13/(m12*m12))  > 1 - m23*m23/(m123*m123))
         passDijet = true;

      if(Rmin_*Rmin_ * (1+m12*m12/(m13*m13))  < 1 - m23*m23/(m123*m123) &&
         Rmax_*Rmax_ * (1+m12*m12/(m13*m13))  > 1 - m23*m23/(m123*m123)) 
         passDijet = true;

      if(m23/m123 > m23OverM123Cut_) 
         passM23OverM123 = true;

      passesDijetCuts.push_back(passDijet);
      passesM23OverM123Cut.push_back(passM23OverM123);    
      bJetsBehave.push_back(bJetPass);
      bJetInRIndices.push_back(tempBJetIndices);

      m123s.push_back(m123);
      m12s.push_back(m12);
      m13s.push_back(m13);
      m23s.push_back(m23);

      hasBJet->push_back(bJetPass);
      passDijetCuts->push_back(passDijet);
      passM23OverM123Cut->push_back(passM23OverM123);
   }

   
   //Choose the selected triplet

   for(unsigned i = 0; i < indices.size(); i++)
   {
      if(!passesDijetCuts[i] || !bJetsBehave[i]) continue;
      
      //Preferentially pick ones that pass the m23/m123 cut
      if(passesM23OverM123Cut[i])
      {
         if(!skipFailedM23OverM123) 
            dTopMin = 999999;
         
         skipFailedM23OverM123 = true;
         
         if(abs(m123s[i] - topMass_) < dTopMin)
         {
            atLeastOnePassed = true;
            selectedIndex = i;
            dTopMin = abs(m123s[i] - topMass_);
         }
      }
      
      //Only consider these if you don't have a triplet passing m23/m123 
      if(!skipFailedM23OverM123)
      {
         if(abs(m123s[i] - topMass_) < dTopMin)
         {
            atLeastOnePassed = true;
            selectedIndex = i;
            dTopMin = abs(m123s[i] - topMass_);
         }
      }
   }
   
   if(atLeastOnePassed)
   {
      
      Jet::Constituents tripletConsts;
      LorentzVector tripletP4(0, 0, 0, 0);
      Candidate::Point tripletVertex = jets[indices[selectedIndex][0]].vertex();
      for( int i = 0; i < int(jets.size()); i++)
      {
         if( i == indices[selectedIndex][0] || 
             i == indices[selectedIndex][1] || 
             i == indices[selectedIndex][2])
         {
            
            tripletP4 = tripletP4 + jets[i].p4();
            for(unsigned j = 0; j < jets[i].numberOfDaughters(); j++)
               tripletConsts.push_back(jets[i].daughterPtr(j));
            
         }
         else
         {
            rSystem.push_back(reco::Jet(jets[i].p4(), jets[i].vertex()));
         }
      }
      
      triplet.push_back(reco::Jet(tripletP4, tripletVertex));
      for(unsigned i = 0; i < bJetInRIndices[selectedIndex].size(); i++)
      {
         bJet.push_back(bJets[bJetInRIndices[selectedIndex][i]]);
      }

      M123 = m123s[selectedIndex];
      M23OverM123 = m23s[selectedIndex]/m123s[selectedIndex];
      

   }
   else
   {
      M123 = -1;
      M23OverM123 = -1;
   }


/*
   iEvent.put(selectedPassDijetCuts, "selectedPassDijetCuts");
   iEvent.put(selectedTriplet, "selectedTriplet");
   iEvent.put(selectedRSystem, "selectedRSystem");
   iEvent.put(selectedBJets, "selectedBJets");
   iEvent.put(selectedTripletIndex, "selectedTripletIndex");
   iEvent.put(selectedBJetsIndex, "selectedBJetsIndex");
   iEvent.put(selectedM12, "selectedM12");
   iEvent.put(selectedM13, "selectedM13");
   iEvent.put(selectedM23, "selectedM23");
   iEvent.put(selectedM123, "selectedM123");
   iEvent.put(selectedM23OverM123, "selectedM23OverM123");
   iEvent.put(hasBJet, "hasBJet");
   iEvent.put(passDijetCuts, "passDijetCuts");
   iEvent.put(passM23OverM123Cut, "passM23OverM123Cut");
*/
}
void RebalanceAndSmear::MTMT2(reco::MET & MET, std::vector<pat::Jet> & triplet,std::vector<pat::Jet> & rSystem, std::vector<pat::Jet> & bJetsInR, double & MT2, double & MTt, double & MTb)
{
   
   unsigned bJetIndex;
   unsigned bJetIndexForDoublet;
   unsigned closestBJetToMETIndex;
   unsigned wJet1Index;
   unsigned wJet2Index;
   int otherJetIndex = -1;
   double dTopMin = 9999999;
   double dRMin = 9999999;
   double dPhiMin = 9999999;
   double dR;
   double dPhi;
   double Et_1;
   double Et_2;

   bool outsideWWindow = false;
   bool foundTriplet = false;

   LorentzVector bJetP4;
   LorentzVector wJet1P4;
   LorentzVector wJet2P4;
   LorentzVector otherTopP4;
   LorentzVector METP4 = MET.p4();
   LorentzVector tripletP4;
   
   Mt2::ChengHanBisect_Mt2_332_Calculator mt2Calculator;   

   
   if( triplet.size() != 0 && rSystem.size() != 0 && bJetsInR.size() != 0)
   {
      tripletP4 = triplet[0].p4();
      for(unsigned i = 0; i < bJetsInR.size(); i++)
      {
         bJetP4 = bJetsInR[i].p4();
         for(unsigned j = 0; j < rSystem.size(); j++)
         {
            wJet1P4 = rSystem[j].p4();
            if((deltaR(wJet1P4.eta(), wJet1P4.phi(), 
                       bJetP4.eta(), bJetP4.phi()) < 0.01)) continue;
            
            for(unsigned k = j + 1; k < rSystem.size(); k++)
            {
               wJet2P4 = rSystem[k].p4();
               
               if((deltaR(wJet2P4.eta(), wJet2P4.phi(), 
                          bJetP4.eta(), bJetP4.phi()) < 0.01) )  continue;
            
               if( fabs((bJetP4 + wJet1P4 + wJet2P4).mass() - mTop_) < dTopMin)
               {

                  dTopMin = fabs((bJetP4 + wJet1P4 + wJet2P4).mass() - mTop_);
                  foundTriplet = true;
                  bJetIndex = i;
                  wJet1Index = j;
                  wJet2Index = k;
               }
            }
         }
      }
      if(foundTriplet)
      {  
         bJetP4 = bJetsInR[bJetIndex].p4();
         wJet1P4 = rSystem[wJet1Index].p4();
         wJet2P4 = rSystem[wJet2Index].p4();
         if ( (wJet1P4 + wJet2P4).mass() < mWMin_ || 
              (wJet1P4 + wJet2P4).mass() > mWMax_) outsideWWindow = true;
      }

      if(!foundTriplet)
      {
         dPhiMin = 999999;
         for(unsigned i = 0; i < bJetsInR.size(); i++)
         {
            dPhi = abs(deltaPhi(bJetsInR[i].phi(), MET.phi()));
            if( dPhi < dPhiMin)
            {
               dPhiMin = dPhi;
               bJetP4 = bJetsInR[i].p4();
               bJetIndexForDoublet = i;
            }
            
         }
      }     
      else
         bJetIndexForDoublet = bJetIndex;
      dRMin = 99999;
      for(unsigned i = 0; i < rSystem.size(); i++)
      {            
         dR = deltaR(bJetP4.eta(), bJetP4.phi(), 
                     rSystem[i].eta(), rSystem[i].phi());
         
         if(dR < 0.01) continue;
         if(dR > 2.0) continue;
         if( (bJetP4 + rSystem[i].p4()).mass() > mTop_) continue;
         
         if(dR < dRMin)
            {            
               otherJetIndex = int(i);
               dRMin = dR;
            }
      }
      
      if(foundTriplet && !outsideWWindow)
      {
         otherTopP4 = bJetP4 + wJet1P4 + wJet2P4;
      }
      else if(outsideWWindow && otherJetIndex != -1)
      {
         otherTopP4 = rSystem[otherJetIndex].p4() + bJetP4;
      }
      else if(otherJetIndex != -1)
      {
         otherTopP4 = rSystem[otherJetIndex].p4() + bJetP4;
      }
      else
      {
         otherTopP4 = bJetP4;            
      }  


      const double massOfSystemA =  tripletP4.M(); 
      const double pxOfSystemA   =  tripletP4.Px(); 
      const double pyOfSystemA   =  tripletP4.Py(); 
      
      const double massOfSystemB =  otherTopP4.M(); 
      const double pxOfSystemB   =  otherTopP4.Px(); 
      const double pyOfSystemB   =  otherTopP4.Py(); 
      
      const double pxMiss        = METP4.Px(); 
      const double pyMiss        = METP4.Py(); 
      
      const double invis_mass    = METP4.M(); 
      
      Mt2::LorentzTransverseVector  vis_A(Mt2::TwoVector(pxOfSystemA, 
                                                         pyOfSystemA), massOfSystemA);
      Mt2::LorentzTransverseVector  vis_B(Mt2::TwoVector(pxOfSystemB, pyOfSystemB), massOfSystemB);
      Mt2::TwoVector                pT_Miss(pxMiss, pyMiss);
   
      MT2 = mt2Calculator.mt2_332(vis_A, vis_B, pT_Miss, invis_mass);
//      cout<<"MT2: "<<MT2<<endl;
   }
   else
      MT2 = -1;
   if(triplet.size() != 0)
   {
      tripletP4 = triplet[0].p4();

      Et_1 = sqrt(tripletP4.mass()*tripletP4.mass() + 
                  tripletP4.pt()*tripletP4.pt());
      Et_2 = sqrt(METP4.mass()*METP4.mass() + METP4.pt()*METP4.pt());

      MTt = sqrt( tripletP4.mass() * tripletP4.mass() + 
                           METP4.mass()*METP4.mass() + 
                           2*(Et_1*Et_2 - tripletP4.px() * METP4.px() - 
                              tripletP4.py()*METP4.py()));
//      cout<<"MTt: "<<MTt<<endl;
   }
   else
   {
      MTt = -1;
   }
   if(bJetsInR.size() != 0)
   {
      dPhiMin = 999999;
      
      for(unsigned i = 0; i < bJetsInR.size(); i++)
      {
         dPhi = abs(deltaPhi( METP4.phi(), bJetsInR[i].phi()));
         if( dPhi < dPhiMin)
         {
            bJetP4 = bJetsInR[i].p4();
            closestBJetToMETIndex = i;
            dPhiMin = dPhi;
         }         

      } 
      otherTopP4 = bJetP4;

      if( otherJetIndex  != -1 && closestBJetToMETIndex == bJetIndexForDoublet)
      {
         otherTopP4 = rSystem[otherJetIndex].p4() + bJetP4;            
      }

      

      Et_1 = sqrt(otherTopP4.mass()*otherTopP4.mass() + otherTopP4.pt() * otherTopP4.pt());
      Et_2 = sqrt(METP4.mass()*METP4.mass() + METP4.pt()*METP4.pt());

      MTb = sqrt( otherTopP4.mass() * otherTopP4.mass() +
                           METP4.mass() * METP4.mass() + 
                           2*(Et_1*Et_2 - otherTopP4.px() * METP4.px() - 
                              otherTopP4.py() * METP4.py()));
//      cout<<"MTb: "<<MTb<<endl;
   }
   else
      MTb = -1;
   

}

//define this as a plug-in
DEFINE_FWK_MODULE(RebalanceAndSmear);

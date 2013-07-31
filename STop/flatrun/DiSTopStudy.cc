// -*- C++ -*-
//
// Package:    DiSTopStudyTree
// Class:      DiSTopStudyTree
// 
/**\class DiSTopStudyTre DiSTopStudyTree.cc SusyAnalysis/DiSTopStudyTree/src/DiSTopStudyTree.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Seema Sharma
//         Created:  Fri Aug 24 12:18:17 CDT 2012
// $Id: DiSTopStudy.cc,v 1.2 2013/04/21 23:13:52 samantha Exp $
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
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

//#include "SandBox/Skims/interface/SAKLooseLepton.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "SusyAnalysis/AnalysisUtils/interface/PrintEventGenInfo.h"
#include "DataFormats/Math/interface/deltaR.h"

// TFile Service
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include "TLorentzVector.h"

//J.. FatJet HPTopTag
#include "DataFormats/JetReco/interface/BasicJet.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Math/interface/deltaR.h"


#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "recipeAUX/OxbridgeMT2/interface/Basic_Mt2_332_Calculator.h"
#include "recipeAUX/OxbridgeMT2/interface/ChengHanBisect_Mt2_332_Calculator.h"

#include "UserCode/TopTagger/interface/Type3TopTagger.h"
//#include "UserCode/TopTagger/src/TopTagger.cc"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "UserCode/TopTagger/interface/combination.h"
#include "UserCode/TopTagger/interface/indexSort.h"
#include "TLorentzVector.h"

using namespace stdindexSort;
using namespace stdcomb;
using namespace std;

typedef unsigned int size;
typedef std::vector<int>::iterator vii;

const std::string defaultOrderingOptArr[] = {"mass", "mass"};
const int defaultMaxIndexForOrderingArr[] = { 1,         1 };


//For sorting by pt
struct GreaterByPtCandPtr 
{
  bool operator()( const edm::Ptr<reco::Candidate> & t1, const edm::Ptr<reco::Candidate> & t2 ) const 
  {return t1->pt() > t2->pt();}
};
//..J

class DiSTopStudyTree : public edm::EDAnalyzer {

  public:

  explicit DiSTopStudyTree(const edm::ParameterSet & iConfig);
  ~DiSTopStudyTree();

  private:

  void beginJob() ;
  void endJob() ;
  void analyze(const edm::Event&, const edm::EventSetup&);
  void         BookHistograms();
  void         clearTreeVectors();  
  edm::ESHandle<ParticleDataTable> pdt_;

  bool          debug_;
  bool          isDimuonSample_;
  bool          isDieleSample_;
  bool          isMuEleSample_;
  bool          isDiLeptonSample_;
  bool          isMC_;
  edm::InputTag vtxSrc_;
  bool          doPUReWeight_;
  //J  edm::InputTag puWeightSrc_;
  edm::InputTag pfMetSrc_;
  edm::InputTag jetAllsrc_;
  std::string   btagname_;
  edm::InputTag mhtSrc_, htSrc_;
  edm::InputTag muonVetoSrc_, eleVetoSrc_;
  edm::InputTag pfCandidateSrc_;

  std::vector<edm::InputTag> varsDoubleTagsV_;//JL
  std::vector<std::string> varsDoubleNamesInTreeV_;//JL
  std::vector<UShort_t> varsDoubleVN_;//JL
  std::vector<Float_t*> varsDoubleV_;//JL

  bool saveAllMuons_;
  edm::InputTag muonSrc_;
  edm::InputTag electronSrc_;
  double        minMuPt_, maxMuEta_, maxMuD0_, maxMuDz_, maxMuRelIso_;
  
  edm::Service<TFileService> fs;
  TTree* tree;
  unsigned int         t_EvtRun, t_EvtLS, t_EvtEvent;
  int                  t_NVertices;
  double               t_PUWeight;
  double               t_TrueNPV;
  double               t_PFMetPx, t_PFMetPy,t_PFMetE,t_PFMetPhi,t_PFMetSignificance;
  double               t_PFOrigMetPx, t_PFOrigMetPy,t_PFOrigMetE,t_PFOrigMetPhi;
  int                  iNJets;
  std::vector<double> *t_PFJetPt,  *t_PFJetEta,  *t_PFJetPhi,  *t_PFJetE,  *t_PFJetBTag, *t_DanTopTaggerIndex;
  double               t_PFht, t_PFmht,t_PFmphi;
  double               t_Genht, t_Genmht;
  int                  t_NVetoMuon, t_NVetoEle;
  int                  t_NJetsPt30Eta2p5, t_NJetsPt30Eta5p0, t_NJetsPt50Eta2p5, t_NJetsPt50Eta5p0;
  int                  t_NJetsPt30Eta2p4, t_NJetsPt50Eta2p4, t_NJetsPt70Eta2p4;

  vector<double> puWeights_;

  // Muon variables
  std::vector<double> *t_PFMuonPt, *t_PFMuonEta, *t_PFMuonPhi, *t_PFMuonE, *t_PFMuonCh;
  std::vector<double> *t_PFMuonChiSq, *t_PFMuonValidMuonHits, *t_PFMuonMatchedStations;
  std::vector<double> *t_PFMuonValidPixelHits, *t_PFMuonTrackerLayerMeasured;
  std::vector<double> *t_PFMuonCharHadEt, *t_PFMuonNeutHadEt, *t_PFMuonPhotEt, *t_PFMuonSumPUPt, *t_PFMuonRelIso;
  std::vector<double> *t_PFMuonCharHadEt03, *t_PFMuonNeutHadEt03, *t_PFMuonPhotEt03, *t_PFMuonSumPUPt03, *t_PFMuonRelIso03;
  std::vector<bool>   *t_PFMuonID, *t_PFMuonVtxAss, *t_PFMuonIsolated, *t_PFMuonIDVtxIso;
  std::vector<double> *t_PFMuonDirIso01, *t_PFMuonDirIso02, *t_PFMuonDirIso03, *t_PFMuonDirIso04, *t_PFMuonDirIso05;
  std::vector<bool>   *t_PFMuonIsWdau, *t_PFMuonIsBdau;
  std::vector<double> *t_PFMuonGenMuPt, *t_PFMuonGenMuEta, *t_PFMuonGenMuPhi;
  double               t_PFMuonZPt, t_PFMuonZEta, t_PFMuonZPhi, t_PFMuonZE, t_PFMuonZMass;

  //std::vector<double>  *t_QCDReweight;

  // Electron variables
  std::vector<double> *t_PFElecPt, *t_PFElecEta, *t_PFElecPhi, *t_PFElecE, *t_PFElecCh; 
  std::vector<bool>   *t_PFEleIsEB, *t_PFEleIsEE;
  std::vector<double> *t_PFEleDEtaIn, *t_PFEleDPhiIn, *t_PFEleSigIEtaIEta;
  std::vector<double> *t_PFEleHOE, *t_PFEleD0Vtx, *t_PFEleZVtx;
  std::vector<double> *t_PFEleIsoCH, *t_PFEleIsoEM, *t_PFEleIsoNH, *t_PFEleIsoRho, *t_PFEleIsoEffA, *t_PFEleRelIso;
  std::vector<bool>   *t_PFEleID, *t_PFEleIsolated, *t_PFEleIDIso;
  std::vector<double> *t_PFEleDirIso01, *t_PFEleDirIso02, *t_PFEleDirIso03, *t_PFEleDirIso04, *t_PFEleDirIso05;
  std::vector<bool>   *t_PFEleIsWdau, *t_PFEleIsBdau;
  std::vector<double> *t_PFEleGenElePt, *t_PFEleGenEleEta, *t_PFEleGenElePhi;
  double               t_PFElecZPt, t_PFElecZEta, t_PFElecZPhi, t_PFElecZE, t_PFElecZMass;

  std::vector<int>     *t_genPartPdgId, *t_genPartStatus;
  std::vector<double>  *t_genPartPt, *t_genPartEta, *t_genPartPhi, *t_genPartE;
  std::vector <int>                   *t_genPartDecayMode;

  //Lepton variables (for general lepton storage)
  std::vector<double> *t_PFLepPt, *t_PFLepEta, *t_PFLepPhi, *t_PFLepE, *t_PFLepCh; 
  std::vector<bool>   *t_PFLepIsMu, *t_PFLepIsEle;
  double               t_PFLepZPt, t_PFLepZEta, t_PFLepZPhi, t_PFLepZE, t_PFLepZMass;



  TLorentzVector lvMuon1;
  TLorentzVector lvMuon2;

  TLorentzVector lvElec1;
  TLorentzVector lvElec2;

  TLorentzVector lvLep1;
  TLorentzVector lvLep2;

  //Dans Tagger
  double               t_DTWMass1,t_DTWMass2;
  double               t_DTWMassChiSqr1,t_DTWMassChiSqr2;
  double               t_DTdPhi1,t_DTdPhi2;
  double               t_DTTopMass1,t_DTTopMass2,t_DTTotalSystemMass;
  double               t_DTTopPt1,t_DTTopPt2,t_DTTotalSystemPt;
  double               t_DTTopMu1,t_DTTopMu2;
  double               t_DTTopEta1,t_DTTopEta2,t_DTTotalSystemEta;
  double               t_DTTopPhi1,t_DTTopPhi2,t_DTTotalSystemPhi;
  int                  indexB1, indexB2, indexW11, indexW12, indexW21, indexW22;


  //J..
  double               t_fatJet1_m123, t_fatJet1_m12, t_fatJet1_m13, t_fatJet1_m23,t_fatJet1_pt123,t_fatJet1_eta123,t_fatJet1_phi123;
  double               t_fatJet2_m123, t_fatJet2_m12, t_fatJet2_m13, t_fatJet2_m23,t_fatJet2_pt123,t_fatJet2_eta123,t_fatJet2_phi123;
  int                  t_nFatJet;
  int                  t_fatJet1_cond1, t_fatJet1_cond2, t_fatJet1_cond3;
  int                  t_fatJet2_cond1, t_fatJet2_cond2, t_fatJet2_cond3;
  std::vector < std::vector<double> > *t_fatJetSubJetPt, *t_fatJetSubJetEta, *t_fatJetSubJetPhi, *t_fatJetSubJetE;
 //..J


  Handle<LHEEventProduct> product;



  //  TopTagger *topTagger;

  double bestTopJetMass;
  bool remainPassCSVS;
  double mTbestTopJet;
  double mTbJet;
  double MT2;
  double mTbestWJet;
  double mTbestbJet;



// doExtraCuts: extra cuts other than the top tagger related, e.g., dphi cuts, HT cut and so on.
    const bool doExtraCuts_;

    const double mTbcut_, mTtcut_, MT2cut_, mTWcut_;
// doMTMT2cuts: numbers are counted. But if enabled, some plots will be filled after all the cuts.
    const bool doMTMT2cuts_;

// PDG values (in GeV)
    const double mW_, mTop_;

    const double mWoverTop_;

// const double lowRatioWoverTop = 0.85, highRatioWoverTop = 1.15;
    const double lowRatioWoverTop_, highRatioWoverTop_;
    const double lowArcTanm13overm12_, highArcTanm13overm12_, lowm23overm123_;

    const double Rmin_, Rmax_;

    const double defaultJetCone_;
    const double simuCAdeltaR_;
// Eta ranges from 0 to 5, phi ranges from 0 to 3.14.
// sqrt(5./2.*5./2. + 3.14/2. * 3.14/2.) ~ 2.95
// const double simuCALargerdeltaR = 3.0;
    const double simuCALargerdeltaR_; // -1 means no deltaR requirement

    const double lowTopCut_, highTopCut_;
    const double lowWCut_, highWCut_;

// Choose CSVM point for now
// --> A good top fat jet might indicate a b-jet already,
// so requiring another tight b-jet is NOT good. (TODO: to be studied)
//    const double CSVL = 0.244, CSVM = 0.679, CSVT = 0.898;
    const double CSVS_;

    const int nSubJetsDiv_;

    const int nJetsSel_;

    const double maxEtaForbJets_;

// mass   : mass ordering --> pick the one with mass closest to the norminal mass
// pt     : pt ordering --> pick the one with larger pt (for the first two fat jets)
// hybrid : pt ordering + mass ordering --> if both of the first two fat jets
//          satisfying criteria, pick the one closer to the norminal mass
    const std::vector<std::string> defaultOrderingOptVec;
//                                                     best   remaining
    const std::vector<std::string>  orderingOptArr_; // (mass, mass) is the best?
    const std::vector<int> defaultMaxIndexForOrderingVec;
    const std::vector<int> maxIndexForOrderingArr_;

    const edm::InputTag metSrc_, jetSrc_;
    const std::string bTagKeyString_;

    edm::Handle<edm::View<reco::Jet> > jets;
    edm::Handle<std::vector<pat::Jet> > patjets;
    size nJets;
    bool isPatJet;
    virtual void loadRecoJets(const edm::Event& iEvent);

    edm::Handle<edm::View<reco::MET> > metHandle;
    double met, metphi;
    virtual void loadMETMHT(const edm::Event& iEvent);

    StringCutObjectSelector<reco::Jet,true> pfJetCutForJetCounting_; // lazy parsing, to allow cutting on variables not in reco::Candidate class
    StringCutObjectSelector<reco::Jet,true> pfJetCutForJetCombining_; // lazy parsing, to allow cutting on variables not in reco::Candidate class
    const edm::InputTag evtWeightInput_;
    edm::Handle<double> evtWeight_;

    unsigned int run, event, lumi; bool isData;
    unsigned int vtxSize;
    edm::Handle<edm::View<reco::Vertex> > vertices;
    void loadEventInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup);

   
    const bool taggingMode_;
 
    const bool dobVetoCS_;
  
    int pickJetsForCombining(std::vector<TLorentzVector>& oriJetsVec, std::vector<double> &recoJetsBtagCSVS);
//    int pickJetsForCombining(std::vector<TLorentzVector>& oriJetsVec, std::vector<double> &recoJetsBtagCSVS, std::vector<std::vector<TLorentzVector> > & oriJetsConstituentsVec, std::vector<std::vector<int> > & oriJetsConstituentsChargeVec, std::vector<std::vector<double> > & oriJetsConstituentsDzVec, std::vector<std::vector<double> > & oriJetsConstituentsD0Vec);
    int countNJets();
  

  
// The counters
    double cntPassnJetsCut, cntTaggedTopEvents, cntTaggedTopEventsWithinMassCuts;
    double cntTaggedAllCutsPlusCSVS;
    double cntPassingMTbestTopJetCut, cntPassingMTclosebJetCut, cntPassingMTbestWJetCut;
    double cntPassingMT2Cut, cntPassingMT2andMTCut;
    double cntTaggedbestFatJetPlusCSVS;

    topTagger::type3TopTagger *type3TopTaggerPtr;
  
  bool kFoundTwoMuons;
  bool kFoundTwoEles;

};


DiSTopStudyTree::DiSTopStudyTree(const edm::ParameterSet & iConfig):
  doExtraCuts_(iConfig.getUntrackedParameter<bool>("doExtraCuts", true)),
   mTbcut_(iConfig.getUntrackedParameter<double>("mTbcut", 500)),
   mTtcut_(iConfig.getUntrackedParameter<double>("mTtcut", 365)),
   MT2cut_(iConfig.getUntrackedParameter<double>("MT2cut", 300)),
   mTWcut_(iConfig.getUntrackedParameter<double>("mTWcut", 600)),
   doMTMT2cuts_(iConfig.getUntrackedParameter<bool>("doMTMT2cuts", true)),
   mW_(iConfig.getUntrackedParameter<double>("mW", 80.385)),
   mTop_(iConfig.getUntrackedParameter<double>("mTop", 173.5)),
   mWoverTop_(mW_/mTop_),
   lowRatioWoverTop_(iConfig.getUntrackedParameter<double>("lowRatioWoverTop", 0.85)),
   highRatioWoverTop_(iConfig.getUntrackedParameter<double>("highRatioWoverTop", 1.25)),
   lowArcTanm13overm12_(iConfig.getUntrackedParameter<double>("lowArcTanm13overm12", 0.2)),
   highArcTanm13overm12_(iConfig.getUntrackedParameter<double>("highArcTanm13overm12", 1.3)),
   lowm23overm123_(iConfig.getUntrackedParameter<double>("lowm23overm123", 0.35)),
   Rmin_(lowRatioWoverTop_*mWoverTop_),
   Rmax_(highRatioWoverTop_*mWoverTop_),
   defaultJetCone_(iConfig.getUntrackedParameter<double>("defaultJetCone", 2.0)),
   simuCAdeltaR_(iConfig.getUntrackedParameter<double>("simuCAdeltaR", 1.5)),
   simuCALargerdeltaR_(iConfig.getUntrackedParameter<double>("simuCALargerdeltaR", -1)),
   lowTopCut_(iConfig.getUntrackedParameter<double>("lowTopCut", 80)),
   highTopCut_(iConfig.getUntrackedParameter<double>("highTopCut", 270)),
   lowWCut_(iConfig.getUntrackedParameter<double>("lowWCut", 50)),
   highWCut_(iConfig.getUntrackedParameter<double>("highWCut", 120)),
   CSVS_(iConfig.getUntrackedParameter<double>("CSVS", 0.679)),
   nSubJetsDiv_(iConfig.getUntrackedParameter<int>("nSubJetsDiv", 3)),
   nJetsSel_(iConfig.getUntrackedParameter<int>("nJetsSel", 5)),
   maxEtaForbJets_(iConfig.getUntrackedParameter<double>("maxEtaForbJets", 2.4)),
   defaultOrderingOptVec(defaultOrderingOptArr, defaultOrderingOptArr+sizeof(defaultOrderingOptArr)/sizeof(defaultOrderingOptArr[0])),
   orderingOptArr_(iConfig.getUntrackedParameter<std::vector<std::string> >("orderingOptArr", defaultOrderingOptVec)),
   defaultMaxIndexForOrderingVec(defaultMaxIndexForOrderingArr, defaultMaxIndexForOrderingArr+sizeof(defaultMaxIndexForOrderingArr)/sizeof(defaultMaxIndexForOrderingArr[0])),
   maxIndexForOrderingArr_(iConfig.getUntrackedParameter<std::vector<int> >("maxIndexForOrderingArr", defaultMaxIndexForOrderingVec)),
   metSrc_(iConfig.getParameter<edm::InputTag>("metSrc")),
   jetSrc_(iConfig.getParameter<edm::InputTag>("jetSrc")),
   bTagKeyString_(iConfig.getUntrackedParameter<std::string>("bTagKeyString", "combinedSecondaryVertexBJetTags")),
   pfJetCutForJetCounting_(iConfig.existsAs<std::string>("pfJetCutForJetCounting") ? iConfig.getParameter<std::string>("pfJetCutForJetCounting") : "", true),
   pfJetCutForJetCombining_(iConfig.existsAs<std::string>("pfJetCutForJetCombining") ? iConfig.getParameter<std::string>("pfJetCutForJetCombining") : "", true),
   evtWeightInput_(iConfig.getParameter<edm::InputTag>("evtWeightInput")),
   vtxSrc_(iConfig.getParameter<edm::InputTag>("vtxSrc")),
   debug_(iConfig.getUntrackedParameter<bool>("debug", true)),
   taggingMode_(iConfig.getUntrackedParameter<bool>("taggingMode", true)),
   dobVetoCS_(iConfig.getUntrackedParameter<bool>("dobVetoCS", false))

 {
  debug_       = iConfig.getParameter<bool>("Debug");
  //isDimuonSample_= iConfig.getParameter<bool>("isDimuonSample");
  //isDieleSample_ = iConfig.getParameter<bool>("isDieleSample");
  //isMuEleSample_ = iConfig.getParameter<bool>("isMuEleSample");
  isDiLeptonSample_ = iConfig.getParameter<bool>("isDiLeptonSample");
  isMC_= iConfig.getParameter<bool>("isMC");
  
  varsDoubleTagsV_ = iConfig.getParameter< std::vector<edm::InputTag> >("VarsDoubleV");//JL
  varsDoubleNamesInTreeV_ = iConfig.getParameter< std::vector<std::string> >("VarsDoubleNamesInTreeV");//JL
  
  //  topTagger = new TopTagger(iConfig);

  vtxSrc_      = iConfig.getParameter<edm::InputTag>("VertexSource");
  doPUReWeight_= iConfig.getParameter<bool>("DoPUReweight");
//J  puWeigthSrc_ = iConfig.getPararmeter<edm::InputTag>("PUWeigthSource");
  pfMetSrc_    = iConfig.getParameter<edm::InputTag>("PFMetSource");
  jetAllsrc_   = iConfig.getParameter<edm::InputTag>("JetAllSource");
  btagname_    = iConfig.getParameter<std::string>  ("bTagName");
  mhtSrc_      = iConfig.getParameter<edm::InputTag>("MHTSource");
  htSrc_       = iConfig.getParameter<edm::InputTag>("HTSource");
  muonVetoSrc_ = iConfig.getParameter<edm::InputTag>("MuonVetoSrc");
  eleVetoSrc_  = iConfig.getParameter<edm::InputTag>("EleVetoSrc");

  saveAllMuons_= iConfig.getParameter<bool>("SaveAllMuons");
  muonSrc_     = iConfig.getParameter<edm::InputTag>("MuonSource");
  electronSrc_ = iConfig.getParameter<edm::InputTag>("ElectronSource");
  minMuPt_     = iConfig.getParameter<double>("MinMuPt");
  maxMuEta_    = iConfig.getParameter<double>("MaxMuEta");
  maxMuD0_     = iConfig.getParameter<double>("MaxMuD0");
  maxMuDz_     = iConfig.getParameter<double>("MaxMuDz");
  maxMuRelIso_ = iConfig.getParameter<double>("MaxMuRelIso");
  pfCandidateSrc_ = iConfig.getParameter<edm::InputTag>("PFCandidateSrc");

   cntPassnJetsCut =0; cntTaggedTopEvents =0; cntTaggedTopEventsWithinMassCuts =0;
   cntTaggedAllCutsPlusCSVS =0;
   cntPassingMTbestTopJetCut =0; cntPassingMTclosebJetCut =0;
   cntPassingMT2Cut =0; cntPassingMT2andMTCut =0; cntPassingMTbestWJetCut =0;
   cntTaggedbestFatJetPlusCSVS =0;

   type3TopTaggerPtr = new topTagger::type3TopTagger();

   type3TopTaggerPtr->setdoExtraCuts(doExtraCuts_);
   type3TopTaggerPtr->setmTbcut(mTbcut_);
   type3TopTaggerPtr->setmTtcut(mTtcut_);
   type3TopTaggerPtr->setMT2cut(MT2cut_);
   type3TopTaggerPtr->setmTWcut(mTWcut_);
   type3TopTaggerPtr->setdoMTMT2cuts(doMTMT2cuts_);
   type3TopTaggerPtr->setPDGmWmTop(mW_, mTop_);
   type3TopTaggerPtr->setlowRatioWoverTop(lowRatioWoverTop_);
   type3TopTaggerPtr->sethighRatioWoverTop(highRatioWoverTop_);
   type3TopTaggerPtr->setlowArcTanm13overm12(lowArcTanm13overm12_);
   type3TopTaggerPtr->sethighArcTanm13overm12(highArcTanm13overm12_);
   type3TopTaggerPtr->setlowm23overm123(lowm23overm123_);
   type3TopTaggerPtr->setdefaultJetCone(defaultJetCone_);
   type3TopTaggerPtr->setsimuCAdeltaR(simuCAdeltaR_);
   type3TopTaggerPtr->setsimuCALargerdeltaR(simuCALargerdeltaR_);
   type3TopTaggerPtr->setlowTopCut(lowTopCut_);
   type3TopTaggerPtr->sethighTopCut(highTopCut_);
   type3TopTaggerPtr->setlowWCut(lowWCut_);
   type3TopTaggerPtr->sethighWCut(highWCut_);
   type3TopTaggerPtr->setCSVS(CSVS_);
   type3TopTaggerPtr->setnSubJetsDiv(nSubJetsDiv_);
   type3TopTaggerPtr->setnJetsSel(nJetsSel_);
   type3TopTaggerPtr->setmaxEtaForbJets(maxEtaForbJets_);
   type3TopTaggerPtr->setorderingOptArr(orderingOptArr_);
   type3TopTaggerPtr->setmaxIndexForOrderingArr(maxIndexForOrderingArr_);
   type3TopTaggerPtr->setdebug(debug_);
   type3TopTaggerPtr->settaggingMode(taggingMode_);
   type3TopTaggerPtr->setdobVetoCS(dobVetoCS_);

   kFoundTwoMuons = false;
   kFoundTwoEles = false;

   varsDoubleVN_ = std::vector<UShort_t>(varsDoubleTagsV_.size(),0);//JL
   for(unsigned int i = 0; i < varsDoubleTagsV_.size(); ++i) //JL
     varsDoubleV_.push_back(new Float_t[120]);//JL
   
}


DiSTopStudyTree::~DiSTopStudyTree() {
  for(unsigned int i = 0; i < varsDoubleTagsV_.size(); ++i)//JL
    delete [] varsDoubleV_.at(i);//JL
}


void DiSTopStudyTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  using namespace edm;

  iSetup.getData( pdt_ );

  clearTreeVectors();

  // fill event ID 
  t_EvtRun   = iEvent.id().run();
  t_EvtLS    = iEvent.luminosityBlock();
  t_EvtEvent = iEvent.id().event();
  
  iEvent.getByLabel("source", product);
   
  // save number of vertices and position of primary vertex
  edm::Handle< std::vector<reco::Vertex> > vertices;
  iEvent.getByLabel(vtxSrc_, vertices);
  if(vertices->size()<1) std::cout << "No vertices are reconstructed - check StdCleaning ?" << std::endl;

  t_NVertices = vertices->size();
  
  reco::Vertex::Point vtxpos = (vertices->size() > 0 ? (*vertices)[0].position() : reco::Vertex::Point());

  // save pileup weight
  /*J
  double pu_event_wt = 1.0;
  edm::Handle<double> puweight;
  if( doPUReWeight_ ) {
    iEvent.getByLabel(puWeigthSrc_, puweight);
    pu_event_wt = *puweight;
  }
  t_PUWeight = pu_event_wt;

  if(debug_)std::cout << "t_NVertices " << t_NVertices << "  t_PUWeight " << t_PUWeight << std::endl;
  J*/
  t_TrueNPV = -1;

  if(isMC_) {
    //save vertex info for PU weight:
    edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
    iEvent.getByLabel( edm::InputTag("addPileupInfo"),PupInfo);
    
    vector<PileupSummaryInfo>::const_iterator PVI;
    
    float Tnpv = -1;
    for(PVI = PupInfo->begin(); PVI !=PupInfo->end(); ++PVI) {
      
      int BX = PVI->getBunchCrossing();
      
      if(BX == 0) {
	Tnpv = PVI->getTrueNumInteractions();
	continue;
      }
    }
    t_TrueNPV = Tnpv;
    
  }
  
  //<JL
  //Save the pdf weights
  //std::cout<<"================"<<std::endl;
  for(unsigned int i = 0; i < varsDoubleTagsV_.size(); ++i) {
    edm::Handle< std::vector<double> > varV;
    iEvent.getByLabel(varsDoubleTagsV_.at(i),varV);
    if( varV.isValid() ) {
      varsDoubleVN_.at(i) = varV->size();
      //std::cout<<"Num="<<varsDoubleVN_.at(i)<<" "<<varV->size()<<std::endl;
      for(unsigned int j = 0; j < varV->size(); ++j) {
	varsDoubleV_.at(i)[j] = varV->at(j);
	//std::cout<< "i: " << i << ", j: " << j << ", weight: " << varsDoubleV_.at(i)[j]<<" "<<varV->at(j)<<std::endl;;
      }
    }
  }
  //JL>

  /* DCH not keeping the puWeights vector, just keep TrueNPV
  //save PU weight:
  if(doPUReWeight_) {
    if(Tnpv < static_cast<int> (puWeights_.size())) {
      t_PUWeight *= puWeights_[Tnpv];
    }
    else {
      cout << "WARNING ... number of PU vertices = " << Tnpv << " out of histogram binning." << endl;
      t_PUWeight = -99;
    }
  }
  else t_PUWeight = -99;
  */

  // save missing transverse energy
  edm::Handle < std::vector<pat::MET> > met;
  iEvent.getByLabel(pfMetSrc_, met);
  //  TLorentzVector *lvMET = new TLorentzVector();
  for (std::vector<pat::MET>::const_iterator it = met->begin(); it != met->end(); ++it) {
    t_PFMetPx = it->px();
    t_PFMetPy = it->py();
    t_PFMetE = it->pt();  
    t_PFMetPhi = it->phi();  
    t_PFMetSignificance =     it->significance();   

    t_PFOrigMetPx = it->px();
    t_PFOrigMetPy = it->py();
    t_PFOrigMetE = it->pt();  
    t_PFOrigMetPhi = it->phi();  

  }
  //  cout<<t_PFMetPx<<" : "<<t_PFMetPy<<" : "<<t_PFMetE<<" : "<<t_PFMetPhi<<" : "<<endl;


  //Build the lorentz vector for MET 
  TLorentzVector lvMET;
  lvMET.SetPx(t_PFMetPx);
  lvMET.SetPy(t_PFMetPy);
  lvMET.SetE(t_PFMetE);
  lvMET.SetPhi(t_PFMetPhi);

  if(isMC_){
    //Handle<std::vector<double> > reweightVector;
    //if(  iEvent.getByLabel("QCDScaleWeightProducer","weight",reweightVector)){
    //for(std::vector<double>::const_iterator dWeight=reweightVector->begin(); dWeight !=reweightVector->end(); ++dWeight){
    //double dTemp = (*dWeight);
    //t_QCDReweight->push_back(dTemp);
    //}
    //}

    // check LHEProduct
    //edm::Handle<LHEEventProduct> product;
    // iEvent.getByLabel("source", product);
    LHEEventProduct::comments_const_iterator c_begin = product->comments_begin();
    LHEEventProduct::comments_const_iterator c_end   = product->comments_end();
    
    /*   for( LHEEventProduct::comments_const_iterator cit=c_begin; cit!=c_end; ++cit) {
	 std::cout << "%%%%%%%%%%%%%%%%%%% "<< *cit << std::endl;
	 }
    */
    
    //get HEPEUP object
    const lhef::HEPEUP hepeup = product->hepeup();
    //const std::vector<lhef::HEPEUP::FiveVector> pup = hepeup.PUP;
    
    int nPartonsLHE = 0;
    double htLHE=0.0, htLepLHE=0.0;
    TLorentzVector mhtLHE   (0.0, 0.0, 0.0, 0.0);
    TLorentzVector mhtLepLHE(0.0, 0.0, 0.0, 0.0);
    TLorentzVector dilepLHE (0.0, 0.0, 0.0, 0.0);
    int            dilepId=-99;
    std::vector<TLorentzVector> lepLHE;
    std::vector<TLorentzVector> partonLHE;
    
    // calculate some global variables based on LHE quantities
    double lheZMass = -1;
    double lheZPt   = -1;
    
    // loop over outgoing particles
    for(int ip = 2; ip < hepeup.NUP; ++ip) {
      
      TLorentzVector vPart(hepeup.PUP[ip][0], hepeup.PUP[ip][1], hepeup.PUP[ip][2], hepeup.PUP[ip][3] );
      
      if( hepeup.IDUP[ip]==23 ) {
	lheZMass = hepeup.PUP[ip][4];
	lheZPt   = vPart.Pt();
      }
      
      // print LHE entries
      if( debug_ ) {
	std::cout << ip << " " << "status " << hepeup.ISTUP[ip]
		  << " PID " << hepeup.IDUP[ip]
		  << " pT " << vPart.Pt()
		  << std::endl;
      }
      
      // check LHE HT, MHT of event using stable particles
      if( hepeup.ISTUP[ip]==1 ) {
	htLepLHE  += vPart.Pt();
	mhtLepLHE -= vPart;
	if( hepeup.IDUP[ip]==21 || (std::abs(hepeup.IDUP[ip])>0 && std::abs(hepeup.IDUP[ip])<7) ) {
	  nPartonsLHE++;
	  htLHE    += vPart.Pt();
	  mhtLHE   -= vPart;
	  partonLHE.push_back(vPart);
	} else if( std::abs(hepeup.IDUP[ip])>10 && std::abs(hepeup.IDUP[ip])<17) {
	  dilepLHE += vPart;
	  lepLHE.push_back(vPart);
	  dilepId = std::abs(hepeup.IDUP[ip]);
	}
      }
    }
    
    //  cout<<"GenLevel HT is : "<<htLHE<<endl;

    t_Genht =  htLHE;
    t_Genmht =  mhtLHE.Pt();

    //========= save gen info
    
    edm::Handle<View<reco::Candidate> > particles;
    iEvent.getByLabel( "genParticles", particles );
    
    for( View<reco::Candidate>::const_iterator p = particles->begin();   p != particles->end(); ++ p ) {
      
      //    const ParticleData * pd = pdt_->particle( p->pdgId() );  
      
      //if(  std::abs(p->pdgId())==13 ){
      if(p->status() == 3) {
	
	t_genPartPdgId  ->push_back(p->pdgId());
	t_genPartStatus ->push_back(p->status());
	t_genPartPt     ->push_back(p->pt());
	t_genPartEta    ->push_back(p->eta());
	t_genPartPhi    ->push_back(p->phi());
	t_genPartE      ->push_back(p->energy());
	
      }//if a status=3 W or b
      
      //delete pd;
    } //loop over gen particles
    
    //================
  }

  if(debug_ && isMC_){
    edm::Handle<View<reco::Candidate> > particles;
    iEvent.getByLabel( "genParticles", particles );
    
    cout << "Level 3 particles: ";
    
    //DCH cout pdg here:
    for( View<reco::Candidate>::const_iterator p2 = particles->begin();   p2 != particles->end(); ++ p2 ) {
      
      if(p2->status() == 3) {
	
	cout << " " << p2->pdgId() << ", ";
	
	//t_genPartPdgId  ->push_back(p->pdgId());
	//t_genPartStatus ->push_back(p->status());
	//t_genPartPt     ->push_back(p->pt());
	//t_genPartEta    ->push_back(p->eta());
	//t_genPartPhi    ->push_back(p->phi());
	//t_genPartE      ->push_back(p->energy());
	
      }//if a status=3 W or b
    }
    cout <<endl;

    edm::Handle<edm::View<pat::Jet> > jetsAllCount;
    iEvent.getByLabel(jetAllsrc_, jetsAllCount);
    
    for(unsigned int ijet=0; ijet<jetsAllCount->size(); ijet++) {
      cout << "DeltaR jet-mu1: " << 
	deltaR((*jetsAllCount)[ijet].eta(), (*jetsAllCount)[ijet].phi(),lvMuon1.Eta(),lvMuon1.Phi());
      
      cout << ", jet-mu2: " << 	deltaR((*jetsAllCount)[ijet].eta(), (*jetsAllCount)[ijet].phi(),lvMuon2.Eta(),lvMuon2.Phi());
      cout << endl;
      
    }
    
  }
  
  edm::Handle<std::vector<reco::GsfElectron> > electrons;
  iEvent.getByLabel(electronSrc_, electrons);
  
  edm::Handle<std::vector<pat::Muon> > muons;
  iEvent.getByLabel(muonSrc_, muons);
    
  //store all lepton info
  if (vertices->size() > 0) {
    for (std::vector<pat::Muon>::const_iterator m = muons->begin(); m != muons->end(); ++m) {
      t_PFLepPt             ->push_back(m->pt() );
      t_PFLepEta            ->push_back(m->eta());    
      t_PFLepPhi            ->push_back(m->phi());
      t_PFLepE              ->push_back(m->energy());
      t_PFLepCh             ->push_back(m->charge());
      t_PFLepIsMu           ->push_back(true);
      t_PFLepIsEle          ->push_back(false);
    } //loop over muons
    
    for (std::vector<reco::GsfElectron>::const_iterator e = electrons->begin(); e != electrons->end(); ++e) {
      t_PFLepPt             ->push_back(e->pt() );
      t_PFLepEta            ->push_back(e->eta());    
      t_PFLepPhi            ->push_back(e->phi());
      t_PFLepE              ->push_back(e->energy());
      t_PFLepCh             ->push_back(e->charge());
      t_PFLepIsMu           ->push_back(false);
      t_PFLepIsEle          ->push_back(true);
    } // loop over eles
    
  } // atleast one vertex, just a protection   
  
  if(isDiLeptonSample_) {
    if((*t_PFLepPt).size() < 2) return;
  }
  else{ 
    if((*t_PFLepPt).size() > 0) return;
  }

  //require at least 2 leptons
  if((*t_PFLepPt).size()>=2) {
    
    //remove these, check OS in the macro (we can then keep SS events or MuMu, EleEle, and MuEle all at once
    //if((*t_PFLepIsMu)[0] && (*t_PFLepIsEle)[1]) {
    //if((*t_PFLepCh)[0]!=(*t_PFLepCh)[1]){
    
    //Find the first Lepton
    lvLep1.SetPtEtaPhiE((*t_PFLepPt)[0],(*t_PFLepEta)[0],(*t_PFLepPhi)[0],(*t_PFLepE)[0]);
    
    //Find the second lepton
    lvLep2.SetPtEtaPhiE((*t_PFLepPt)[1],(*t_PFLepEta)[1],(*t_PFLepPhi)[1],(*t_PFLepE)[1]);
    
    TLorentzVector lvZ =  lvLep1 + lvLep2;
    
    t_PFLepZPt    = lvZ.Pt();
    t_PFLepZEta   = lvZ.Eta();
    t_PFLepZPhi   = lvZ.Phi();
    t_PFLepZE     = lvZ.E();
    t_PFLepZMass  = lvZ.M();
    
    double metPx = lvMET.Px() + lvLep1.Px() + lvLep2.Px();
    double metPy = lvMET.Py() + lvLep1.Py() + lvLep2.Py();
    
    //lvMET = lvMET +  lvMuon1 + lvMuon2;
    
    TLorentzVector newMet(0.,0.,0.,0.);
    newMet.SetPxPyPzE(metPx,metPy,0.0, sqrt(pow(metPx,2)+pow(metPy,2)));
    lvMET = newMet;
    
    t_PFMetPx = lvMET.Px();
    t_PFMetPy = lvMET.Py();
    t_PFMetE = lvMET.Pt();  
    t_PFMetPhi = lvMET.Phi();  
  }  
  
  // save ra2 ht/mht mainly for debugging only
  edm::Handle<double> ht;
  iEvent.getByLabel(htSrc_, ht);
  t_PFht  = *ht;
  
  edm::Handle<edm::View<reco::MET> > mht;
  iEvent.getByLabel(mhtSrc_, mht);
  t_PFmht = (*mht)[0].pt() ;

  // std::cout << "ht " << t_PFht << " mht " << t_PFmht << std::endl;

  // save all jets along with btag information
  edm::Handle<edm::View<pat::Jet> > jetsAll;
  iEvent.getByLabel(jetAllsrc_, jetsAll);

  t_NJetsPt30Eta2p5 = 0;
  t_NJetsPt30Eta5p0 = 0; 
  t_NJetsPt50Eta2p5 = 0;
  t_NJetsPt50Eta5p0 = 0;

  t_NJetsPt30Eta2p4 = 0;
  t_NJetsPt50Eta2p4 = 0;
  t_NJetsPt70Eta2p4 = 0;

  iNJets = 0;
  Int_t iBTagTight = 0;
  Int_t iBTagLoose = 0;
 
  for(unsigned int ijet=0; ijet<jetsAll->size(); ijet++) {

    //check if there are 6 jets with pt > 30 
    if( (*jetsAll)[ijet].pt() > 30.0 ) {
   
      t_PFJetPt  ->push_back((*jetsAll)[ijet].pt() );
      t_PFJetEta ->push_back((*jetsAll)[ijet].eta());    
      t_PFJetPhi ->push_back((*jetsAll)[ijet].phi());
      t_PFJetE   ->push_back((*jetsAll)[ijet].energy()  );
      t_PFJetBTag->push_back((*jetsAll)[ijet].bDiscriminator(btagname_.c_str()) );
      
      if( (*jetsAll)[ijet].pt()>30.0 && std::abs((*jetsAll)[ijet].eta())<5.0 )
	t_NJetsPt30Eta5p0++;
      if( (*jetsAll)[ijet].pt()>30.0 && std::abs((*jetsAll)[ijet].eta())<2.5 )
	t_NJetsPt30Eta2p5++;
      if( (*jetsAll)[ijet].pt()>50.0 && std::abs((*jetsAll)[ijet].eta())<5.0 )
	t_NJetsPt50Eta5p0++;
      if( (*jetsAll)[ijet].pt()>50.0 && std::abs((*jetsAll)[ijet].eta())<2.5 )
	t_NJetsPt50Eta2p5++;
      
      if( (*jetsAll)[ijet].pt()>30.0 && std::abs((*jetsAll)[ijet].eta())<2.4 )
	t_NJetsPt30Eta2p4++;
      if( (*jetsAll)[ijet].pt()>50.0 && std::abs((*jetsAll)[ijet].eta())<2.4 )
	t_NJetsPt50Eta2p4++;
      if( (*jetsAll)[ijet].pt()>70.0 && std::abs((*jetsAll)[ijet].eta())<2.4 )
	t_NJetsPt70Eta2p4++;
    }
    
  }
  
  // save no. of isolated/identified muons for veto
  edm::Handle<edm::View<reco::Muon> > muonVeto;
  iEvent.getByLabel(muonVetoSrc_, muonVeto);
  t_NVetoMuon = muonVeto->size();
  
  // save no. of isolated/identified electrons for veto
  edm::Handle<std::vector<reco::GsfElectron> > eleVeto;
  iEvent.getByLabel(eleVetoSrc_, eleVeto);
  t_NVetoEle  = eleVeto->size();
  
  edm::Handle<std::vector<reco::PFCandidate> > pfCandidates;
  iEvent.getByLabel(pfCandidateSrc_, pfCandidates);
  const reco::PFCandidateCollection & pfCands = *pfCandidates;
  
 
  // **************************** Code for Type 3 Tagger *********************************************************//
  loadEventInfo(iEvent, iSetup);
  
  loadRecoJets(iEvent);
  //   if( !isPatJet ){ std::cout<<"Not a pat::jet input! Need it for bTag information..."<<std::endl; return false; }
  
  double evtWeight = 1.0;
  if( evtWeight_.isValid() ) evtWeight = (*evtWeight_);
  
  cntPassnJetsCut += evtWeight;
  
  // The index has to be the same between jets & btags!
  vector<TLorentzVector> oriJetsVec; vector<double> recoJetsBtagCSVS;
  iNJets = pickJetsForCombining(oriJetsVec, recoJetsBtagCSVS);
  
  bool pass = type3TopTaggerPtr->processEvent(oriJetsVec, recoJetsBtagCSVS, lvMET);
  
  bestTopJetMass = type3TopTaggerPtr->bestTopJetMass;
  remainPassCSVS = type3TopTaggerPtr->remainPassCSVS ;
  mTbestTopJet = type3TopTaggerPtr->mTbestTopJet;
  mTbJet = type3TopTaggerPtr->mTbJet;
  MT2 = type3TopTaggerPtr->MT2;
  mTbestWJet  = type3TopTaggerPtr->mTbestWJet;
  mTbestbJet = type3TopTaggerPtr->mTbestbJet;
  
  /*   if( type3TopTaggerPtr->bestTopJetMass > 80 && type3TopTaggerPtr->bestTopJetMass < 270 )
       cout<<type3TopTaggerPtr->bestTopJetMass<<endl;
       else
       return;
  */
  
  
  if( type3TopTaggerPtr->isTopEvent ) cntTaggedTopEvents += evtWeight;
  /*
  // Must have a fat top jet
  if( !taggingMode_ && type3TopTaggerPtr->bestTopJetIdx == -1 ) return false;
  
  // Some mass cuts on the fat top jet
  if( type3TopTaggerPtr->bestTopJetMass > lowTopCut_ && type3TopTaggerPtr->bestTopJetMass < highTopCut_) cntTaggedTopEventsWithinMassCuts += evtWeight;
  if( !taggingMode_ && !(type3TopTaggerPtr->bestTopJetMass > lowTopCut_ && type3TopTaggerPtr->bestTopJetMass < highTopCut_) ) return false;
  
  if( type3TopTaggerPtr->remainPassCSVS ) cntTaggedAllCutsPlusCSVS += evtWeight;
  if( type3TopTaggerPtr->bestPassCSVS ) cntTaggedbestFatJetPlusCSVS += evtWeight;
  // Must have a jet b-tagged!
  if( !taggingMode_ && !type3TopTaggerPtr->remainPassCSVS ) return false;
  
  if( !taggingMode_ && type3TopTaggerPtr->pickedRemainingCombfatJetIdx == -1 && oriJetsVec.size()>=6 ) return false; 
  
  if( type3TopTaggerPtr->mTbestTopJet > mTtcut_ ) cntPassingMTbestTopJetCut += evtWeight;
  if( (type3TopTaggerPtr->mTbJet + 0.5*type3TopTaggerPtr->mTbestTopJet) > mTbcut_ ) cntPassingMTclosebJetCut += evtWeight;
  
  if( type3TopTaggerPtr->MT2 > MT2cut_ ) cntPassingMT2Cut += evtWeight;
  
  if( (type3TopTaggerPtr->mTbestbJet + type3TopTaggerPtr->mTbestWJet) > mTWcut_ ) cntPassingMTbestWJetCut += evtWeight;
  
  if( type3TopTaggerPtr->MT2 > MT2cut_ && (type3TopTaggerPtr->mTbJet + 0.5*type3TopTaggerPtr->mTbestTopJet) > mTbcut_ ){
  cntPassingMT2andMTCut += evtWeight;
  }
  
  if( !taggingMode_ && doMTMT2cuts_ && !(type3TopTaggerPtr->MT2 > MT2cut_ && (type3TopTaggerPtr->mTbJet + 0.5*type3TopTaggerPtr->mTbestTopJet) > mTbcut_) ) return false;
  
  return true;
  */
  
  tree->Fill();
  
}

// ------------ method called once each job just before starting event loop  ------------
void DiSTopStudyTree::beginJob() {
  
  BookHistograms();
  
  puWeights_.clear();
  
  /* DCH
    //REMOVE THIS:
    //Saving TrueNPV from MC block.  Can then use this reweighting scheme on the fly:
  if(doPUReWeight) {

    TFile* puDataFile = new TFile("/uscms_data/d3/dhare/CMSSW_5_3_5/src/DiSTopStudy/DataPileupHistogram_RA2Summer12_190456-208686_ABCD.root");
    TH1D* data_npu_estimated = (TH1D*) puDataFile->Get("pileup");
    
    // Distribution used for Summer2012 MC.
    Double_t Summer2012[60] = {
      2.344E-05,
      2.344E-05,
      2.344E-05,
      2.344E-05,
      4.687E-04,
      4.687E-04,
      7.032E-04,
      9.414E-04,
      1.234E-03,
      1.603E-03,
      2.464E-03,
      3.250E-03,
      5.021E-03,
      6.644E-03,
      8.502E-03,
      1.121E-02,
      1.518E-02,
      2.033E-02,
      2.608E-02,
      3.171E-02,
      3.667E-02,
      4.060E-02,
      4.338E-02,
      4.520E-02,
      4.641E-02,
      4.735E-02,
      4.816E-02,
      4.881E-02,
      4.917E-02,
      4.909E-02,
      4.842E-02,
      4.707E-02,
      4.501E-02,
      4.228E-02,
      3.896E-02,
      3.521E-02,
      3.118E-02,
      2.702E-02,
      2.287E-02,
      1.885E-02,
      1.508E-02,
      1.166E-02,
      8.673E-03,
      6.190E-03,
      4.222E-03,
      2.746E-03,
      1.698E-03,
      9.971E-04,
      5.549E-04,
      2.924E-04,
      1.457E-04,
      6.864E-05,
      3.054E-05,
      1.282E-05,
      5.081E-06,
      1.898E-06,
      6.688E-07,
      2.221E-07,
      6.947E-08,
      2.047E-08
    };
    
    puWeights_.resize(60);
    double s = 0.0;
    for (int npu = 0; npu < 60; ++npu) {
      double npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));
      puWeights_[npu] = npu_estimated/ Summer2012[npu];
      s += npu_estimated;
     
    }
   
    //normalized weights such that the total sum of weights over whole sample is 1.0
    for(int npu = 0; npu < 60; ++npu) {
      puWeights_[npu] /= s;
    } 
    
  }
  */
  
}

// ------------ method called once each job just after ending the event loop  ------------
void DiSTopStudyTree::endJob() {

}

void DiSTopStudyTree::BookHistograms() {

  // book tree here
  tree = fs->make<TTree>("SAT", "tree");
  tree->SetAutoSave(10000);

  tree->Branch("t_EvtRun",    &t_EvtRun,   "t_EvtRun/i");
  tree->Branch("t_EvtLS",     &t_EvtLS,    "t_EvtLS/i");
  tree->Branch("t_EvtEvent",  &t_EvtEvent, "t_EvtEvent/i");
  tree->Branch("t_NVertices", &t_NVertices,"t_NVertices/I");
  tree->Branch("t_PUWeight",  &t_PUWeight, "t_PUWeight/D");
  tree->Branch("t_TrueNPV",   &t_TrueNPV,  "t_TrueNPV/D");
  tree->Branch("t_PFMetPx",   &t_PFMetPx,  "t_PFMetPx/D");
  tree->Branch("t_PFMetPy",   &t_PFMetPy,  "t_PFMetPy/D");
  tree->Branch("t_PFMetE",   &t_PFMetE,  "t_PFMetE/D");
  tree->Branch("t_PFMetPhi",   &t_PFMetPhi,  "t_PFMetPhi/D");
  tree->Branch("t_PFMetSignificance",   &t_PFMetSignificance,  "t_PFMetSignificance/D");

  tree->Branch("t_PFOrigMetPx",   &t_PFOrigMetPx,  "t_PFOrigMetPx/D");
  tree->Branch("t_PFOrigMetPy",   &t_PFOrigMetPy,  "t_PFOrigMetPy/D");
  tree->Branch("t_PFOrigMetE",    &t_PFOrigMetE,   "t_PFOrigMetE/D");
  tree->Branch("t_PFOrigMetPhi",  &t_PFOrigMetPhi, "t_PFOrigMetPhi/D");

  t_PFJetPt  = new std::vector<double>();
  t_PFJetEta = new std::vector<double>();    
  t_PFJetPhi = new std::vector<double>();
  t_PFJetE   = new std::vector<double>();
  t_PFJetBTag= new std::vector<double>();
  tree->Branch("NJets",   &iNJets,  "NJets/I");
  tree->Branch("t_PFJetPt",  "vector<double>", &t_PFJetPt );
  tree->Branch("t_PFJetEta", "vector<double>", &t_PFJetEta);
  tree->Branch("t_PFJetPhi", "vector<double>", &t_PFJetPhi);
  tree->Branch("t_PFJetE",   "vector<double>", &t_PFJetE);
  tree->Branch("t_PFJetBTag","vector<double>", &t_PFJetBTag);
  tree->Branch("t_NJetsPt30Eta2p5", &t_NJetsPt30Eta2p5, "t_NJetsPt30Eta2p5/I");
  tree->Branch("t_NJetsPt30Eta5p0", &t_NJetsPt30Eta5p0, "t_NJetsPt30Eta5p0/I"); 
  tree->Branch("t_NJetsPt50Eta2p5", &t_NJetsPt50Eta2p5, "t_NJetsPt50Eta2p5/I"); 
  tree->Branch("t_NJetsPt50Eta5p0", &t_NJetsPt50Eta5p0, "t_NJetsPt50Eta5p0/I");
  tree->Branch("t_NJetsPt30Eta2p4", &t_NJetsPt30Eta2p4, "t_NJetsPt30Eta2p4/I");
  tree->Branch("t_NJetsPt50Eta2p4", &t_NJetsPt50Eta2p4, "t_NJetsPt50Eta2p4/I");
  tree->Branch("t_NJetsPt70Eta2p4", &t_NJetsPt70Eta2p4, "t_NJetsPt70Eta2p4/I");

  tree->Branch("t_PFht",      &t_PFht,      "t_PFht/D");
  tree->Branch("t_PFmht",     &t_PFmht,     "t_PFmht/D");
  tree->Branch("t_PFmphi",     &t_PFmphi,     "t_PFmphi/D");
  
  tree->Branch("t_NVetoMuon", &t_NVetoMuon,"t_NVetoMuon/I");
  tree->Branch("t_NVetoEle",  &t_NVetoEle, "t_NVetoEle/I");
 
  //t_QCDReweight = new std::vector<double>();
  //tree->Branch("t_QCDReweight", "vector<double>", &t_QCDReweight);
    
  t_PFLepPt  = new std::vector<double>();
  t_PFLepEta = new std::vector<double>();    
  t_PFLepPhi = new std::vector<double>();
  t_PFLepE   = new std::vector<double>();
  t_PFLepCh  = new std::vector<double>();
  tree->Branch("t_PFLepPt",  "vector<double>", &t_PFLepPt );
  tree->Branch("t_PFLepEta", "vector<double>", &t_PFLepEta);
  tree->Branch("t_PFLepPhi", "vector<double>", &t_PFLepPhi);
  tree->Branch("t_PFLepE",   "vector<double>", &t_PFLepE  );
  tree->Branch("t_PFLepCh",  "vector<double>", &t_PFLepCh);
    
  t_PFLepIsMu     = new std::vector<bool>();
  t_PFLepIsEle    = new std::vector<bool>();
  tree->Branch("t_PFLepIsMu",   "vector<bool>",   &t_PFLepIsMu  );
  tree->Branch("t_PFLepIsEle",   "vector<bool>",   &t_PFLepIsEle  );
  
  tree->Branch("t_DiLepZPt"  , &t_PFLepZPt  , "t_DiLepZPt/D");
  tree->Branch("t_DiLepZEta" , &t_PFLepZEta , "t_DiLepZEta/D");
  tree->Branch("t_DiLepZPhi" , &t_PFLepZPhi , "t_DiLepZPhi/D");
  tree->Branch("t_DiLepZE"   , &t_PFLepZE   , "t_DiLepZE/D");
  tree->Branch("t_DiLepZMass", &t_PFLepZMass, "t_DiLepZMass/D");
  
  if(isMC_){
    t_genPartPdgId     = new std::vector<int>();
    t_genPartStatus    = new std::vector<int>();
    t_genPartPt        = new std::vector<double>();
    t_genPartEta       = new std::vector<double>();
    t_genPartPhi       = new std::vector<double>();
    t_genPartE         = new std::vector<double>();
    t_genPartDecayMode = new std::vector<int>();
    tree->Branch("t_genPartPdgId",     "vector<int>",    &t_genPartPdgId    );
    tree->Branch("t_genPartStatus",    "vector<int>",    &t_genPartStatus   );
    tree->Branch("t_genPartPt",        "vector<double>", &t_genPartPt       );
    tree->Branch("t_genPartEta",       "vector<double>", &t_genPartEta      );
    tree->Branch("t_genPartPhi",       "vector<double>", &t_genPartPhi      );
    tree->Branch("t_genPartE",         "vector<double>", &t_genPartE        );
    tree->Branch("t_genPartDecayMode", "vector<int>",              &t_genPartDecayMode);


    tree->Branch("t_Genht",      &t_Genht,      "t_Genht/D");
    tree->Branch("t_Genmht",     &t_Genmht,     "t_Genmht/D");
  }

  //  tree->Branch("LHEEventProduct", "LHEEventProduct", &product );

  bestTopJetMass = -1;
  remainPassCSVS =-10;
  mTbestTopJet = -1;
  mTbJet = -1;
  MT2 =-1;
  mTbestWJet = -1;
  mTbestbJet = -1;
  
  tree->Branch("t_bestTopJetMass", &bestTopJetMass, "t_bestTopJetMass/D");
  tree->Branch("t_remainPassCSVS", &remainPassCSVS, "t_remainPassCSVS/B");
  tree->Branch("t_mTbestTopJet", &mTbestTopJet, "t_mTbestTopJet/D");
  tree->Branch("t_mTbJet", &mTbJet, "t_mTbJet/D");
  tree->Branch("t_MT2", &MT2, "t_MT2/D");
  tree->Branch("t_mTbestWJet", &mTbestWJet, "t_mTbestWJet/D");
  tree->Branch("t_mTbestbJet", &mTbestbJet, "t_mTbestbJet/D");

  //<JL
  for(unsigned int i = 0; i < varsDoubleV_.size(); ++i) {
    std::string name = varsDoubleTagsV_.at(i).label();
    if( i < varsDoubleNamesInTreeV_.size() ) {
      name = varsDoubleNamesInTreeV_.at(i);
    }
    tree->Branch((name+"Num").c_str(),&(varsDoubleVN_.at(i)),(name+"Num/s").c_str());
    tree->Branch(name.c_str(),varsDoubleV_.at(i),(name+"["+name+"Num]/F").c_str());

  }
  //JL>

}

void DiSTopStudyTree::clearTreeVectors() {

  t_PFJetPt  ->clear();
  t_PFJetEta ->clear();    
  t_PFJetPhi ->clear();
  t_PFJetE   ->clear();
  t_PFJetBTag->clear();

  //t_QCDReweight->clear();

  /*
    if(isDimuonSample_){
    t_PFMuonPt  ->clear();
    t_PFMuonEta ->clear();    
    t_PFMuonPhi ->clear();
    t_PFMuonE   ->clear();
    t_PFMuonCh  ->clear();
    t_PFMuonChiSq                ->clear(); 
    t_PFMuonValidMuonHits        ->clear(); 
    t_PFMuonMatchedStations      ->clear(); 
    t_PFMuonValidPixelHits       ->clear(); 
    t_PFMuonTrackerLayerMeasured ->clear();
    t_PFMuonCharHadEt            ->clear(); 
    t_PFMuonNeutHadEt            ->clear(); 
    t_PFMuonPhotEt               ->clear(); 
    t_PFMuonSumPUPt              ->clear(); 
    t_PFMuonRelIso               ->clear();
    t_PFMuonCharHadEt03          ->clear(); 
    t_PFMuonNeutHadEt03          ->clear(); 
    t_PFMuonPhotEt03             ->clear(); 
    t_PFMuonSumPUPt03            ->clear(); 
    t_PFMuonRelIso03             ->clear();
    t_PFMuonID       ->clear();
    t_PFMuonVtxAss   ->clear();
    t_PFMuonIsolated ->clear();
    t_PFMuonIDVtxIso ->clear();
    t_PFMuonIsWdau   ->clear();
    t_PFMuonIsBdau   ->clear();
    t_PFMuonGenMuPt  ->clear();
    t_PFMuonGenMuEta ->clear();
    t_PFMuonGenMuPhi ->clear();
    
    t_PFMuonDirIso01->clear();
    t_PFMuonDirIso02->clear();
    t_PFMuonDirIso03->clear();
    t_PFMuonDirIso04->clear();
    t_PFMuonDirIso05->clear();

    t_PFMuonZPt=-1; 
    t_PFMuonZEta=999; 
    t_PFMuonZPhi=999; 
    t_PFMuonZE=-1;
    t_PFMuonZMass=-1;
  
    }
  */

  /*
    if(isDieleSample_){
    t_PFElecPt  ->clear();
    t_PFElecEta ->clear();    
    t_PFElecPhi ->clear();
    t_PFElecE   ->clear();
    t_PFElecCh  ->clear();
    
    t_PFEleIsEB        ->clear();
    t_PFEleIsEE        ->clear();
    t_PFEleDEtaIn      ->clear();
    t_PFEleDPhiIn      ->clear();
    t_PFEleSigIEtaIEta ->clear();
    t_PFEleHOE         ->clear();
    t_PFEleD0Vtx       ->clear();
    t_PFEleZVtx        ->clear();
    t_PFEleIsoCH       ->clear();
    t_PFEleIsoEM       ->clear();
    t_PFEleIsoNH       ->clear();
    t_PFEleIsoRho      ->clear();
    t_PFEleIsoEffA     ->clear();
    t_PFEleRelIso      ->clear();
    t_PFEleID          ->clear();
    t_PFEleIsolated    ->clear();
    t_PFEleIDIso       ->clear();
    
    t_PFEleDirIso01->clear();
    t_PFEleDirIso02->clear();
    t_PFEleDirIso03->clear();
    t_PFEleDirIso04->clear();
    t_PFEleDirIso05->clear();
    
    t_PFEleIsWdau   ->clear();
    t_PFEleIsBdau   ->clear();
    t_PFEleGenElePt  ->clear();
    t_PFEleGenEleEta ->clear();
    t_PFEleGenElePhi ->clear();

    t_PFElecZPt=-1; 
    t_PFElecZEta=999; 
    t_PFElecZPhi=999; 
    t_PFElecZE=-1;
    t_PFElecZMass=-1;
    }
  */

  t_PFLepPt  ->clear();
  t_PFLepEta ->clear();    
  t_PFLepPhi ->clear();
  t_PFLepE   ->clear();
  t_PFLepCh  ->clear();
  
  t_PFLepIsMu   ->clear();
  t_PFLepIsEle  ->clear();
  
  t_PFLepZPt=-1; 
  t_PFLepZEta=999; 
  t_PFLepZPhi=999; 
  t_PFLepZE=-1;
  t_PFLepZMass=-1;

  if(isMC_){
    t_genPartPdgId     ->clear();
    t_genPartStatus    ->clear();
    t_genPartPt        ->clear();
    t_genPartEta       ->clear();
    t_genPartPhi       ->clear();
    t_genPartE         ->clear();
    t_genPartDecayMode ->clear();
    t_Genht = -1 ;
    t_Genmht = -1;
  }

  
  //J..
  //t_fatJetSubJetPt ->clear();
  //t_fatJetSubJetEta->clear();
  //t_fatJetSubJetPhi->clear();
  //t_fatJetSubJetE  ->clear();
  //..J

  bestTopJetMass = -1;
  remainPassCSVS =-10;
  mTbestTopJet = -1;
  mTbJet = -1;
  MT2 =-1;
  mTbestWJet = -1;
  mTbestbJet = -1;

 for(unsigned int i = 0; i < varsDoubleV_.size(); ++i) {
    varsDoubleVN_.at(i) = 0.;
    for(unsigned int j = 0; j < 120; ++j) {
      varsDoubleV_.at(i)[j] = -9999.;
    }
  }
  //JL>
}




int DiSTopStudyTree::pickJetsForCombining(std::vector<TLorentzVector>& oriJetsVec, std::vector<double> &recoJetsBtagCSVS){

   oriJetsVec.clear(); recoJetsBtagCSVS.clear();

//   reco::Vertex::Point vtxpos = (vertices->size() > 0 ? (*vertices)[0].position() : reco::Vertex::Point());


   int cntNJets =0;

   int cntBJets =0;
   double topCSV = 0;

   for(size ij=0; ij<nJets; ij++){
     const reco::Jet & jet = (*jets)[ij];
     if( jet.pt() < 30.0 ) continue;
     if( !pfJetCutForJetCombining_(jet) ) continue;

     //if(isDimuonSample_&&kFoundTwoMuons){
     //if(deltaR(jet.eta(), jet.phi(),lvMuon1.Eta(),lvMuon1.Phi())<0.3) {
     //  cout << "WARNING: Found a muon matched to jet, warningB" << endl;
     // continue;
     //}
     //if(deltaR(jet.eta(), jet.phi(),lvMuon2.Eta(),lvMuon2.Phi())<0.3) {
     // cout << "WARNING: Found a muon matched to jet, warningB" << endl;
     // continue;
     //}
     //}

     TLorentzVector jetLVec;
     jetLVec.SetPtEtaPhiE(jet.pt(), jet.eta(), jet.phi(), jet.energy());
     oriJetsVec.push_back(jetLVec);
     const pat::Jet & patjet = (*patjets)[ij];
     double btag = patjet.bDiscriminator(bTagKeyString_.c_str());
     recoJetsBtagCSVS.push_back(btag);
     //if(btag > 0.679 && jet.eta() < 2.4) cntBJets++;
     if(btag > topCSV) topCSV = btag;
     
     cntNJets ++;
   }
   
   if(topCSV >= 0 && topCSV <= 0.679) {
   //DCH, hack to set highest csv jet as the b
     for(int q= 0; q < cntNJets; q++ ) {
       if(recoJetsBtagCSVS[q] == topCSV) recoJetsBtagCSVS[q] = 1;
       else recoJetsBtagCSVS[q] = 0;

     }
   }

   return cntNJets;
}

int DiSTopStudyTree::countNJets(){
   int cntNJets =0;
   for(size ij=0; ij<nJets; ij++){
      const reco::Jet & jet = (*jets)[ij];
      if( !pfJetCutForJetCounting_(jet) ) continue;
      if( jet.pt() < 30.0 ) continue;
      cntNJets ++;
   }
   return cntNJets;
}

void DiSTopStudyTree::loadRecoJets(const edm::Event& iEvent){
   iEvent.getByLabel(jetSrc_, jets); nJets = jets->size();

   isPatJet = false;
   edm::ProductID jetProdID = jets.id();
   const edm::Provenance & jetProv = iEvent.getProvenance(jetProdID);
   const std::string jetclassName = jetProv.className();
   TString jetclassNameT(jetclassName);
   if( jetclassNameT.Contains("pat::Jet") ) isPatJet = true;

   if( isPatJet ) iEvent.getByLabel(jetSrc_, patjets);
}

void DiSTopStudyTree::loadMETMHT(const edm::Event& iEvent){
   iEvent.getByLabel(metSrc_, metHandle);
   met = (*metHandle)[0].pt(); metphi = (*metHandle)[0].phi();
}

void DiSTopStudyTree::loadEventInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup){
// Determine if it's data
   if( !iEvent.isRealData() ) isData = false;

// Get run, event, lumi info
   run = iEvent.id().run();
   event = iEvent.id().event();
   lumi = iEvent.luminosityBlock();

// Get vertices
   iEvent.getByLabel(vtxSrc_, vertices); vtxSize = vertices->size();

// Get event weight
   iEvent.getByLabel(evtWeightInput_, evtWeight_);
}



#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(DiSTopStudyTree);

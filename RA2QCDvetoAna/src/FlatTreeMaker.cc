#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/GenMET.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "DataFormats/Common/interface/Ptr.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TRandom3.h"

typedef unsigned int size;
static const int sgnfnDof = 2;

using namespace edm;
using namespace std;

class FlatTreeMaker : public edm::EDFilter{

  public:

    explicit FlatTreeMaker(const edm::ParameterSet & iConfig);
    ~FlatTreeMaker();
    virtual bool filter(edm::Event&, const edm::EventSetup&);
    virtual void beginJob();
    virtual void endJob();
    virtual bool beginRun(edm::Run&, const edm::EventSetup&);
    virtual bool endRun(edm::Run&, const edm::EventSetup&);

	static const UInt_t Njet_max = 50;

    size run, event, ls; bool isdata;
    edm::InputTag vtxSrc_;
    edm::Handle<edm::View<reco::Vertex> > vertices;
    size nVtxPUcut_, vtxSize;
//    edm::InputTag evtWeightInput_;
 //   edm::Handle<double> evtWeight_;
    void loadEventInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup);

    edm::InputTag jetSrc_;
    edm::Handle<edm::View<reco::Jet> > jets;
    size nJets;
    reco::Jet jet1, jet2, jet3, jetother;
    virtual void loadRecoJets(const edm::Event& iEvent);

    edm::InputTag muonSrc_;
    edm::Handle<edm::View<reco::Muon > > muons;
    edm::InputTag eleSrc_;
    edm::Handle<edm::View<pat::Electron > > electrons;
    edm::InputTag tauSrc_;
    edm::Handle<edm::View<pat::Tau > > taus;
    size nMuons, nElectrons, nTaus;
    edm::InputTag forVetoMuonSrc_;
    edm::Handle<edm::View<reco::Muon > > forVetoMuons;
    edm::InputTag forVetoElectronSrc_;
    edm::Handle<edm::View<pat::Electron > > forVetoElectrons;
    size nMuonsForVeto, nElectronsForVeto;
    reco::Muon muon1, muon2; int muon1Charge, muon2Charge;
    pat::Electron ele1, ele2; int ele1Charge, ele2Charge;
    virtual void loadLeptons(const edm::Event& iEvent);

    edm::InputTag photonSrc_;
    edm::Handle<edm::View<pat::Photon> > photons;
    size nPhotons;
    virtual void loadPhotons(const edm::Event& iEvent);

    edm::InputTag metSrc_;
    edm::Handle<edm::View<reco::MET> > met;
    edm::InputTag mhtSrc_;
    edm::Handle<edm::View<reco::MET> > mht;
    virtual void loadMETMHT(const edm::Event& iEvent);

    edm::InputTag htSrc_;
    edm::Handle<double> ht;
    virtual void loadHT(const edm::Event& iEvent);

    double pthat, scalePDF;
    edm::InputTag genJetSrc_;
    edm::Handle<edm::View<reco::GenJet > > genJets;
    size nGenJets;
    reco::GenJet genJet1, genJet2, genJet3, genJetother;
    edm::InputTag genParticleSrc_;
    edm::Handle<edm::View<reco::GenParticle > > genParticles;
    edm::InputTag genMETSrc_;
    edm::Handle<edm::View<reco::GenMET > > genMET;
    Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    int npv; double avg_npv;
    virtual void loadGenInfo(const edm::Event& iEvent);

    edm::InputTag dPhis_CUT_vec_Src_;
    edm::Handle<std::vector<double> > dPhis_CUT_vec_;
    edm::InputTag nJets_CUT_Src_;
    edm::Handle<int> nJets_CUT_;
//    edm::InputTag externalBitToTree_Src_;
//    edm::Handle<int> externalBitToTree_;
//    virtual void loadAUX(const edm::Event& iEvent);

// Indices ==> 0 : for muon  1 : for electron 
// Combination is the following:
// 00 : veto all muons and electrons ==> standard RA2
// 10 : muon control sample ==> for lost-lepton, hadronic tau, W+Jets
// 20 : muon control sample ==> for Z+Jets
// others : currently NOT used
    std::vector<int> nLeptonsSels_;

    bool debug_;

    bool doFillTree_;

    bool isData;

    bool printOnce_;
    vector<std::string> availJECLevels_;

    std::string outRootName_;
    TFile *outRootFile;
    TTree *outTree;
    double evtWeight_TR;
    double mht_TR, ht_TR, met_TR, mt_TR, meff;
    double mhtphi_TR, metphi_TR;
    double metSgnf_TR, metSgnfProb_TR;
    double jet1pt_TR, jet1eta_TR, jet1phi_TR;
    double jet2pt_TR, jet2eta_TR, jet2phi_TR;
    double jet3pt_TR, jet3eta_TR, jet3phi_TR;
    std::vector<double> *otherJetspt_TR, *otherJetseta_TR, *otherJetsphi_TR;
    double mu1pt_TR, mu1eta_TR, mu1phi_TR;
    double mu2pt_TR, mu2eta_TR, mu2phi_TR;
    std::vector<double> *otherMuspt_TR, *otherMuseta_TR, *otherMusphi_TR;
    double ele1pt_TR, ele1eta_TR, ele1phi_TR;
    double ele2pt_TR, ele2eta_TR, ele2phi_TR;
    std::vector<double> *otherElespt_TR, *otherEleseta_TR, *otherElesphi_TR;
    double dPhi0_CUT, dPhi1_CUT, dPhi2_CUT;
    int nJets_CUT;
    double genJet1pt_TR, genJet1eta_TR, genJet1phi_TR;
    double genJet2pt_TR, genJet2eta_TR, genJet2phi_TR;
    double genJet3pt_TR, genJet3eta_TR, genJet3phi_TR;
    std::vector<double> *otherGenJetspt_TR, *otherGenJetseta_TR, *otherGenJetsphi_TR;
    double genMET_TR;
    double jet1Res_TR, jet2Res_TR, jet3Res_TR;
    std::vector<double> *otherJetsRes_TR;
    int externalBitToTree_TR;
	 UInt_t jet_num;
	 Float_t alljets_E[Njet_max], alljets_Px[Njet_max], alljets_Py[Njet_max], alljets_Pz[Njet_max];
	 Float_t alljets_Pt[Njet_max], alljets_Eta[Njet_max];

    void setTreeDefaultVars();

    TH1D *h1_cutFlow;

    TH1D *h1_evtWeight;

    TH1D *h1_metSgnf, *h1_lowPU_metSgnf, *h1_highPU_metSgnf;
    TH1D *h1_metSgnf_Prob, *h1_lowPU_metSgnf_Prob, *h1_highPU_metSgnf_Prob;
    TH1D *h1_nVtx, *h1_nJets, *h1_nMuons, *h1_nEles;

    TH1D *h1_mhtphi, *h1_metphi;
    TH1D *h1_mht, *h1_ht, *h1_met, *h1_mt;
    TH1D *h1_lowPU_mht, *h1_lowPU_ht, *h1_lowPU_met, *h1_lowPU_mt;
    TH1D *h1_highPU_mht, *h1_highPU_ht, *h1_highPU_met, *h1_highPU_mt;

    TH1D *h1_alljetspt, *h1_jet1pt, *h1_jet2pt, *h1_jet3pt;
    TH1D *h1_alljetseta, *h1_jet1eta, *h1_jet2eta, *h1_jet3eta;
    TH1D *h1_alljetsphi, *h1_jet1phi, *h1_jet2phi, *h1_jet3phi;

    TH1D *h1_allelept, *h1_alleleeta, *h1_allelephi;
    TH1D *h1_allmupt, *h1_allmueta, *h1_allmuphi;
    TH1D *h1_ele1pt, *h1_ele2pt;
    TH1D *h1_ele1eta, *h1_ele2eta;
    TH1D *h1_ele1phi, *h1_ele2phi;
    TH1D *h1_mu1pt, *h1_mu2pt;
    TH1D *h1_mu1eta, *h1_mu2eta;
    TH1D *h1_mu1phi, *h1_mu2phi;

    TH1D *h1_allgenJetspt, *h1_genJet1pt, *h1_genJet2pt, *h1_genJet3pt;
    TH1D *h1_allgenJetseta, *h1_genJet1eta, *h1_genJet2eta, *h1_genJet3eta;
    TH1D *h1_allgenJetsphi, *h1_genJet1phi, *h1_genJet2phi, *h1_genJet3phi;

    TH1D *h1_genMET;
 
    TH1D *h1_allJetsRes, *h1_jet1Res, *h1_jet2Res, *h1_jet3Res;
 
    double recoGenJetsDR_;

    vector<TH1D*> h1_jecLevAllJetsPtVec, h1_jecLevAllJetsEtaVec, h1_jecLevAllJetsPhiVec;

    bool doSgnf_;

	 edm::LumiReWeighting LumiWeights_;
	 std::vector<float> vMCNvtxDist, vDATANvtxDist;
	 double sumLumiWeights; //sum of lumi weights for cross checks
	 float storedWeight; //for flat MC sample weights
	 float PUWeight;
	 int mcFlag;   //0=DATA, 1= MC
	 unsigned int uProcessed; // number of processed events
	 unsigned int uPassed; //number of events passed the filter
};

FlatTreeMaker::FlatTreeMaker(const edm::ParameterSet & iConfig) {

   isData = true;
 
   jetSrc_ = iConfig.getParameter<edm::InputTag>("jetSrc");
   muonSrc_ = iConfig.getParameter<edm::InputTag>("muonSrc");
   eleSrc_ = iConfig.getParameter<edm::InputTag>("eleSrc");
//   tauSrc_ = iConfig.getParameter<edm::InputTag>("tauSrc");
  
   forVetoMuonSrc_ = iConfig.getParameter<edm::InputTag>("forVetoMuonSrc");
   forVetoElectronSrc_ = iConfig.getParameter<edm::InputTag>("forVetoElectronSrc");
 
   genJetSrc_ = iConfig.getParameter<edm::InputTag>("genJetSrc");
   genMETSrc_ = iConfig.getParameter<edm::InputTag>("genMETSrc");
   genParticleSrc_ = iConfig.getParameter<edm::InputTag>("genParticleSrc");

   recoGenJetsDR_ = iConfig.getParameter<double>("recoGenJetsDR");
 
   vtxSrc_ = iConfig.getParameter<edm::InputTag>("vtxSrc");
   //evtWeightInput_ = iConfig.getParameter<edm::InputTag>("evtWeightInput");
   nVtxPUcut_ = iConfig.getParameter<unsigned int>("nVtxPUcut");
 
   metSrc_ = iConfig.getParameter<edm::InputTag>("metSrc");
   mhtSrc_ = iConfig.getParameter<edm::InputTag>("mhtSrc");
   htSrc_ = iConfig.getParameter<edm::InputTag>("htSrc");

   //dPhis_CUT_vec_Src_ = iConfig.getParameter<edm::InputTag>("dPhis_CUT_vec_Src");
   //nJets_CUT_Src_ = iConfig.getParameter<edm::InputTag>("nJets_CUT_Src");

   //nLeptonsSels_ = iConfig.getParameter<std::vector<int> >("nLeptonsSels");
 
   debug_ = iConfig.getParameter<bool>("debug");

   doFillTree_ = iConfig.getParameter<bool>("doFillTree");

   //externalBitToTree_Src_ = iConfig.getParameter<edm::InputTag>("externalBitToTree_Src");

   doSgnf_ = iConfig.getParameter<bool>("doSgnf");
 
   outRootName_ = iConfig.getParameter<std::string>("outRootName");
 
   std::cout<<"producing root file : "<<outRootName_.c_str()<<std::endl;
 
   outRootFile = new TFile( outRootName_.c_str(), "RECREATE" );

   otherJetspt_TR = new std::vector<double>; otherJetseta_TR = new std::vector<double>; otherJetsphi_TR = new std::vector<double>;
   otherMuspt_TR = new std::vector<double>; otherMuseta_TR = new std::vector<double>; otherMusphi_TR = new std::vector<double>;
   otherElespt_TR = new std::vector<double>; otherEleseta_TR = new std::vector<double>; otherElesphi_TR = new std::vector<double>;
   otherGenJetspt_TR = new std::vector<double>; otherGenJetseta_TR = new std::vector<double>; otherGenJetsphi_TR = new std::vector<double>;
   otherJetsRes_TR = new std::vector<double>;

   npv = -1; avg_npv = -1;

   setTreeDefaultVars();

   if( doFillTree_ ){
      outTree = new TTree("AUX", "aux info");
      outTree->Branch("run", &run, "run/I");
      outTree->Branch("event", &event, "event/I");
      outTree->Branch("lumi", &ls, "lumi/I");
      outTree->Branch("npv", &npv, "npv/I");
      outTree->Branch("avg_npv", &avg_npv, "avg_npv/D");
      outTree->Branch("vtxSize", &vtxSize, "vtxSize/I");
      outTree->Branch("nJets", &nJets, "nJets/I");
      outTree->Branch("nMuons", &nMuons, "nMuons/I");
      outTree->Branch("nMuons_CUT", &nMuonsForVeto, "nMuonsForVeto/I");
      outTree->Branch("nElectrons", &nElectrons, "nElectrons/I");
      outTree->Branch("nElectrons_CUT", &nElectronsForVeto, "nElectronsForVeto/I");
      outTree->Branch("evtWeight", &evtWeight_TR, "evtWeight/D");
      outTree->Branch("meff", &meff, "meff/D");
      outTree->Branch("mht", &mht_TR, "mht/D");
      outTree->Branch("ht", &ht_TR, "ht/D");
      outTree->Branch("met", &met_TR, "met/D");
      outTree->Branch("mt", &mt_TR, "mt/D");
      outTree->Branch("mhtphi", &mhtphi_TR, "mhtphi/D");
      outTree->Branch("metphi", &metphi_TR, "metphi/D");
      outTree->Branch("metSgnf", &metSgnf_TR, "metSgnf/D");
      outTree->Branch("metSgnfProb", &metSgnfProb_TR, "metSgnfProb/D");
      outTree->Branch("jet1pt", &jet1pt_TR, "jet1pt/D");
      outTree->Branch("jet1eta", &jet1eta_TR, "jet1eta/D");
      outTree->Branch("jet1phi", &jet1phi_TR, "jet1phi/D");
      outTree->Branch("jet2pt", &jet2pt_TR, "jet2pt/D");
      outTree->Branch("jet2eta", &jet2eta_TR, "jet2eta/D");
      outTree->Branch("jet2phi", &jet2phi_TR, "jet2phi/D");
      outTree->Branch("jet3pt", &jet3pt_TR, "jet3pt/D");
      outTree->Branch("jet3eta", &jet3eta_TR, "jet3eta/D");
      outTree->Branch("jet3phi", &jet3phi_TR, "jet3phi/D");
      outTree->Branch("otherJetspt", "std::vector<double>", &otherJetspt_TR, 32000, 0);
      outTree->Branch("otherJetseta", "std::vector<double>", &otherJetseta_TR, 32000, 0);
      outTree->Branch("otherJetsphi", "std::vector<double>", &otherJetsphi_TR, 32000, 0);
      outTree->Branch("mu1pt", &mu1pt_TR, "mu1pt/D");
      outTree->Branch("mu1eta", &mu1eta_TR, "mu1eta/D");
      outTree->Branch("mu1phi", &mu1phi_TR, "mu1phi/D");
      outTree->Branch("mu2pt", &mu2pt_TR, "mu2pt/D");
      outTree->Branch("mu2eta", &mu2eta_TR, "mu2eta/D");
      outTree->Branch("mu2phi", &mu2phi_TR, "mu2phi/D");
      outTree->Branch("otherMuspt", "std::vector<double>", &otherMuspt_TR, 32000, 0);
      outTree->Branch("otherMuseta", "std::vector<double>", &otherMuseta_TR, 32000, 0);
      outTree->Branch("otherMusphi", "std::vector<double>", &otherMusphi_TR, 32000, 0);
      outTree->Branch("ele1pt", &ele1pt_TR, "ele1pt/D");
      outTree->Branch("ele1eta", &ele1eta_TR, "ele1eta/D");
      outTree->Branch("ele1phi", &ele1phi_TR, "ele1phi/D");
      outTree->Branch("ele2pt", &ele2pt_TR, "ele2pt/D");
      outTree->Branch("ele2eta", &ele2eta_TR, "ele2eta/D");
      outTree->Branch("ele2phi", &ele2phi_TR, "ele2phi/D");
      outTree->Branch("otherElespt", "std::vector<double>", &otherElespt_TR, 32000, 0);
      outTree->Branch("otherEleseta", "std::vector<double>", &otherEleseta_TR, 32000, 0);
      outTree->Branch("otherElesphi", "std::vector<double>", &otherElesphi_TR, 32000, 0);
      outTree->Branch("dPhi0_CUT", &dPhi0_CUT, "dPhi0_CUT/D"); 
      outTree->Branch("dPhi1_CUT", &dPhi1_CUT, "dPhi1_CUT/D"); 
      outTree->Branch("dPhi2_CUT", &dPhi2_CUT, "dPhi2_CUT/D");
      outTree->Branch("nJets_CUT", &nJets_CUT, "nJets_CUT/I");
      outTree->Branch("genJet1pt", &genJet1pt_TR, "genJet1pt/D");
      outTree->Branch("genJet1eta", &genJet1eta_TR, "genJet1eta/D");
      outTree->Branch("genJet1phi", &genJet1phi_TR, "genJet1phi/D");
      outTree->Branch("genJet2pt", &genJet2pt_TR, "genJet2pt/D");
      outTree->Branch("genJet2eta", &genJet2eta_TR, "genJet2eta/D");
      outTree->Branch("genJet2phi", &genJet2phi_TR, "genJet2phi/D");
      outTree->Branch("genJet3pt", &genJet3pt_TR, "genJet3pt/D");
      outTree->Branch("genJet3eta", &genJet3eta_TR, "genJet3eta/D");
      outTree->Branch("genJet3phi", &genJet3phi_TR, "genJet3phi/D");
      outTree->Branch("otherGenJetspt", "std::vector<double>", &otherGenJetspt_TR, 32000, 0);
      outTree->Branch("otherGenJetseta", "std::vector<double>", &otherGenJetseta_TR, 32000, 0);
      outTree->Branch("otherGenJetsphi", "std::vector<double>", &otherGenJetsphi_TR, 32000, 0);
      outTree->Branch("genMET", &genMET_TR, "genMET/D");
      outTree->Branch("jet1Res", &jet1Res_TR, "jet1Res/D");
      outTree->Branch("jet2Res", &jet2Res_TR, "jet2Res/D");
      outTree->Branch("jet3Res", &jet3Res_TR, "jet3Res/D");
      outTree->Branch("otherJetsRes", "std::vector<double>", &otherJetsRes_TR, 32000, 0);
		outTree->Branch("jet_num", &jet_num, "jet_num/i");
		outTree->Branch("alljets_E", alljets_E, "alljets_E[jet_num]/F");
		outTree->Branch("alljets_Px", alljets_Px, "alljets_Px[jet_num]/F");
		outTree->Branch("alljets_Py", alljets_Py, "alljets_Py[jet_num]/F");
		outTree->Branch("alljets_Pz", alljets_Pz, "alljets_Pz[jet_num]/F");
		outTree->Branch("alljets_Pt", alljets_Pt, "alljets_Pt[jet_num]/F");
		outTree->Branch("alljets_Eta", alljets_Eta, "alljets_Eta[jet_num]/F");
		outTree->Branch("mcFlag", &mcFlag, "mcFlag/i");
		outTree->Branch("PUWeight", &PUWeight, "PUWeight/F");
		outTree->Branch("storedWeight", &storedWeight, "storedWeight/F");
      
      outTree->Branch("externalBit", &externalBitToTree_TR, "externalBitToTree_TR/I");
   }
 
//   char names[200], dispt[200];
   h1_cutFlow = new TH1D("cutFlow", "cut flow;cut", 50, 0, 50);
   int binNumCnt =0;
   binNumCnt ++; h1_cutFlow->GetXaxis()->SetBinLabel(binNumCnt,"All"); // 0
   binNumCnt ++; h1_cutFlow->GetXaxis()->SetBinLabel(binNumCnt,"AftLepReq"); // 1
 
   h1_evtWeight = new TH1D("evtWeight", "evtWeight;weight", 1000, 0, 10);
 
   h1_nVtx = new TH1D("nVtx", "nVtx;n_{Vtx}", 50, 0, 50);
   h1_nJets = new TH1D("nJets", "nJets;n_{Jets}", 20, 0, 20);
   h1_nMuons = new TH1D("nMuons", "nMuons;n_{Muons}", 20, 0, 20);
   h1_nEles = new TH1D("nEles", "nEles;n_{Electrons}", 20, 0, 20);
 
   h1_metSgnf = new TH1D("METsgnf", "MET significance;S_{MET}", 200, 0, 20);     
   h1_lowPU_metSgnf = new TH1D("lowPU_METsgnf", "MET significance;S_{MET}", 200, 0, 20);
   h1_highPU_metSgnf = new TH1D("highPU_METsgnf", "MET significance;S_{MET}", 200, 0, 20);
     
   h1_metSgnf_Prob = new TH1D("METsgnf_Prob", "MET significance probability;P(#chi^{2})", 200, 0, 1);    
   h1_lowPU_metSgnf_Prob = new TH1D("lowPU_METsgnf_Prob", "MET significance probability;P(#chi^{2})", 200, 0, 1);
   h1_highPU_metSgnf_Prob = new TH1D("highPU_METsgnf_Prob", "MET significance probability;P(#chi^{2})", 200, 0, 1);
     
   h1_mhtphi = new TH1D("MHTphi", "MHT #phi;MHT #phi", 200, -3.2, 3.2);
   h1_metphi = new TH1D("METphi", "MET #phi;MET #phi", 200, -3.2, 3.2);
     
   h1_mht = new TH1D("MHT", "MHT;MHT(GeV)", 100, 0, 800);
   h1_ht = new TH1D("HT", "HT;HT(GeV)", 100, 0, 1500); 
   h1_met = new TH1D("MET", "PF MET;MET(GeV)", 100, 0, 800);
   h1_mt = new TH1D("MT", "M_{T};M_{T}(GeV)", 100, 0, 150);
     
   h1_genMET = new TH1D("GenMET", "Gen MET;MET(GeV)", 100, 0, 800);

   h1_lowPU_mht = new TH1D("lowPU_MHT", "MHT;MHT(GeV)", 100, 0, 800);
   h1_lowPU_ht = new TH1D("lowPU_HT", "HT;HT(GeV)", 100, 0, 1500); 
   h1_lowPU_met = new TH1D("lowPU_MET", "PF MET;MET(GeV)", 100, 0, 800);
   h1_lowPU_mt = new TH1D("lowPU_MT", "M_{T};M_{T}(GeV)", 100, 0, 150);
 
   h1_highPU_mht = new TH1D("highPU_MHT", "MHT;MHT(GeV)", 100, 0, 800);
   h1_highPU_ht = new TH1D("highPU_HT", "HT;HT(GeV)", 100, 0, 1500);
   h1_highPU_met = new TH1D("highPU_MET", "PF MET;MET(GeV)", 100, 0, 800);
   h1_highPU_mt = new TH1D("highPU_MT", "M_{T};M_{T}(GeV)", 100, 0, 150);
 
   h1_alljetspt = new TH1D("AllJetsPt", "P_{T} of all jets;P_{T}(GeV)", 200, 0, 400);
   h1_jet1pt = new TH1D("Jet1Pt", "P_{T} of jet_{1};P_{T}(GeV)", 200, 0, 400);
   h1_jet2pt = new TH1D("Jet2Pt", "P_{T} of jet_{2};P_{T}(GeV)", 200, 0, 400);
   h1_jet3pt = new TH1D("Jet3Pt", "P_{T} of jet_{3};P_{T}(GeV)", 200, 0, 400);
   h1_alljetseta = new TH1D("AllJetsEta", "#eta of all jets;#eta", 200, -3.2, 3.2);
   h1_jet1eta = new TH1D("Jet1Eta", "#eta of jet_{1};#eta", 200, -3.2, 3.2);
   h1_jet2eta = new TH1D("Jet2Eta", "#eta of jet_{2};#eta", 200, -3.2, 3.2);
   h1_jet3eta = new TH1D("Jet3Eta", "#eta of jet_{3};#eta", 200, -3.2, 3.2);
   h1_alljetsphi = new TH1D("AllJetsPhi", "#phi of all jets;#phi", 200, -3.2, 3.2);
   h1_jet1phi = new TH1D("Jet1Phi", "#phi of jet_{1};#phi", 200, -3.2, 3.2);
   h1_jet2phi = new TH1D("Jet2Phi", "#phi of jet_{2};#phi", 200, -3.2, 3.2);
   h1_jet3phi = new TH1D("Jet3Phi", "#phi of jet_{3};#phi", 200, -3.2, 3.2);
 
   h1_allmupt = new TH1D("AllMuonsPt", "P_{T} of all muons;P_{T}(GeV)", 200, 0, 300);
   h1_mu1pt = new TH1D("Mu1Pt", "P_{T} of #mu_{1};P_{T}(GeV)", 200, 0, 300);
   h1_mu2pt = new TH1D("Mu2Pt", "P_{T} of #mu_{2};P_{T}(GeV)", 200, 0, 300);
   h1_allmueta = new TH1D("AllMuonsEta", "#eta of all muons;#eta", 200, -3.2, 3.2);
   h1_mu1eta = new TH1D("Mu1Eta", "#eta of #mu_{1};#eta", 200, -3.2, 3.2);
   h1_mu2eta = new TH1D("Mu2Eta", "#eta of #mu_{2};#eta", 200, -3.2, 3.2);
   h1_allmuphi = new TH1D("AllMuonsPhi", "#phi of all muons;#phi", 200, -3.2, 3.2);
   h1_mu1phi = new TH1D("Mu1Phi", "#phi of #mu_{1};#phi", 200, -3.2, 3.2);
   h1_mu2phi = new TH1D("Mu2Phi", "#phi of #mu_{2};#phi", 200, -3.2, 3.2);
 
   h1_allelept = new TH1D("AllElesPt", "P_{T} of all eles;P_{T}(GeV)", 200, 0, 300);
   h1_ele1pt = new TH1D("Ele1Pt", "P_{T} of ele_{1};P_{T}(GeV)", 200, 0, 300);
   h1_ele2pt = new TH1D("Ele2Pt", "P_{T} of ele_{2};P_{T}(GeV)", 200, 0, 300);
   h1_alleleeta = new TH1D("AllElesEta", "#eta of all eles;#eta", 200, -3.2, 3.2);
   h1_ele1eta = new TH1D("Ele1Eta", "#eta of ele_{1};#eta", 200, -3.2, 3.2);
   h1_ele2eta = new TH1D("Ele2Eta", "#eta of ele_{2};#eta", 200, -3.2, 3.2);
   h1_allelephi = new TH1D("AllElesPhi", "#phi of all eles;#phi", 200, -3.2, 3.2);
   h1_ele1phi = new TH1D("Ele1Phi", "#phi of ele_{1};#phi", 200, -3.2, 3.2);
   h1_ele2phi = new TH1D("Ele2Phi", "#phi of ele_{2};#phi", 200, -3.2, 3.2);

   h1_allgenJetspt = new TH1D("AllGenJetsPt", "P_{T} of all genJets;P_{T}(GeV)", 200, 0, 400);
   h1_genJet1pt = new TH1D("GenJet1Pt", "P_{T} of genJet_{1};P_{T}(GeV)", 200, 0, 400);
   h1_genJet2pt = new TH1D("GenJet2Pt", "P_{T} of genJet_{2};P_{T}(GeV)", 200, 0, 400);
   h1_genJet3pt = new TH1D("GenJet3Pt", "P_{T} of genJet_{3};P_{T}(GeV)", 200, 0, 400);
   h1_allgenJetseta = new TH1D("AllGenJetsEta", "#eta of all genJets;#eta", 200, -3.2, 3.2);
   h1_genJet1eta = new TH1D("GenJet1Eta", "#eta of genJet_{1};#eta", 200, -3.2, 3.2);
   h1_genJet2eta = new TH1D("GenJet2Eta", "#eta of genJet_{2};#eta", 200, -3.2, 3.2);
   h1_genJet3eta = new TH1D("GenJet3Eta", "#eta of genJet_{3};#eta", 200, -3.2, 3.2);
   h1_allgenJetsphi = new TH1D("AllGenJetsPhi", "#phi of all genJets;#phi", 200, -3.2, 3.2);
   h1_genJet1phi = new TH1D("GenJet1Phi", "#phi of genJet_{1};#phi", 200, -3.2, 3.2);
   h1_genJet2phi = new TH1D("GenJet2Phi", "#phi of genJet_{2};#phi", 200, -3.2, 3.2);
   h1_genJet3phi = new TH1D("GenJet3Phi", "#phi of genJet_{3};#phi", 200, -3.2, 3.2);

   h1_allJetsRes = new TH1D("AllJetsRes", "Response of all jets;jet_{P_{T}}^{RECO}/jet_{P_{T}}^{GEN}", 200, 0, 2.5);
   h1_jet1Res = new TH1D("Jet1Res", "Response of jet_{1};jet_{P_{T}}^{RECO}/jet_{P_{T}}^{GEN}", 200, 0, 2.5);
   h1_jet2Res = new TH1D("Jet2Res", "Response of jet_{2};jet_{P_{T}}^{RECO}/jet_{P_{T}}^{GEN}", 200, 0, 2.5);
   h1_jet3Res = new TH1D("Jet3Res", "Response of jet_{3};jet_{P_{T}}^{RECO}/jet_{P_{T}}^{GEN}", 200, 0, 2.5);
 
}

FlatTreeMaker::~FlatTreeMaker() {

   outRootFile->cd();

   h1_cutFlow->Write();

   h1_evtWeight->Write();
    
   h1_nVtx->Write(); h1_nJets->Write(); h1_nMuons->Write(); h1_nEles->Write();
    
   h1_metSgnf->Write(); h1_lowPU_metSgnf->Write(); h1_highPU_metSgnf->Write();
   h1_metSgnf_Prob->Write(); h1_lowPU_metSgnf_Prob->Write(); h1_highPU_metSgnf_Prob->Write();

   h1_mhtphi->Write(); h1_metphi->Write();
   h1_mht->Write(); h1_ht->Write(); h1_met->Write(); h1_mt->Write();
   h1_lowPU_mht->Write(); h1_lowPU_ht->Write(); h1_lowPU_met->Write(); h1_lowPU_mt->Write();
   h1_highPU_mht->Write(); h1_highPU_ht->Write(); h1_highPU_met->Write(); h1_highPU_mt->Write();
    
   h1_alljetspt->Write(); h1_jet1pt->Write(); h1_jet2pt->Write(); h1_jet3pt->Write();
   h1_alljetseta->Write(); h1_jet1eta->Write(); h1_jet2eta->Write(); h1_jet3eta->Write();
   h1_alljetsphi->Write(); h1_jet1phi->Write(); h1_jet2phi->Write(); h1_jet3phi->Write();
    
   h1_allelept->Write(); h1_alleleeta->Write(); h1_allelephi->Write();
   h1_allmupt->Write(); h1_allmueta->Write(); h1_allmuphi->Write();
   h1_ele1pt->Write(); h1_ele2pt->Write();
   h1_ele1eta->Write(); h1_ele2eta->Write();
   h1_ele1phi->Write(); h1_ele2phi->Write();
   h1_mu1pt->Write(); h1_mu2pt->Write();
   h1_mu1eta->Write(); h1_mu2eta->Write();
   h1_mu1phi->Write(); h1_mu2phi->Write();

   h1_allgenJetspt->Write(); h1_genJet1pt->Write(); h1_genJet2pt->Write(); h1_genJet3pt->Write();
   h1_allgenJetseta->Write(); h1_genJet1eta->Write(); h1_genJet2eta->Write(); h1_genJet3eta->Write();
   h1_allgenJetsphi->Write(); h1_genJet1phi->Write(); h1_genJet2phi->Write(); h1_genJet3phi->Write();

   h1_genMET->Write();

   h1_allJetsRes->Write(); h1_jet1Res->Write(); h1_jet2Res->Write(); h1_jet3Res->Write();

   if( doFillTree_ ){ outTree->Write(); delete outTree; }
   outRootFile->Write(); outRootFile->Close();
   if( outRootFile ) delete outRootFile;

}

void FlatTreeMaker::setTreeDefaultVars(){

   evtWeight_TR = 1.0;
   mht_TR= -99, ht_TR= -99, met_TR= -99, mt_TR= -99;
   mhtphi_TR= -99, metphi_TR= -99;
   metSgnf_TR= -99, metSgnfProb_TR= -99;
   jet1pt_TR= -99, jet1eta_TR= -99, jet1phi_TR= -99;
   jet2pt_TR= -99, jet2eta_TR= -99, jet2phi_TR= -99;
   jet3pt_TR= -99, jet3eta_TR= -99, jet3phi_TR= -99;
   otherJetspt_TR->clear(); otherJetseta_TR->clear(); otherJetsphi_TR->clear();
   mu1pt_TR= -99, mu1eta_TR= -99, mu1phi_TR= -99;
   mu2pt_TR= -99, mu2eta_TR= -99, mu2phi_TR= -99;
   otherMuspt_TR->clear(); otherMuseta_TR->clear(); otherMusphi_TR->clear();
   ele1pt_TR= -99, ele1eta_TR= -99, ele1phi_TR= -99;
   ele2pt_TR= -99, ele2eta_TR= -99, ele2phi_TR= -99;
   otherElespt_TR->clear(); otherEleseta_TR->clear(); otherElesphi_TR->clear();

   dPhi0_CUT = -99, dPhi1_CUT = -99, dPhi2_CUT = -99;
   nJets_CUT = -99;

   genJet1pt_TR= -99, genJet1eta_TR= -99, genJet1phi_TR= -99;
   genJet2pt_TR= -99, genJet2eta_TR= -99, genJet2phi_TR= -99;
   genJet3pt_TR= -99, genJet3eta_TR= -99, genJet3phi_TR= -99;
   otherGenJetspt_TR->clear(); otherGenJetseta_TR->clear(); otherGenJetsphi_TR->clear();
   genMET_TR = -99;
   jet1Res_TR = -99, jet2Res_TR = -99, jet3Res_TR = -99;
   otherJetsRes_TR->clear();
   externalBitToTree_TR = -99;

	jet_num = 0;
	meff = -99;
	mcFlag = -1;
	storedWeight = 1;

}

// ------------ method called on each new Event  ------------
bool FlatTreeMaker::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {

	++uProcessed;

   setTreeDefaultVars();

// Event setup
   loadEventInfo(iEvent, iSetup);
	mcFlag = (int) (! isData);

   loadGenInfo(iEvent);
   loadRecoJets(iEvent);
   loadLeptons(iEvent);
   loadMETMHT(iEvent);
   loadHT(iEvent);
//   loadAUX(iEvent); if( externalBitToTree_.isValid() ) externalBitToTree_TR = (*externalBitToTree_);

//   if( nLeptonsSels_.size() <2 ){ std::cout<<"\nERROR ... input nLeptonsSels_.size : "<<nLeptonsSels_.size()<<"  less than 2?!"<<std::endl; } 

// If muons or electrons are not the size we require, don't go further
//   if( (nLeptonsSels_[0] !=-1 && (int)nMuonsForVeto != nLeptonsSels_[0]) || (nLeptonsSels_[1] !=-1 && (int)nElectronsForVeto != nLeptonsSels_[1]) ){
//      if( debug_ ){
//         std::cout<<"nMuonsForVeto : "<<nMuonsForVeto<<"  nElectronsForVeto : "<<nElectronsForVeto<<std::endl;
//         std::cout<<"nLeptonsSels_[0] : "<<nLeptonsSels_[0]<<"  nLeptonsSels_[1] : "<<nLeptonsSels_[1]<<std::endl;
//      }
//      return true;
//   }


	/* PU S3 reweighting for QCD MC sample
	 */
	PUWeight = 1;
	if (! isData)
	{
		PUWeight = LumiWeights_.weight( iEvent );
		//std::cout << "lum wgt = " << PUWeight << std::endl;
		sumLumiWeights += PUWeight;
	}

   double sgnfMET =0, probMET =0;
   if( doSgnf_ ){ sgnfMET = (*met)[0].significance(); probMET = TMath::Prob(sgnfMET, sgnfnDof); }

   double mtw = -1;
   if( nMuons == 1 ){
      mtw = sqrt( 2*( (*met)[0].pt()*muon1.pt() -( (*met)[0].px()*muon1.px() + (*met)[0].py()*muon1.py() ) ) ); 
   }
   if( nMuons ==2 && muon1Charge*muon2Charge == -1 ){
      mtw = sqrt( 2*( muon2.pt()*muon1.pt() -( muon2.px()*muon1.px() + muon2.py()*muon1.py() ) ) );
   }

   mht_TR = (*mht)[0].pt(); ht_TR = (*ht); met_TR = (*met)[0].pt(); mt_TR = mtw;
   mhtphi_TR = (*mht)[0].phi(); metphi_TR = (*met)[0].phi();
   metSgnf_TR = sgnfMET; metSgnfProb_TR = probMET;
	meff = mht_TR + ht_TR;


   if( !isData ){
      genMET_TR = (*genMET)[0].pt();
      for(size ig=0; ig<nGenJets; ig++){
         if( ig ==0 ){
            genJet1 = (*genJets)[ig];
            genJet1pt_TR = genJet1.pt(); genJet1eta_TR = genJet1.eta(); genJet1phi_TR = genJet1.phi();
         }
         if( ig ==1 ){
            genJet2 = (*genJets)[ig];
            genJet2pt_TR = genJet2.pt(); genJet2eta_TR = genJet2.eta(); genJet2phi_TR = genJet2.phi();
         }
         if( ig ==2 ){
            genJet3 = (*genJets)[ig];
            genJet3pt_TR = genJet3.pt(); genJet3eta_TR = genJet3.eta(); genJet3phi_TR = genJet3.phi();
         }
         if( ig >2 ){
            genJetother = (*genJets)[ig];
            otherGenJetspt_TR->push_back(genJetother.pt()); otherGenJetseta_TR->push_back(genJetother.eta()); otherGenJetsphi_TR->push_back(genJetother.phi()); 
         }
      }   
   }

	if (nJets > Njet_max) 
	{ 
		std::cout << __FILE__ << ":" << __LINE__ << ":Max njet limi of "<< Njet_max << " is exceeded!" << std::endl;
		assert (false); 
	}
	jet_num = nJets;
   for(size ij=0; ij<nJets; ij++){
		alljets_E[ij] = (*jets)[ij].energy();	
		alljets_Px[ij] = (*jets)[ij].px();	
		alljets_Py[ij] = (*jets)[ij].py();	
		alljets_Pz[ij] = (*jets)[ij].pz();	
		alljets_Pt[ij] = (*jets)[ij].pt();	
		alljets_Eta[ij] = (*jets)[ij].eta();	
	
      if( !isData ){
         double jeteta = (*jets)[ij].eta(), jetphi = (*jets)[ij].phi(), jetpt = (*jets)[ij].pt();
         double minDR = 999.0, jetRes = -999.0;
         for(size ig=0; ig<nGenJets; ig++){
            double genjeteta = (*genJets)[ig].eta(), genjetphi = (*genJets)[ig].phi(), genjetpt = (*genJets)[ig].pt();
            const double dR = reco::deltaR(jeteta, jetphi, genjeteta, genjetphi);
            if( minDR > dR ){ minDR = dR; jetRes = jetpt/genjetpt; }
         }
         if( minDR < recoGenJetsDR_ ){
            if( ij == 0 ){
               jet1Res_TR = jetRes;
            }
            if( ij == 1 ){
               jet2Res_TR = jetRes;
            }
            if( ij == 2 ){
               jet3Res_TR = jetRes;
            }
            if( ij > 2 ){
               otherJetsRes_TR->push_back(jetRes);
            }
         }
      }
      if( ij ==0 ){
         jet1 = (*jets)[ij];
         jet1pt_TR = jet1.pt(); jet1eta_TR = jet1.eta(); jet1phi_TR = jet1.phi();
      }
      if( ij ==1 ){
         jet2 = (*jets)[ij];
         jet2pt_TR = jet2.pt(); jet2eta_TR = jet2.eta(); jet2phi_TR = jet2.phi();
      }
      if( ij ==2 ){
         jet3 = (*jets)[ij]; 
         jet3pt_TR = jet3.pt(); jet3eta_TR = jet3.eta(); jet3phi_TR = jet3.phi();
      }
      if( ij > 2 ){
         jetother = (*jets)[ij];
         otherJetspt_TR->push_back(jetother.pt()); otherJetseta_TR->push_back(jetother.eta()); otherJetsphi_TR->push_back(jetother.phi()); 
      }
   }

   for(size im=0; im<nMuons; im++){
      if( im == 0){
         mu1pt_TR = muon1.pt(); mu1eta_TR = muon1.eta(); mu1phi_TR = muon1.phi();
      }
      if( im == 1){
         mu2pt_TR = muon2.pt(); mu2eta_TR = muon2.eta(); mu2phi_TR = muon2.phi();
      }
      if( im >1 ){
         otherMuspt_TR->push_back((*muons)[im].pt()); otherMuseta_TR->push_back((*muons)[im].eta()); otherMusphi_TR->push_back((*muons)[im].phi());
      }
   }

   for(size ie=0; ie<nElectrons; ie++){
      if( ie == 0){
         ele1pt_TR = ele1.pt(); ele1eta_TR = ele1.eta(); ele1phi_TR = ele1.phi();
      }
      if( ie == 1){
         ele2pt_TR = ele2.pt(); ele2eta_TR = ele2.eta(); ele2phi_TR = ele2.phi();
      }
      if( ie >1 ){
         otherElespt_TR->push_back((*electrons)[ie].pt()); otherEleseta_TR->push_back((*electrons)[ie].eta()); otherElesphi_TR->push_back((*electrons)[ie].phi());
      }
   }

	if( doFillTree_ ){
		//      dPhi0_CUT = (*dPhis_CUT_vec_)[0]; dPhi1_CUT = (*dPhis_CUT_vec_)[1]; dPhi2_CUT = (*dPhis_CUT_vec_)[2];
		//      nJets_CUT = (*nJets_CUT_);
		outTree->Fill(); 
		++uPassed;
	}

   return true;

}

// ------------ method called once each job just before starting event loop  ------------
void FlatTreeMaker::beginJob() {

	sumLumiWeights = 0;
	uProcessed = 0;
	uPassed = 0;

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
void FlatTreeMaker::endJob() {
	std::cout << ">>>>>" <<  __FILE__ << "::" << __FUNCTION__  << "\n" << std::endl;
	std::cout << "[ATM:00] Job settings " << std::endl;
	std::cout << "[ATM:01] Events Processed --- = " << uProcessed << std::endl;
	std::cout << "[ATM:02] Events Passed ------ = " << uPassed << std::endl;
	std::cout << "[ATM:03] LumiWeights Avg ---- = " << sumLumiWeights/(double)uPassed << std::endl;

}

// ------------ method called once each run just before starting event loop  ------------
bool FlatTreeMaker::beginRun(edm::Run &run, const edm::EventSetup& iSetup) {
  return true;
}

// ------------ method called once each run just after starting event loop  ------------
bool FlatTreeMaker::endRun(edm::Run &run, const edm::EventSetup& iSetup) {
  return true;
}


void FlatTreeMaker::loadEventInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup){

// Determine if it's data
   if( !iEvent.isRealData() ) isData = false;

// Get run, event, lumi info
   run = iEvent.id().run();
   event = iEvent.id().event();
   ls = iEvent.luminosityBlock();

// Get vertices
   iEvent.getByLabel(vtxSrc_, vertices); vtxSize = vertices->size();
   
// Get event weight
//   iEvent.getByLabel(evtWeightInput_, evtWeight_);
	if (isData)
	{
		storedWeight = 1;
	} else {
		edm::Handle<GenEventInfoProduct> genEvtInfoHandle;
		iEvent.getByLabel("generator", genEvtInfoHandle);
		storedWeight = genEvtInfoHandle->weight();
		//std::cout << "storedWeight = " << storedWeight << std::endl;
	}

}

void FlatTreeMaker::loadGenInfo(const edm::Event& iEvent){

// MC generate level related info
   scalePDF = -1; pthat = -1;
   if (!isData) {
      edm::Handle<edm::HepMCProduct> evt;
      edm::Handle< GenEventInfoProduct > GenInfoHandle;

      iEvent.getByLabel("generator", evt);
      if (evt.isValid()) {
         HepMC::GenEvent * myGenEvent = new HepMC::GenEvent(*(evt->GetEvent()));
         scalePDF = myGenEvent->event_scale();
         if( myGenEvent ) delete myGenEvent;
      }

      iEvent.getByLabel( "generator", GenInfoHandle );
      if (GenInfoHandle.isValid()) { pthat = ( GenInfoHandle->hasBinningValues() ? (GenInfoHandle->binningValues())[0] : 0.0); }

      iEvent.getByLabel(genParticleSrc_, genParticles);
      iEvent.getByLabel(genJetSrc_, genJets); nGenJets = genJets->size();
      iEvent.getByLabel(genMETSrc_, genMET); 

      iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);
      std::vector<PileupSummaryInfo>::const_iterator PVI;

      npv = -1; avg_npv = 0;
      for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {

         int BX = PVI->getBunchCrossing();

         avg_npv += double(PVI->getPU_NumInteractions());

         if(BX == 0) {
            npv = PVI->getPU_NumInteractions();
            continue;
         }
      }
      avg_npv /= 3.0;
   }
}

void FlatTreeMaker::loadRecoJets(const edm::Event& iEvent){
   iEvent.getByLabel(jetSrc_, jets); nJets = jets->size();
}

void FlatTreeMaker::loadLeptons(const edm::Event& iEvent){
   iEvent.getByLabel(muonSrc_, muons); nMuons = muons->size();
   iEvent.getByLabel(eleSrc_, electrons); nElectrons = electrons->size();
//   iEvent.getByLabel(tauSrc_, taus); nTaus = taus->size();
   iEvent.getByLabel(forVetoMuonSrc_, forVetoMuons); nMuonsForVeto = forVetoMuons->size();
   iEvent.getByLabel(forVetoElectronSrc_, forVetoElectrons); nElectronsForVeto = forVetoElectrons->size();

   for( size im=0; im<nMuons; im++){
      if( im ==0 ){ muon1 = (*muons)[im]; muon2 = (*muons)[im]; }
      if( im ==1 ) muon2 = (*muons)[im];
   }

   muon1Charge = muon1.charge(); muon2Charge = muon2.charge();

   for( size ie=0; ie<nElectrons; ie++){
      if( ie ==0 ){ ele1 = (*electrons)[ie]; ele2 = (*electrons)[ie]; }
      if( ie ==1 ) ele2 = (*electrons)[ie];
   }

   ele1Charge = ele1.charge(); ele2Charge = ele2.charge();

}

void FlatTreeMaker::loadPhotons(const edm::Event& iEvent){
   iEvent.getByLabel(photonSrc_, photons); nPhotons = photons->size();
}

void FlatTreeMaker::loadMETMHT(const edm::Event& iEvent){
   iEvent.getByLabel(metSrc_, met);
   iEvent.getByLabel(mhtSrc_, mht);
}

void FlatTreeMaker::loadHT(const edm::Event& iEvent){
   iEvent.getByLabel(htSrc_, ht);
}

//void FlatTreeMaker::loadAUX(const edm::Event& iEvent){
//   iEvent.getByLabel(dPhis_CUT_vec_Src_, dPhis_CUT_vec_);
//   iEvent.getByLabel(nJets_CUT_Src_, nJets_CUT_);

//   iEvent.getByLabel(externalBitToTree_Src_, externalBitToTree_);
//}

#include "FWCore/Framework/interface/MakerMacros.h"
//define this as a plug-in
DEFINE_FWK_MODULE(FlatTreeMaker);

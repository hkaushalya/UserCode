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
#include "TLorentzVector.h"

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

	 static const UInt_t Njet_max = 100;

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


    bool debug_;

    bool doFillTree_;

    bool isData;

    bool printOnce_;
    vector<std::string> availJECLevels_;

    std::string outRootName_;
    TFile *outRootFile;
    TTree *outTree;
    double evtWeight_TR;
    double mht_TR, ht_TR, met_TR, meff;
    double mhtphi_TR, metphi_TR;
    double metSgnf_TR, metSgnfProb_TR;
    double genMET_TR;
    double jet1Res_TR, jet2Res_TR, jet3Res_TR;
    std::vector<double> *otherJetsRes_TR;
    int externalBitToTree_TR;
	 UInt_t jet_num;

    void setTreeDefaultVars();

    double recoGenJetsDR_;
    bool doSgnf_;

	 edm::LumiReWeighting LumiWeights_;
	 std::vector<float> vMCNvtxDist, vDATANvtxDist;
	 double sumLumiWeights; //sum of lumi weights for cross checks
	 float storedWeight; //for flat MC sample weights
	 float PUWeight;
	 int mcFlag;   //0=DATA, 1= MC
	 std::vector<TLorentzVector> ak5pt50eta25jets;
	 std::vector<TLorentzVector> ak5pt30eta5jets;
	 std::vector<TLorentzVector> ak5jets;
    std::vector<TLorentzVector> GenJets;
	 unsigned int uProcessed; // number of processed events
	 unsigned int uPassed; //number of events passed the filter
};

FlatTreeMaker::FlatTreeMaker(const edm::ParameterSet & iConfig) {

   isData = true;
 
   jetSrc_ = iConfig.getParameter<edm::InputTag>("jetSrc");
 
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

 
   debug_ = iConfig.getParameter<bool>("debug");

   doFillTree_ = iConfig.getParameter<bool>("doFillTree");

   //externalBitToTree_Src_ = iConfig.getParameter<edm::InputTag>("externalBitToTree_Src");

   doSgnf_ = iConfig.getParameter<bool>("doSgnf");
 
   outRootName_ = iConfig.getParameter<std::string>("outRootName");
 
   std::cout<<"producing root file : "<<outRootName_.c_str()<<std::endl;
 
   outRootFile = new TFile( outRootName_.c_str(), "RECREATE" );

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
      outTree->Branch("evtWeight", &evtWeight_TR, "evtWeight/D");
      outTree->Branch("meff", &meff, "meff/D");
      outTree->Branch("mht", &mht_TR, "mht/D");
      outTree->Branch("ht", &ht_TR, "ht/D");
      outTree->Branch("met", &met_TR, "met/D");
      outTree->Branch("mhtphi", &mhtphi_TR, "mhtphi/D");
      outTree->Branch("metphi", &metphi_TR, "metphi/D");
      outTree->Branch("metSgnf", &metSgnf_TR, "metSgnf/D");
      outTree->Branch("metSgnfProb", &metSgnfProb_TR, "metSgnfProb/D");
      outTree->Branch("jet1Res", &jet1Res_TR, "jet1Res/D");
      outTree->Branch("jet2Res", &jet2Res_TR, "jet2Res/D");
      outTree->Branch("jet3Res", &jet3Res_TR, "jet3Res/D");
      outTree->Branch("otherJetsRes", "std::vector<double>", &otherJetsRes_TR, 32000, 0);
		outTree->Branch("jet_num", &jet_num, "jet_num/i");
		outTree->Branch("mcFlag", &mcFlag, "mcFlag/i");
		outTree->Branch("PUWeight", &PUWeight, "PUWeight/F");
		outTree->Branch("storedWeight", &storedWeight, "storedWeight/F");
      
      outTree->Branch("externalBit", &externalBitToTree_TR, "externalBitToTree_TR/I");
   }
 
}

FlatTreeMaker::~FlatTreeMaker() {

   outRootFile->cd();

   if( doFillTree_ ){ outTree->Write(); delete outTree; }
   outRootFile->Write(); outRootFile->Close();
   if( outRootFile ) delete outRootFile;

}

void FlatTreeMaker::setTreeDefaultVars(){

   evtWeight_TR = 1.0;
   mht_TR= -99, ht_TR= -99, met_TR= -99;
   mhtphi_TR= -99, metphi_TR= -99;
   metSgnf_TR= -99, metSgnfProb_TR= -99;

   genMET_TR = -99;

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
   loadMETMHT(iEvent);
   loadHT(iEvent);

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

   mht_TR = (*mht)[0].pt(); ht_TR = (*ht); met_TR = (*met)[0].pt(); 
   mhtphi_TR = (*mht)[0].phi(); metphi_TR = (*met)[0].phi();
   metSgnf_TR = sgnfMET; metSgnfProb_TR = probMET;
	meff = mht_TR + ht_TR;


   if( !isData ){
      genMET_TR = (*genMET)[0].pt();
      for(size ig=0; ig<nGenJets; ig++){
      //      otherGenJetspt_TR->push_back(genJetother.pt()); otherGenJetseta_TR->push_back(genJetother.eta()); otherGenJetsphi_TR->push_back(genJetother.phi()); 
      }   
   }

	if (nJets > Njet_max) 
	{ 
		std::cout << __FILE__ << ":" << __LINE__ << ":Max njet limi of "<< Njet_max << " is exceeded!" << std::endl;
		assert (false); 
	}
	jet_num = nJets;
   for(size ij=0; ij<nJets; ij++){
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

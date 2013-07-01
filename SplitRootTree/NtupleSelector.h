//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Apr 12 09:49:03 2013 by ROOT version 5.32/00
// from TTree tree/tree
// found on file: 74907_13.root
//////////////////////////////////////////////////////////

#ifndef NtupleSelector_h
#define NtupleSelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <string>

// Fixed size dimensions of array or collections stored in the TTree if any.
using namespace std;

class NtupleSelector : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          t_EvtRun;
   UInt_t          t_EvtLS;
   UInt_t          t_EvtEvent;
   Int_t           t_NVertices;
   Int_t           t_tru_Npv;
   Int_t           t_avg_Npv;
   Int_t           t_MCflag;
   Double_t        t_PUWeight;
   Double_t        t_PUWeightAB;
   Double_t        t_PUWeightABC;
   Double_t        t_PUWeightRA2;
   Double_t        t_EvtWeight;
   Double_t        t_PFMetPx;
   Double_t        t_PFMetPy;
   vector<double>  *t_PFJetPt;
   vector<double>  *t_PFJetEta;
   vector<double>  *t_PFJetPhi;
   vector<double>  *t_PFJetE;
   vector<double>  *t_PFJetBTag;
   vector<double>  *t_PFJetNHF;
   vector<double>  *t_PFJetEMF;
   vector<double>  *t_PFJetHOEne;
   Int_t           t_NJetsPt30Eta2p5;
   Int_t           t_NJetsPt30Eta5p0;
   Int_t           t_NJetsPt50Eta2p5;
   Int_t           t_NJetsPt50Eta5p0;
   Double_t        t_PFht;
   Double_t        t_PFmht;
   vector<double>  *t_genJetPt;
   vector<double>  *t_genJetEta;
   vector<double>  *t_genJetPhi;
   vector<double>  *t_genJetE;
   Int_t           t_beamHaloFilter;
   Int_t           t_eeBadScFilter;
   Int_t           t_eeNoiseFilter;
   Int_t           t_greedyMuons;
   Int_t           t_hcalLaserEventFilter;
   Int_t           t_inconsistentMuons;
   Int_t           t_ra2EcalBEFilter;
   Int_t           t_ra2EcalTPFilter;
   Int_t           t_trackingFailureFilter;
   Int_t           t_HBHENoiseFilterRA2;
   Int_t           t_ecalLaserCorrFilter;
   vector<string>  *t_firedTrigs;
   vector<double>  *t_firedTrigsPrescale;
   vector<double>  *t_genParPt;
   vector<double>  *t_genParEta;
   vector<double>  *t_genParPhi;
   vector<double>  *t_genParE;
   vector<double>  *t_genParStatus;
   vector<double>  *t_genParID;

   // List of branches
   TBranch        *b_t_EvtRun;   //!
   TBranch        *b_t_EvtLS;   //!
   TBranch        *b_t_EvtEvent;   //!
   TBranch        *b_t_NVertices;   //!
   TBranch        *b_t_tru_Npv;   //!
   TBranch        *b_t_avg_Npv;   //!
   TBranch        *b_t_MCflag;   //!
   TBranch        *b_t_PUWeight;   //!
   TBranch        *b_t_PUWeightAB;   //!
   TBranch        *b_t_PUWeightABC;   //!
   TBranch        *b_t_PUWeightRA2;   //!
   TBranch        *b_t_EvtWeight;   //!
   TBranch        *b_t_PFMetPx;   //!
   TBranch        *b_t_PFMetPy;   //!
   TBranch        *b_t_PFJetPt;   //!
   TBranch        *b_t_PFJetEta;   //!
   TBranch        *b_t_PFJetPhi;   //!
   TBranch        *b_t_PFJetE;   //!
   TBranch        *b_t_PFJetBTag;   //!
   TBranch        *b_t_PFJetNHF;   //!
   TBranch        *b_t_PFJetEMF;   //!
   TBranch        *b_t_PFJetHOEne;   //!
   TBranch        *b_t_NJetsPt30Eta2p5;   //!
   TBranch        *b_t_NJetsPt30Eta5p0;   //!
   TBranch        *b_t_NJetsPt50Eta2p5;   //!
   TBranch        *b_t_NJetsPt50Eta5p0;   //!
   TBranch        *b_t_PFht;   //!
   TBranch        *b_t_PFmht;   //!
   TBranch        *b_t_genJetPt;   //!
   TBranch        *b_t_genJetEta;   //!
   TBranch        *b_t_genJetPhi;   //!
   TBranch        *b_t_genJetE;   //!
   TBranch        *b_t_beamHaloFilter;   //!
   TBranch        *b_t_eeBadScFilter;   //!
   TBranch        *b_t_eeNoiseFilter;   //!
   TBranch        *b_t_greedyMuons;   //!
   TBranch        *b_t_hcalLaserEventFilter;   //!
   TBranch        *b_t_inconsistentMuons;   //!
   TBranch        *b_t_ra2EcalBEFilter;   //!
   TBranch        *b_t_ra2EcalTPFilter;   //!
   TBranch        *b_t_trackingFailureFilter;   //!
   TBranch        *b_t_HBHENoiseFilterRA2;   //!
   TBranch        *b_t_ecalLaserCorrFilter;   //!
   TBranch        *b_t_firedTrigs;   //!
   TBranch        *b_t_firedTrigsPrescale;   //!

   NtupleSelector(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~NtupleSelector() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

	//When this class is created using TTree->MakeSelector(), following line is enables.
	//But it seems to cause trouble with the current Makefile setup. 
   //ClassDef(NtupleSelector,0);
};

#endif

//This must match with what is in the cc file.
#ifdef NtupleSelector_cxx
void NtupleSelector::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   t_PFJetPt = 0;
   t_PFJetEta = 0;
   t_PFJetPhi = 0;
   t_PFJetE = 0;
   t_PFJetBTag = 0;
   t_PFJetNHF = 0;
   t_PFJetEMF = 0;
   t_PFJetHOEne = 0;
   t_genJetPt = 0;
   t_genJetEta = 0;
   t_genJetPhi = 0;
   t_genJetE = 0;
   t_firedTrigs = 0;
   t_firedTrigsPrescale = 0;
	t_genParPt  = 0; 
	t_genParEta = 0;
	t_genParPhi = 0;
	t_genParE   = 0;
	t_genParStatus  = 0;
	t_genParID  = 0;

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("t_EvtRun", &t_EvtRun, &b_t_EvtRun);
   fChain->SetBranchAddress("t_EvtLS", &t_EvtLS, &b_t_EvtLS);
   fChain->SetBranchAddress("t_EvtEvent", &t_EvtEvent, &b_t_EvtEvent);
   fChain->SetBranchAddress("t_NVertices", &t_NVertices, &b_t_NVertices);
   fChain->SetBranchAddress("t_tru_Npv", &t_tru_Npv, &b_t_tru_Npv);
   fChain->SetBranchAddress("t_avg_Npv", &t_avg_Npv, &b_t_avg_Npv);
   fChain->SetBranchAddress("t_MCflag", &t_MCflag, &b_t_MCflag);
   fChain->SetBranchAddress("t_PUWeight", &t_PUWeight, &b_t_PUWeight);
   fChain->SetBranchAddress("t_PUWeightAB", &t_PUWeightAB, &b_t_PUWeightAB);
   fChain->SetBranchAddress("t_PUWeightABC", &t_PUWeightABC, &b_t_PUWeightABC);
   fChain->SetBranchAddress("t_PUWeightRA2", &t_PUWeightRA2, &b_t_PUWeightRA2);
   fChain->SetBranchAddress("t_EvtWeight", &t_EvtWeight, &b_t_EvtWeight);
   fChain->SetBranchAddress("t_PFMetPx", &t_PFMetPx, &b_t_PFMetPx);
   fChain->SetBranchAddress("t_PFMetPy", &t_PFMetPy, &b_t_PFMetPy);
   fChain->SetBranchAddress("t_PFJetPt", &t_PFJetPt, &b_t_PFJetPt);
   fChain->SetBranchAddress("t_PFJetEta", &t_PFJetEta, &b_t_PFJetEta);
   fChain->SetBranchAddress("t_PFJetPhi", &t_PFJetPhi, &b_t_PFJetPhi);
   fChain->SetBranchAddress("t_PFJetE", &t_PFJetE, &b_t_PFJetE);
   fChain->SetBranchAddress("t_PFJetBTag", &t_PFJetBTag, &b_t_PFJetBTag);
   fChain->SetBranchAddress("t_PFJetNHF", &t_PFJetNHF, &b_t_PFJetNHF);
   fChain->SetBranchAddress("t_PFJetEMF", &t_PFJetEMF, &b_t_PFJetEMF);
   fChain->SetBranchAddress("t_PFJetHOEne", &t_PFJetHOEne, &b_t_PFJetHOEne);
   fChain->SetBranchAddress("t_NJetsPt30Eta2p5", &t_NJetsPt30Eta2p5, &b_t_NJetsPt30Eta2p5);
   fChain->SetBranchAddress("t_NJetsPt30Eta5p0", &t_NJetsPt30Eta5p0, &b_t_NJetsPt30Eta5p0);
   fChain->SetBranchAddress("t_NJetsPt50Eta2p5", &t_NJetsPt50Eta2p5, &b_t_NJetsPt50Eta2p5);
   fChain->SetBranchAddress("t_NJetsPt50Eta5p0", &t_NJetsPt50Eta5p0, &b_t_NJetsPt50Eta5p0);
   fChain->SetBranchAddress("t_PFht", &t_PFht, &b_t_PFht);
   fChain->SetBranchAddress("t_PFmht", &t_PFmht, &b_t_PFmht);
   fChain->SetBranchAddress("t_genJetPt", &t_genJetPt, &b_t_genJetPt);
   fChain->SetBranchAddress("t_genJetEta", &t_genJetEta, &b_t_genJetEta);
   fChain->SetBranchAddress("t_genJetPhi", &t_genJetPhi, &b_t_genJetPhi);
   fChain->SetBranchAddress("t_genJetE", &t_genJetE, &b_t_genJetE);
   fChain->SetBranchAddress("t_beamHaloFilter", &t_beamHaloFilter, &b_t_beamHaloFilter);
   fChain->SetBranchAddress("t_eeBadScFilter", &t_eeBadScFilter, &b_t_eeBadScFilter);
   fChain->SetBranchAddress("t_eeNoiseFilter", &t_eeNoiseFilter, &b_t_eeNoiseFilter);
   fChain->SetBranchAddress("t_greedyMuons", &t_greedyMuons, &b_t_greedyMuons);
   fChain->SetBranchAddress("t_hcalLaserEventFilter", &t_hcalLaserEventFilter, &b_t_hcalLaserEventFilter);
   fChain->SetBranchAddress("t_inconsistentMuons", &t_inconsistentMuons, &b_t_inconsistentMuons);
   fChain->SetBranchAddress("t_ra2EcalBEFilter", &t_ra2EcalBEFilter, &b_t_ra2EcalBEFilter);
   fChain->SetBranchAddress("t_ra2EcalTPFilter", &t_ra2EcalTPFilter, &b_t_ra2EcalTPFilter);
   fChain->SetBranchAddress("t_trackingFailureFilter", &t_trackingFailureFilter, &b_t_trackingFailureFilter);
   fChain->SetBranchAddress("t_HBHENoiseFilterRA2", &t_HBHENoiseFilterRA2, &b_t_HBHENoiseFilterRA2);
	fChain->SetBranchAddress("t_ecalLaserCorrFilter", &t_ecalLaserCorrFilter, &b_t_ecalLaserCorrFilter);
   fChain->SetBranchAddress("t_firedTrigs", &t_firedTrigs, &b_t_firedTrigs);
   fChain->SetBranchAddress("t_firedTrigsPrescale", &t_firedTrigsPrescale, &b_t_firedTrigsPrescale);
}

Bool_t NtupleSelector::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef NtupleSelector_cxx

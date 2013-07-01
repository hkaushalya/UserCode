#define SplitTree_cxx

#include <iostream>
#include <vector>
#include <cstring>
#include <string>
#include <map>
#include "SplitTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "Utils.hh"
#include <algorithm>
#include <iomanip>
#include "TStyle.h"
#include <cctype>
#include "TColor.h"
#include "TSystem.h"
#include <fstream>
#include <assert.h>

using namespace std;

/***************************************************************
 * Description: This is to split my LostLeptonTree flat ntuples
 * into much smaller files to decrease time use in the farm
 * when smearing.
 * Author: Sam Hewamanage
 * Institution: Florida Internationa University, USA.
 **************************************************************/

void PrintVector(const vector<double>& v)
{
	for (unsigned i=0; i < v.size(); ++i)
	{
		cout << v.at(i);
		if (i+1<v.size()) cout << ",";
		else cout << endl;
	}
}


int main(int argc, char* argv[])
{
	//cout << __FUNCTION__ << ": number of arguments = " << argc << endl;
	
	if (argc <= 3) {
		cerr <<"Please give 3 arguments: inputFileList outputFileName evtsToProcess" << endl;
		cerr <<" Valid configurations are " << std::endl;
		cerr << " ./runsplit inputFileList outputFileName N p" << std::endl;
		cerr << "  inputFileList : a text file with ROOT files with complete paths" << endl;
		cerr << "  outputFileName: file name to be use for output files (withou the extension). '.root' extension will be added automatically" << endl;
		cerr << "       N        : Total number of events to process. -1 to process all events in the inputFileList" << endl;
		cerr << "       p        : Splitting value. Number of events per output file to be produced. Default is 1000 events per file." << endl;
		cerr << "Eg:  ./runsplit filelist.txt out 10000 100" << std::endl;
		return -1;
	}

	const char *inputFileList = argv[1];
	const char *outFileName   = argv[2];
	const char *evt2Proc      = argv[3];
		
	int evts    = -1;

	if (isdigit(evt2Proc[0]))
	{
		int num = atoi(argv[3]);
		if (num>=1) evts = num;
		//cout << __FUNCTION__ << ": evts = " << evts << endl;
	} else {
			cout << "argument 4 is not a number. using default value for evts = " << evts << endl;
	}
	
	int evtsPerFile = 1000;
	if (argc>4)
	{
		//cout << __FUNCTION__ << ": processing arg4 ..." << endl;
		const char *g4 = argv[4];
		if (isdigit(g4[0])) evtsPerFile = atoi(g4);
		else {
			cout << "argument 4 is not a number. using default value for evtsPerFile = " << evtsPerFile << endl;
		}
	}

	SplitTree smear(inputFileList);
	smear.EventLoop(outFileName, evts, evtsPerFile);

	return 0;
}

void SplitTree::EventLoop( 
									const string outFileName,
									const int evts2Process,
									const unsigned evtsPerFile
									)
{
	if (fChain == 0) return;


	cout << "Inputs:" << endl;
	cout << "outFileName  = " << outFileName <<endl;
	cout << "evts2Process = " << evts2Process << endl;
	cout << "evtsPerFile  = " << evtsPerFile << endl;

	/**** events ti process ***********************/
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t entriesFound = nentries;

	if (evts2Process>=1 && evts2Process < nentries) 
	{	
		nentries = evts2Process;
		cout << "Requested events to process = " << nentries << endl;
	}



	Long64_t nbytes = 0, nb = 0;
	int decade = 0;
	bool bDEBUG = 0;
	unsigned nProcessed = 0;
	
		
	TFile *curOutFile = 0;
	TTree *curOutTree = 0;


	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		if (jentry ==0 || jentry%evtsPerFile==0) 
		{
			cout << "j = " << jentry << endl;
			if (jentry!=0) 
			{
				CloseFile(curOutFile);
			}
			//NewOutputFile(outFileName, curOutFile, curOutTree);

			++fileCounter;
			stringstream name;
			name << outFileName << "_" << fileCounter << ".root";
			curOutFile = new TFile (name.str().c_str(),"RECREATE");
			//	fnew->Print();
			TDirectory *dir = curOutFile->mkdir("treeMaker");
			dir->cd();

			curOutTree = new TTree("tree", "tree");
			curOutTree->SetAutoSave(10000);

			curOutTree->Branch("t_EvtRun",    &t_EvtRun,   "t_EvtRun/i");
			curOutTree->Branch("t_EvtLS",     &t_EvtLS,    "t_EvtLS/i");
			curOutTree->Branch("t_EvtEvent",  &t_EvtEvent, "t_EvtEvent/i");
			curOutTree->Branch("t_NVertices", &t_NVertices,"t_NVertices/I");
			curOutTree->Branch("t_tru_Npv",   &t_tru_Npv,  "t_tru_Npv/I");
			curOutTree->Branch("t_avg_Npv",   &t_avg_Npv,  "t_avg_Npv/I");
			curOutTree->Branch("t_MCflag",    &t_MCflag,   "t_MCflag/I");
			curOutTree->Branch("t_PUWeight",  &t_PUWeight, "t_PUWeight/D");
			curOutTree->Branch("t_PUWeightAB",   &t_PUWeightAB,  "t_PUWeightAB/D");
			curOutTree->Branch("t_PUWeightABC",  &t_PUWeightABC, "t_PUWeightABC/D");
			curOutTree->Branch("t_PUWeightRA2",  &t_PUWeightRA2, "t_PUWeightRA2/D");
			curOutTree->Branch("t_EvtWeight", &t_EvtWeight, "t_EvtWeight/D");
			curOutTree->Branch("t_PFMetPx",   &t_PFMetPx,  "t_PFMetPx/D");
			curOutTree->Branch("t_PFMetPy",   &t_PFMetPy,  "t_PFMetPy/D");

			/*t_PFJetPt  = new std::vector<double>();
			t_PFJetEta = new std::vector<double>();    
			t_PFJetPhi = new std::vector<double>();
			t_PFJetE   = new std::vector<double>();
			t_PFJetBTag  = new std::vector<double>();
			t_PFJetNHF   = new std::vector<double>();
			t_PFJetEMF   = new std::vector<double>();
			t_PFJetHOEne = new std::vector<double>();
			*/
			curOutTree->Branch("t_PFJetPt",  "vector<double>", &t_PFJetPt );
			curOutTree->Branch("t_PFJetEta", "vector<double>", &t_PFJetEta);
			curOutTree->Branch("t_PFJetPhi", "vector<double>", &t_PFJetPhi);
			curOutTree->Branch("t_PFJetE",   "vector<double>", &t_PFJetE);
			curOutTree->Branch("t_PFJetBTag",  "vector<double>", &t_PFJetBTag);
			curOutTree->Branch("t_PFJetNHF",   "vector<double>", &t_PFJetNHF);
			curOutTree->Branch("t_PFJetEMF",   "vector<double>", &t_PFJetEMF);
			curOutTree->Branch("t_PFJetHOEne", "vector<double>", &t_PFJetHOEne);
			curOutTree->Branch("t_NJetsPt30Eta2p5", &t_NJetsPt30Eta2p5, "t_NJetsPt30Eta2p5/I");
			curOutTree->Branch("t_NJetsPt30Eta5p0", &t_NJetsPt30Eta5p0, "t_NJetsPt30Eta5p0/I"); 
			curOutTree->Branch("t_NJetsPt50Eta2p5", &t_NJetsPt50Eta2p5, "t_NJetsPt50Eta2p5/I"); 
			curOutTree->Branch("t_NJetsPt50Eta5p0", &t_NJetsPt50Eta5p0, "t_NJetsPt50Eta5p0/I");

			curOutTree->Branch("t_PFht",      &t_PFht,      "t_PFht/D");
			curOutTree->Branch("t_PFmht",     &t_PFmht,     "t_PFmht/D");

			/*t_genJetPt  = new std::vector<double>();
			t_genJetEta = new std::vector<double>();    
			t_genJetPhi = new std::vector<double>();
			t_genJetE   = new std::vector<double>();
			*/
			curOutTree->Branch("t_genJetPt",  "vector<double>", &t_genJetPt );
			curOutTree->Branch("t_genJetEta", "vector<double>", &t_genJetEta);
			curOutTree->Branch("t_genJetPhi", "vector<double>", &t_genJetPhi);
			curOutTree->Branch("t_genJetE",   "vector<double>", &t_genJetE  );

			/*	t_genParPt  = new std::vector<double>();
				t_genParEta = new std::vector<double>();    
				t_genParPhi = new std::vector<double>();
				t_genParE   = new std::vector<double>();
				t_genParStatus  = new std::vector<double>();
				t_genParID  = new std::vector<double>();
				*/
			/*	curOutTree->Branch("t_genParPt",  "vector<double>", &t_genParPt );
				curOutTree->Branch("t_genParEta", "vector<double>", &t_genParEta);
				curOutTree->Branch("t_genParPhi", "vector<double>", &t_genParPhi);
				curOutTree->Branch("t_genParE",   "vector<double>", &t_genParE  );
				curOutTree->Branch("t_genParStatus",  "vector<double>", &t_genParStatus );
				curOutTree->Branch("t_genParID",  "vector<double>", &t_genParID );
				*/

			curOutTree->Branch("t_beamHaloFilter", &t_beamHaloFilter, "t_beamHaloFilter/I");
			curOutTree->Branch("t_eeBadScFilter", &t_eeBadScFilter, "t_eeBadScFilter/I");
			curOutTree->Branch("t_eeNoiseFilter", &t_eeNoiseFilter, "t_eeNoiseFilter/I");
			curOutTree->Branch("t_greedyMuons", &t_greedyMuons, "t_greedyMuons/I");
			curOutTree->Branch("t_hcalLaserEventFilter", &t_hcalLaserEventFilter, "t_hcalLaserEventFilter/I");
			curOutTree->Branch("t_inconsistentMuons", &t_inconsistentMuons, "t_inconsistentMuons/I");
			curOutTree->Branch("t_ra2EcalBEFilter", &t_ra2EcalBEFilter, "t_ra2EcalBEFilter/I");
			curOutTree->Branch("t_ra2EcalTPFilter", &t_ra2EcalTPFilter, "t_ra2EcalTPFilter/I");
			curOutTree->Branch("t_trackingFailureFilter", &t_trackingFailureFilter, "t_trackingFailureFilter/I");
			curOutTree->Branch("t_HBHENoiseFilterRA2", &t_HBHENoiseFilterRA2, "t_HBHENoiseFilterRA2/I");
			curOutTree->Branch("t_ecalLaserCorrFilter", &t_ecalLaserCorrFilter, "t_ecalLaserCorrFilter/I");

			//t_firedTrigs  = new std::vector<string>();
			curOutTree->Branch("t_firedTrigs",  "vector<string>", &t_firedTrigs );
			//t_firedTrigsPrescale   = new std::vector<double>();
			curOutTree->Branch("t_firedTrigsPrescale",  "vector<double>", &t_firedTrigsPrescale );



			bool stop = false;
			if (curOutTree == 0)
			{
				cout << "new tree not found!" << endl;
				stop = true;
			}
			if (curOutFile == 0)
			{
				cout << "new file not found!" << endl;
				stop = true;
			}
			if (stop)  assert(false);
			cout << "Opened new output file " << curOutFile->GetName() << endl;
		}


		// ==============print number of events done == == == == == == == =
		double progress = 10.0 * jentry / (1.0 * nentries);
		int k = int (progress);
		if (k > decade)
		{
			cout << green <<  10 * k << " % (" << jentry << ")" << clearatt << endl;
			//printProgBar( (int)100.0 * jentry/ (double) nentries);
		}
		decade = k;
		
		// ===============read this entry == == == == == == == == == == == 
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;

		nb = fChain->GetEntry(jentry);   nbytes += nb;

		++nProcessed;
		curOutTree->Fill();


	} // event loop

	CloseFile (curOutFile);


	/***********************************************
	 *  End job summary
	 **********************************************/

	cout << ">>>>>>> " << __FILE__ << ":" << __FUNCTION__ << ": End Job " << endl;
	if (bDEBUG) cout << red << " ------ DEBUG MODE ---------- " << clearatt << endl;
	cout << "Entries found/proces'd = " << entriesFound << " / " << nProcessed << endl;
}


void SplitTree::DumpJets(const vector<TLorentzVector>& jets, const double& minPt) const
{
	cout << setw(3) << "#" 
			<< setw(15) << "Pt>" << minPt 
			<< setw(15) << "Eta"  
			<< setw(15) << "Phi" 
			<< setw(15) << "M" 
			<< endl;  
	for (unsigned i=0; i < jets.size(); i++)
	{
		if (jets.at(i).Pt()>minPt)
		{
			cout << setprecision(1)  << fixed
				<< setw(3) << i 
				<< setw(15) << jets.at(i).Pt() 
				<< setw(15) << jets.at(i).Eta() 
				<< setw(15) << jets.at(i).Phi() 
				<< setw(15) << jets.at(i).M() 
				<< endl;  
		}
	}
}


void SplitTree::DumpJetsAndCSV(const vector<TLorentzVector>& jets, const vector<double>& csv, 
							const TLorentzVector& met) const
{
	cout << setw(3) << "#" 
			<< setw(15) << "Pt>0" 
			<< setw(15) << "Eta"  
			<< setw(15) << "Phi" 
			<< setw(15) << "M" 
			<< setw(15) << "CSV" 
			<< setw(15) << "dphi(met)" 
			<< endl;  
	for (unsigned i=0; i < jets.size(); i++)
	{
		if (jets.at(i).Pt()>0.0)
		{
			cout << setprecision(1)  << fixed
				<< setw(3) << i 
				<< setw(15) << jets.at(i).Pt() 
				<< setw(15) << jets.at(i).Eta() 
				<< setw(15) << jets.at(i).Phi() 
				<< setw(15) << jets.at(i).M() 
				<< setw(15) << csv.at(i) 
				<< setw(15) << fabs(jets.at(i).DeltaPhi(met))
				<< endl;  
		}
	}
}

void SplitTree::DumpJet(const TLorentzVector& jet) const
{
	cout << " jet pt/eta/phi/m=" <<setprecision(1)  << fixed
		<< setw(15) << jet.Pt() 
		<< setw(15) << jet.Eta() 
		<< setw(15) << jet.Phi() 
		<< setw(15) << jet.M() 
		<< endl;  
}

void SplitTree::PrintEventNumber() const
{
	cout << "Run/Ls/Evt = " << t_EvtRun << " / " 
		<<  t_EvtLS << " / " <<  t_EvtEvent << endl;
}

/********************************************************************
 *				C O N S T R U C T O R
 *******************************************************************/
SplitTree::SplitTree(
				const TString &inputFileList 
				) {

	TChain *tree = new TChain("treeMaker/tree");  

	if ( ! FillChain(tree, inputFileList) ) {
		std::cerr << "Cannot get the tree " << std::endl;
		assert(false);
	}

	Init(tree);
}

void SplitTree::CloseFile(TFile* f)
{
	if (f != NULL)
	{
		cout << "Closing file "<< f->GetName() << endl;
		f->cd();
		f->Write();
		f->Close();
	}
	delete f;
}


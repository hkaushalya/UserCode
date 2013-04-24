#define FactorizationBySmearing_cxx

#include <iostream>
#include <vector>
#include <cstring>
#include <string>
#include <map>
#include "FactorizationBySmearing.h"
#include "TGraph.h"
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
//#include "Basic_Mt2_332_Calculator.h"
#include "ChengHanBisect_Mt2_332_Calculator.h"
//#include "combination.h"
//#include "indexSort.h"

using namespace std;

/***************************************************************
 * This class is created for Factorization to utilize gen-jet
 * smearing to better quantify constant 'c'. This works only 
 * flat ntuples created using a modifed version of
 * lostLeptonTree.
 * Author: Sam Hewamanage
 * Institution: Florida Internationa University, USA.
 **************************************************************/

/*sorts two related vector*/
void domysort(vector<TLorentzVector>& vjets, vector<double>& bDisc)
{
	bool done = false;
	unsigned i =0;
	if (vjets.size()>1)
	{
		unsigned loop  = 1;
		do {
			done = true;
			for (i=0; i < (vjets.size()-1); ++i)
			{
				if (vjets.at(i).Pt()<vjets.at(i+1).Pt())
				{
					std::swap(vjets.at(i), vjets.at(i+1));
					std::swap(bDisc.at(i), bDisc.at(i+1));
					done = false & done;
				} 
			}
			++loop;
		} while ( ! done);

		for (i=0; i<(vjets.size()-1); ++i)
		{
			if (vjets.at(i).Pt()<vjets.at(i+1).Pt())
			{
				cout << __FUNCTION__ << ": sort failed! i/pti/pti+1= "<< i << "/" << vjets.at(i).Pt() << "/" << vjets.at(i+1).Pt() << endl;
				assert(false);
			} else {
				//	cout << __FUNCTION__ << ": i/pt = " << i << "-> " << vjets.at(i).Pt() << endl;
			}
		}
	}
}


static bool sort_using_less_than(double u, double v)
{
	   return u < v;
}

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
		cerr <<"Please give 3 arguments " << "runList " << " " << "outputFileName" << " " << "# of evts  nJet50Min  nJet50Max" << endl;
		cerr <<" Valid configurations are " << std::endl;
		cerr << " ./optimize filelist.txt out.root nevts2process[-1 = all] smearingSyst[0=mean, 1-6 systs]" << std::endl;
		cerr << "Eg:  ./optimize filelist.txt out.root 100" << std::endl;
		return -1;
	}

	const char *inputFileList = argv[1];
	const char *outFileName   = argv[2];
	const char *evt2Proc      = argv[3];
		
	int evts    = -1;
	int systematic_var = 0;
	unsigned cutmask = 127;  //enable all cuts

	if (isdigit(evt2Proc[0]))
	{
		int num = atoi(argv[3]);
		if (num>=1) evts = num;
		//cout << __FUNCTION__ << ": evts = " << evts << endl;
	} else {
			cout << "argument 4 is not a number. using default value for evts = " << evts << endl;
	}
	
	if (argc>4)
	{
		//cout << __FUNCTION__ << ": processing arg4 ..." << endl;
		const char *g4 = argv[4];
		if (isdigit(g4[0])) systematic_var = atoi(g4);
		else {
			cout << "argument 4 is not a number. using default value for systematic_var = " << systematic_var << endl;
		}
	}

	if (argc>5)
	{
		const char *g5 = argv[5];
		if (isdigit(g5[0])) cutmask = atoi(g5);
		else {
			cout << "argument 5 is not a number. using default value for cutmask = " << cutmask << endl;
		}
	}

	//cout << "systematic_var = " << systematic_var << endl;
	FactorizationBySmearing smear(inputFileList, outFileName);

	smear.EventLoop(inputFileList, evts, cutmask, systematic_var);

	return 0;
}

void FactorizationBySmearing::EventLoop(const char *datasetname, 
									const int evts2Process,
									const unsigned cutmask,
									const int systematicVarition
									)
{
	if (fChain == 0) return;


	bAPPLY_NJET_CUT        =  1 & cutmask;
	bAPPLY_MET_CUT         =  2 & cutmask;
	bAPPLY_TRIPLET_CUT     =  4 & cutmask;
	bAPPLY_TOPMASS_CUT     =  8 & cutmask;
	bAPPLY_TOPPLUSBJET_CUT = 16 & cutmask;
	bAPPLY_MT2_CUT         = 32 & cutmask;
	bAPPLY_DPHI_CUT        = 64 & cutmask;
	
	uMinNjet70Eta2p4_ = 2; uMinNjet50Eta2p4_ = 4; uMinNjet30Eta2p4_ = 5;
	dMinMet_ = 175.0;
	uMinTriplets_ = 1;
	dMinTopMass_ = 80.0; dMaxTopMass_ = 270.0;
	dMinTopPlusBjetMass_ = 500.0;
	dMinMt2_ = 300.0;


	/*****************************************************
	 * Smearing function constants and variables
	 ******************************************************/
	smearedJetPt_        = 13.0;


	/**** events ti process ***********************/
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t entriesFound = nentries;

	if (evts2Process>=1 && evts2Process < nentries) 
	{	
		nentries = evts2Process;
		cout << "Requested events to process = " << nentries << endl;
	}

	const double dDATA_LUMI = 19458.0; //used in MC lumi weighing
	bNON_STD_MODE        = 0;  //
	sNON_STD_MODE_EXPLAIN = "Excluding trig wgt>1000 by setting it to 1!";
	bRUNNING_ON_MC       = 1; 
	bDO_TRIG_SELECTION   = 0;   //DATA
	bDO_TRIG_PRESCALING  = 0;   //DATA
	//bAPPLY_DPHI_CUT      = 0; 
	bDO_PU_WEIGHING      = false;  //MC
	bDO_LUMI_WEIGHING    = false;  //MC 
	bDO_GENJET_SMEARING  = 1;  //MC
	uNTRIES              = 1; //number of pseudo experiments per event
	bDEBUG               = false;


	//* overide above with for std run settings for ease
	const bool bUseStdDataSettings = 0;
	const bool bUseStdMCSettings   = 0;
 	
	/*if (bUseStdMCSettings) { bRUNNING_ON_MC = 1; bDO_PU_WEIGHING = 1; bDO_LUMI_WEIGHING = 0; bDO_GENJET_SMEARING = 0; }
	else if (bUseStdDataSettings) { bRUNNING_ON_MC = 0; bDO_TRIG_SELECTION = 1; bDO_TRIG_PRESCALING = 0; }
	*/

	//sanity check
	if (bRUNNING_ON_MC && (bDO_TRIG_PRESCALING || bDO_TRIG_SELECTION) )
	{
		cout << __FUNCTION__ << ": Contradicting settings. Can't run on MC with TriggerSelection or TriggerPrescaling or vise-versa!! (bRUNNING_ON_MC = " 
						<< bRUNNING_ON_MC << ", bDO_TRIG_PRESCALING = " << bDO_TRIG_PRESCALING << ", bDO_TRIG_SELECTION = " << bDO_TRIG_SELECTION << endl;
		return;
	}

	if (bDO_TRIG_PRESCALING && ! bDO_TRIG_SELECTION)
	{
		cout << __FUNCTION__ << ": Contradicting settings. Require trigger selection inorder to get trigger prescaling! " 
						<< " bDO_TRIG_PRESCALING = " << bDO_TRIG_PRESCALING <<  " / bDO_TRIG_SELECTION = " << bDO_TRIG_SELECTION << endl;
		return;
	}

	if (! bRUNNING_ON_MC && bDO_PU_WEIGHING)
	{
		cout << __FUNCTION__ << ": Contradicting settings. PU Weighing can be applied only with MC.! " << endl;
		return;
	}

	if (! bRUNNING_ON_MC && bDO_LUMI_WEIGHING)
	{
		cout << __FUNCTION__ << ": Contradicting settings. Lumi Weighing can be applied only with MC.! " << endl;
		return;
	}


	double dMC_lumiWgt = 1.0; 
	if (bRUNNING_ON_MC && bDO_LUMI_WEIGHING)
	{
		dMC_lumiWgt = GetLumiWgt(datasetname, dDATA_LUMI);
	}
	const double dMC_LUMI_WGT = dMC_lumiWgt; 

	if (bDEBUG) uNTRIES  = 1;

	const double smearingWgt = 1.0/(double)uNTRIES;

	std::cout << red << "Dataset = " << datasetname << clearatt << endl;
	if (bRUNNING_ON_MC) std::cout << red << "Using dMC_LUMI_WGT / smearingWgt = " << dMC_LUMI_WGT << " / " 
				<< smearingWgt << clearatt << std::endl;
	
	//book histograms
	BookHistogram(outRootFile, bRUNNING_ON_MC);

	Long64_t nbytes = 0, nb = 0;
	int decade = 0;
	
	if (bRUNNING_ON_MC && bDO_GENJET_SMEARING)
	{
		smearFunc_ = new SmearFunction();
		//Pythia
		//smearFunc_->SetSmearingFile("/share/store/users/samantha/CMSSW_5_2_5/src/UserCode/RA2QCDvetoAna/flatrun/input_files/MCJetResolutions_Summer12_DR53X_QCD_Pt_15to3000_TuneZ2star_Flat_8TeV_pythia6_withCHS_withoutPUReweighting.root");
		//MG
		smearFunc_->SetSmearingFile("/share/store/users/samantha/CMSSW_5_2_5/src/UserCode/RA2QCDvetoAna/flatrun/input_files/MCJetResolutions_Summer12_DR53X_QCD_HT_100ToInf_TuneZ2star_8TeV_madgraph_withCHS_withoutPUReweighting.root");

		//according to Kristin's config file this is always false
//		if (systematicVarition>0)
//		{
			absoluteTailScaling_ = false;  //false for systematics
//		} else {
//			absoluteTailScaling_ = true; 
//		}
		smearFunc_->SetAbsoluteTailScaling(absoluteTailScaling_);

//		const float sysVarScale = 2.0; //this is the default 
		//const float sysVarScale = 5.0;  //this is what is used in 2010 note
//		const float systvarLowFactor = 1.0 / sysVarScale;
//		const float systvarHiFactor  = 1.0 * sysVarScale;


		//these are the dafault (or mean)
		/*const double ptBinEdges[] = {0, 220, 270, 300, 350, 500, 2200};
		const double etaBinEdges[]= {0, 0.5, 1.1,  1.7, 2.3, 5.0};
		const double tailScales[] = {0.953,1.418,1.156,1.305,1.342,1.353,1.096,1.083,1.083,1.195,1.248,1.248,0.965,1.035,1.035,1.035,1.358,1.358,0.938,1.196,1.196,1.196,1.196,1.196,1.069,1.069,1.069,1.069,1.069,1.069};
		const double additionalSmearing[] = {1.052,1.052,1.052,1.052,1.052,1.052,1.057,1.057,1.057,1.057,1.057,1.057,1.096,1.096,1.096,1.096,1.096,1.096,1.134,1.134,1.134,1.134,1.134,1.134,1.288,1.288,1.288,1.288,1.288,1.288};
*/
		//debugging
		const double ptBinEdges[] = {0,7000};
		const double etaBinEdges[] = {0.0, 5.0};
		const double tailScales[] = {1.0};
		const double additionalSmearing[] = {1.0};


		vector<double> vPtBinEdges(ptBinEdges, ptBinEdges + sizeof(ptBinEdges) / sizeof(double));
		vector<double> vEtaBinEdges(etaBinEdges, etaBinEdges + sizeof(etaBinEdges) / sizeof(double));
		vector<double> vTailScales(tailScales, tailScales + sizeof(tailScales) / sizeof(double));
		vector<double> vAdditionalSmearing(additionalSmearing, additionalSmearing + sizeof(additionalSmearing) / sizeof(double));

		smearFunc_->SetPtBinEdges_scaling(vPtBinEdges);
		smearFunc_->SetEtaBinEdges_scaling(vEtaBinEdges);
		smearFunc_->SetLowerTail_scaling(vTailScales);
		smearFunc_->SetUpperTail_scaling(vTailScales);
		smearFunc_->SetAdditionalSmear_scaling(vAdditionalSmearing);
		smearFunc_->UncertaintyName("Mean");


		//no vary scaling to do the systemtatics by resetting some
		//of the scalings.

		switch (systematicVarition)
		{
			case 1: //tailUP
			{
				//const double tailScales_1[] = {1.236,1.92,1.505,1.655,1.735,1.703,1.47,1.455,1.455,1.52,1.672,1.672,1.298,1.33,1.33,1.33,1.685,1.685,1.224,1.621,1.621,1.621,1.621,1.621,1.839,1.839,1.839,1.839,1.839,1.839};
				const double tailScales_1[] = {2.0};
				vector<double> vTailScales_1(tailScales_1, tailScales_1 + sizeof(tailScales_1) / sizeof(double));

				//smearFunc_->SetAbsoluteTailScaling(true);
				smearFunc_->SetLowerTail_scaling(vTailScales_1);
				//smearFunc_->SetUpperTail_scaling(vTailScales_1);
				//smearFunc_->UncertaintyName("tailUP");
				break;
			}

			case 2: //tailUP
			{
				//const double tailScales_1[] = {1.236,1.92,1.505,1.655,1.735,1.703,1.47,1.455,1.455,1.52,1.672,1.672,1.298,1.33,1.33,1.33,1.685,1.685,1.224,1.621,1.621,1.621,1.621,1.621,1.839,1.839,1.839,1.839,1.839,1.839};
				const double tailScales_1[] = {0.5};
				vector<double> vTailScales_1(tailScales_1, tailScales_1 + sizeof(tailScales_1) / sizeof(double));

				//smearFunc_->SetAbsoluteTailScaling(true);
				smearFunc_->SetLowerTail_scaling(vTailScales_1);
				//smearFunc_->SetUpperTail_scaling(vTailScales_1);
				//smearFunc_->UncertaintyName("tailUP");
				break;
			}


			case 3: //tailDOWN
			{
				//const double tailScales_1[] = {0.67,0.916,0.807,0.955,0.949,1.003,0.722,0.711,0.711,0.87,0.824,0.824,0.632,0.74,0.74,0.74,1.031,1.031,0.652,0.771,0.771,0.771,0.771,0.771,0.299,0.299,0.299,0.299,0.299,0.299};
				const double tailScales_1[] = {2.0};
				vector<double> vTailScales_1(tailScales_1, tailScales_1 + sizeof(tailScales_1) / sizeof(double));
				//smearFunc_->SetLowerTail_scaling(vTailScales_1);
				smearFunc_->SetUpperTail_scaling(vTailScales_1);
				//smearFunc_->UncertaintyName("tailDOWN");
				break;
			}

			case 4: //tailDOWN
			{
				//const double tailScales_1[] = {0.67,0.916,0.807,0.955,0.949,1.003,0.722,0.711,0.711,0.87,0.824,0.824,0.632,0.74,0.74,0.74,1.031,1.031,0.652,0.771,0.771,0.771,0.771,0.771,0.299,0.299,0.299,0.299,0.299,0.299};
				const double tailScales_1[] = {0.5};
				vector<double> vTailScales_1(tailScales_1, tailScales_1 + sizeof(tailScales_1) / sizeof(double));
				//smearFunc_->SetLowerTail_scaling(vTailScales_1);
				smearFunc_->SetUpperTail_scaling(vTailScales_1);
				//smearFunc_->UncertaintyName("tailDOWN");
				break;
			}

			case 5:  //coreUP
			{
				//const double additionalSmearing_1[] = {1.116,1.116,1.116,1.116,1.116,1.116,1.117,1.117,1.117,1.117,1.117,1.117,1.166,1.166,1.166,1.166,1.166,1.166,1.237,1.237,1.237,1.237,1.237,1.237,1.511,1.511,1.511,1.511,1.511,1.511}; 
				const double additionalSmearing_1[] = {1.10};
				vector<double> vAdditionalSmearing_1(additionalSmearing_1, additionalSmearing_1 + sizeof(additionalSmearing_1) / sizeof(double));
				smearFunc_->SetAdditionalSmear_scaling(vAdditionalSmearing_1);
				smearFunc_->UncertaintyName("coreUP");
				break;
			}

			case 6:  //coreDOWN
			{
				//const double additionalSmearing_1[] = {0.988,0.988,0.988,0.988,0.988,0.988,0.999,0.999,0.999,0.999,0.999,0.999,0.998,0.998,0.998,0.998,0.998,0.998,1.033,1.033,1.033,1.033,1.033,1.033,1.067,1.067,1.067,1.067,1.067,1.067};
				const double additionalSmearing_1[] = {0.9};
				vector<double> vAdditionalSmearing_1(additionalSmearing_1, additionalSmearing_1 + sizeof(additionalSmearing_1) / sizeof(double));
				smearFunc_->SetAdditionalSmear_scaling(vAdditionalSmearing_1);
				smearFunc_->UncertaintyName("coreDOWN");
				break;
			}
		}

		smearFunc_->Init();

		BookJerDebugHists();
	}

	if (! bRUNNING_ON_MC) LoadBadHcalLaserEvents();

	vector<TLorentzVector> recoJets, genJets, smearedGenJets;
	vector<double> bDiscrminators_reco, bDiscrminators_gen, bDiscrminators_smearedgen;

	unsigned nProcessed = 0, nCleaningFailed = 0, nBadHcalLaserEvts = 0;
	unsigned nHOfilterFailed = 0;
	nBadEcalLaserEvts = 0;
	nSmearedJetEvts = 0.0; nRecoJetEvts = 0.0; nGenJetEvts = 0.0;
	unsigned int topsFound = 0;


	for (Long64_t jentry=0; jentry<nentries;jentry++) {
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

		//to study the seed sample effects
		//if (bNON_STD_MODE && t_EvtEvent%2 == 0) continue;

		++nProcessed;

		if (bDEBUG) 
			cout << red << "===================== run/ls/evt= " << t_EvtRun << " / " 
				<<  t_EvtLS << " / " <<  t_EvtEvent << clearatt  << endl;

		/* preselection cuts */
		if (t_NVertices < 1) continue;  //at least one good vertex

		
		CreateRecoJetVec(recoJets, bDiscrminators_reco);
		if (! PassHOfilter())
		{
			++nHOfilterFailed;
			continue;
		}
		
		if (! PassCleaning()) 
		{
			++nCleaningFailed;
			if (bDEBUG) cout << "Failed cleaning." << endl;;
			continue;
		}
		
		//reject hcal laser bad events in data
		if (! bRUNNING_ON_MC) 
		{
			stringstream strevt;
			strevt << t_EvtRun << ":" << t_EvtLS << ":" << t_EvtEvent;
			if (find(vBadHcalLaserEvts.begin(), vBadHcalLaserEvts.end(), strevt.str()) != vBadHcalLaserEvts.end())
			{
				++nBadHcalLaserEvts;
				continue;
			}
		}

		const double evtWeight = t_EvtWeight; //==1 except for Flat QCD sample
		const double puWeight  = t_PUWeight;  //PU weight for MC
	//	cout << "PU weight = " << t_PUWeight << endl;

		/******************************************************
		 * For MC;   recoEvtTotWgt = dMC_LUMI_WGT * evtWeight
		 * For DATA; recoEvtTotWgt = trigPrescaleWeight
		 *****************************************************/
		double recoEvtTotWgt = 1;
		if (bRUNNING_ON_MC) 
		{
			recoEvtTotWgt = dMC_LUMI_WGT * evtWeight;
			if (bDO_PU_WEIGHING) recoEvtTotWgt *= puWeight;
		} else 
		{
			bool passTrigger = true;

			/* select on trigger */
			if (bDO_TRIG_SELECTION)
			{
				double prescaleWgt = 1.0;
				TrigPrescaleWeight(passTrigger, prescaleWgt); 
				if (! passTrigger) 
				{
					if (bDEBUG) cout << "Failed trigger " << endl;
					continue;   //when running on data accept event passing trigger only
				} else {
					if (bDEBUG) cout << "Passed trigger " << endl;
					//cout << "Passed Trigger:"; 
					PrintEventNumber();
				}

				/* apply trigger prescale weight */
				if (bDO_TRIG_PRESCALING) 
				{
					//NON STD MODE TO CHECK EFFECT OF THIS LARGE PRESCALED EVENTS ON THE FINAL PREDICTIONS 03262013
					if (prescaleWgt>500.0) 
					{
						prescaleWgt ==1.0;
					}
					recoEvtTotWgt = prescaleWgt;
				}
			}
		}

		/* This is to calculate uncl MET = pfMET - MHT */
		TVector3 pfmetv(t_PFMetPx, t_PFMetPy,0);
		TLorentzVector  pfmetvec(0,0,0,0);
		pfmetvec.SetVectMag(pfmetv,0);
		TLorentzVector reco_mhtvec(MHT(recoJets));	
		const TLorentzVector reco_uncl_met_vec = pfmetvec - reco_mhtvec;  

		const bool accept_reco_evt = FillHistogram(recoJets, bDiscrminators_reco, reco_uncl_met_vec, 0, recoEvtTotWgt);
		if (accept_reco_evt) 
		{
			nRecoJetEvts += recoEvtTotWgt;
		}

		if (bRUNNING_ON_MC)
		{
			CreateGenJetVec(genJets);
			SetGenJetBdiscriminators(recoJets, bDiscrminators_reco, genJets, bDiscrminators_gen); 
			const bool accept_gen_evt = FillHistogram(genJets, bDiscrminators_gen, reco_uncl_met_vec, 1, recoEvtTotWgt);
			if (accept_gen_evt) ++nGenJetEvts;

			if (bDO_GENJET_SMEARING)
			{
				for (unsigned n = 0; n < uNTRIES; ++n)
				{
					SmearingGenJets(genJets, bDiscrminators_gen, smearedGenJets, bDiscrminators_smearedgen);
					
					double totWeight = smearingWgt * dMC_LUMI_WGT * evtWeight;
					if (bDO_PU_WEIGHING) totWeight *= puWeight;
					const bool accept_smear_evt = FillHistogram(smearedGenJets, bDiscrminators_smearedgen, reco_uncl_met_vec, 2, totWeight);
					if (accept_smear_evt) {
						nSmearedJetEvts += totWeight;
					}
					if (accept_reco_evt && !accept_reco_evt)
					{
						cout << red << __LINE__ << ": Smeared event thrown out at ntries= " << n << clearatt <<endl;
					}
				}
			}
		}

	} // event loop


	/***** NORMALIZE VARIABLE BINNED HIST BY BIN WIDTH *****/
	/*
	TCanvas *pf = new TCanvas("pf");
	gPad->Print("passfail.eps[");
	for (unsigned jetbin=0; jetbin < JetBins_.size(); ++jetbin)
	{
		for (unsigned htbin=0; htbin < HtBins_.size() -1 ; ++htbin)
		{
			for (unsigned mhtbin=0; mhtbin < MhtBins_.size() -1 ; ++mhtbin)
			{
				DivideByBinWidth(Hist.at(jetbin).at(htbin).at(mhtbin).hv_SmearedEvt.signal);
				Hist.at(jetbin).at(htbin).at(mhtbin).hv_SmearedEvt.signal->SetMarkerStyle(20);
				Hist.at(jetbin).at(htbin).at(mhtbin).hv_SmearedEvt.signal->SetMarkerColor(kRed);
				TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
				leg->AddEntry(Hist.at(jetbin).at(htbin).at(mhtbin).hv_SmearedEvt.signal,"Pass RA2 #Delta #Phi cuts");
				TCanvas *c1 = new TCanvas("c1");
				gPad->SetLogy();
				Hist.at(jetbin).at(htbin).at(mhtbin).hv_SmearedEvt.signal->Draw();
				TCanvas *c2 = new TCanvas("c2");
				gPad->SetLogy();

				for (int i=0; i< vDphiVariations.size(); ++i)
				{
					DivideByBinWidth(Hist.at(jetbin).at(htbin).at(mhtbin).hv_SmearedEvt.pass.at(i));
					DivideByBinWidth(Hist.at(jetbin).at(htbin).at(mhtbin).hv_SmearedEvt.fail.at(i));
					string drawOption("");
					if (i!=0) drawOption +="same";
					c1->cd();
					Hist.at(jetbin).at(htbin).at(mhtbin).hv_SmearedEvt.fail.at(i)->SetMarkerStyle(22+i);
					Hist.at(jetbin).at(htbin).at(mhtbin).hv_SmearedEvt.fail.at(i)->SetMarkerColor(13+i);
					Hist.at(jetbin).at(htbin).at(mhtbin).hv_SmearedEvt.fail.at(i)->Draw("same");
					stringstream fail_leg_text;
					fail_leg_text << "#Delta #Phi_{min}<" << vDphiVariations.at(i);
					leg->AddEntry(Hist.at(jetbin).at(htbin).at(mhtbin).hv_SmearedEvt.fail.at(i), fail_leg_text.str().c_str());

					c2->cd();
					TH1* temp = dynamic_cast<TH1*> (Hist.at(jetbin).at(htbin).at(mhtbin).hv_SmearedEvt.signal->Clone("copy"));
					temp->SetDirectory(0);
					temp->SetTitle("Signal Region/ Control Region");
					temp->Divide(Hist.at(jetbin).at(htbin).at(mhtbin).hv_SmearedEvt.fail.at(i));
					temp->SetMarkerStyle(22+i);
					temp->SetMarkerColor(13+i);
					temp->DrawCopy(drawOption.c_str());

				}

				c1->cd();
				leg->Draw();
				gPad->Print("passfail.eps");
				c2->cd();
				leg->Draw();
				gPad->Print("passfail.eps");
				delete leg;
				delete c1;
				delete c2;

			}
		}
	}

	gPad->Print("passfail.eps]");
	delete pf;
	*/


	/****************************************************************************
	 * FOR QUICK DEBUGGING ONLY
	 ***************************************************************************/
	TCanvas *c = new TCanvas();
	gStyle->SetOptStat(0);
	gErrorIgnoreLevel = kWarning; //supress print messages
	gPad->Print("samples.eps[");

	for (unsigned jetbin=0; jetbin < JetBins_.size(); ++jetbin)
	{
		for (unsigned htbin=0; htbin < HtBins_.size() -1 ; ++htbin)
		{
			for (unsigned mhtbin=0; mhtbin < MhtBins_.size() -1 ; ++mhtbin)
			{
				TH1* hrecomht = (TH1*) Hist.at(jetbin).at(htbin).at(mhtbin).hv_RecoEvt.h_Mht->Clone("reco_mht");
				hrecomht->SetDirectory(0);
				TH1* hsmearmht = (TH1*) Hist.at(jetbin).at(htbin).at(mhtbin).hv_SmearedEvt.h_Mht->Clone("smear_mht");
				hsmearmht->SetDirectory(0);
				TH1* hrecodphimin = (TH1*) Hist.at(jetbin).at(htbin).at(mhtbin).hv_RecoEvt.h_DphiMin->Clone("reco_dphimin");
				hrecodphimin->SetDirectory(0);
				TH1* hsmeardphimin = (TH1*) Hist.at(jetbin).at(htbin).at(mhtbin).hv_SmearedEvt.h_DphiMin->Clone("smear_dphimin");
				hsmeardphimin->SetDirectory(0);

				hsmearmht->SetLineColor(kRed);
				hsmearmht->SetMarkerColor(kRed);
				hsmearmht->SetMarkerStyle(24);
				hsmearmht->SetLineWidth(2);
				hrecomht->SetLineWidth(2);
				hrecomht->SetMarkerStyle(kPlus);
				hrecomht->GetXaxis()->SetRangeUser(0,300);
				hsmearmht->GetXaxis()->SetRangeUser(0,300);

				stringstream recomhtleg,smearmhtleg;
				const double sum_recomht = hrecomht->Integral(1, hrecomht->GetNbinsX()+1);
				const double sum_smearmht = hsmearmht->Integral(1, hsmearmht->GetNbinsX()+1);
				recomhtleg << "reco (" << sum_recomht << ")";
				smearmhtleg << "smeared (" << sum_smearmht << ")";

				TLegend *l2 = new TLegend(0.6,0.6,0.9,0.9);
				l2->AddEntry(hrecomht, recomhtleg.str().c_str());
				l2->AddEntry(hsmearmht, smearmhtleg.str().c_str());

				hrecomht->DrawCopy();
				hsmearmht->DrawCopy("same");
				l2->Draw();
				gPad->Print("samples.eps");

				//njets
				TH1* hreconjet50 = (TH1*) Hist.at(jetbin).at(htbin).at(mhtbin).hv_RecoEvt.h_Njet50eta2p5->Clone("reco_njet50");
				hreconjet50->SetDirectory(0);
				TH1* hsmearnjet50 = (TH1*) Hist.at(jetbin).at(htbin).at(mhtbin).hv_SmearedEvt.h_Njet50eta2p5->Clone("smear_njet50");
				hsmearnjet50->SetDirectory(0);
				stringstream reconjet50leg, smnjet50leg;
				const double sum_reconj = hreconjet50->Integral(1,hreconjet50->GetNbinsX()+1);
				const double sum_smearnj= hsmearnjet50->Integral(1, hsmearnjet50->GetNbinsX()+1);
				reconjet50leg << "reco (" << sum_reconj << ")";
				smnjet50leg << "smeared (" << sum_smearnj << ")";

				//cout << red << "mht reco/smear  = " << sum_recomht << " / " << sum_smearmht << endl; 
				//cout << "nj50 reco/smear = " << sum_reconj << " / " << sum_smearnj  << clearatt << endl; 

				hreconjet50->SetLineWidth(2);
				hsmearnjet50->SetLineWidth(2);
				hsmearnjet50->SetLineColor(kRed);
				hreconjet50->SetMarkerStyle(kPlus);
				hsmearnjet50->SetMarkerColor(kRed);
				hsmearnjet50->SetMarkerStyle(24);

				TLegend *l3 = new TLegend(0.6,0.6,0.9,0.9);
				l3->AddEntry(hreconjet50, reconjet50leg.str().c_str());
				l3->AddEntry(hsmearnjet50, smnjet50leg.str().c_str());

				gPad->SetLogy(1);
				hreconjet50->DrawCopy();
				hsmearnjet50->DrawCopy("same");
				l3->Draw();
				gPad->Print("samples.eps");
				gPad->SetLogy(0);

				for (unsigned i=0; i < 3 ; ++i)
				{
					TH1* hrecojetpt = (TH1*) Hist.at(jetbin).at(htbin).at(mhtbin).hv_RecoJets.at(i).h_Jet_pt->Clone("reco_jet");
					hrecojetpt->SetDirectory(0);
					TH1* hgenjetpt = (TH1*) Hist.at(jetbin).at(htbin).at(mhtbin).hv_GenJets.at(i).h_Jet_pt->Clone("gen_jet");
					hgenjetpt->SetDirectory(0);
					TH1* hsmearjetpt = (TH1*) Hist.at(jetbin).at(htbin).at(mhtbin).hv_SmearedJets.at(i).h_Jet_pt->Clone("smear_jet");
					hsmearjetpt->SetDirectory(0);

					hrecojetpt->SetLineWidth(2);
					hgenjetpt->SetLineWidth(2);
					hsmearjetpt->SetLineWidth(2);
					hgenjetpt->SetLineColor(kBlue);
					hsmearjetpt->SetLineColor(kRed);
					hrecojetpt->SetMarkerStyle(kPlus);
					hgenjetpt->SetMarkerStyle(22);
					hsmearjetpt->SetMarkerStyle(24);
					hrecojetpt->SetMarkerColor(kBlack);
					hgenjetpt->SetMarkerColor(kBlue);
					hsmearjetpt->SetMarkerColor(kRed);

					hrecojetpt->DrawCopy();
					//hgenjetpt->DrawCopy("same");
					hsmearjetpt->DrawCopy("same");

					stringstream recoleg, genleg, smearleg;
					recoleg << "reco (" << hrecojetpt->Integral(1, hrecojetpt->GetNbinsX()+1) << ")";
					//genleg << "gen (" << hgenjetpt->Integral() << ")";
					smearleg << "smeared (" << hsmearjetpt->Integral(1, hsmearjetpt->GetNbinsX()+1) << ")";

					TLegend *l = new TLegend(0.6,0.6,0.9,0.9);
					l->AddEntry(hrecojetpt, recoleg.str().c_str());
					//l->AddEntry(hgenjetpt, genleg.str().c_str());
					l->AddEntry(hsmearjetpt, smearleg.str().c_str());
					l->Draw();
					gPad->Print("samples.eps");

					TH1* hsmearjetpt_copy = (TH1*) hsmearjetpt->Clone("smear_jet_copy");
					hsmearjetpt_copy->SetDirectory(0);
					hsmearjetpt_copy->Divide(hrecojetpt);
					hsmearjetpt_copy->SetTitle("Smeared Gen/Reco");
					hsmearjetpt_copy->Draw();
					//gPad->Print("samples.eps");

					delete l;
				}
			}
		}
	}
	gPad->Print("samples.eps]");


	vector<double> ptBins;
	vector<double> etaBins;
	if (bRUNNING_ON_MC && bDO_GENJET_SMEARING) 
	{
		ptBins  = *(smearFunc_->GetPtBinEdges());
		etaBins = *(smearFunc_->GetEtaBinEdges());

		/*TCanvas *c1 = new TCanvas("c1");
		  gStyle->SetOptStat(11111);
		  gPad->Print("JERs.eps[");

		  for (unsigned ptbin=0; ptbin < ptBins.size() -1; ++ptbin)
		  {
		  for (unsigned etabin=0; etabin < etaBins.size() -1; ++etabin)
		  {
		  jerHist.at(ptbin).at(etabin)->Draw();
		  gPad->Print("JERs.eps");
		  }
		  }
		  gPad->Print("JERs.eps]");
		  */
	}


	/***********************************************
	 *  End job summary
	 **********************************************/

	cout << ">>>>>>> " << __FILE__ << ":" << __FUNCTION__ << ": End Job " << endl;
	if (bDEBUG) cout << red << " ------ DEBUG MODE ---------- " << clearatt << endl;
	if (bNON_STD_MODE) 
	{
		cout << red << "      !!!!!!! NOT STD MODE !!!!!! : " << sNON_STD_MODE_EXPLAIN << clearatt << endl;
	}
	cout << red <<  "MC Flag  ------------- = " << bRUNNING_ON_MC << clearatt << endl;
	cout << "Entries found/proces'd = " << entriesFound << " / " << nProcessed << endl;
	cout << "Cleaning Failed Evts   = " << nCleaningFailed << endl;
	cout << "HO filter Failed Evts  = " << nHOfilterFailed << endl;
	cout << "Bad Hcal Laser Evts    = " << nBadHcalLaserEvts << endl;
	cout << "Bad Ecal Laser Evts    = " << nBadEcalLaserEvts << endl;
	cout << "TOP candidates Found   = " << topsFound << endl;
	cout << "---------- Settings ----------------" << endl;
	cout << "Cut Mask               = " << cutmask << endl;
	cout << "Cuts Enabled           = "; 
	stringstream cutlist;
	if (bAPPLY_NJET_CUT) cutlist << "NJET(70/50/30>=" << uMinNjet70Eta2p4_ << "/" << uMinNjet50Eta2p4_ << "/" << uMinNjet30Eta2p4_ << "), ";
	if (bAPPLY_MET_CUT) cutlist << "MET(>" << dMinMet_ << "), ";
	if (bAPPLY_TRIPLET_CUT) cutlist << "TRIPLETS(>" << uMinTriplets_ << "), ";
	if (bAPPLY_TOPMASS_CUT) cutlist << "TOPMASS(" << dMinTopMass_ << "," << dMaxTopMass_ << "), ";
	if (bAPPLY_TOPPLUSBJET_CUT) cutlist << "TOP+0.5*BJET(>" << dMinTopPlusBjetMass_ << "), ";
	if (bAPPLY_MT2_CUT) cutlist << "MT2(>" << dMinMt2_ << "), ";
	if (bAPPLY_DPHI_CUT) cutlist << "DPHI(.5,.5,.3)";
	cutlist << endl;

	cout << cutlist.str();

	if (bRUNNING_ON_MC) {
		cout << red <<  "PU weight applied?     = " << bDO_PU_WEIGHING << clearatt << endl;
		if (dMC_LUMI_WGT != 1) 
		{
			cout << "Data Lumi scaled to    = " << dDATA_LUMI << endl;
			cout << "Lumi wgt               = " << dMC_LUMI_WGT << endl;
		}
	} else {
		cout << "DO_TRIG_SELECTION      = " << bDO_TRIG_SELECTION << endl;
		cout << "DO_TRIG_PRESCALING     = " << bDO_TRIG_PRESCALING << endl;
		if (bDO_TRIG_SELECTION)
		{
			cout << "Trigger selection      = ";
			for (unsigned int i =0; i < vTriggersToUse.size(); ++i)
			{
				cout << vTriggersToUse.at(i);
				if (i+1<vTriggersToUse.size()) cout << "/";
				else cout << endl;
			}
		}
	}
	cout << red << "bAPPLY_DPHI_CUT        = " << bAPPLY_DPHI_CUT  << clearatt << endl;
	cout << "Jetbins                = "; for (int bin=0; bin < JetBins_.size(); ++bin) { cout << JetBins_.at(bin).first << "-" << JetBins_.at(bin).second << ", "; }; cout << endl; 
	cout << "Htbins                 = "; for (int bin=0; bin < HtBins_.size(); ++bin) { cout << HtBins_.at(bin) << ", "; }; cout << endl; 
	cout << "Mhtbins                = "; for (int bin=0; bin < MhtBins_.size(); ++bin) { cout << MhtBins_.at(bin) << ", "; }; cout << endl; 
	cout << "nReco/nGen/nSm Evts    = " << nRecoJetEvts << "/" << nGenJetEvts << "/" << nSmearedJetEvts << endl;
	
	if (bRUNNING_ON_MC && bDO_GENJET_SMEARING)
	{
		cout << red << "---- Smearing Function Settings ----" << clearatt << endl;
		cout << "Sampling (nTries_)     = " << uNTRIES << endl;
		cout << "Smear wgt              = " << smearingWgt << endl;
		cout << "smearedJetPt_          = " << smearedJetPt_ << endl;

		const bool  absTailScalingFact_val   = smearFunc_->GetAbsoluteTailScaling();
		const float lowerTailScalingFact_val = smearFunc_->GetLowerTailScalingVariation();
		const float upperTailScalingFact_val = smearFunc_->GetUpperTailScalingVariation();
		const float additionalSmearingFact_val   = smearFunc_->GetAdditionalSmearingVariation();

		const bool bSystematicMode = (lowerTailScalingFact_val != 1 
				|| upperTailScalingFact_val != 1 || additionalSmearingFact_val != 1);

		if (bSystematicMode) {	
			cout << red << " >>>> SYSTEMATIC MODE <<<< " << endl;
		}
		cout << "Uncertainty Mode       = " << smearFunc_->GetUncertaintyName() << endl;
		vector<double> values;
		values = smearFunc_->GetPtBinEdges_scaling();
		cout << "PtBinEdges_scaling     = "; PrintVector(values);
		values = smearFunc_->GetEtaBinEdges_scaling();
		cout << "EtaBinEdges_scaling    = "; PrintVector(values);
		values = smearFunc_->GetLowerTail_scaling();
		cout << "LowerTail_scaling      = "; PrintVector(values);
		values = smearFunc_->GetUpperTail_scaling();
		cout << "UpperTail_scaling      = "; PrintVector(values);
		values = smearFunc_->GetAdditionalSmear_scaling();
		cout << "Additional Smearing    = "; PrintVector(values);

		cout << "AbsoluteTailScaling_val= " << absTailScalingFact_val << " (false = systematic mode)" << endl;
		cout << "LowerTailScaling_val   = " << lowerTailScalingFact_val << endl;
		cout << "UpperTailScaling_val   = " << upperTailScalingFact_val << endl;
		cout << "AdditionalSmearing_val = " << additionalSmearingFact_val << endl;

		if (bSystematicMode) {	
			cout << clearatt << endl;
		}
		cout << "Smearing file          = " << smearFunc_->SmearingFile() << endl;
		cout << "Jet Res. Collection    = " << smearFunc_->GetResFuncCollType() << endl;
	}

	cout << "---- Actual Event Counts in Reco (Smeared) MHT hist ---- " << endl;
	cout << setw(8) << "jet bin" << setw(10) << "ht bin" << setw(12) << "mht bin" << setw(10) << "Reco" << setw(10) << "Smear"  << setw(15) << "Smear/Reco" << endl; 
	for (unsigned jetbin=0; jetbin < JetBins_.size(); ++jetbin)
	{
		stringstream jetbin_range;
		jetbin_range << JetBins_.at(jetbin).first << "-" << JetBins_.at(jetbin).second;
		for (unsigned htbin=0; htbin < HtBins_.size() -1 ; ++htbin)
		{
			stringstream htbin_range;
			htbin_range << HtBins_.at(htbin) << "-" << setw(5) << HtBins_.at(htbin+1); 
			for (unsigned mhtbin=0; mhtbin < MhtBins_.size() -1 ; ++mhtbin)
			{
				stringstream mhtbin_range;
				mhtbin_range << setw(6) << MhtBins_.at(mhtbin) << "-" << setw(5) << MhtBins_.at(mhtbin+1); 
				double int_smear = 0; 
				if (bRUNNING_ON_MC) int_smear = (Hist.at(jetbin).at(htbin).at(mhtbin).hv_SmearedEvt.h_Mht)->Integral();
				const double int_reco = (Hist.at(jetbin).at(htbin).at(mhtbin).hv_RecoEvt.h_Mht)->Integral();
				double ratio = 0;
				if (bRUNNING_ON_MC) ratio = int_reco>0 ? int_smear/int_reco : 0 ;
				
				cout << setprecision(1);
				cout << fixed << setw(6) << jetbin_range.str() << setw(12) << htbin_range.str() 
						<< setw(12) << mhtbin_range.str() << setw(10) << int_reco << setw(10) << int_smear  << setw(15) << ratio << endl; 
			}
		}
	}
	if (nVectorInexWarnings>0)
	{
		cout << "---- ERROR/WARNING SUMMARY ---------" << endl;
		cout << "Values out of range errors = " << nVectorInexWarnings << endl;
		cout << red << "Type-1 MET Phi Corrections Not applied!!! as they are not in the NTuples!!" << clearatt << endl;
	}

}

double FactorizationBySmearing::HT (const std::vector<TLorentzVector>& vjets) {
	double ht = 0.0;
	for( unsigned int ijet=0; ijet<vjets.size(); ijet++) {    
		if( vjets[ijet].Pt()>50.0 && std::abs(vjets[ijet].Eta())<2.5 ) 
			ht += vjets[ijet].Pt();
	}
	return ht;
}

TLorentzVector FactorizationBySmearing::GetMET(const vector<TLorentzVector>& jets)
{
	TVector3 MHT;
	MHT.SetXYZ(0.0,0.0,0.0);
	for (int i = 0; i < jets.size(); i++)
	{
		if (jets.at(i).Pt()<10.0 || fabs(jets.at(i).Eta())>5.0) continue;
		MHT -= jets[i].Vect();
	}
	MHT.SetZ(0.0);
	
	TLorentzVector metvec(0,0,0,0);
	metvec.SetVectMag(MHT,0);
	
	return metvec;
}

TLorentzVector FactorizationBySmearing::MHT(const std::vector<TLorentzVector>& vjets) 
{
	TLorentzVector mht(0.0, 0.0, 0.0, 0.0);
	for ( unsigned int ijet=0; ijet<vjets.size(); ijet++) 
	{    
		if ( vjets[ijet].Pt()>30.0 && std::fabs(vjets[ijet].Eta())<5.0 ) 
		{
			mht -= vjets[ijet];
		}
	}
	return mht;
}

void FactorizationBySmearing::CreateRecoJetVec(std::vector<TLorentzVector>& vjets, std::vector<double>& bDisc)
{
	vjets.clear();
	bDisc.clear();
	if (t_PFJetPt->size() != t_PFJetBTag->size())
	{
		cout << __FUNCTION__ << ":: jet size (" << t_PFJetPt->size() << ") and bDisc (" 
			<< t_PFJetBTag->size() << ") vec sizes do not match!!!" << endl;
		assert(false);
	}
	for (unsigned i=0; i < t_PFJetPt->size(); i++)
	{
		 //this is the correct way to set vector
		 //with the info I have save in the tree
		TLorentzVector tl(0,0,0,0); 
		tl.SetPtEtaPhiE( (*t_PFJetPt)[i], (*t_PFJetEta)[i],
							  (*t_PFJetPhi)[i], (*t_PFJetE)[i]   );
		vjets.push_back(tl);
		bDisc.push_back( (*t_PFJetBTag)[i]);
	}

	//std::sort(vjets.begin(), vjets.end(), PtAComparator);
	domysort(vjets, bDisc);
}

void FactorizationBySmearing::CreateGenJetVec(std::vector<TLorentzVector>& vjets)
{
	vjets.clear();
	for (unsigned i=0; i < t_genJetPt->size(); i++)
	{
		TLorentzVector tl(0,0,0,0);
		tl.SetPtEtaPhiE( t_genJetPt->at(i), t_genJetEta->at(i),
							 t_genJetPhi->at(i), t_genJetE->at(i));
		vjets.push_back(tl);
	}

	std::sort(vjets.begin(), vjets.end(), PtAComparator);
}

//--------------------------------------------------------------------------
void FactorizationBySmearing::SmearingGenJets(const vector<TLorentzVector>& jets_gen, 
		const vector<double>& bDisc_gen, std::vector<TLorentzVector>& genJets_smeared, 
		std::vector<double>& bDisc_smeared)
{
	if (bDEBUG) cout << __FUNCTION__ << ":" << __LINE__ << endl;

	genJets_smeared.clear();
	bDisc_smeared.clear();
	
	int i_jet = 0;
	unsigned i =0;

	for (vector<TLorentzVector>::const_iterator it = jets_gen.begin(); it != jets_gen.end(); ++it) {

		if (it->Pt() > smearedJetPt_) {
			const double scale  = JetResolutionHist_Pt_Smear(it->Pt(), it->Eta(), i_jet); //what is this i_jet purpose??? jet id index

			const double newPt  = it->Pt() * scale;
			const double newEta = it->Eta();
			const double newPhi = it->Phi();
			const double newM   = it->M();
			if (bDEBUG) cout << "old/new pt = " << it->Pt() << "/" << newPt << endl; 

			//for JER reconstruction for debugging
			int i_Pt  = GetIndex(it->Pt(), smearFunc_->GetPtBinEdges());
			int i_eta = GetIndex(it->Eta(), smearFunc_->GetEtaBinEdges());
			jerHist.at(i_Pt).at(i_eta)->Fill(scale);

			TLorentzVector newP4(0,0,0,0);
			newP4.SetPtEtaPhiM(newPt, newEta, newPhi, it->M());
			TLorentzVector smearedJet(0,0,0,0);
			smearedJet.SetPxPyPzE(newP4.Px(), newP4.Py(), newP4.Pz(), newP4.E());
			genJets_smeared.push_back(smearedJet);
			bDisc_smeared.push_back(bDisc_gen.at(i));
			++i_jet;
		} else {
			TLorentzVector smearedJet(*it);
			genJets_smeared.push_back(smearedJet);
			bDisc_smeared.push_back(bDisc_gen.at(i));
		}
		++i;

	}

	if (genJets_smeared.size() != bDisc_smeared.size())
	{
		cout << __LINE__ << ": jet/disc mismatch = " << genJets_smeared.size() << "/ " << bDisc_smeared.size() << endl;
		assert(false);
	}
//	std::sort(genJets_smeared.begin(), genJets_smeared.end(), PtAComparator);
	domysort(genJets_smeared, bDisc_smeared);

	if (bDEBUG) cout << __FUNCTION__ << ":" << __LINE__ << endl;
}

//--------------------------------------------------------------------------
// pt resolution for smearing
double FactorizationBySmearing::JetResolutionHist_Pt_Smear(const double& pt, const double& eta, const int& i)
{
	if (bDEBUG) cout << __FUNCTION__ << ":" << __LINE__ << endl;
   int i_jet;
   i < 2 ? i_jet = i : i_jet = 2;
   int i_Pt = GetIndex(pt, smearFunc_->GetPtBinEdges());
   int i_eta = GetIndex(eta, smearFunc_->GetEtaBinEdges());

	//const double res = smearFunc_->getSmearFunc(1,i_jet,i_eta,i_Pt)->GetRandom();
	double res = 1.0;
	//int tries = 0;
	//for (tries; tries< 10; ++tries)
	///{
		const double res_temp = smearFunc_->getSmearFunc(0,i_jet,i_eta,i_Pt)->GetRandom();
		//if (res_temp <0.00001) 
		//{
		//	cout << i_jet << ":" << i_eta << ":" << i_Pt << ": tries=" << tries << " = " << res_temp << endl;
		//	continue;
		//}
		res = res_temp;
		//break;
	//}
	//if (tries == 10) cout << i_jet << ":" << i_eta << ":" << i_Pt << ": MAX TRIES REACHED!!!!" << endl;
	if (bDEBUG) cout << __FUNCTION__ << ":" << __LINE__ << ": res = " << res << endl;

   return res;
}

//--------------------------------------------------------------------------
int FactorizationBySmearing::GetIndex(const double& x, const std::vector<double>* vec)
{
   int index = -1;
   // this is a check
   //int index = 0;
	if (bDEBUG) cout << __FUNCTION__ << ":" << __LINE__ << ": x / vec size = " << x << " / " << vec->size() << endl; 
   for (std::vector<double>::const_iterator it = vec->begin(); it != vec->end(); ++it) {
      if ((*it) > fabs(x))
         break;
      ++index;
   }
   if (index < 0)
      index = 0;
   if (index > (int) vec->size() - 2)
      index = vec->size() - 2;


	if (index < 0) 
	{
		cout << __FUNCTION__ <<
					":" << __LINE__ << ":: WARN !!! INDEX = " << index << " < 0!!! "<< endl; 
	}

   return index;
}

//--------------------------------------------------------------------------
// HF probability for smearing
/*double FactorizationBySmearing::GetHFProb(const double& pt, const double& eta, const int& i_jet) {
  
   int i_bin;
   if( i_jet == 0 ) {
      i_bin = h_bProb_jet1->FindBin(pt, eta);
      return h_bProb_jet1->GetBinContent(i_bin);
   }
   if( i_jet == 1 ) {
      i_bin = h_bProb_jet2->FindBin(pt, eta);
      return h_bProb_jet2->GetBinContent(i_bin);
   }
   else {
      i_bin = h_bProb_jet3p->FindBin(pt, eta);
      return h_bProb_jet3p->GetBinContent(i_bin);
   }
}
*/

void FactorizationBySmearing::DumpJets(const vector<TLorentzVector>& jets)
{
	cout << setw(3) << "#" 
			<< setw(15) << "Pt" 
			<< setw(15) << "eta"  
			<< setw(15) << "phi" << endl;  
	for (unsigned i=0; i < jets.size(); i++)
	{
		cout << setw(3) << i 
			<< setw(15) << jets.at(i).Pt() 
			<< setw(15) << jets.at(i).Eta() 
			<< setw(15) << jets.at(i).Phi() << endl;  
	}
}


bool FactorizationBySmearing::FillHistogram(const vector<TLorentzVector>& jets, const vector<double>& bDisc, 
				const TLorentzVector& reco_uncl_met_vec, const int& jetcoll, const double& wgt)
{
	//0=reco
	//1=gen
	//2=smeared
	
//	const TLorentzVector mhtvec(MHT(jets)); 
	const TLorentzVector clus_metvec(GetMET(jets)); 

	TVector3 mv3 = reco_uncl_met_vec.Vect() + clus_metvec.Vect(); // NEED TO APPLY Type-1 MET phi corrections
	mv3.SetZ(0);
	TLorentzVector tot_metvec(0,0,0,0);
	tot_metvec.SetVectMag(mv3,0);

//	const double mht        = mhtvec.Pt(); 
	const double met        = tot_metvec.Pt(); 
	const double mht        = mht;  //this is just lazy! so I do not have to rename mht everywhere
	const double ht         = HT(jets);
	const int njet50eta2p5  = CountJets(jets, 50.0, 2.5); 
	const int njet30eta5p0  = CountJets(jets, 30.0, 5.0); 
	const int nVtx          = t_NVertices;
		
	unsigned njet_pt70 =0, njet_pt50=0, njet_pt30=0;
	for (unsigned i=0; i < jets.size(); ++i)
	{
		if (fabs(jets.at(i).Eta()) > 2.4) continue;
		const float pt = jets.at(i).Pt();
		if (pt>30.0) ++njet_pt30;
		if (pt>50.0) ++njet_pt50;
		if (pt>70.0) ++njet_pt70;
	}

	//top tagger stuff

	const double tripletJetPtCut_ = 5.0; //fake for now
	vector<TLorentzVector> Jets_thisTry_PtCut;
	for(unsigned jet_i = 0; jet_i < jets.size(); jet_i++)
	{
		if(jets.at(jet_i).Pt() > tripletJetPtCut_)
		{
			Jets_thisTry_PtCut.push_back(jets.at(jet_i));
		}
	}

	std::vector<TLorentzVector> triplet;
	std::vector<TLorentzVector> rSystem;
	std::vector<TLorentzVector> bJet;

	double M23OverM123 = -9999.99;
	double M123        = -9999.99; 

	TripletSelector(Jets_thisTry_PtCut, bDisc, triplet, rSystem, bJet, M23OverM123, M123);

//	..M123 mass of the top
   //MTMT2(vMET,triplet,rSystem,bJet,MT2,MTt,MTb);
	double MT2 = 0, MTt =0 , MTb=0; 
   MTMT2(tot_metvec,triplet,rSystem,bJet,MT2,MTt,MTb);
	//end top tagger stuff
	bool pass_topselection = true;
	const double MTb_p_MTt = MTb + 0.5 * MTt;


	//now apply cuts
	const bool pass_njetcut       = (njet_pt70>=2 && njet_pt50>=4 && njet_pt30>=5) ? true : false;           // 0000001 = 1
	const bool pass_metcut        = (tot_metvec.Pt()>175.0 ) ? true : false;                                 // 0000010 = 2
	const bool pass_tripletcut    = (triplet.size()>=1) ? true : false;                                      // 0000100 = 4
	const bool pass_topmasscut    = ( M123>80.0 && M123<270.0 ) ? true : false;                              // 0001000 = 8
	const bool pass_top_plus_bjet = (MTb_p_MTt> 500.0) ? true : false;                                       // 0010000 = 16
	const bool pass_mt2cut        = (MT2>300.0) ? true : false;                                              // 0100000 = 32
 	// this is not obsolete as the jet cut is more complex now.
	// Need to fix the code for this
	//const bool pass_dphicut    = PassDphiCut(jets, tot_metvec, JetBins_.at(i_JetBin).first, 0.5, 0.5, 0.3);  // 1000000 = 64
	const bool pass_dphicut    = PassDphiCut(jets, tot_metvec, 3, 0.5, 0.5, 0.3);  // 1000000 = 64

	if (bAPPLY_NJET_CUT && ! pass_njetcut ) return 0;
	if (bAPPLY_MET_CUT  && ! pass_metcut ) return 0;
	if (bAPPLY_TRIPLET_CUT && ! pass_tripletcut ) return 0;
	if (bAPPLY_TOPMASS_CUT && ! pass_topmasscut ) return 0;
	if (bAPPLY_TOPPLUSBJET_CUT && ! pass_top_plus_bjet ) return 0;
	if (bAPPLY_MT2_CUT && ! pass_mt2cut ) return 0;
	if (bAPPLY_DPHI_CUT && ! pass_dphicut ) return 0;


	const unsigned nbjets = bJet.size();
	double bjetmass = 0;
	double bjetpt   = 0;
	if (bJet.size()>0)
	{
		bjetmass = bJet.at(0).M();
		bjetpt   = bJet.at(0).Pt();
	}

	const unsigned i_JetBin = GetVectorIndex(JetBins_, (unsigned) njet50eta2p5);
	const unsigned i_HtBin  = GetVectorIndex(HtBins_, ht);
	const unsigned i_MhtBin = GetVectorIndex(MhtBins_, mht);


	//This is primarily for gen-jet smearing
	//discard event that do not fall into given njet/ht/mht bin
	//this is primarily for counting when low bound of HT =500
	//the upper bounds are set much higher as not to throw events
	bool discard = false;
	if ( i_JetBin > JetBins_.size() ) { discard = discard || true; /*if ((unsigned) nVectorInexWarnings % 500 == 0) cout << "iJetBin = " << i_JetBin << " is out of range. Discarding event! " << endl; */} 
	if ( i_HtBin  > HtBins_.size()  ) { discard = discard || true; /*if ((unsigned) nVectorInexWarnings % 500 == 0) cout << "iHtBin  = " << i_HtBin  << " is out of range. Discarding event! " << endl; */} 
	if ( i_MhtBin > MhtBins_.size() ) { discard = discard || true; /*if ((unsigned) nVectorInexWarnings % 500 == 0) cout << "iMhtBin = " << i_MhtBin << " is out of range. Discarding event! " << endl;*/ } 
	if (discard) return discard;

	//const float dPhiMin = DelPhiMin(jets,	tot_metvec, JetBins_.at(i_JetBin).first); // this is not obsolete as the jet cut is more complex now.
	//Need to fix the code for this
	const float dPhiMin = DelPhiMin(jets,	tot_metvec, 3);

	/*****************************************************************
	*                 for systematics
	*  Try dphi cuts slightly different from RA2
	*  to probe how the control region would change
	*  with such selections.
	****************************************************************/
	const bool bPassDphiSystVariation_1 = PassDphiCut(jets, tot_metvec, JetBins_.at(i_JetBin).first, 0.5, 0.5, 0.1);
	const bool bPassDphiSystVariation_2 = PassDphiCut(jets, tot_metvec, JetBins_.at(i_JetBin).first, 0.4, 0.4, 0.2);

	if (bDEBUG) 
	{
		cout << "[jetcoll = " << jetcoll << "] ht/mht/mymht = " 
				<< ht << "/" << mht 
				//<< "/" << mymht 
				<< endl;
	}

	const vector<TLorentzVector> jets50 = GetPt50Eta2p5Jets(jets);

	if (jetcoll == 0) //reco hists
	{ 
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_RecoEvt.h_Mht->Fill(mht, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_RecoEvt.h_Ht->Fill(ht, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_RecoEvt.h_Njet50eta2p5->Fill(njet50eta2p5, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_RecoEvt.h_Njet30eta5p0->Fill(njet30eta5p0, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_RecoEvt.h_DphiMin->Fill(dPhiMin, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_RecoEvt.h_DphiMinVsMht->Fill(mht, dPhiMin, wgt);
		//the wgt = trigPrescale when running over data
		//wgt = puweight when running over MC OR
		//wgt = smear (and/or) lumi weight when running over MC 
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_RecoEvt.h_evtWeight->Fill(wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_RecoEvt.h_nVtx->Fill(nVtx, wgt);

		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_RecoEvt.h_nbjets->Fill(nbjets, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_RecoEvt.h_bjetMass->Fill(bjetmass, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_RecoEvt.h_bjetPt->Fill(bjetpt, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_RecoEvt.h_M123->Fill(M123, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_RecoEvt.h_M23OverM123->Fill(M23OverM123, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_RecoEvt.h_MT2->Fill(MT2, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_RecoEvt.h_MTb->Fill(MTb, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_RecoEvt.h_MTt->Fill(MTt, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_RecoEvt.h_MTb_p_MTt->Fill(MTb_p_MTt, wgt);


		FillJetHistogram(jets50, Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_RecoJets, jetcoll, tot_metvec, wgt);

		if (pass_dphicut) 
		{
			Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_RecoEvt.signal->Fill(mht,wgt);
			Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_RecoEvt.signalFineBin->Fill(mht,wgt);
			//in case of data, the only weight comes from trigPrescale
			Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_RecoEvt.signal_trigPrescales->Fill(wgt);
		}

		for (unsigned i=0; i < vDphiVariations.size(); ++i)
		{
			if (dPhiMin > vDphiVariations.at(i)) 
			{
				Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_RecoEvt.pass.at(i)->Fill(mht, wgt);
				Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_RecoEvt.passFineBin.at(i)->Fill(mht, wgt);
				//in case of data, the only weight comes from trigPrescale
				Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_RecoEvt.pass_trigPrescales.at(i)->Fill(wgt);
			} else 
			{
				Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_RecoEvt.fail.at(i)->Fill(mht, wgt);
				Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_RecoEvt.failFineBin.at(i)->Fill(mht, wgt);
				//in case of data, the only weight comes from trigPrescale
				Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_RecoEvt.fail_trigPrescales.at(i)->Fill(wgt);
			}

		}

		if (! bPassDphiSystVariation_1) //we want the denominator (or failed events)
		{
			Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_RecoEvt.sidebandSyst[0]->Fill(mht, wgt);
			Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_RecoEvt.sidebandSystFineBin[0]->Fill(mht, wgt);
		}
		if (! bPassDphiSystVariation_2)
		{
			Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_RecoEvt.sidebandSyst[1]->Fill(mht, wgt);
			Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_RecoEvt.sidebandSystFineBin[1]->Fill(mht, wgt);
		}


	} else if (jetcoll == 1)  //gen hists
	{
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_GenEvt.h_Mht->Fill(mht, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_GenEvt.h_Ht->Fill(ht, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_GenEvt.h_Njet50eta2p5->Fill(njet50eta2p5, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_GenEvt.h_Njet30eta5p0->Fill(njet30eta5p0, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_GenEvt.h_DphiMin->Fill(dPhiMin, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_GenEvt.h_DphiMinVsMht->Fill(mht, dPhiMin, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_GenEvt.h_evtWeight->Fill(wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_GenEvt.h_nVtx->Fill(nVtx, wgt);

		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_GenEvt.h_nbjets->Fill(nbjets, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_GenEvt.h_bjetMass->Fill(bjetmass, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_GenEvt.h_bjetPt->Fill(bjetpt, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_GenEvt.h_M123->Fill(M123, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_GenEvt.h_M23OverM123->Fill(M23OverM123, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_GenEvt.h_MT2->Fill(MT2, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_GenEvt.h_MTb->Fill(MTb, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_GenEvt.h_MTt->Fill(MTt, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_GenEvt.h_MTb_p_MTt->Fill(MTb_p_MTt, wgt);
		FillJetHistogram(jets50, Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_GenJets, jetcoll, tot_metvec, wgt);
	} else if (jetcoll == 2)  //smeared hists
	{
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.h_Mht->Fill(mht, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.h_Ht->Fill(ht, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.h_Njet50eta2p5->Fill(njet50eta2p5, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.h_Njet30eta5p0->Fill(njet30eta5p0, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.h_DphiMin->Fill(dPhiMin, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.h_DphiMinVsMht->Fill(mht, dPhiMin, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.h_evtWeight->Fill(wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.h_nVtx->Fill(nVtx, wgt);

		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.h_nbjets->Fill(nbjets, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.h_bjetMass->Fill(bjetmass, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.h_bjetPt->Fill(bjetpt, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.h_M123->Fill(M123, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.h_M23OverM123->Fill(M23OverM123, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.h_MT2->Fill(MT2, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.h_MTb->Fill(MTb, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.h_MTt->Fill(MTt, wgt);
		Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.h_MTb_p_MTt->Fill(MTb_p_MTt, wgt);
		FillJetHistogram(jets50, Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_SmearedJets, jetcoll, tot_metvec, wgt);


		if (pass_dphicut) 
		{
			Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.signal->Fill(mht,wgt);
			Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.signalFineBin->Fill(mht,wgt);
		}

		for (unsigned i=0; i < vDphiVariations.size(); ++i)
		{
			if (dPhiMin > vDphiVariations.at(i)) 
			{
				Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.pass.at(i)->Fill(mht, wgt);
				Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.passFineBin.at(i)->Fill(mht, wgt);
			} else 
			{
				Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.fail.at(i)->Fill(mht, wgt);
				Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.failFineBin.at(i)->Fill(mht, wgt);
			}
		}


		if (! bPassDphiSystVariation_1) //we want the denominator (or failed events)
		{
			Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.sidebandSyst[0]->Fill(mht, wgt);
			Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.sidebandSystFineBin[0]->Fill(mht, wgt);
		}
		if (! bPassDphiSystVariation_2)
		{
			Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.sidebandSyst[1]->Fill(mht, wgt);
			Hist.at(i_JetBin).at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.sidebandSystFineBin[1]->Fill(mht, wgt);
		}

	} else {
		cout << __FUNCTION__ << ": Invalid jet colelctions!" << endl;
		assert (false);
	}
	
	if (jetcoll == 0 && mht>1000)
	{
		cout << "njet/ht/mht=" << njet50eta2p5 << "/" << ht << "/" << mht;
		PrintEventNumber();
	}
	return true;
}
void FactorizationBySmearing::FillJetHistogram(const vector<TLorentzVector>& jets, 
				vector<JetHist_t> hist, const int jetcoll, const TLorentzVector& mhtvec, 
				const double& wgt) 
{
	const unsigned njetsExist = jets.size();
	const unsigned njetsIwant = 10;
	unsigned njet2loop = njetsIwant;
	if (njetsExist<njetsIwant) njet2loop = njetsExist;

	int njetsfound = 0;
	for (unsigned i=0; i < njetsExist; ++i)
	{
		if (njetsfound>= njetsIwant) break;

		const float pt   = jets.at(i).Pt();
		const float eta  = jets.at(i).Eta();
		const float phi  = jets.at(i).Phi();
		const float dphi = fabs(jets.at(i).DeltaPhi(mhtvec));

		++njetsfound;
		const int j = njetsfound - 1;
		hist.at(j).h_Jet_pt->Fill(pt, wgt);
		hist.at(j).h_Jet_eta->Fill(eta, wgt);
		hist.at(j).h_Jet_phi->Fill(phi, wgt);
		hist.at(j).h_Jet_dphi->Fill(dphi, wgt);
	}
}

int FactorizationBySmearing::CountJets(const std::vector<TLorentzVector>& vjets, 
					const double minPt, const double maxEta)
{
	int njet = 0;
	for (unsigned int ijet=0; ijet<vjets.size(); ijet++) {    
		if ( vjets[ijet].Pt()>minPt && std::fabs(vjets[ijet].Eta())<maxEta ) 
		{
			njet++;
		}
	}
	return njet;
}

float FactorizationBySmearing::DelPhiMin(const vector<TLorentzVector>& jets, 
					const TLorentzVector& mhtVec, const unsigned& nJet50min)
{
	std::vector<float> vDelPhi_jetmht;
	const unsigned loopTo = std::min(3, (int) nJet50min); 

	for (unsigned i = 0 ; i < jets.size() ; ++i)
	{	//use only three leading jets
		//or in the case >=2 min jets, use only 2 jets
		const float delphi_jetmht = fabs(jets.at(i).DeltaPhi(mhtVec));
		vDelPhi_jetmht.push_back(delphi_jetmht);
		if ( vDelPhi_jetmht.size() >= loopTo ) break; 
	}

	std::sort(vDelPhi_jetmht.begin(), vDelPhi_jetmht.end(), sort_using_less_than);	

	
	const float dPhiMin = vDelPhi_jetmht.at(0);
//	cout << __FUNCTION__ << ": delphimin = " <<  dPhiMin << endl;
	return dPhiMin;
}

bool FactorizationBySmearing::PassDphiCut(const vector<TLorentzVector>& jets, 
					const TLorentzVector& mht, const unsigned nJet50min, 
					const float& cut1, const float& cut2, const float& cut3)
{
	bool pass = true;
	//use only three leading jets or in the case >=2 min jets, use only 2 jets
	const unsigned upTo = std::min(3, (int) nJet50min); 

	if (jets.size() >=1 && upTo>=1) { pass = pass && (fabs(jets.at(0).DeltaPhi(mht)) > cut1); }
	if (jets.size() >=2 && upTo>=2) { pass = pass && (fabs(jets.at(1).DeltaPhi(mht)) > cut2); }
	if (jets.size() >=3 && upTo>=3) { pass = pass && (fabs(jets.at(2).DeltaPhi(mht)) > cut3); }
	return pass;
}
void FactorizationBySmearing::BookJerDebugHists()
{
	vector<double> ptBins = *(smearFunc_->GetPtBinEdges());
	vector<double> etaBins = *(smearFunc_->GetEtaBinEdges());
	
	TDirectory *dir = (TDirectory*) outRootFile->FindObject("Hist");
	if (dir == NULL)
	{
		cout << __FUNCTION__ << ": dir named Hist not found!" << endl;
		assert(false);
	}
	TDirectory *subdir = dir->mkdir("JERS");
	subdir->cd();
	
	for (unsigned ptbin=0; ptbin < ptBins.size() -1; ++ptbin)
	{
		vector<TH1*> h;
		for (unsigned etabin=0; etabin < etaBins.size() -1; ++etabin)
		{
			const float ptmin  = ptBins.at(ptbin);
			const float ptmax  = ptBins.at(ptbin+1);
			const float etamin = etaBins.at(etabin);
			const float etamax = etaBins.at(etabin+1);
			stringstream name, title;
			name << "jer_pt" << ptmin << "to" << ptmax << "_eta" << etamin << "to" << etamax;  
			title << "Pt" << ptbin << "_Eta"<< etabin 
				   << ":[" <<  ptmin << "<PT<" << ptmax << ", " 
				   << etamin << "< #eta <" << etamax << "]";  
			TH1* hist = new TH1D(name.str().c_str(), title.str().c_str(), 2000,0,2); 
			//hist->Sumw2();
			h.push_back(hist);
		}
		jerHist.push_back(h);
	}
}

unsigned FactorizationBySmearing::GetVectorIndex(const vector<double>& binEdges, const double& val)
{
	for (unsigned bin=0; bin < binEdges.size() -1; ++bin)
	{
		const double min = binEdges.at(bin);
		const double max = binEdges.at(bin+1);
		if (val>=min && val<max) return bin;
	}

	++nVectorInexWarnings;
	/* //commenting out the warnings
	 stringstream msg;
	msg << __FUNCTION__ << ": WARNNING! val = " << val << " is out of bin ranges[ ";

	for (unsigned bin=0; bin < binEdges.size(); ++bin)
	{
		msg << binEdges.at(bin);
		if (bin< (binEdges.size()-1)) msg << ", ";
	}
	msg << "]" << endl;
	msg << "   Assigning extreme value!" << endl;
	if ( (unsigned)nVectorInexWarnings % 10000  == 0) cout << msg.str();
	*/
	return 99999;
}
unsigned FactorizationBySmearing::GetVectorIndex(const vector< pair<unsigned, unsigned> >& binEdges, const unsigned& val)
{
	for (unsigned bin=0; bin < binEdges.size(); ++bin)
	{
		const double min = binEdges.at(bin).first;
		const double max = binEdges.at(bin).second;
		if (val>=min && val<=max) return bin;
	}

	++nVectorInexWarnings;

	/* //commenting out the warnings
	stringstream msg;
	msg << __FUNCTION__ << ": WARNNING! val = " << val << " is out of bin ranges[ ";

	for (unsigned bin=0; bin < binEdges.size(); ++bin)
	{
	 	msg  << binEdges.at(bin).first << "-" << binEdges.at(bin).second;
		if (bin< (binEdges.size()-1)) msg << ", ";
	}
	msg << "]" << endl;
	msg << "   Assigning extreme value!" << endl;
	if ( (unsigned) nVectorInexWarnings % 10000 ==0 ) cout << msg.str();
	*/
	return 99999;
}


double FactorizationBySmearing::GetLumiWgt(const string& datasetname, const double& dataLumi)
{
	/*************************************************
	 * X sec of Summer 12 Pythis MC samples
	 * **********************************************/

	//cross section for each sample in pt order
	const float xSec_Pythia[] = {
		1759.549,  //pt 300-470
		113.8791,  //470-600
		26.9921,  //600-800
		3.550036, //8000-1000
		0.737844, //1000-1400
		0.03352235, //1400-1800
		0.001829005 //1800
	};

	const float nEvts_Pythia[] = {
		5927300, // 300-470  #numbers from DBS, PREP page numbers are approximate
		3994848, // 470-600
		3992760, // 600-800
		3998563, //800-1000
		1964088, //1000-1400
		2000062, //1400-1800
		977586 //1800
	};

	const float xSec_MG[] = {
		276000, //250-HT-500
		8426,   //500-HT-1000
		204 	  //HT>1000
	};

	const float nEvts_MG[] = {
		27002490,
		30599239,
		13829995
	};

	double lumiWgt = 1;

	if (datasetname.find("pythia") != string::npos && datasetname.find("MG") == string::npos)
	{
		if (datasetname.find("QCD_Pt_300to470") != string::npos) 
		{
			cout << "Lum for Pythia qcd 1 used." << endl; 
			lumiWgt = dataLumi/(nEvts_Pythia[0] / xSec_Pythia[0]);
		} else if (datasetname.find("QCD_Pt_470to600") != string::npos) 
		{
			cout << "Lum for Pythia qcd 2 used." << endl; 
			lumiWgt = dataLumi/(nEvts_Pythia[1] / xSec_Pythia[1]);
		} else if (datasetname.find("QCD_Pt_600to800") != string::npos) 
		{
			cout << "Lum for Pythia qcd 3 used." << endl; 
			lumiWgt = dataLumi/(nEvts_Pythia[2] / xSec_Pythia[2]);
		} else if (datasetname.find("QCD_Pt_800to1000") != string::npos) 
		{
			cout << "Lum for Pythia qcd 4 used." << endl; 
			lumiWgt = dataLumi/(nEvts_Pythia[3] / xSec_Pythia[3]);
		} else if (datasetname.find("QCD_Pt_1000to1400") != string::npos) 
		{
			cout << "Lum for Pythia qcd 5 used." << endl; 
			lumiWgt = dataLumi/(nEvts_Pythia[4] / xSec_Pythia[4]);
		} else if (datasetname.find("QCD_Pt_1400to1800") != string::npos) 
		{ 
			cout << "Lum for Pythia qcd 6 used." << endl; 
			lumiWgt = dataLumi/(nEvts_Pythia[5] / xSec_Pythia[5]);
		} else if (datasetname.find("QCD_Pt_1800") != string::npos) 
		{
			cout << "Lum for Pythia qcd 7 used." << endl; 
			lumiWgt = dataLumi/(nEvts_Pythia[6] / xSec_Pythia[6]);
		} else {
			cout << red << "WARNING! UNKNOWN PYTHIA dataset name, " << datasetname 
				<< " given. returning lumiWgt = 1 !!!" << clearatt << endl;  
		}
		
	//QCD MG Samples
	} else if (datasetname.find("MGPythia") != string::npos )
	{
		if (datasetname.find("QCD_HT_250To500") != string::npos)
		{
			cout << "Lum for MG qcd 1 used." << endl; 
			lumiWgt = dataLumi/(nEvts_MG[0] / xSec_MG[0]);
		} else if (datasetname.find("QCD_HT_500To1000") != string::npos)
		{
			cout << "Lum for MG qcd 2 used." << endl; 
			lumiWgt = dataLumi/(nEvts_MG[1] / xSec_MG[1]);
		} else if (datasetname.find("QCD_HT_1000ToInf") != string::npos)
		{
			cout << "Lum for MG qcd 3 used." << endl; 
			lumiWgt = dataLumi/(nEvts_MG[2] / xSec_MG[2]);
		} else {
			cout << red << "WARNING! UNKNOWN MadGraph dataset name, " << datasetname 
				<< " given. returning lumiWgt = 1 !!!" << clearatt << endl;  
		}
	
	} else if (datasetname.find("ZJets_400_HT_inf") != string::npos )
	{
		lumiWgt = dataLumi/ (1006922.0/6.26); //NNLO
		cout << "Lumi weight of " << lumiWgt << " for ZJets_400_HT_inf is using for dataset " << datasetname << endl;

	} else if (datasetname.find("TT_CT10") != string::npos )
	{
		lumiWgt = dataLumi/ (21693253/234.0); //NNLO
		cout << "Lumi weight of " << lumiWgt << " for TT_CT10 is using for dataset " << datasetname << endl;

	} else if (datasetname.find("WJetsToLNu_HT-400ToInf") != string::npos )
	{
		lumiWgt = dataLumi/ (4971837/30.08); //NNLO
		cout << "Lumi weight of " << lumiWgt << " for WJetsToLNu_HT-400ToInf is using for dataset " << datasetname << endl;

	} else if (datasetname.find("WJetsToLNu_HT-300To400") != string::npos )
	{
		lumiWgt = dataLumi/ (1028198/30.08); //NNLO
		cout << "Lumi weight of " << lumiWgt << " for WJetsToLNu_HT-300To400 is using for dataset " << datasetname << endl;
	} else 
	{
		cout << red << "WARNING! UNKNOWN dataset name, " << datasetname 
				<< " given. returning lumiWgt = 1 !!!" << clearatt << endl;  
	}

	return lumiWgt;
}

/********************************************************************
 *  BOOK HISTOGRAMS
 *******************************************************************/
void FactorizationBySmearing::BookHistogram(TFile *oFile, const bool mcFlag)
{
	//oFile = new TFile(outFileName, "recreate");
	TDirectory *dir = oFile->mkdir("Hist");

	//if (bDEBUG) cout << __FUNCTION__ << " dir = " << dir->GetName() << endl;

	for (unsigned jetbin=0; jetbin < JetBins_.size(); ++jetbin)
	{
		const unsigned njetmin = JetBins_.at(jetbin).first;
		const unsigned njetmax = JetBins_.at(jetbin).second;
		vector<vector<Hist_t> > htcoll;

		for (unsigned htbin=0; htbin < HtBins_.size() -1; ++htbin)
		{
			vector<Hist_t> mhtcoll;
			for (unsigned mhtbin=0; mhtbin < MhtBins_.size() -1; ++mhtbin)
			{
				const pair<float, float> htrange (HtBins_.at(htbin), HtBins_.at(htbin+1));
				const pair<float, float> mhtrange (MhtBins_.at(mhtbin), MhtBins_.at(mhtbin+1));
				const float htmin  = HtBins_.at(htbin);
				const float htmax  = HtBins_.at(htbin+1);
				const float mhtmin = MhtBins_.at(mhtbin);
				const float mhtmax = MhtBins_.at(mhtbin+1);
				stringstream folder;
				folder << "Njet" << njetmin << "to" << njetmax << "HT" << htmin << "to" << htmax
					<<  "MHT" << mhtmin << "to" << mhtmax;
				TDirectory *subdir = dir->mkdir(folder.str().c_str());
				subdir->cd();

				Hist_t hist;
				GetHist(subdir, hist, JetBins_.at(jetbin), htrange, mhtrange, mcFlag);
				mhtcoll.push_back(hist);
			}
			htcoll.push_back(mhtcoll);
		}
		Hist.push_back(htcoll);
	}

}

void FactorizationBySmearing::GetHist(TDirectory *dir, Hist_t& hist, 
					const pair<unsigned, unsigned> njetrange,
					const pair<unsigned, unsigned> htrange,
					const pair<unsigned, unsigned> mhtrange,
					const bool mcFlag
					)
{
	stringstream htmhtrange;
	htmhtrange << njetrange.first << "-" << njetrange.second << "Jets, " 
			<< htrange.first  << "<HT<"  << htrange.second << " GeV, " 
			<< mhtrange.first << "<MHT<" << mhtrange.second << " GeV";
	
	vector<string> jetcoll;
	jetcoll.push_back("reco");
	jetcoll.push_back("gen");
	jetcoll.push_back("smeared");


	const int nBins_mht = 200; const double min_mht = 0, max_mht = 1500;
	const int nBins_ht  = 100; const double min_ht  = 0, max_ht  = 5000;
	//for ratio plot
	//const double evt_mht_max = 1100, evt_mht_bins = 1100;
	const double evt_mht_max = 1100, evt_mht_bins = 550;
	const double evt_ht_max  = 4000, evt_ht_bins  = 800;
	const float npassFailHistBins  = 14;
	const float passFailHistBins[] = {50,60,70,80,90,100,110,120,140,160,200,350,500,800,1000};
	//const float npassFail_min = 50, npassFail_max = 1050, npassFail_nbins = 100;
	const float npassFail_min = 50, npassFail_max = 1550, npassFail_nbins = 150;
	const int passFailBinOption = 1;

	//for evtWeight hist 
	int evtWgt_nbins = 100;
	double evtWgt_min = 0, evtWgt_max = 10;
	if (! mcFlag) {
		evtWgt_nbins = 1100;
		evtWgt_min = 0; evtWgt_max = 1100;
	}


	for (unsigned i = 0;i < jetcoll.size(); ++i)
	{
		stringstream njet50eta2p5title, njet30eta5p0title, httitle, mhttitle;
		njet50eta2p5title << htmhtrange.str().c_str() << ";Njets [ET>50 GeV,| #eta |<2.5];Events;";
		njet30eta5p0title << htmhtrange.str().c_str() << ";Njets [ET>30 GeV,| #eta |<5.0];Events;";
		httitle  << htmhtrange.str().c_str() << ";HT [3 Jets, PT>50 GeV, | #eta |<2.5];Events;";
		mhttitle << htmhtrange.str().c_str() << ";MET [PT>30 GeV, | #eta |<5.0];Events;";

		stringstream njet50eta2p5name, njet30eta5p0name, htname, mhtname;
		njet50eta2p5name << jetcoll.at(i) << "_njet50eta2p5";
		njet30eta5p0name << jetcoll.at(i) << "_njet30eta5p0";
		htname  << jetcoll.at(i) << "_ht";
		mhtname << jetcoll.at(i) << "_met";

		stringstream dphiminname, dphiminvsmhtname;
		stringstream dphimintitle, dphiminvsmhttitle;
		dphimintitle << htmhtrange.str().c_str() << ";#Delta #Phi_{min};Events;";
		dphiminname  << jetcoll.at(i) << "_dphimin";
		dphiminvsmhttitle << htmhtrange.str().c_str() << ";MET;#Delta #Phi_{min};";
		dphiminvsmhtname  << jetcoll.at(i) << "_dphiminVsMet";
		
		stringstream evtWeight_name, nvtx_name;	
		stringstream evtWeight_title, nvtx_title;	
		
		evtWeight_name << jetcoll.at(i) << "_evtWeight";
		nvtx_name << jetcoll.at(i) << "_nvtx";

		evtWeight_title << "Total Event Weight;;;";
		nvtx_title << "Nvtx;;;";

		stringstream nbjets_name, nbjets_title, bjetmass_name, bjetmass_title, bjetpt_name, bjetpt_title;
		stringstream m123_name, m123_title, m23overm123_name, m23overm123_title, mt2_name, mt2_title;
		stringstream mtb_name, mtb_title, mtt_name, mtt_title, mtb_p_mtt_name, mtb_p_mtt_title;

		nbjets_name << jetcoll.at(i) << "_nbjets";
		bjetmass_name << jetcoll.at(i) << "_bjetMass";
		bjetpt_name << jetcoll.at(i) << "_bjetPt";
		m123_name << jetcoll.at(i) << "_M123";
		m23overm123_name << jetcoll.at(i) << "_M23overM123";
		mt2_name << jetcoll.at(i) << "_MT2";
		mtb_name << jetcoll.at(i) << "_MTb";
		mtt_name << jetcoll.at(i) << "_MTt";
		mtb_p_mtt_name << jetcoll.at(i) << "_MTb_p_MTt";

		nbjets_title << ";Number of b-jets;Events;";
		bjetmass_title << ";Mass [b-jet];Events;";
		bjetpt_title << "; P_{T} of b-jet;Events;";
		m123_title << ";M123 [Top Mass];Events;";
		m23overm123_title << ";M23/M123;Events;";
		mt2_title << ";MT2;Events;";
		mtb_title << ";MTb;Events;";
		mtt_title << ";MTt;Events;";
		mtb_p_mtt_title << ";MTb+1/2 * MTt;Events;";


		if (i==0)
		{
			hist.hv_RecoEvt.h_Njet50eta2p5 = new TH1D(njet50eta2p5name.str().c_str(), njet50eta2p5title.str().c_str(), 20, 0, 20);
			hist.hv_RecoEvt.h_Njet30eta5p0 = new TH1D(njet30eta5p0name.str().c_str(), njet30eta5p0title.str().c_str(), 20, 0, 20);
			hist.hv_RecoEvt.h_Mht          = new TH1D(mhtname.str().c_str()         , mhttitle.str().c_str()         , nBins_mht, min_mht, max_mht);
			hist.hv_RecoEvt.h_Ht           = new TH1D(htname.str().c_str()          , httitle.str().c_str()          , nBins_ht , min_ht ,  max_ht);
			hist.hv_RecoEvt.h_DphiMin      = new TH1D(dphiminname.str().c_str()     , dphimintitle.str().c_str()     , 125, 0, 2.5);
			hist.hv_RecoEvt.h_DphiMinVsMht = new TH2D(dphiminvsmhtname.str().c_str(), dphiminvsmhttitle.str().c_str(), 1400, 0, 350, 250, 0, 2.5);
			hist.hv_RecoEvt.h_evtWeight    = new TH1D(evtWeight_name.str().c_str()  , evtWeight_title.str().c_str()  , evtWgt_nbins, evtWgt_min, evtWgt_max); 
			hist.hv_RecoEvt.h_nVtx         = new TH1D(nvtx_name.str().c_str()       , nvtx_title.str().c_str()       , 40, 0, 40);
			hist.hv_RecoEvt.h_nbjets       = new TH1D(nbjets_name.str().c_str()     , nbjets_title.str().c_str()     , 10, 0, 10);
			hist.hv_RecoEvt.h_bjetMass     = new TH1D(bjetmass_name.str().c_str()   , bjetmass_title.str().c_str()   , 120, 0, 600);
			hist.hv_RecoEvt.h_bjetPt       = new TH1D(bjetpt_name.str().c_str()     , bjetpt_title.str().c_str()     , 60, 0, 300);
			hist.hv_RecoEvt.h_M123         = new TH1D(m123_name.str().c_str()       , m123_title.str().c_str()       , 120, 0, 600);
			hist.hv_RecoEvt.h_M23OverM123  = new TH1D(m23overm123_name.str().c_str(), m23overm123_title.str().c_str()       , 100, 0, 2);
			hist.hv_RecoEvt.h_MT2          = new TH1D(mt2_name.str().c_str()        , mt2_title.str().c_str()        , 100, 0, 1000);
			hist.hv_RecoEvt.h_MTb          = new TH1D(mtb_name.str().c_str()        , mtb_title.str().c_str()        , 40, 0, 400);
			hist.hv_RecoEvt.h_MTt          = new TH1D(mtt_name.str().c_str()        , mtt_title.str().c_str()        , 40, 0, 400);
			hist.hv_RecoEvt.h_MTb_p_MTt    = new TH1D(mtb_p_mtt_name.str().c_str()  , mtb_p_mtt_title.str().c_str()  , 100, 0, 1000);

			hist.hv_RecoEvt.h_Njet50eta2p5->Sumw2();
			hist.hv_RecoEvt.h_Njet30eta5p0->Sumw2();
			hist.hv_RecoEvt.h_Mht->Sumw2();
			hist.hv_RecoEvt.h_Ht->Sumw2();
			hist.hv_RecoEvt.h_DphiMin->Sumw2();
			hist.hv_RecoEvt.h_DphiMinVsMht->Sumw2();
			hist.hv_RecoEvt.h_evtWeight->Sumw2();
			hist.hv_RecoEvt.h_nVtx->Sumw2();

			hist.hv_RecoEvt.h_nbjets->Sumw2();
			hist.hv_RecoEvt.h_bjetMass->Sumw2();
			hist.hv_RecoEvt.h_bjetPt->Sumw2();
			hist.hv_RecoEvt.h_M123->Sumw2();
			hist.hv_RecoEvt.h_M23OverM123->Sumw2();
			hist.hv_RecoEvt.h_MT2->Sumw2();
			hist.hv_RecoEvt.h_MTb->Sumw2();
			hist.hv_RecoEvt.h_MTt->Sumw2();
			hist.hv_RecoEvt.h_MTb_p_MTt->Sumw2();

			GetJetHist(hist.hv_RecoJets,jetcoll.at(i), htmhtrange.str());

			//pass variations
			//
			for (unsigned j=0; j < vDphiVariations.size(); ++j)
			{
				const float dphival = vDphiVariations.at(j); 

				stringstream pass_name, fail_name, passFineBin_name, failFineBin_name, pass_trigPres_name, fail_trigPres_name;
				stringstream pass_title, fail_title, passFineBin_title, failFineBin_title, pass_trigPres_title, fail_trigPres_title;

				pass_name << jetcoll.at(i) << "_pass" << j;
				fail_name << jetcoll.at(i) << "_fail" << j;
				passFineBin_name << jetcoll.at(i) << "_passFineBin" << j;
				failFineBin_name << jetcoll.at(i) << "_failFineBin" << j;
				pass_trigPres_name << jetcoll.at(i) << "_pass_trigPrescales" << j;
				fail_trigPres_name << jetcoll.at(i) << "_fail_trigPrescales" << j;

				pass_title <<"Events with #Delta#Phi_{min}>" << dphival << ";MHT;Events;";
				fail_title <<"Events with #Delta#Phi_{min}<" << dphival << ";MHT;Events;";
				passFineBin_title <<"Events with #Delta#Phi_{min}>" << dphival << ";MHT;Events;";
				failFineBin_title <<"Events with #Delta#Phi_{min}<" << dphival << ";MHT;Events;";

				pass_trigPres_title <<"Trig Prescales of Events with #Delta#Phi_{min}>" << dphival << ";MHT;Events;";
				fail_trigPres_title <<"Trig Prescales of Events with #Delta#Phi_{min}<" << dphival << ";MHT;Events;";

				if (passFailBinOption == 1)
				{
					hist.hv_RecoEvt.pass.push_back( new TH1D(pass_name.str().c_str(), pass_title.str().c_str(), npassFail_nbins, npassFail_min, npassFail_max) ); 
					hist.hv_RecoEvt.fail.push_back( new TH1D(fail_name.str().c_str(), fail_title.str().c_str(), npassFail_nbins, npassFail_min, npassFail_max) );
				} else {
					hist.hv_RecoEvt.pass.push_back( new TH1D(pass_name.str().c_str(), pass_title.str().c_str(), npassFailHistBins, passFailHistBins) ); 
					hist.hv_RecoEvt.fail.push_back( new TH1D(fail_name.str().c_str(), fail_title.str().c_str(), npassFailHistBins, passFailHistBins) );
				}
				hist.hv_RecoEvt.passFineBin.push_back( new TH1D(passFineBin_name.str().c_str(), passFineBin_title.str().c_str(), evt_mht_bins, 0, evt_mht_max) ); 
				hist.hv_RecoEvt.failFineBin.push_back( new TH1D(failFineBin_name.str().c_str(), failFineBin_title.str().c_str(), evt_mht_bins, 0, evt_mht_max) );

				hist.hv_RecoEvt.pass_trigPrescales.push_back( new TH1D(pass_trigPres_name.str().c_str(), pass_trigPres_title.str().c_str(),1000,0,1000) );
				hist.hv_RecoEvt.fail_trigPrescales.push_back( new TH1D(fail_trigPres_name.str().c_str(), fail_trigPres_title.str().c_str(),1000,0,1000) );

				hist.hv_RecoEvt.pass.at(j)->Sumw2(); 
				hist.hv_RecoEvt.fail.at(j)->Sumw2(); 
				hist.hv_RecoEvt.passFineBin.at(j)->Sumw2(); 
				hist.hv_RecoEvt.failFineBin.at(j)->Sumw2(); 
				hist.hv_RecoEvt.pass_trigPrescales.at(j)->Sumw2(); 
				hist.hv_RecoEvt.fail_trigPrescales.at(j)->Sumw2(); 

			}

			if (passFailBinOption == 1)
			{
				hist.hv_RecoEvt.sidebandSyst[0] = new TH1D("reco_sidebandSyst1","Reco Events: Sideband Syst1: FAIL from dphi (j1 & j2 , mht) >0.5 and dphi (j3, mht)>0.1 ", npassFail_nbins, npassFail_min, npassFail_max);
				hist.hv_RecoEvt.sidebandSyst[1] = new TH1D("reco_sidebandSyst2","Reco Events: Sideband Syst2: FAIL from dphi (j1 & j2 , mht) >0.5 and dphi (j3, mht)>0.1 ", npassFail_nbins, npassFail_min, npassFail_max);
			} else {
				hist.hv_RecoEvt.sidebandSyst[0] = new TH1D("reco_sidebandSyst1","Reco Events: Sideband Syst1: FAIL from dphi (j1 & j2 , mht) >0.5 and dphi (j3, mht)>0.1 ", npassFailHistBins, passFailHistBins);
				hist.hv_RecoEvt.sidebandSyst[1] = new TH1D("reco_sidebandSyst2","Reco Events: Sideband Syst2: FAIL from dphi (j1 & j2 , mht) >0.5 and dphi (j3, mht)>0.1 ", npassFailHistBins, passFailHistBins);
			}
			hist.hv_RecoEvt.sidebandSystFineBin[0] = new TH1D("reco_sidebandSyst1_fineBin","Reco Events: Sideband Syst1 (FineBinned) : FAIL from dphi (j1 & j2 , mht) >0.5 and dphi (j3, mht)>0.1 ", evt_mht_bins, 0, evt_mht_max);
			hist.hv_RecoEvt.sidebandSystFineBin[1] = new TH1D("reco_sidebandSyst2_fineBin","Reco Events: Sideband Syst2 (FineBinned) : FAIL from dphi (j1 & j2 , mht) >0.4 and dphi (j3, mht)>0.2 ", evt_mht_bins, 0, evt_mht_max);

			if (passFailBinOption == 1)
			{
				hist.hv_RecoEvt.signal = new TH1D("reco_signal" ,"Reco Events: Signal Region", npassFail_nbins, npassFail_min, npassFail_max);
			} else {
				hist.hv_RecoEvt.signal = new TH1D("reco_signal" ,"Reco Events: Signal Region", npassFailHistBins, passFailHistBins);
			}
			hist.hv_RecoEvt.signalFineBin        = new TH1D("reco_signalFineBin" ," Reco Events: Signal Region", evt_mht_bins, 0, evt_mht_max);
			hist.hv_RecoEvt.signal_trigPrescales = new TH1D("reco_signal_trigPrescales", "Reco Events: Prescales of Signal Region Events",1000,0,1000);

			hist.hv_RecoEvt.sidebandSyst[0]->Sumw2();
			hist.hv_RecoEvt.sidebandSyst[1]->Sumw2();
			hist.hv_RecoEvt.sidebandSystFineBin[0]->Sumw2();
			hist.hv_RecoEvt.sidebandSystFineBin[1]->Sumw2();

			hist.hv_RecoEvt.signal->Sumw2();
			hist.hv_RecoEvt.signalFineBin->Sumw2();
			hist.hv_RecoEvt.signal_trigPrescales->Sumw2();

		} else if (i==1)
		{
			hist.hv_GenEvt.h_Njet50eta2p5 = new TH1D(njet50eta2p5name.str().c_str(), njet50eta2p5title.str().c_str(), 20, 0, 20);
			hist.hv_GenEvt.h_Njet30eta5p0 = new TH1D(njet30eta5p0name.str().c_str(), njet30eta5p0title.str().c_str(), 20, 0, 20);
			hist.hv_GenEvt.h_Mht          = new TH1D(mhtname.str().c_str()         , mhttitle.str().c_str()         , nBins_mht, min_mht, max_mht);
			hist.hv_GenEvt.h_Ht           = new TH1D(htname.str().c_str()          , httitle.str().c_str()          , nBins_ht, min_ht, max_ht);
			hist.hv_GenEvt.h_DphiMin      = new TH1D(dphiminname.str().c_str()     , dphimintitle.str().c_str()     , 120, 0, 2.5);
			hist.hv_GenEvt.h_DphiMinVsMht = new TH2D(dphiminvsmhtname.str().c_str(), dphiminvsmhttitle.str().c_str(), 1400, 0, 350, 250, 0, 2.5);
			hist.hv_GenEvt.h_evtWeight    = new TH1D(evtWeight_name.str().c_str()  , evtWeight_title.str().c_str()  , evtWgt_nbins,evtWgt_min,evtWgt_max);
			hist.hv_GenEvt.h_nVtx         = new TH1D(nvtx_name.str().c_str()       , nvtx_title.str().c_str()       , 40,0,40);

			hist.hv_GenEvt.h_nbjets       = new TH1D(nbjets_name.str().c_str()     , nbjets_title.str().c_str()     , 10, 0, 10);
			hist.hv_GenEvt.h_bjetMass     = new TH1D(bjetmass_name.str().c_str()   , bjetmass_title.str().c_str()   , 120, 0, 600);
			hist.hv_GenEvt.h_bjetPt       = new TH1D(bjetpt_name.str().c_str()     , bjetpt_title.str().c_str()     , 60, 0, 300);
			hist.hv_GenEvt.h_M123         = new TH1D(m123_name.str().c_str()       , m123_title.str().c_str()       , 120, 0, 600);
			hist.hv_GenEvt.h_M23OverM123  = new TH1D(m23overm123_name.str().c_str(), m23overm123_title.str().c_str()       , 100, 0, 2);
			hist.hv_GenEvt.h_MT2          = new TH1D(mt2_name.str().c_str()        , mt2_title.str().c_str()        , 100, 0, 1000);
			hist.hv_GenEvt.h_MTb          = new TH1D(mtb_name.str().c_str()        , mtb_title.str().c_str()        , 40, 0, 400);
			hist.hv_GenEvt.h_MTt          = new TH1D(mtt_name.str().c_str()        , mtt_title.str().c_str()        , 40, 0, 400);
			hist.hv_GenEvt.h_MTb_p_MTt    = new TH1D(mtb_p_mtt_name.str().c_str()  , mtb_p_mtt_title.str().c_str()  , 100, 0, 1000);

			hist.hv_GenEvt.h_Njet50eta2p5->Sumw2();
			hist.hv_GenEvt.h_Njet30eta5p0->Sumw2();
			hist.hv_GenEvt.h_Mht->Sumw2();
			hist.hv_GenEvt.h_Ht->Sumw2();
			hist.hv_GenEvt.h_DphiMin->Sumw2();
			hist.hv_GenEvt.h_DphiMinVsMht->Sumw2();
			hist.hv_GenEvt.h_evtWeight->Sumw2();
			hist.hv_GenEvt.h_nVtx->Sumw2();
			
			hist.hv_GenEvt.h_nbjets->Sumw2();
			hist.hv_GenEvt.h_bjetMass->Sumw2();
			hist.hv_GenEvt.h_bjetPt->Sumw2();
			hist.hv_GenEvt.h_M123->Sumw2();
			hist.hv_GenEvt.h_M23OverM123->Sumw2();
			hist.hv_GenEvt.h_MT2->Sumw2();
			hist.hv_GenEvt.h_MTb->Sumw2();
			hist.hv_GenEvt.h_MTt->Sumw2();
			hist.hv_GenEvt.h_MTb_p_MTt->Sumw2();
			GetJetHist(hist.hv_GenJets,jetcoll.at(i), htmhtrange.str());

		} else if (i==2)
		{
			hist.hv_SmearedEvt.h_Njet50eta2p5 = new TH1D(njet50eta2p5name.str().c_str(), njet50eta2p5title.str().c_str(), 20, 0, 20);
			hist.hv_SmearedEvt.h_Njet30eta5p0 = new TH1D(njet30eta5p0name.str().c_str(), njet30eta5p0title.str().c_str(), 20, 0, 20);
			hist.hv_SmearedEvt.h_Mht          = new TH1D(mhtname.str().c_str()         , mhttitle.str().c_str()         , nBins_mht, min_mht, max_mht);
			hist.hv_SmearedEvt.h_Ht           = new TH1D(htname.str().c_str()          , httitle.str().c_str()          ,  nBins_ht,  min_ht,  max_ht);
			hist.hv_SmearedEvt.h_DphiMin      = new TH1D(dphiminname.str().c_str()     , dphimintitle.str().c_str()     , 125, 0, 2.5);
			hist.hv_SmearedEvt.h_DphiMinVsMht = new TH2D(dphiminvsmhtname.str().c_str(), dphiminvsmhttitle.str().c_str(), 1400, 0, 350, 250, 0, 2.5);
			hist.hv_SmearedEvt.h_evtWeight    = new TH1D(evtWeight_name.str().c_str()  , evtWeight_title.str().c_str()  , evtWgt_nbins, evtWgt_min, evtWgt_max);
			hist.hv_SmearedEvt.h_nVtx         = new TH1D(nvtx_name.str().c_str()       , nvtx_title.str().c_str()       , 40, 0, 40);

			hist.hv_SmearedEvt.h_nbjets       = new TH1D(nbjets_name.str().c_str()     , nbjets_title.str().c_str()     , 10, 0, 10);
			hist.hv_SmearedEvt.h_bjetMass     = new TH1D(bjetmass_name.str().c_str()   , bjetmass_title.str().c_str()   , 120, 0, 600);
			hist.hv_SmearedEvt.h_bjetPt       = new TH1D(bjetpt_name.str().c_str()     , bjetpt_title.str().c_str()     , 60, 0, 300);
			hist.hv_SmearedEvt.h_M123         = new TH1D(m123_name.str().c_str()       , m123_title.str().c_str()       , 120, 0, 600);
			hist.hv_SmearedEvt.h_M23OverM123  = new TH1D(m23overm123_name.str().c_str(), m23overm123_title.str().c_str()       , 100, 0, 2);
			hist.hv_SmearedEvt.h_MT2          = new TH1D(mt2_name.str().c_str()        , mt2_title.str().c_str()        , 100, 0, 1000);
			hist.hv_SmearedEvt.h_MTb          = new TH1D(mtb_name.str().c_str()        , mtb_title.str().c_str()        , 40, 0, 400);
			hist.hv_SmearedEvt.h_MTt          = new TH1D(mtt_name.str().c_str()        , mtt_title.str().c_str()        , 40, 0, 400);
			hist.hv_SmearedEvt.h_MTb_p_MTt    = new TH1D(mtb_p_mtt_name.str().c_str()  , mtb_p_mtt_title.str().c_str()  , 100, 0, 1000);


			hist.hv_SmearedEvt.h_Njet50eta2p5->Sumw2();
			hist.hv_SmearedEvt.h_Njet30eta5p0->Sumw2();
			hist.hv_SmearedEvt.h_Mht->Sumw2();
			hist.hv_SmearedEvt.h_Ht->Sumw2();
			hist.hv_SmearedEvt.h_DphiMin->Sumw2();
			hist.hv_SmearedEvt.h_DphiMinVsMht->Sumw2();
			hist.hv_SmearedEvt.h_evtWeight->Sumw2();
			hist.hv_SmearedEvt.h_nVtx->Sumw2();

			hist.hv_SmearedEvt.h_nbjets->Sumw2();
			hist.hv_SmearedEvt.h_bjetMass->Sumw2();
			hist.hv_SmearedEvt.h_bjetPt->Sumw2();
			hist.hv_SmearedEvt.h_M123->Sumw2();
			hist.hv_SmearedEvt.h_M23OverM123->Sumw2();
			hist.hv_SmearedEvt.h_MT2->Sumw2();
			hist.hv_SmearedEvt.h_MTb->Sumw2();
			hist.hv_SmearedEvt.h_MTt->Sumw2();
			hist.hv_SmearedEvt.h_MTb_p_MTt->Sumw2();
			GetJetHist(hist.hv_SmearedJets, jetcoll.at(i), htmhtrange.str());

			//pass variations
			//
			for (unsigned j=0; j < vDphiVariations.size(); ++j)
			{
				const float dphival = vDphiVariations.at(j); 

				stringstream pass_name, fail_name, passFineBin_name, failFineBin_name;
				stringstream pass_title, fail_title, passFineBin_title, failFineBin_title;

				pass_name << jetcoll.at(i) << "_pass" << j;
				fail_name << jetcoll.at(i) << "_fail" << j;
				passFineBin_name << jetcoll.at(i) << "_passFineBin" << j;
				failFineBin_name << jetcoll.at(i) << "_failFineBin" << j;
				pass_title <<"Smeared events with #Delta#Phi_{min}>" << dphival << ";MHT;Events;";
				fail_title <<"Smeared events with #Delta#Phi_{min}<" << dphival << ";MHT;Events;";
				passFineBin_title <<"Smeared events with #Delta#Phi_{min}>" << dphival << ";MHT;Events;";
				failFineBin_title <<"Smeared events with #Delta#Phi_{min}<" << dphival << ";MHT;Events;";

				if (passFailBinOption == 1)
				{
					hist.hv_SmearedEvt.pass.push_back( new TH1D(pass_name.str().c_str(), pass_title.str().c_str(), npassFail_nbins, npassFail_min, npassFail_max) ); 
					hist.hv_SmearedEvt.fail.push_back( new TH1D(fail_name.str().c_str(), fail_title.str().c_str(), npassFail_nbins, npassFail_min, npassFail_max) );
				} else {
					hist.hv_SmearedEvt.pass.push_back( new TH1D(pass_name.str().c_str(), pass_title.str().c_str(), npassFailHistBins, passFailHistBins) ); 
					hist.hv_SmearedEvt.fail.push_back( new TH1D(fail_name.str().c_str(), fail_title.str().c_str(), npassFailHistBins, passFailHistBins) );
				}

				hist.hv_SmearedEvt.passFineBin.push_back( new TH1D(passFineBin_name.str().c_str(), passFineBin_title.str().c_str(), evt_mht_bins, 0, evt_mht_max) ); 
				hist.hv_SmearedEvt.failFineBin.push_back( new TH1D(failFineBin_name.str().c_str(), failFineBin_title.str().c_str(), evt_mht_bins, 0, evt_mht_max) );

				hist.hv_SmearedEvt.pass.at(j)->Sumw2(); 
				hist.hv_SmearedEvt.fail.at(j)->Sumw2(); 
				hist.hv_SmearedEvt.passFineBin.at(j)->Sumw2(); 
				hist.hv_SmearedEvt.failFineBin.at(j)->Sumw2(); 
			}

			if (passFailBinOption == 1)
			{
				hist.hv_SmearedEvt.sidebandSyst[0] = new TH1D("smeared_sidebandSyst1"," Smeared Events: Sideband Syst1: FAIL from dphi (j1 & j2 , mht) >0.5 and dphi (j3, mht)>0.1 ", npassFail_nbins, npassFail_min, npassFail_max);
				hist.hv_SmearedEvt.sidebandSyst[1] = new TH1D("smeared_sidebandSyst2"," Smeared Events: Sideband Syst2: FAIL from dphi (j1 & j2 , mht) >0.5 and dphi (j3, mht)>0.1 ", npassFail_nbins, npassFail_min, npassFail_max);
			} else {
				hist.hv_SmearedEvt.sidebandSyst[0] = new TH1D("smeared_sidebandSyst1"," Smeared Events: Sideband Syst1: FAIL from dphi (j1 & j2 , mht) >0.5 and dphi (j3, mht)>0.1 ", npassFailHistBins, passFailHistBins);
				hist.hv_SmearedEvt.sidebandSyst[1] = new TH1D("smeared_sidebandSyst2"," Smeared Events: Sideband Syst2: FAIL from dphi (j1 & j2 , mht) >0.5 and dphi (j3, mht)>0.1 ", npassFailHistBins, passFailHistBins);
			}
			hist.hv_SmearedEvt.sidebandSystFineBin[0] = new TH1D("smeared_sidebandSyst1_fineBin"," Smeared Events:Sideband Syst1 (FineBinned) : FAIL from dphi (j1 & j2 , mht) >0.5 and dphi (j3, mht)>0.1 ", evt_mht_bins, 0, evt_mht_max);
			hist.hv_SmearedEvt.sidebandSystFineBin[1] = new TH1D("smeared_sidebandSyst2_fineBin"," Smeared Events:Sideband Syst2 (FineBinned) : FAIL from dphi (j1 & j2 , mht) >0.4 and dphi (j3, mht)>0.2 ", evt_mht_bins, 0, evt_mht_max);

			if (passFailBinOption == 1)
			{
				hist.hv_SmearedEvt.signal = new TH1D("smeared_signal", " Smeared Events:Signal Region", npassFail_nbins, npassFail_min, npassFail_max);
			} else {
				hist.hv_SmearedEvt.signal = new TH1D("smeared_signal", " Smeared Events:Signal Region", npassFailHistBins, passFailHistBins);
			}
			hist.hv_SmearedEvt.signalFineBin = new TH1D("smeared_signalFineBin", " Smeared Events:Signal Region", evt_mht_bins,0, evt_mht_max);

			hist.hv_SmearedEvt.sidebandSyst[0]->Sumw2();
			hist.hv_SmearedEvt.sidebandSyst[1]->Sumw2();
			hist.hv_SmearedEvt.sidebandSystFineBin[0]->Sumw2();
			hist.hv_SmearedEvt.sidebandSystFineBin[1]->Sumw2();

			hist.hv_SmearedEvt.signal->Sumw2();
			hist.hv_SmearedEvt.signalFineBin->Sumw2();

		}
	}


}

void FactorizationBySmearing::GetJetHist(vector<JetHist_t>& Hist, const string jetcoll,
					const string htmhtrange)
{
	//save up to 10 leading jets
	const int njets = 10;
	for (int i=0; i<njets; ++i)
	{
		stringstream ptname, etaname, phiname, dphiname; 
		stringstream pttitle, etatitle, phititle, dphititle; 
		const int index = i + 1;
		ptname   << jetcoll << "jet" << index << "_pt";
		etaname  << jetcoll << "jet" << index << "_eta";
		phiname  << jetcoll << "jet" << index << "_phi";
		dphiname << jetcoll << "jet" << index << "_dphi";
		pttitle  << htmhtrange << ";Jet-" << index << " P_{T} [GeV];Events;";
		etatitle << htmhtrange << ";Jet-" << index << " #eta;Events;";
		pttitle  << htmhtrange << ";Jet-" << index << " #phi ;Events;";
		pttitle  << htmhtrange << ";#delta #phi [Jet-" << index << ", MHT];Events;";

		JetHist_t hist;
		hist.h_Jet_pt   = new TH1D(ptname.str().c_str()  , pttitle.str().c_str()  , 150,  0.0, 1500.0);  
		hist.h_Jet_eta  = new TH1D(etaname.str().c_str() , etatitle.str().c_str() , 120, -6.0,    6.0);  
		hist.h_Jet_phi  = new TH1D(phiname.str().c_str() , phititle.str().c_str() ,  70, -3.5,    3.5);  
		hist.h_Jet_dphi = new TH1D(dphiname.str().c_str(), dphititle.str().c_str(),  70,    0,    3.5);  
		hist.h_Jet_pt->Sumw2();
		hist.h_Jet_eta->Sumw2();
		hist.h_Jet_phi->Sumw2();
		hist.h_Jet_dphi->Sumw2();
		Hist.push_back(hist);
	}

}

void FactorizationBySmearing::DivideByBinWidth(TH1* h)
{
	if (h == NULL)
	{
		cout << __FUNCTION__ << ": null pointer passed! retuning!." << endl; 
		return;
	}
	
	//normalization of underflow and overflow bins
	//does not seem accurate as the bin width seems undefined
	//for these two bins.
	for (int bin=0; bin<=h->GetNbinsX(); ++bin)
	{
		const double v  = h->GetBinContent(bin);
		const double e  = h->GetBinError(bin);
		const double w  = h->GetBinWidth(bin);
		const double nv = v/w;
		const double ne = e/w;
		h->SetBinContent(bin, nv);
		h->SetBinError(bin, ne);
	}
}
vector<TLorentzVector> FactorizationBySmearing::GetPt50Eta2p5Jets(const vector<TLorentzVector>& jets)
{
	vector<TLorentzVector> newjets;
	for (unsigned i=0; i< jets.size(); ++i)
	{
		TLorentzVector jetvec(jets.at(i));
		if (jetvec.Pt()<50 || fabs(jetvec.Eta())>2.5) continue;
		newjets.push_back(jetvec);
	}
	return newjets;
}
bool FactorizationBySmearing::PassCleaning()
{
	 //52x recipe recomments not using eeNoiseFilter due to overtagging 
	const bool pass = ( 
					(bool) t_beamHaloFilter 
				&& (bool) t_eeBadScFilter 
				//&& ( ! ((bool) t_eeNoiseFilter))
				&& (bool) t_greedyMuons
				&& (bool) t_hcalLaserEventFilter
				&& (bool) t_inconsistentMuons
				&& (bool) t_ra2EcalBEFilter
				&& (bool) t_ra2EcalTPFilter
				&& (bool) t_trackingFailureFilter
				&& (bool) t_HBHENoiseFilterRA2
				);

			stringstream strevt;
			strevt << t_EvtRun << ":" << t_EvtLS << ":" << t_EvtEvent;
			bool badEcalLaser = false;
			if (find(vBadEcalLaserEvts.begin(), vBadEcalLaserEvts.end(), strevt.str()) != vBadEcalLaserEvts.end())
			{
				badEcalLaser = true;
				++nBadEcalLaserEvts;
			}

/*	if (! pass) 
		cout << __FUNCTION__ << ":"
			<< t_beamHaloFilter << "/"
			<< t_eeBadScFilter << "/"
			<< t_eeNoiseFilter << "/"
			<< t_greedyMuons << "/"
			<< t_hcalLaserEventFilter << "/"
			<< t_inconsistentMuons << "/"
			<< t_ra2EcalBEFilter << "/"
			<< t_ra2EcalTPFilter << "/"
			<< t_trackingFailureFilter << "/"
			<< t_HBHENoiseFilterRA2
			<< endl;
*/

	return (pass && ! badEcalLaser);
}

void FactorizationBySmearing::TrigPrescaleWeight(bool &passTrig, double &weight) const
{
	/*****************************************************************
	 *  For DATA only. Need to find the highest prescaled HT trigger
	 *  from the given list of triggers
	 *  Assume the trigger names are given without wildcards!!!
	 ****************************************************************/

	unsigned highestTrigPrescale = 999999;
	bool fired = false;
	int nFiredTrigs = 0;
	string highestTrigPrescale_name("");
	for (unsigned j =0; j < vTriggersToUse.size(); ++j)
	{
		//cout << __LINE__ << ":: trigtouse[" << j << "] = " << vTriggersToUse.at(j) << endl; 
		for (unsigned i =0; i < t_firedTrigs->size(); ++i)
		{
			if (t_firedTrigs->at(i).find(vTriggersToUse.at(j)) != string::npos)
			{	//found fired trigger
				//cout << __LINE__ << ":: fired[" << t_firedTrigs->at(i) << "->" << t_firedTrigsPrescale->at(i) << "]" << endl; 
				if (t_firedTrigsPrescale->at(i) < highestTrigPrescale)
				{
					fired = true;
					++nFiredTrigs;
					highestTrigPrescale = t_firedTrigsPrescale->at(i);
					highestTrigPrescale_name = t_firedTrigs->at(i);
					//cout << __LINE__ << ":: fired[" << i << "] = " << t_firedTrigs->at(i) << "->" << t_firedTrigsPrescale->at(i) << endl; 
				}
			}
		}
	}

	passTrig = fired;
	if (fired) 
	{
		weight = (double) highestTrigPrescale;
		//if (nFiredTrigs>1) cout << ">>>>>>>>>>>> LOOK HERE <<<<<<<<<<<<" << endl;
		//if (weight>1) 
		//{
			//cout << ">>>>>>>>>>>> LOOK HERE <<<<<<<<<<<<" << endl;
			//cout << "\t" << __LINE__ << ": selected trig/prescale = " << highestTrigPrescale_name << "/" << highestTrigPrescale << ":";
			//PrintEventNumber();
		//}
	}

}

void FactorizationBySmearing::LoadBadHcalLaserEvents()
{
	ifstream hcalfile;
	hcalfile.open("/share/store/users/samantha/CMSSW_5_2_5/src/UserCode/RA2QCDvetoAna/flatrun/AllBadHCALLaser.txt");
	if (hcalfile.is_open())
	{
		string line;
		unsigned nLines = 0;
		while (! hcalfile.eof())
		{
			++nLines;
			getline(hcalfile,line);
			if (line.find(" ") != string::npos)
			{
				cout << "AllBadHCALLaser.txt file has spaces in line# " << nLines << "(" << line <<  ")! Removed them first.!" << endl;
				assert(false);
			} else
			{
				//cout << line  << endl;
				vBadHcalLaserEvts.push_back(line);
			}
		}
		cout << __FUNCTION__ << ": Read " << nLines << " lines of data from AllBadHCALLaser.txt" << endl;
		//cout << vBadHcalLaserEvts.size() << endl;
	} else {
		cout << __FUNCTION__ << ":AllBadHCALLaser.txt file not found!!!" << endl;
		assert (false);
	}

}

void FactorizationBySmearing::LoadBadEcalLaserEvents()
{
	//these are several bad events found to fail ecalLaserFilter
	vBadEcalLaserEvts.push_back("196453:734:630436469");
	vBadEcalLaserEvts.push_back("195552:192:318100959");
	vBadEcalLaserEvts.push_back("195552:721:975903831");
	vBadEcalLaserEvts.push_back("195552:1253:1444299746");
	vBadEcalLaserEvts.push_back("194912:1147:1612156042");
	vBadEcalLaserEvts.push_back("194912:743:1160203093");
	vBadEcalLaserEvts.push_back("196495:118:27004272");
	vBadEcalLaserEvts.push_back("194424:199:266603387");
	vBadEcalLaserEvts.push_back("194428:240:250267732");
	vBadEcalLaserEvts.push_back("196364:533:513663352");
	vBadEcalLaserEvts.push_back("195398:1345:1069637397");
	vBadEcalLaserEvts.push_back("205158:462:644937029");
	vBadEcalLaserEvts.push_back("205339:226:301238193");
	vBadEcalLaserEvts.push_back("203994:112:115834698");
	vBadEcalLaserEvts.push_back("207515:1012:1413733443");
	vBadEcalLaserEvts.push_back("208487:331:537336833");
	vBadEcalLaserEvts.push_back("190895:80:49855369");
	vBadEcalLaserEvts.push_back("191721:6:2368465");
}


void FactorizationBySmearing::PrintEventNumber() const
{
	cout << "Run/Ls/Evt = " << t_EvtRun << " / " 
		<<  t_EvtLS << " / " <<  t_EvtEvent << endl;
}

/********************************************************************
 *				C O N S T R U C T O R
 *******************************************************************/
FactorizationBySmearing::FactorizationBySmearing(
				const TString &inputFileList, 
				const char *outFileName
				) {

	TChain *tree = new TChain("treeMaker/tree");  

	if ( ! FillChain(tree, inputFileList) ) {
		std::cerr << "Cannot get the tree " << std::endl;
		assert(false);
	}

	Init(tree);

	smearFunc_ = 0;
	HtBins_.push_back(0);
//	HtBins_.push_back(500);
//	HtBins_.push_back(800);
//	HtBins_.push_back(1000);
//	HtBins_.push_back(1250);
//	HtBins_.push_back(1500);
	HtBins_.push_back(8000);

	MhtBins_.push_back(0);
/*	MhtBins_.push_back(200);
	MhtBins_.push_back(300);
	MhtBins_.push_back(450);
	MhtBins_.push_back(600);
*/	//MhtBins_.push_back(0);
	//MhtBins_.push_back(200);
	//MhtBins_.push_back(500);
	MhtBins_.push_back(8000);

	//jet bins: 2, 3-5, 6-7,>=8
	//JetBins_.push_back(make_pair(2,2));	
//	JetBins_.push_back(make_pair(3,5));	
//	JetBins_.push_back(make_pair(6,7));	
//	JetBins_.push_back(make_pair(8,1000));	
	JetBins_.push_back(make_pair(0,1000));	

	bNON_STD_MODE = false;
	nRecoJetEvts  = 0;
	nGenJetEvts   = 0;
	nSmearedJetEvts     = 0;
	nVectorInexWarnings = 0;
	uNTRIES = 1;  //default for testing

	//sanity check to have at least 1 bin in njet/ht/mht
	bool ready = true;
	if (JetBins_.size() < 1) { ready = ready && false; cout << __FUNCTION__ << ": Require at least one Jet bin!" << endl; }
	if (HtBins_.size()  < 2) { ready = ready && false; cout << __FUNCTION__ << ": Require at least one Ht bin!" << endl;  }
	if (MhtBins_.size() < 2) { ready = ready && false; cout << __FUNCTION__ << ": Require at least one Mht bin!" << endl; }

	outRootFile = new TFile(outFileName, "recreate");
	if (outRootFile->IsZombie())
	{
		cout << __FUNCTION__ << ": Unable to create output root file!" << endl;
		ready = false;
	}

	
	//difference variation of the dPhiMin selections.
	//make sure the book the correct number of histograms 
	//when these are changed!!!
	vDphiVariations.push_back(0.15);
	vDphiVariations.push_back(0.20);
	//vDphiVariations.push_back(0.25);
	//vDphiVariations.push_back(0.30);
	//vDphiVariations.push_back(0.35);
	//vDphiVariations.push_back(0.40);

	//List all the triggers to be used for data WITHOUT wildcards (i.e. * )


//2012AJuly13 rereco
//2012Aaug6
//2012bjULY13reco
//2012Cprompteco
//2012C_rereco
//HLTPathsByName_[0] = HLT_HT*");
//vTriggersToUse.push_back("HLT_HT200_v");
//vTriggersToUse.push_back("HLT_HT250_v");
//vTriggersToUse.push_back("HLT_HT300_v");
//vTriggersToUse.push_back("HLT_HT350_v");
//vTriggersToUse.push_back("HLT_HT400_v");
//vTriggersToUse.push_back("HLT_HT450_v");
//vTriggersToUse.push_back("HLT_HT500_v");
//vTriggersToUse.push_back("HLT_HT550_v");
//vTriggersToUse.push_back("HLT_HT650_v");
//vTriggersToUse.push_back("HLT_HT750_v");
//HLTPathsByName_[1] = HLT_PFHT*");
//Kristin confirmed that she is only using HLT_PFHT triggers.
//No jet triggers are used either. Jan 25th, 2013
vTriggersToUse.push_back("HLT_PFHT350_v");
vTriggersToUse.push_back("HLT_PFHT650_v");
vTriggersToUse.push_back("HLT_PFHT700_v");
vTriggersToUse.push_back("HLT_PFHT750_v");
vTriggersToUse.push_back("HLT_PFNoPUHT350_v");
vTriggersToUse.push_back("HLT_PFNoPUHT650_v");
vTriggersToUse.push_back("HLT_PFNoPUHT700_v");
vTriggersToUse.push_back("HLT_PFNoPUHT750_v");

	if (! ready) 
	{ 
		cout << __FUNCTION__ << ": Minimum run conditions failed. returning!!" << endl;
		assert (false);
	} 

	//BookHistogram(outFileName);
}

bool FactorizationBySmearing::PassHOfilter()
{

	for (int i=0; i < t_PFJetE->size(); ++i)
	{
		const float hofrac = t_PFJetHOEne->at(i)/t_PFJetE->at(i);
		if (hofrac>0.4) return false;  //fail
	}
	return true;

}


void FactorizationBySmearing::TripletSelector(const std::vector<TLorentzVector> & jets, 
			const vector<double> bDiscriminator, std::vector<TLorentzVector> & triplet, 
			std::vector<TLorentzVector> & rSystem, std::vector<TLorentzVector> & bJet, 
			double & M23OverM123, double & M123)
{

	/* one thing, whenever you calculate a MET in TLorentzVector
	 * set the mass of it to be Zero (force it)
	 * because the mass of MET is used for MT2 calculation.
	 * and you know the default MET in CMS has mass zero
	 * when you modify it, you probably can end up with a non-zero met mass
	 * so you have to force it to be zero manually
	 */

	const double Rmin_ = 0.85 * 80.385/173.5;
	const double Rmax_ = 1.25 * 80.385/173.5;
	const double arctanmin_ = 0.2;
	const double arctanmax_ = 1.3;
	const double m23OverM123Cut_ = 0.35;
	const double topMass_   = 173.5;

	const double mTop_  = 173.5;
	const double mWMin_ =  50.0;
	const double mWMax_ = 120.0;

	const double tripletJetPtCut_ = 30.0;
	const double tripletDRCut_    = 1.5;

	const double bJetPtCut_   = 30;
	const double bJetEtaCut_  = 2.4;
	//    bJetDisc = "combinedSecondaryVertexBJetTags";
	const double bJET_DISC_CUT_ = 0.679;

	double highPt  = 0;
	double highCSV = -1;
	unsigned highPtIndex  = -1;
	unsigned highCSVIndex = -1;

	vector<TLorentzVector> bJets;

	for (unsigned i = 0; i<jets.size(); i++)
	{
		const double disc = bDiscriminator.at(i);
		if ( jets[i].Pt() < bJetPtCut_ 
			|| abs(jets[i].Eta()) > bJetEtaCut_ 
			|| disc < bJET_DISC_CUT_) 
		{
			continue;
		}
		bJets.push_back(jets[i]);
	}
	// cou << __LINE__ << ": bJets = " << bJets.size() << endl;

	//if no b-jet is found based on bDiscrminator cut
	//use the jet with highest CSV value at the b-jet.
	if (bJets.size() == 0)
	{
		for (unsigned i = 0; i < jets.size(); i++)
		{
			/*         csv = jets[i].bDiscriminator(bJetDisc_.c_str());
						  cout<<"ETA: "<<jets[i].Eta()<<" PHI: "<<jets[i].phi()<<" PT: "<<jets[i].Pt()<<" CSV: "<<csv<<endl;
						  */
			if (abs(jets[i].Eta()) > bJetEtaCut_) continue;

//			const double pt  = jets[i].Pt();
			const double csv = bDiscriminator.at(i);   //????

/*			if (pt > highPt)
			{
				highPt = pt;
				highPtIndex = i;
			}
*/			if (csv > highCSV)
			{
				highCSV = csv;
				highCSVIndex = i;
			}
		}
	} //if no bjet is found

	if (highCSV != -1 && highCSV != -10)
	{
		bJets.push_back(jets[highCSVIndex]);
	}

	TLorentzVector p41(0,0,0,0);
	TLorentzVector p42(0,0,0,0);
	TLorentzVector p43(0,0,0,0);

	bool bJetPass              = false;
	bool bJetInTrip            = false;
	bool bJetOutTrip           = false;
	bool passDijet             = false;
	bool passM23OverM123       = false;
	bool atLeastOnePassed      = false;
	bool skipFailedM23OverM123 = false;

	unsigned selectedIndex = 0;
	unsigned rank          = 5;
	const double LN = -9999.99; 
	double bJetEta  = LN;
	double bJetPhi  = LN;
	double m123     = LN ;
	double m12      = LN;
	double m13      = LN ;
	double m23      = LN;
	double dTopMin  = 999999;


	vector<double> m123s;
	vector<double> m12s;
	vector<double> m13s;
	vector<double> m23s;
	vector<vector<int> > bJetInRIndices;

	vector<bool> bJetsBehave;
	vector<bool> passesDijetCuts;
	vector<bool> passesM23OverM123Cut;

	auto_ptr<vector<bool> > selectedPassDijetCuts(new vector<bool>());
	auto_ptr<vector<TLorentzVector> > selectedTriplet(new vector<TLorentzVector>());
	auto_ptr<vector<TLorentzVector> > selectedRSystem(new vector<TLorentzVector>());
	auto_ptr<vector<TLorentzVector> > selectedBJets(new vector<TLorentzVector>());
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

	for (int i = 0; i < int(jets.size()); i++)
	{
		for (int j = i+1; j < int(jets.size()); j++)
		{
			for (int k = j+1; k < int(jets.size()); k++)
			{

				TLorentzVector tempLor = jets[i] + jets[j] + jets[k];
				/*            double dR1 = deltaR( jets[i].Eta(), jets[i].phi(),
								  tempLor.Eta(), tempLor.phi());
								  double dR2 = deltaR( jets[j].Eta(), jets[j].phi(),
								  tempLor.Eta(), tempLor.phi());
								  double dR3 = deltaR( jets[k].Eta(), jets[k].phi(),
								  tempLor.Eta(), tempLor.phi());
								  */				
				const double dR1 = tempLor.DeltaR(jets[i]);
				const double dR2 = tempLor.DeltaR(jets[j]);
				const double dR3 = tempLor.DeltaR(jets[k]);

				if (dR1 > tripletDRCut_) continue;
				if (dR2 > tripletDRCut_) continue;
				if (dR3 > tripletDRCut_) continue;

				vector<int> tempVec;     
				tempVec.push_back(i);
				tempVec.push_back(j);
				tempVec.push_back(k);
				indices.push_back(tempVec);
			}
		}
	}


	for (unsigned i = 0; i < indices.size(); i++)
	{
		//      cout<<indices[i][0]<<" "<<indices[i][1]<<" "<<indices[i][2]<<" "<<endl;
		bJetPass  = true;
		passDijet = false;
		passM23OverM123 = false;

		vector<int> tempBJetIndices;
		if (indices[i].size() != 3) 
			cout<<indices[i].size()<<endl;

		p41 = jets[indices[i][0]];
		p42 = jets[indices[i][1]];
		p43 = jets[indices[i][2]];

		m123 = (p41 + p42 + p43).M();
		m12 = (p41 + p42).M();
		m13 = (p41 + p43).M();
		m23 = (p42 + p43).M();

		//Check if the bjet(s) is in the right place
		if ( bJets.size() == 0)
		{
			bJetPass = false;
		} else if ( bJets.size() == 1)
		{         
			bJetEta = bJets[0].Eta();
			bJetPhi = bJets[0].Phi();
			//if( (deltaR(p41.Eta(), p41.phi(), bJetEta, bJetPhi) < 0.01)   || 
			//  (deltaR(p42.Eta(), p42.phi(), bJetEta, bJetPhi) < 0.01)   || 
			// (deltaR(p43.Eta(), p43.phi(), bJetEta, bJetPhi) < 0.01) )
			if (  p41.DeltaR(bJets.at(0)) < 0.01
				|| p42.DeltaR(bJets.at(0)) < 0.01
				|| p43.DeltaR(bJets.at(0)) < 0.01 )
			{
				bJetPass = false;
			} else
			{
				tempBJetIndices.push_back(0); //why is this 0?
			}
		} else if( bJets.size() >= 2)
		{
			bJetInTrip = false;
			bJetOutTrip = false;
			for (unsigned j = 0; j < bJets.size(); j++)
			{
				//bJetEta = bJets[j].Eta();
				//bJetPhi = bJets[j].Phi();
				//           if( (deltaR(p41.Eta(), p41.phi(), bJetEta, bJetPhi) < 0.01)  || 
				//               (deltaR(p42.Eta(), p42.phi(), bJetEta, bJetPhi) < 0.01)  || 
				//              (deltaR(p43.Eta(), p43.phi(), bJetEta, bJetPhi) < 0.01) )
				if (  p41.DeltaR(bJets.at(j)) < 0.01 
					|| p42.DeltaR(bJets.at(j)) < 0.01 
					|| p43.DeltaR(bJets.at(j)) < 0.01 )
				{
					bJetInTrip = true;
				} else
				{
					tempBJetIndices.push_back(j);
					bJetOutTrip = true;
				}
			}

			
			if ( !bJetInTrip || !bJetOutTrip)
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

		TLorentzVector tripletP4(0, 0, 0, 0);
		for( int i = 0; i < int(jets.size()); i++)
		{
			if( i == indices[selectedIndex][0] || 
					i == indices[selectedIndex][1] || 
					i == indices[selectedIndex][2])
			{

				tripletP4 = tripletP4 + jets[i];

			}
			else
			{
				//rSystem.push_back(reco::Jet(jets[i].p4(), jets[i].vertex()));
				rSystem.push_back(jets.at(i));
			}
		}

		//so how would I creat a jet that asscoaited with the original vertex as the 
		//3 jets used to create it??
		//NEED TO CHECK THIS! use sum of p4 for now
		//triplet.push_back(reco::Jet(tripletP4, tripletVertex));
		triplet.push_back(tripletP4);
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

}
/* This will mathc the RECO jets to gen jets and assign the CSV of RECO jet
 * to matching GenJet
 */
void FactorizationBySmearing::SetGenJetBdiscriminators(const vector<TLorentzVector>& jets_reco, std::vector<double>& bDisc_reco,
					const std::vector<TLorentzVector>& jets_gen, std::vector<double>& bDisc_gen)
{
	bDisc_gen.clear();

	for (unsigned j =0; j < jets_gen.size(); ++j)
	{
		bool found = false;
		for (unsigned i =0; i < jets_reco.size(); ++i)
		{
			if (jets_reco.at(i).DeltaR(jets_gen.at(j))<0.5) //match found
			{
				found = true;
				bDisc_gen.push_back(bDisc_reco.at(i));
//				Print4vec(jets_reco.at(i), jets_gen.at(j));
				break;
			}
		}

		if (!found) 
		{
			bDisc_gen.push_back(-1); //to keep the vector sizes same
		}
	}
}

void FactorizationBySmearing::Print4vec(const TLorentzVector& tl1, const TLorentzVector& tl2) const 
{
	cout << setw(10) << tl1.Pt() << " / " << tl1.Eta() << " / " << tl1.Phi() << " / " << tl1.E() << endl;  
	cout << setw(10) << tl2.Pt() << " / " << tl2.Eta() << " / " << tl2.Phi() << " / " << tl2.E() << endl;  
}

/* apply stop baseline selection including top tagging */
/*This asumes the lepton veto is already applied at the ntuple making!*/
bool FactorizationBySmearing::PassBaselineSelection(const vector<TLorentzVector>& jets, const TLorentzVector met_vec)
{
		const double max_eta = 2.4;
		unsigned njet_pt70 =0, njet_pt50=0, njet_pt30=0;
		vector<TLorentzVector> jets_eta2p5;
		for (unsigned i=0; i < jets.size(); ++i)
		{
			if (fabs(jets.at(i).Eta()) > max_eta) continue;
			jets_eta2p5.push_back(jets.at(i));
		}

		for (unsigned i=0; i < jets_eta2p5.size(); ++i)
		{
			const float pt = jets_eta2p5.at(i).Pt();

			if (pt>30.0) ++njet_pt30;
			if (pt>50.0) ++njet_pt50;
			if (pt>70.0) ++njet_pt70;
		}

		bool pass_njetcut = false;
		if (njet_pt70>=2 && njet_pt50>=4 && njet_pt30>=5) pass_njetcut = true;
	
		const bool pass_dphicut = PassDphiCut(jets_eta2p5, met_vec,3,0.5,0.5,0.3); 
		const bool pass_metcut  = (met_vec.Pt()>175.0 ) ? true : false;
		//mass
		//1 btag outside triplet
		//MT2>300M
		//MTb+1/MTtop >500

	return (pass_njetcut && pass_dphicut && pass_metcut);
}

void FactorizationBySmearing::MTMT2(const TLorentzVector& MetVec, 
			const std::vector<TLorentzVector> & triplet,
			const std::vector<TLorentzVector> & rSystem, 
			const std::vector<TLorentzVector> & bJetsInR, 
			double & MT2, double & MTt, double & MTb)
{

	//predefined inputs
   //mTop_ = iConfig.getParameter<double>("mTop");
   //mWMin_ = iConfig.getParameter<double>("mWMin");
   //mWMax_ = iConfig.getParameter<double>("mWMax");

	const double mTop_ = 173.5;
	const double mWMin_ = 50.0;
	const double mWMax_ = 120.0;

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

	TLorentzVector bJetP4;
	TLorentzVector wJet1P4;
	TLorentzVector wJet2P4;
	TLorentzVector otherTopP4;
	TLorentzVector tripletP4;

	Mt2::ChengHanBisect_Mt2_332_Calculator mt2Calculator;   


	if( triplet.size() != 0 && rSystem.size() != 0 && bJetsInR.size() != 0)
	{
		tripletP4 = triplet[0];
		for(unsigned i = 0; i < bJetsInR.size(); i++)
		{
			bJetP4 = bJetsInR[i];
			for(unsigned j = 0; j < rSystem.size(); j++)
			{
				wJet1P4 = rSystem[j];
				//if((deltaR(wJet1P4.eta(), wJet1P4.phi(), 
				//           bJetP4.eta(), bJetP4.phi()) < 0.01)) continue;
				if(wJet1P4.DeltaR(bJetP4) < 0.01) continue;

				for(unsigned k = j + 1; k < rSystem.size(); k++)
				{
					wJet2P4 = rSystem[k];

					//if((deltaR(wJet2P4.eta(), wJet2P4.phi(), 
					//           bJetP4.eta(), bJetP4.phi()) < 0.01) )  continue;
					if(wJet2P4.DeltaR(bJetP4) < 0.01) continue;

					if( fabs((bJetP4 + wJet1P4 + wJet2P4).M() - mTop_) < dTopMin)
					{

						dTopMin = fabs((bJetP4 + wJet1P4 + wJet2P4).M() - mTop_);
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
			bJetP4 = bJetsInR[bJetIndex];
			wJet1P4 = rSystem[wJet1Index];
			wJet2P4 = rSystem[wJet2Index];
			if ( (wJet1P4 + wJet2P4).M() < mWMin_ || 
					(wJet1P4 + wJet2P4).M() > mWMax_) 
			{
				outsideWWindow = true;
			}
		}

		if(!foundTriplet)
		{
			dPhiMin = 999999;
			for(unsigned i = 0; i < bJetsInR.size(); i++)
			{
				//dPhi = abs(deltaPhi(bJetsInR[i].phi(), MET.phi()));
				dPhi = fabs(bJetsInR[i].DeltaPhi(MetVec));
				if( dPhi < dPhiMin)
				{
					dPhiMin = dPhi;
					bJetP4 = bJetsInR[i];
					bJetIndexForDoublet = i;
				}
			}
		} else
		{
			bJetIndexForDoublet = bJetIndex;
		}
		dRMin = 99999;
		for(unsigned i = 0; i < rSystem.size(); i++)
		{            
			//dR = deltaR(bJetP4.eta(), bJetP4.phi(), 
			//            rSystem[i].eta(), rSystem[i].phi());
			dR = bJetP4.DeltaR(rSystem[i]);

			if(dR < 0.01) continue;
			if(dR > 2.0) continue;
			if( (bJetP4 + rSystem[i]).M() > mTop_) continue;

			if(dR < dRMin)
			{            
				otherJetIndex = int(i);
				dRMin = dR;
			}
		}

		if(foundTriplet && !outsideWWindow)
		{
			otherTopP4 = bJetP4 + wJet1P4 + wJet2P4;
		} else if(outsideWWindow && otherJetIndex != -1)
		{
			otherTopP4 = rSystem[otherJetIndex] + bJetP4;
		} else if(otherJetIndex != -1)
		{
			otherTopP4 = rSystem[otherJetIndex] + bJetP4;
		} else {
			otherTopP4 = bJetP4;            
		}  


		const double massOfSystemA =  tripletP4.M(); 
		const double pxOfSystemA   =  tripletP4.Px(); 
		const double pyOfSystemA   =  tripletP4.Py(); 

		const double massOfSystemB =  otherTopP4.M(); 
		const double pxOfSystemB   =  otherTopP4.Px(); 
		const double pyOfSystemB   =  otherTopP4.Py(); 

		const double pxMiss        = MetVec.Px(); 
		const double pyMiss        = MetVec.Py(); 

		const double invis_mass    = MetVec.M(); 

		Mt2::LorentzTransverseVector  vis_A(Mt2::TwoVector(pxOfSystemA, 
					pyOfSystemA), massOfSystemA);
		Mt2::LorentzTransverseVector  vis_B(Mt2::TwoVector(pxOfSystemB, pyOfSystemB), massOfSystemB);
		Mt2::TwoVector                pT_Miss(pxMiss, pyMiss);

		MT2 = mt2Calculator.mt2_332(vis_A, vis_B, pT_Miss, invis_mass);
		//      cout<<"MT2: "<<MT2<<endl;
	} else
	{
		MT2 = -1;
	}

	if(triplet.size() != 0)
	{
		tripletP4 = triplet[0];

		Et_1 = sqrt(tripletP4.M()*tripletP4.M() + 
				tripletP4.Pt()*tripletP4.Pt());
		Et_2 = sqrt(MetVec.M()*MetVec.M() + MetVec.Pt()*MetVec.Pt());

		MTt = sqrt( tripletP4.M() * tripletP4.M() + 
				MetVec.M()*MetVec.M() + 
				2*(Et_1*Et_2 - tripletP4.Px() * MetVec.Px() - 
					tripletP4.Py()*MetVec.Py()));
		//      cout<<"MTt: "<<MTt<<endl;
	} else
	{
		MTt = -1;
	}

	if(bJetsInR.size() != 0)
	{
		dPhiMin = 999999;

		for(unsigned i = 0; i < bJetsInR.size(); i++)
		{
			//dPhi = abs(deltaPhi( MetVec.Phi(), bJetsInR[i].phi()));
			dPhi = fabs(MetVec.DeltaPhi(bJetsInR[i]));
			if( dPhi < dPhiMin)
			{
				bJetP4 = bJetsInR[i];
				closestBJetToMETIndex = i;
				dPhiMin = dPhi;
			}         
		} 

		otherTopP4 = bJetP4;

		if( otherJetIndex  != -1 && closestBJetToMETIndex == bJetIndexForDoublet)
		{
			otherTopP4 = rSystem[otherJetIndex] + bJetP4;            
		}

		Et_1 = sqrt(otherTopP4.M()*otherTopP4.M() + otherTopP4.Pt() * otherTopP4.Pt());
		Et_2 = sqrt(MetVec.M()*MetVec.M() + MetVec.Pt()*MetVec.Pt());

		MTb = sqrt( otherTopP4.M() * otherTopP4.M() +
				MetVec.M() * MetVec.M() + 
				2*(Et_1*Et_2 - otherTopP4.Px() * MetVec.Px() - 
					otherTopP4.Py() * MetVec.Py()));
		//      cout<<"MTb: "<<MTb<<endl;
	} else {
		MTb = -1;
	}


}

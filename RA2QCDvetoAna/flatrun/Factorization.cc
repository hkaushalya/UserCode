//
// Original Author:  samantha hewamanage
// $Id$
//
//


// system include files
#include <memory>
#include <vector>
#include <algorithm>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

//ROOT
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "assert.h"
#include "TVector3.h"
#include <iomanip>
#include "TPad.h"
#include "TBenchmark.h"
#include "TChain.h"
#include "TFile.h"
#include "TROOT.h"
#include "Factorization.hh"
using namespace std;

#define M_PI   3.1415926535897932385
static bool sort_using_less_than(double u, double v)
{
	   return u < v;
}

double deltaPhi(double phi1, double phi2)
{ 
 double result = phi1 - phi2;
 while (result > M_PI) result -= 2*M_PI;
 while (result <= -M_PI) result += 2*M_PI;
 return result;
}



//
// constructors and destructor
//
Factorization::Factorization():
sHistFileName("Default.root"),
iProgressBy(100000),
iVerbose(0),
dMinHT(0),
dMinMHT(0)
{
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
int Factorization::Run(TChain* myChain, int iRunEvents)
{

	std::cout << __FILE__ << ": Events Requested = " << iRunEvents << std::endl;
	TBenchmark timer;
	timer.Start("phoana_time");
	
	if (myChain == NULL) {
		std::cout << "NULL chain!" << std::endl;
		assert(false);
	}
	if (myChain->GetEntries()<1)
	{
		std::cout << "No entries found!" << std::endl;
		return 1;
	}

	gROOT->Reset();
	
	BookHistograms();
	
	Init(myChain);
	unsigned int iENTRIES = (unsigned int) myChain->GetEntries();
	unsigned int iEvt2Process = 0; 	//_____________________________________________ events to be processed
	if (iRunEvents < 0 ) {
		iRunEvents = iENTRIES;			//______________________________________________ requested all entries
		iEvt2Process = iENTRIES;
	} else if (iRunEvents > 0 && iRunEvents <= iENTRIES) iEvt2Process = iRunEvents;
	else iEvt2Process = iENTRIES;		//______________________________________________ default

	std::cout << " ==> Entries found :: " << iENTRIES << std::endl;

	//do a check for the MC. this is a common mistake I make
	
	  myChain->GetEntry(1);
/*	  if (stuple.evt_McFlag != GetMcFlag()) 
	  {
		  std::cout << red << "MC Flag," << GetMcFlag()
			  << ", setting does not match what is in Stuple, " 
			  << stuple.evt_McFlag << ". pleas check." 
			  << " returning!" << clearatt << std::endl;
		  assert (false);
	  }
*/
	unsigned int iEvtProc = 0;			//______________________________________________ number of events processed
	std::cout << "Init done. Ready, Set , GO... " << std::endl; 
	
	timer.Start("looptimer");
	
	int iCurrTree = -1;
	for ( ; iEvtProc < iEvt2Process; ++iEvtProc)
	{
		if (iCurrTree != myChain->GetTreeNumber())
		{
			std::cout << "Opening file " << myChain->GetFile()->GetName() << std::endl;
			iCurrTree = myChain->GetTreeNumber();
		}
		myChain->GetEntry(iEvtProc);

		if (iEvtProc > 0 && (iEvtProc%iProgressBy == 0)) {
			std::cout << iEvtProc << "\t events processed [" << (int)( iEvtProc/(double)iEvt2Process * 100)<< "%] ";
			timer.Show("looptimer");
			timer.Start("looptimer");
		}


		++uProcessed;
		kRunLumiEvent.run  = run;
		kRunLumiEvent.lumi = event; 
		kRunLumiEvent.evt  = ls;

		if (iVerbose) 
		{
			std::cout << __LINE__<< ":: Processing event: "
			<<  kRunLumiEvent.run << ":" << kRunLumiEvent.lumi 
				<< ":" << kRunLumiEvent.evt << std::endl;
			std::cout << "jet_num = " <<  jet_num  << "/" << ht_TR << std::endl;
		}

		Weight = 1;

		/* PU S3 reweighting for QCD MC sample
		*/
		double lumiWeight = 1;
		if ( doLumiWeighing )
		{
			lumiWeight = PUWeight;
			//std::cout << "lum wgt = " << lumiWeight << std::endl;
			sumLumiWeights += lumiWeight;
			Weight *= lumiWeight;
		}

		//event weights for flat QCD samples
		if ( doEventWeighing )
		{
			Weight *= storedWeight;
		}

		//APPLY RA2 cuts

		if ( jet_num <3 ) continue;
		if ( ht_TR < dMinHT) { ++uFailMinHTCut; continue;}
		if ( mht_TR < dMinMHT ) {++uFailMinPFMHTCut; continue; }

		/*
		 * Fill hists after RA2b base selection
		 */
		evtHist.njet->Fill(jet_num, Weight);
		evtHist.mht->Fill( mht_TR, Weight);
		evtHist.ht->Fill(ht_TR, Weight);
		evtHist.meff->Fill(meff, Weight);

		pf30_jet1Hist.pt->Fill(jet1pt_TR, Weight);
		pf30_jet2Hist.pt->Fill(jet2pt_TR, Weight);
		pf30_jet3Hist.pt->Fill(jet3pt_TR, Weight);
		pf30_jet1Hist.eta->Fill(jet1eta_TR, Weight);
		pf30_jet2Hist.eta->Fill(jet2eta_TR, Weight);
		pf30_jet3Hist.eta->Fill(jet3eta_TR, Weight);

		DoDelMinStudy();
	}


	histFile->Write();
	histFile->Close();
	endJob();


	++uPassed;
   return true;
}

// ------------ method called once each job just before starting event loop  ------------
/*void 
Factorization::PUWeighing()
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
*/
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
	std::cout << "[ATM:50] LumiWeights Avg ---- = " << sumLumiWeights/(double)uPassed << std::endl;

}

void Factorization::DoDelMinStudy(
				)
{
	//PrintHeader();
	std::vector<float> vDelPhi_jetmht;
	const float mht = mht_TR;

	for (unsigned i = 0 ; i < jet_num ; ++i)
	{
		TLorentzVector vec(alljets_Px[i], alljets_Py[i], alljets_Pz[i], alljets_E[i]);

		const float delphi_jetmht = fabs(TVector2::Phi_mpi_pi(vec.Phi() - mhtphi_TR));
		vDelPhi_jetmht.push_back(delphi_jetmht);
		if (i>2) break; //use only three leading jets
		//std::cout << "i = " <<  i << std::endl;
	}

	std::sort(vDelPhi_jetmht.begin(), vDelPhi_jetmht.end(), sort_using_less_than);	

	const float dPhiMin = vDelPhi_jetmht.at(0);

	dphiVsMHT_2D->Fill(mht, dPhiMin);

	// dphimin in slices of MHT
	if ( mht_TR >= 60 && mht_TR< 80 ) hDelPhiMin_mht[0]->Fill(dPhiMin, Weight);
	else if ( mht_TR >= 80 && mht_TR< 100 ) hDelPhiMin_mht[1]->Fill(dPhiMin, Weight);
	else if ( mht_TR >= 100 && mht_TR< 120 ) hDelPhiMin_mht[2]->Fill(dPhiMin, Weight);
	else if ( mht_TR >= 120 && mht_TR< 140 ) hDelPhiMin_mht[3]->Fill(dPhiMin, Weight);
	else if ( mht_TR >= 140 && mht_TR< 170 ) hDelPhiMin_mht[4]->Fill(dPhiMin, Weight);
	else if ( mht_TR >= 170 && mht_TR< 200 ) hDelPhiMin_mht[5]->Fill(dPhiMin, Weight);
	else if ( mht_TR >= 200 && mht_TR< 250 ) hDelPhiMin_mht[6]->Fill(dPhiMin, Weight);
	else if ( mht_TR >= 250 ) hDelPhiMin_mht[7]->Fill(dPhiMin, Weight);

	// mht in slices of dphimin
	if ( dPhiMin < 0.1 ) MHT_by_phislice[0]->Fill(mht_TR);
	else if ( dPhiMin >= 0.1 && dPhiMin < 0.2 ) MHT_by_phislice[1]->Fill(mht_TR, Weight);
	else if ( dPhiMin >= 0.2 && dPhiMin < 0.3 ) MHT_by_phislice[2]->Fill(mht_TR, Weight);
	else if ( dPhiMin >= 0.3 && dPhiMin < 0.5 ) MHT_by_phislice[3]->Fill(mht_TR, Weight);
	else if ( dPhiMin >= 0.5 && dPhiMin < 0.8 ) MHT_by_phislice[4]->Fill(mht_TR, Weight);
	else if ( dPhiMin > 0.8 ) MHT_by_phislice[5]->Fill(mht_TR, Weight);


	//make PASS/FAIL plots with RA2 cut

	if (dPhiMin>0.15) hPass[0]->Fill(mht_TR, Weight);
	else 	hFail[0]->Fill(mht_TR, Weight);

	if (dPhiMin>0.2)
	{
		hPass[1]->Fill(mht_TR, Weight);
		hFail[6]->Fill(mht_TR, Weight);
	} else 	hFail[1]->Fill(mht_TR, Weight);

	if (dPhiMin>0.25) hPass[2]->Fill(mht_TR, Weight);
	else 	hFail[2]->Fill(mht_TR, Weight);

	if (dPhiMin>=0.3)
	{
		hPass[3]->Fill(mht_TR, Weight);
		hFail[7]->Fill(mht_TR, Weight);
	} else 	hFail[3]->Fill(mht_TR, Weight);

	if (dPhiMin>=0.35) hPass[4]->Fill(mht_TR, Weight);
	else 	hFail[4]->Fill(mht_TR, Weight);

	if (dPhiMin>=0.4) hPass[5]->Fill(mht_TR, Weight);
	else 	hFail[5]->Fill(mht_TR, Weight);

	//Pass selection with RA2 dphi cuts
	bool passed = true;
	if (jet_num >= 1) passed = passed && (fabs(deltaPhi(jet1phi_TR, mhtphi_TR)) > 0.5);
	if (jet_num >= 2) passed = passed && (fabs(deltaPhi(jet2phi_TR, mhtphi_TR)) > 0.5);
	if (jet_num >= 3) passed = passed && (fabs(deltaPhi(jet3phi_TR, mhtphi_TR)) > 0.3);

	if (passed) 
	{
		hPass[6]->Fill(mht_TR, Weight);
		if ( ht_TR > 500)
		{
			hPass[7]->Fill(mht_TR, Weight);
			if (mht_TR>200) hSignalRegion[0]->Fill(mht_TR,Weight);
			if (mht_TR>350) hSignalRegion[1]->Fill(mht_TR,Weight);
			if (mht_TR>500) hSignalRegion[2]->Fill(mht_TR,Weight);

		}
		if ( ht_TR > 800) 
		{
			hPass[8]->Fill(mht_TR, Weight);
		}

	}


	if (dPhiMin<0.2)
	{
		if ( ht_TR > 500) hFail[8]->Fill(mht_TR, Weight);
		if ( ht_TR > 800) hFail[9]->Fill(mht_TR, Weight);
		if ( ht_TR > 1000) hFail[10]->Fill(mht_TR, Weight);
		if ( ht_TR > 1200) hFail[11]->Fill(mht_TR, Weight);
	}

	//for systematics
	passed = true;
	if (jet_num >= 1) passed = passed && (fabs(deltaPhi(jet1phi_TR, mhtphi_TR)) > 0.5);
	if (jet_num >= 2) passed = passed && (fabs(deltaPhi(jet1phi_TR, mhtphi_TR)) > 0.5);
	if (jet_num >= 3) passed = passed && (fabs(deltaPhi(jet1phi_TR, mhtphi_TR)) > 0.1);
	
	if (! passed)
	{
		if ( ht_TR > 500) 
		{
			hFail[12]->Fill(mht_TR, Weight);
			hFail[13]->Fill(mht_TR, Weight);
		}
		if ( ht_TR > 800)
		{
			hFail[14]->Fill(mht_TR, Weight);
			hFail[15]->Fill(mht_TR, Weight);
		}
	}


	passed = true;
	if (jet_num >= 1) passed = passed && (fabs(deltaPhi(jet1phi_TR, mhtphi_TR)) > 0.4);
	if (jet_num >= 2) passed = passed && (fabs(deltaPhi(jet1phi_TR, mhtphi_TR)) > 0.4);
	if (jet_num >= 3) passed = passed && (fabs(deltaPhi(jet1phi_TR, mhtphi_TR)) > 0.2);
	
	if (! passed)
	{
		if ( ht_TR > 500) 
		{
			hFail[16]->Fill(mht_TR, Weight);
			hFail[17]->Fill(mht_TR, Weight);
		}
		if ( ht_TR > 800)
		{
			hFail[18]->Fill(mht_TR, Weight);
			hFail[19]->Fill(mht_TR, Weight);
		}
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

void Factorization::Init(TChain* outTree)
{

		outTree->SetBranchStatus("*", 1);
      outTree->SetBranchAddress("run", &run);
      outTree->SetBranchAddress("event", &event);
      outTree->SetBranchAddress("lumi", &ls);
      outTree->SetBranchAddress("npv", &npv);
      outTree->SetBranchAddress("vtxSize", &vtxSize);
      outTree->SetBranchAddress("nJets", &nJets);
      outTree->SetBranchAddress("evtWeight", &evtWeight_TR);
      outTree->SetBranchAddress("meff", &meff);
      outTree->SetBranchAddress("mht", &mht_TR);
      outTree->SetBranchAddress("ht", &ht_TR);
      outTree->SetBranchAddress("met", &met_TR);
      outTree->SetBranchAddress("mt", &mt_TR);
      outTree->SetBranchAddress("mhtphi", &mhtphi_TR);
      outTree->SetBranchAddress("metphi", &metphi_TR);
      outTree->SetBranchAddress("metSgnf", &metSgnf_TR);
      outTree->SetBranchAddress("metSgnfProb", &metSgnfProb_TR);
      outTree->SetBranchAddress("jet1pt", &jet1pt_TR);
      outTree->SetBranchAddress("jet1eta", &jet1eta_TR);
      outTree->SetBranchAddress("jet1phi", &jet1phi_TR);
      outTree->SetBranchAddress("jet2pt", &jet2pt_TR);
      outTree->SetBranchAddress("jet2eta", &jet2eta_TR);
      outTree->SetBranchAddress("jet2phi", &jet2phi_TR);
      outTree->SetBranchAddress("jet3pt", &jet3pt_TR);
      outTree->SetBranchAddress("jet3eta", &jet3eta_TR);
      outTree->SetBranchAddress("jet3phi", &jet3phi_TR);
      //outTree->SetBranchAddress("nJets_CUT", &nJets_CUT);
/*      outTree->SetBranchAddress("jet1Res", &jet1Res_TR);
      outTree->SetBranchAddress("jet2Res", &jet2Res_TR);
      outTree->SetBranchAddress("jet3Res", &jet3Res_TR);*/
		outTree->SetBranchAddress("jet_num", &jet_num);
		outTree->SetBranchAddress("alljets_E", alljets_E);
		outTree->SetBranchAddress("alljets_Px", alljets_Px);
		outTree->SetBranchAddress("alljets_Py", alljets_Py);
		outTree->SetBranchAddress("alljets_Pz", alljets_Pz);
		outTree->SetBranchAddress("alljets_Pt", alljets_Pt);
		outTree->SetBranchAddress("alljets_Eta", alljets_Eta);
		outTree->SetBranchAddress("mcFlag", &mcFlag);
		outTree->SetBranchAddress("PUWeight", &PUWeight);
		outTree->SetBranchAddress("storedWeight", &storedWeight);


	uProcessed = 0;
	uPassed = 0;

} 



void Factorization::BookHistograms()
{
	//generate hists

	histFile = new TFile(sHistFileName.c_str(),"RECREATE");

   dphiVsMHT_2D = new TH2F ("dphiMinVsMHT" ,";#slash{H}_{T} [GeV];#Delta#Phi_{min}(jet_{1,2,3},#slash{H}_{T});", 800,0,800,300,0,3);


	const float met_min = 0, met_max=500, met_bins=100;
	MHT_by_phislice[0] = new TH1F ("mht_phislice_lt0.1" ,"MHT (#Delta#Phi_{min}<0.1);MHT [GeV];Events;", met_bins, met_min, met_max);
	MHT_by_phislice[1] = new TH1F ("mht_phislice_lt0.2" ,"MHT (0.1<#Delta#Phi_{min}<0.2);MHT [GeV];Events;", met_bins, met_min, met_max);
	MHT_by_phislice[2] = new TH1F ("mht_phislice_lt0.3" ,"MHT (0.2<#Delta#Phi_{min}<0.3);MHT [GeV];Events;", met_bins, met_min, met_max);
	MHT_by_phislice[3] = new TH1F ("mht_phislice_lt0.5" ,"MHT (0.3<#Delta#Phi_{min}<0.5);MHT [GeV];Events;", met_bins, met_min, met_max);
	MHT_by_phislice[4] = new TH1F ("mht_phislice_lt0.8" ,"MHT (0.5<#Delta#Phi_{min}<0.8);MHT [GeV];Events;", met_bins, met_min, met_max);
	MHT_by_phislice[5] = new TH1F ("mht_phislice_gt0.8" ,"MHT (#Delta#Phi_{min}>0.8);MHT [GeV];Events;", met_bins, met_min, met_max);

	//these are general event hist to check the PAT tuples cuts
	const double evt_met_max = 800, evt_met_bins = 400;
	const double evt_ht_max = 4000, evt_ht_bins = 80;

	evtHist.mht = new TH1F ("mht" ,"RA2: (MHT from PFmetHandle);MHT [GeV];Events;", evt_met_bins, 0, evt_met_max);
	evtHist.ht = new TH1F ("ht" ,"RA2: HT from Jets ET>50 && |#Eta|<2.4;HT [GeV];Events;", evt_ht_bins, 0, evt_ht_max);
	evtHist.njet = new TH1F ("njet_et50eta24" ,"RA2: Njets (Et>50 && |#Eta|<2.4;NJETS;Events;", 10, 0, 10);
	evtHist.meff = new TH1F ("meteff" ,"RA2:;MEff;Events;", 50, 0, 5000);

	const double pt_bins = 200, pt_max = 2000;
	pf30_jet1Hist.pt = new TH1F ("pf30_jet1_pt" ,"RA2: PF30-Jet1 pt;pt [GeV];Events;", pt_bins, 0, pt_max);
	pf30_jet1Hist.eta = new TH1F ("pf30_jet1_eta" ,"RA2: PF30-Jet1 eta;eta;Events;", 100, -5, 5);
	pf30_jet1Hist.phi = new TH1F ("pf30_jet1_phi" ,"RA2: PF30-Jet1 phi;phi;Events;", 160, -8, 8);
//	pf30_jet1Hist.delphi = new TH1F ("pf30_jet1_delphi" ,"RA2: PF30-Jet1: delphi;delphi;Events;", 160, -8, 8);

	pf30_jet2Hist.pt = new TH1F ("pf30_jet2_pt" ,"RA2: PF30-Jet2 pt;pt [GeV];Events;", pt_bins, 0, pt_max);
	pf30_jet2Hist.eta = new TH1F ("pf30_jet2_eta" ,"RA2: PF30-Jet2 eta;eta;Events;", 100, -5, 5);
	pf30_jet2Hist.phi = new TH1F ("pf30_jet2_phi" ,"RA2: PF30-Jet2 phi;phi;Events;", 160, -8, 8);
//	pf30_jet2Hist.delphi = new TH1F ("pf30_jet2_delphi" ,"RA2: PF30-Jet2: delphi;delphi;Events;", 160, -8, 8);

	pf30_jet3Hist.pt = new TH1F ("pf30_jet3_pt" ,"RA2: PF30-Jet3 pt;pt [GeV];Events;", pt_bins, 0, pt_max);
	pf30_jet3Hist.eta = new TH1F ("pf30_jet3_eta" ,"RA2: PF30-Jet3 eta;eta;Events;", 100, -5, 5);
	pf30_jet3Hist.phi = new TH1F ("pf30_jet3_phi" ,"RA2: PF30-Jet3 phi;phi;Events;", 160, -8, 8);
//	pf30_jet3Hist.delphi = new TH1F ("pf30_jet3_delphi" ,"RA2: PF30-Jet3: delphi;delphi;Events;", 160, -8, 8);

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

	hPass[0] = new TH1F ("Pass_0","PASS from #Delta#Phi_{min} cut >0.15", npassFailHistBins, passFailHistBins);
	hPass[1] = new TH1F ("Pass_1","PASS from #Delta#Phi_{min} cui >0.2", npassFailHistBins, passFailHistBins);
	hPass[2] = new TH1F ("Pass_2","PASS from #Delta#Phi_{min} cut >0.25", npassFailHistBins, passFailHistBins);
	hPass[3] = new TH1F ("Pass_3","PASS from #Delta#Phi_{min} cut >0.3", npassFailHistBins, passFailHistBins);
	hPass[4] = new TH1F ("Pass_4","PASS from #Delta#Phi_{min} cut >0.35", npassFailHistBins, passFailHistBins);
	hPass[5] = new TH1F ("Pass_5","PASS from #Delta#Phi_{min} cut >0.4", npassFailHistBins, passFailHistBins);
	hPass[6] = new TH1F ("Pass_RA2dphi","PASS from RA2 dPhi Selection: dphi (j1 & j1 , mht) >0.5 and dphi (j3, mht)>0.3)", npassFailHistBins, passFailHistBins);
	hPass[7] = new TH1F ("Pass_RA2dphi_HT500","PASS from RA2 dPhi Selection: dphi (j1 & j1 , mht) >0.5 and dphi (j3, mht)>0.3) && HT>500 GeV", npassFailHistBins, passFailHistBins);
	hPass[8] = new TH1F ("Pass_RA2dphi_HT800","PASS from RA2 dPhi Selection: dphi (j1 & j1 , mht) >0.5 and dphi (j3, mht)>0.3) && HT>800 GeV", npassFailHistBins, passFailHistBins);

	hFail[0] = new TH1F ("Fail_0","FAIL from #Delta#Phi_{min} cut <0.15", npassFailHistBins, passFailHistBins);
	hFail[1] = new TH1F ("Fail_1","FAIL from #Delta#Phi_{min} cut <0.2", npassFailHistBins, passFailHistBins);
	hFail[2] = new TH1F ("Fail_2","FAIL from #Delta#Phi_{min} cut <0.25", npassFailHistBins, passFailHistBins);
	hFail[3] = new TH1F ("Fail_3","FAIL from #Delta#Phi_{min} cut <0.3", npassFailHistBins, passFailHistBins);
	hFail[4] = new TH1F ("Fail_4","FAIL from #Delta#Phi_{min} cut <0.35", npassFailHistBins, passFailHistBins);
	hFail[5] = new TH1F ("Fail_5","FAIL from #Delta#Phi_{min} cut <0.4", npassFailHistBins, passFailHistBins);
	hFail[6] = new TH1F ("Fail_lt_point2","Fail from #Delta#Phi_{min} cut <0.2", 1500, 0, 1500);
	hFail[7] = new TH1F ("Fail_lt_point3","Fail from #Delta#Phi_{min} cut <0.3", 1500, 0, 1500);
	hFail[8] = new TH1F ("Fail_lt_point2_HT500","Fail from #Delta#Phi_{min} cut <0.2 && HT>500 GeV", 1500, 0, 1500);
	hFail[9] = new TH1F ("Fail_lt_point2_HT800","Fail from #Delta#Phi_{min} cut <0.2 && HT>800 GeV", 1500, 0, 1500);
	hFail[10] = new TH1F ("Fail_lt_point2_HT1000","Fail from #Delta#Phi_{min} cut <0.2 && HT>1000 GeV", 1500, 0, 1500);
	hFail[11] = new TH1F ("Fail_lt_point2_HT1200","Fail from #Delta#Phi_{min} cut <0.2 && HT>1200 GeV", 1500, 0, 1500);

	hFail[12] = new TH1F ("Syst1_Fail_ht500","HT>500: FAIL from dphi (j1 & j1 , mht) >0.5 and dphi (j3, mht)>0.1 ", npassFailHistBins, passFailHistBins);
	hFail[13] = new TH1F ("Syst1_Fail_ht500_fineBin","FineBin:HT>500: FAIL from dphi (j1 & j1 , mht) >0.5 and dphi (j3, mht)>0.1 ", 1500, 0, 1500);
	hFail[14] = new TH1F ("Syst1_Fail_ht800","HT>800: FAIL from dphi (j1 & j1 , mht) >0.5 and dphi (j3, mht)>0.1 ", npassFailHistBins, passFailHistBins);
	hFail[15] = new TH1F ("Syst1_Fail_ht800_fineBin","FineBin: HT>800: FAIL from dphi (j1 & j1 , mht) >0.5 and dphi (j3, mht)>0.1 ", 1500, 0, 1500);
	hFail[16] = new TH1F ("Syst2_Fail_ht500","HT>500: FAIL from dphi (j1 & j1 , mht) >0.4 and dphi (j3, mht)>0.2", npassFailHistBins, passFailHistBins);
	hFail[17] = new TH1F ("Syst2_Fail_ht500_fineBin","FineBin: HT>500: FAIL from dphi (j1 & j1 , mht) >0.4 and dphi (j3, mht)>0.2 ", 1500, 0, 1500);
	hFail[18] = new TH1F ("Syst2_Fail_ht800","HT>800: FAIL from dphi (j1 & j1 , mht) >0.4 and dphi (j3, mht)>0.2", npassFailHistBins, passFailHistBins);
	hFail[19] = new TH1F ("Syst2_Fail_ht800_fineBin","FineBin: HT>800: FAIL from dphi (j1 & j1 , mht) >0.4 and dphi (j3, mht)>0.2 ", 1500, 0, 1500);

	for (int i = 0; i <=8; ++i ) hPass[i]->Sumw2();
	for (int i = 0; i <=20; ++i ) hFail[i]->Sumw2();

	hSignalRegion[0] = new TH1F ("Signal_0","Signal Region: MHT>200 GeV && pass RA2 #Delta#Phi_{min}", npassFailHistBins, passFailHistBins);
	hSignalRegion[1] = new TH1F ("Signal_1","Signal Region: MHT>350 GeV && pass RA2 #Delta#Phi_{min}", npassFailHistBins, passFailHistBins);
	hSignalRegion[2] = new TH1F ("Signal_2","Signal Region: MHT>500 GeV && pass RA2 #Delta#Phi_{min}", npassFailHistBins, passFailHistBins);

	hSignalRegion[0]->Sumw2();
	hSignalRegion[1]->Sumw2();
	hSignalRegion[2]->Sumw2();

	hDelPhiMin_mht[0] = new TH1F ("ra2_delPhiMin_MHT_60TO80","RA2: #Delta#Phi_{min} for 60 <#slash{H}_{T}<80 GeV", 150, 0, 3);
	hDelPhiMin_mht[1] = new TH1F ("ra2_delPhiMin_MHT_80to100","RA2: #Delta#Phi_{min} for 80 <#slash{H}_{T}<100 GeV", 150, 0, 3);
	hDelPhiMin_mht[2] = new TH1F ("ra2_delPhiMin_MHT_100to120","RA2: #Delta#Phi_{min} for 100<#slash{H}_{T}<120 GeV", 150, 0, 3);
	hDelPhiMin_mht[3] = new TH1F ("ra2_delPhiMin_MHT_120to140","RA2: #Delta#Phi_{min} for 120<#slash{H}_{T}<140 GeV", 150, 0, 3);
	hDelPhiMin_mht[4] = new TH1F ("ra2_delPhiMin_MHT_140to170","RA2: #Delta#Phi_{min} for 140<#slash{H}_{T}<170 GeV", 150, 0, 3);
	hDelPhiMin_mht[5] = new TH1F ("ra2_delPhiMin_MHT_170to200","RA2: #Delta#Phi_{min} for 170<#slash{H}_{T}<200 GeV", 150, 0, 3);
	hDelPhiMin_mht[6] = new TH1F ("ra2_delPhiMin_MHT_200to250","RA2: #Delta#Phi_{min} for 200<#slash{H}_{T}<250 GeV", 150, 0, 3);
	hDelPhiMin_mht[7] = new TH1F ("ra2_delPhiMin_MHT_250up","RA2: #Delta#Phi_{min} for #slash{H}_{T}>250 GeV", 150, 0, 3);

	hDelPhiMin_mht[1]->Sumw2();
	hDelPhiMin_mht[2]->Sumw2();
	hDelPhiMin_mht[3]->Sumw2();
	hDelPhiMin_mht[4]->Sumw2();
	hDelPhiMin_mht[5]->Sumw2();
	hDelPhiMin_mht[6]->Sumw2();
	hDelPhiMin_mht[7]->Sumw2();


}

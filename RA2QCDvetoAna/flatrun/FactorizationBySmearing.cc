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

using namespace std;


/***************************************************************
 * This class is created for Factorization to utilize gen-jet
 * smearing to better quantify constant 'c'. This works only 
 * flat ntuples created using a modifed version of
 * lostLeptonTree.
 * Author: Sam Hewamanage
 * Institution: Florida Internationa University, USA.
 **************************************************************/



static bool sort_using_less_than(double u, double v)
{
	   return u < v;
}

int main(int argc, char* argv[])
{

	if (argc <= 3) {
		cerr <<"Please give 5 arguments " << "runList " << " " << "outputFileName" << " " << "# of evts  nJet50Min  nJet50Max" << endl;
		cerr <<" Valid configurations are " << std::endl;
		cerr << " ./optimize runlist_ttjets.txt isoplots.root qcd.files 100 2 5" << std::endl;
		return -1;
	}

	const char *inputFileList = argv[1];
	const char *outFileName   = argv[2];
	const char *arg3          = argv[3];
	//const char *arg4          = argv[4];
	//const char *arg5          = argv[5];
	int evts    = -1;

	if (isdigit(arg3[0]))
	{
		int num = atoi(argv[3]);
		//evts = 1000;
		if (num>=1) evts = num;
		cout << __FUNCTION__ << ": evts = " << evts << endl;
	} else 
	{
		cout << "argument 3 is not a number. using default value for evts = " << evts << endl;
	}
	

	FactorizationBySmearing smear(inputFileList, outFileName);
	//cout << "dataset " << dataset << " " << endl;

	smear.EventLoop(inputFileList, evts);

	return 0;
}

void FactorizationBySmearing::EventLoop(const char *datasetname, const int evts2Process)
{
	if (fChain == 0) return;

	//settings for the day
/*	PtBinEdges_scaling_  = NumStringToVec("0.,5000.");
	EtaBinEdges_scaling_ = NumStringToVec("0.,5.");
	LowerTailScaling_    = NumStringToVec("1.0");
	UpperTailScaling_    = NumStringToVec("1.0");
	AdditionalSmearing_  = NumStringToVec("1.0");
*/
	/*****************************************************
	 * Smearing function constants and variables
	 ******************************************************/
	absoluteTailScaling_ = true;  //false for systematics
	smearedJetPt_        = 13.0;
	MHTcut_low_          = 200.;
	MHTcut_medium_       = 350.;
	MHTcut_high_         = 500.;
	HTcut_low_           = 500.;
	HTcut_medium_        = 800.;
	HTcut_high_          = 1000.;
	HTcut_veryhigh_      = 1200.;
	HTcut_extremehigh_   = 1400.;

	PtBinEdges_.push_back(0);
	PtBinEdges_.push_back(20);
	PtBinEdges_.push_back(30);
	PtBinEdges_.push_back(50);
	PtBinEdges_.push_back(80);
	PtBinEdges_.push_back(120);
	PtBinEdges_.push_back(170);
	PtBinEdges_.push_back(230);
	PtBinEdges_.push_back(300);
	PtBinEdges_.push_back(380);
	PtBinEdges_.push_back(470);
	PtBinEdges_.push_back(570);
	PtBinEdges_.push_back(680);
	PtBinEdges_.push_back(800);
	PtBinEdges_.push_back(1000);
	PtBinEdges_.push_back(1300);
	PtBinEdges_.push_back(1700);
	PtBinEdges_.push_back(2200);
	PtBinEdges_.push_back(2800);
	PtBinEdges_.push_back(3500);

	EtaBinEdges_.push_back(0); 
	EtaBinEdges_.push_back(0.3); 
	EtaBinEdges_.push_back(0.5); 
	EtaBinEdges_.push_back(0.8); 
	EtaBinEdges_.push_back(1.1); 
	EtaBinEdges_.push_back(1.4); 
	EtaBinEdges_.push_back(1.7); 
	EtaBinEdges_.push_back(2.0); 
	EtaBinEdges_.push_back(2.3); 
	EtaBinEdges_.push_back(2.8); 
	EtaBinEdges_.push_back(3.2); 
	EtaBinEdges_.push_back(4.1); 
	EtaBinEdges_.push_back(5.0);


	//difference variation of the dPhiMin selections.
	//make sure the book the correct number of histograms 
	//when these are changed!!!
	vDphiVariations.push_back(0.15);
	vDphiVariations.push_back(0.20);
	vDphiVariations.push_back(0.25);
	vDphiVariations.push_back(0.30);
	vDphiVariations.push_back(0.35);
	vDphiVariations.push_back(0.40);

	/**** events ti process ***********************/
	Long64_t nentries = fChain->GetEntriesFast();
	cout << "nentries = " << nentries << endl;
	if (evts2Process>=1 && evts2Process < nentries) 
	{	
		nentries = evts2Process;
		cout << "Requested events to process = " << nentries << endl;
	}

	bDEBUG = false;
	const double dDATA_LUMI = 10000.0; // 10 fb-1
	const double lumiWgt = GetLumiWgt(datasetname, dDATA_LUMI); 
	//const double lumiWgt = 1; 
	//jet bins: 2, 3-5, 6-7,>=8
	NJet50_min_          = 8;
	NJet50_max_          = 1000;
	unsigned nTries_     = 1000; //number of pseudo experiments per event
	if (bDEBUG) nTries_  = 1;
	const double smearingWgt = 1.0/(double)nTries_;

	std::cout << red << "Dataset = " << datasetname << endl;
	std::cout << "Using lumiWgt / smearingWgt = " << lumiWgt << " / " 
				<< smearingWgt << clearatt << std::endl;

	Long64_t nbytes = 0, nb = 0;
	int decade = 0;
	
   smearFunc_ = new SmearFunction();
		smearFunc_->SetAbsoluteTailScaling(absoluteTailScaling_);

	vector<TLorentzVector> recoJets, genJets, smearedGenJets;
	BookJerDebugHists();

	unsigned nProcessed = 0;


	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		++nProcessed;
		// ==============print number of events done == == == == == == == =
		double progress = 10.0 * jentry / (1.0 * nentries);
		int k = int (progress);
		if (k > decade)
		{
			cout << green <<  10 * k << " % (" << jentry << ")" << clearatt << endl;
		}
		decade = k;

		// ===============read this entry == == == == == == == == == == == 
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;

		if (bDEBUG) cout << red << "run/ls/evt= " << t_EvtRun << " / " 
				<<  t_EvtLS << " / " <<  t_EvtEvent << clearatt  << endl;
		//cout << red << "run/ls/evt= " << t_EvtRun << " / " 
		//		<<  t_EvtLS << " / " <<  t_EvtEvent << clearatt  << endl;

		CreateRecoJetVec(recoJets);
		CreateGenJetVec(genJets);
		const int njet50_reco = CountJets(recoJets, 50, 2.5);
		const int njet50_gen = CountJets(genJets, 50, 2.5);
		if (njet50_reco>= NJet50_min_ && njet50_reco<= NJet50_max_) FillHistogram(recoJets,0, lumiWgt);
		if (njet50_gen >= NJet50_min_ && njet50_gen <= NJet50_max_) FillHistogram(genJets,1, lumiWgt);
		const double recoHT = HT(recoJets);
		const double genHT  = HT(genJets);
		const double recoMHT= (MHT(recoJets)).Pt();
		const double genMHT = (MHT(genJets)).Pt();

		for (unsigned n = 0; n < nTries_; ++n)
		{
			SmearingGenJets(genJets, smearedGenJets);
			const double smHT  = HT(smearedGenJets);
			const double smMHT = (MHT(smearedGenJets)).Pt();

			if (bDEBUG) 
			{
				cout << "entry/recoht/genht/smht=" << jentry 
					<< "/" << recoHT << "/" << genHT << "/" << smHT << endl;
				cout << "entry/recomht/genmht/smmht=" << jentry 
					<< "/" << recoMHT << "/" << genMHT << "/" << smMHT << endl;
				DumpJets(genJets);
				DumpJets(smearedGenJets);
			}

			const int njet50_smeared = CountJets(smearedGenJets, 50.0, 2.5);
			if (njet50_smeared>= NJet50_min_ && njet50_smeared<= NJet50_max_) 
			{
				FillHistogram(smearedGenJets,2, smearingWgt * lumiWgt);
			}
		}
		
		if (bDEBUG) cout << "size reco/gen/smeared = " << recoJets.size() << "/" << genJets.size() 
					<< "/" << smearedGenJets.size() << endl;

	} // event loop


	/***** NORMALIZE VARIABLE BINNED HIST BY BIN WIDTH *****/

	TCanvas *pf = new TCanvas("pf");
	gPad->Print("passfail.eps[");
	for (unsigned htbin=0; htbin < HtBins_.size() -1 ; ++htbin)
	{
		for (unsigned mhtbin=0; mhtbin < MhtBins_.size() -1 ; ++mhtbin)
		{
			for (int i=0; i<6; ++i)
			{
				DivideByBinWidth(Hist.at(htbin).at(mhtbin).hv_SmearedEvt.pass[i]);
				DivideByBinWidth(Hist.at(htbin).at(mhtbin).hv_SmearedEvt.fail[i]);
					
				Hist.at(htbin).at(mhtbin).hv_SmearedEvt.pass[i]->Draw();
				gPad->Print("passfail.eps");
				Hist.at(htbin).at(mhtbin).hv_SmearedEvt.fail[i]->Draw();
				gPad->Print("passfail.eps");
			}
		}
	}

	gPad->Print("passfail.eps]");
	delete pf;




	/****************************************************************************
	 * FOR QUICK DEBUGGING ONLY
	 ***************************************************************************/
   //if (smearFunc_) delete smearFunc_;
	TCanvas *c = new TCanvas();
	gStyle->SetOptStat(0);
	gPad->Print("samples.eps[");

	for (unsigned htbin=0; htbin < HtBins_.size() -1 ; ++htbin)
	{
		for (unsigned mhtbin=0; mhtbin < MhtBins_.size() -1 ; ++mhtbin)
		{
			TH1* hrecomht = (TH1*) Hist.at(htbin).at(mhtbin).hv_RecoEvt.h_Mht->Clone("reco_mht");
			hrecomht->SetDirectory(0);
			TH1* hsmearmht = (TH1*) Hist.at(htbin).at(mhtbin).hv_SmearedEvt.h_Mht->Clone("smear_mht");
			hsmearmht->SetDirectory(0);
			TH1* hrecodphimin = (TH1*) Hist.at(htbin).at(mhtbin).hv_RecoEvt.h_DphiMin->Clone("reco_dphimin");
			hrecodphimin->SetDirectory(0);
			TH1* hsmeardphimin = (TH1*) Hist.at(htbin).at(mhtbin).hv_SmearedEvt.h_DphiMin->Clone("smear_dphimin");
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
			//

			TH1* hreconjet50 = (TH1*) Hist.at(htbin).at(mhtbin).hv_RecoEvt.h_Njet50eta2p5->Clone("reco_njet50");
			hreconjet50->SetDirectory(0);
			TH1* hsmearnjet50 = (TH1*) Hist.at(htbin).at(mhtbin).hv_SmearedEvt.h_Njet50eta2p5->Clone("smear_njet50");
			hsmearnjet50->SetDirectory(0);
			stringstream reconjet50leg, smnjet50leg;
			const double sum_reconj = hreconjet50->Integral(1,hreconjet50->GetNbinsX()+1);
			const double sum_smearnj= hsmearnjet50->Integral(1, hsmearnjet50->GetNbinsX()+1);
			reconjet50leg << "reco (" << sum_reconj << ")";
			smnjet50leg << "smeared (" << sum_smearnj << ")";


			cout << red << "mht reco/smear  = " << sum_recomht << " / " << sum_smearmht << endl; 
			cout << "nj50 reco/smear = " << sum_reconj << " / " << sum_smearnj  << clearatt << endl; 

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
				TH1* hrecojetpt = (TH1*) Hist.at(htbin).at(mhtbin).hv_RecoJets.at(i).h_Jet_pt->Clone("reco_jet");
				hrecojetpt->SetDirectory(0);
				TH1* hgenjetpt = (TH1*) Hist.at(htbin).at(mhtbin).hv_GenJets.at(i).h_Jet_pt->Clone("gen_jet");
				hgenjetpt->SetDirectory(0);
				TH1* hsmearjetpt = (TH1*) Hist.at(htbin).at(mhtbin).hv_SmearedJets.at(i).h_Jet_pt->Clone("smear_jet");
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
				hsmearjetpt_copy->Divide(hrecojetpt);
				hsmearjetpt_copy->SetTitle("Smeared Gen/Reco");
				hsmearjetpt_copy->Draw();
				//gPad->Print("samples.eps");

				delete l;
			}
		}
	}
	gPad->Print("samples.eps]");


	vector<double> ptBins  = smearFunc_->GetPtBinEdges();
	vector<double> etaBins = smearFunc_->GetEtaBinEdges();

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
	//end job summary
	cout << ">>>>>>> " << __FILE__ << ":" << __FUNCTION__ << ": End Job " << endl;
	if (bDEBUG) cout << red << " ------ DEBUG MODE ---------- " << clearatt << endl;
	cout << "Entries found/processed = " << nentries << " / " << nProcessed << endl;
	cout << "---------- Settings --------------- " << endl;
	cout << "Lumi wgt             = " << lumiWgt << endl;
	cout << "Smear wgt            = " << smearingWgt << endl;
	cout << "nTries_              = " << nTries_ << endl;
	cout << "smearedJetPt_        = " << smearedJetPt_ << endl;
	cout << "absoluteTailScaling_ = " << smearFunc_->GetAbsoluteTailScaling() << endl;
	cout << "MHTcut_low_          = " << MHTcut_low_  << endl;
	cout << "MHTcut_medium_       = " << MHTcut_medium_  << endl;
	cout << "MHTcut_high_         = " << MHTcut_high_  << endl;
	cout << "HTcut_low_           = " << HTcut_low_  << endl;
	cout << "HTcut_medium_        = " << HTcut_medium_  << endl;
	cout << "HTcut_high_          = " << HTcut_high_  << endl;
	cout << "HTcut_veryhigh_      = " << HTcut_veryhigh_  << endl;
	cout << "HTcut_extremehigh_   = " << HTcut_extremehigh_  << endl;

}

double FactorizationBySmearing::HT (const std::vector<TLorentzVector>& vjets) {
	double ht = 0.0;
	for( unsigned int ijet=0; ijet<vjets.size(); ijet++) {    
		if( vjets[ijet].Pt()>50.0 && std::abs(vjets[ijet].Eta())<2.5 ) 
			ht += vjets[ijet].Pt();
	}
	return ht;
}

TVector3 FactorizationBySmearing::GetMHT(vector<TLorentzVector> jets){
	TVector3 MHT;
	MHT.SetXYZ(0.0,0.0,0.0);
	for(int i = 0; i < jets.size(); i++){
		MHT -= jets[i].Vect();
	}
	MHT.SetZ(0.0);
	return MHT;
}

TLorentzVector FactorizationBySmearing::MHT(const std::vector<TLorentzVector>& vjets) {;
	TLorentzVector mht(0.0, 0.0, 0.0, 0.0);
	for( unsigned int ijet=0; ijet<vjets.size(); ijet++) {    
		if( vjets[ijet].Pt()>30.0 && std::abs(vjets[ijet].Eta())<5.0 ) 
			mht -= vjets[ijet];
	}
	//cout << __FUNCTION__ << "mht pt/phi = " << mht.Pt() << "/" << mht.Phi() << endl;

	return mht;
}

void FactorizationBySmearing::CreateRecoJetVec(std::vector<TLorentzVector>& vjets) {

	vjets.clear();
	for (unsigned i=0; i < t_PFJetPt->size(); i++)
	{
		 //this is the correct way to set vector
		 //with the info I have save in the tree
		TLorentzVector tl(0,0,0,0); 
		tl.SetPtEtaPhiE( (*t_PFJetPt)[i], (*t_PFJetEta)[i],
							  (*t_PFJetPhi)[i], (*t_PFJetE)[i]   );
		vjets.push_back(tl);

		//if (bDEBUG) cout << __FUNCTION__ << "::jet [" << i 
		//		<< "] pt/eta/phi = " << tl.Pt() << "/" 
		//		<< tl.Eta() << "/" << tl.Phi() <<endl;
	}

	std::sort(vjets.begin(), vjets.end(), PtAComparator);

}

void FactorizationBySmearing::CreateGenJetVec(std::vector<TLorentzVector>& vjets) {

	vjets.clear();
	for (unsigned i=0; i < t_genJetPt->size(); i++)
	{
		TLorentzVector tl(0,0,0,0);
		tl.SetPtEtaPhiE( t_genJetPt->at(i), t_genJetEta->at(i),
							 t_genJetPhi->at(i), t_genJetE->at(i));
		vjets.push_back(tl);
		//if (bDEBUG) cout << __FUNCTION__ << "::jet [" << i 
		//		<< "] pt/eta/phi = " << tl.Pt() << "/" 
		//		<< tl.Eta() << "/" << tl.Phi() <<endl;
	}
	std::sort(vjets.begin(), vjets.end(), PtAComparator);
}

//--------------------------------------------------------------------------
void FactorizationBySmearing::SmearingGenJets(const vector<TLorentzVector>& jets_gen, 
		std::vector<TLorentzVector>& genJets_smeared) {

	if (bDEBUG) cout << __FUNCTION__ << ":" << __LINE__ << endl;

	genJets_smeared.clear();
	int i_jet = 0;

	for (vector<TLorentzVector>::const_iterator it = jets_gen.begin(); it != jets_gen.end(); ++it) {

		if (it->Pt() > smearedJetPt_) {
			const double scale  = JetResolutionHist_Pt_Smear(it->Pt(), it->Eta(), i_jet); //what is this i_jet purpose??? jet id index

			const double newPt  = it->Pt() * scale;
			const double newEta = it->Eta();
			const double newPhi = it->Phi();
			const double newM   = it->M();
			if (bDEBUG) cout << "old/new pt = " << it->Pt() << "/" << newPt << endl; 

			//for JER reconstruction for debugging
			int i_Pt  = GetIndex(it->Pt(), &PtBinEdges_);
			int i_eta = GetIndex(it->Eta(), &EtaBinEdges_);
			jerHist.at(i_Pt).at(i_eta)->Fill(scale);

          //pat::Jet::PolarLorentzVector newP4(newPt, newEta, newPhi, it->mass());
			 //typedef math::PtEtaPhiMLorentzVector PolarLorentzVector;
			TLorentzVector newP4(0,0,0,0);
			//cout << __LINE__ << "it->E / it->M = " << it->E() << " / " << it->M() << endl;
			newP4.SetPtEtaPhiM(newPt, newEta, newPhi, it->M());
			TLorentzVector smearedJet(0,0,0,0);
			smearedJet.SetPxPyPzE(newP4.Px(), newP4.Py(), newP4.Pz(), newP4.E());
			//smearedJet.SetPxPyPzE(newPt, newEta, newPhi, newM);
			genJets_smeared.push_back(smearedJet);
			++i_jet;
		} else {
			TLorentzVector smearedJet(*it);
			genJets_smeared.push_back(smearedJet);
		}
	}
	std::sort(genJets_smeared.begin(), genJets_smeared.end(), PtAComparator);
	if (bDEBUG) cout << __FUNCTION__ << ":" << __LINE__ << endl;
}

//--------------------------------------------------------------------------
// pt resolution for smearing
double FactorizationBySmearing::JetResolutionHist_Pt_Smear(const double& pt, const double& eta, const int& i) {
	if (bDEBUG) cout << __FUNCTION__ << ":" << __LINE__ << endl;
   int i_jet;
   i < 2 ? i_jet = i : i_jet = 2;
   int i_Pt = GetIndex(pt, &PtBinEdges_);
   int i_eta = GetIndex(eta, &EtaBinEdges_);

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
int FactorizationBySmearing::GetIndex(const double& x, const std::vector<double>* vec) {
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

void FactorizationBySmearing::PtBinEdges_scaling(const string s)
{
	//PtBinEdges_scaling_ = NumStringToVec(s); 
}  

void FactorizationBySmearing::EtaBinEdges_scaling(const string s)
{
	//EtaBinEdges_scaling_ = NumStringToVec(s); 
}  

void FactorizationBySmearing::AdditionalSmearing(const string s)
{
	//AdditionalSmearing_ = NumStringToVec(s); 
}  

void FactorizationBySmearing::LowerTailScaling(const string s)
{
	//LowerTailScaling_ = NumStringToVec(s); 
}  

void FactorizationBySmearing::UpperTailScaling(const string s)
{
	//UpperTailScaling_ = NumStringToVec(s); 
}  
void FactorizationBySmearing::PtBinEdges(const string s)
{
   //PtBinEdges_ = NumStringToVec(s);
}
void FactorizationBySmearing::EtaBinEdges(const string s)
{
   //EtaBinEdges_ = NumStringToVec(s);
}

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


void FactorizationBySmearing::FillHistogram(const vector<TLorentzVector>& jets, 
				const int jetcoll, const double& wgt)
{
	//0=reco
	//1=gen
	//2=smeared
	
	const TLorentzVector mhtvec(MHT(jets)); 
	const double mht    = mhtvec.Pt(); 
	const double ht     = HT(jets);
	const float dPhiMin = DelPhiMin(jets,	mhtvec);
	const bool bPassRA2dphiCut = PassRA2dphiCut(jets,	mhtvec);

	/**************************************************/
	//                 for systematics
	/**************************************************/
	bool bSidebandSyst1 = true;
	if (jets.size() >= 1) bSidebandSyst1 = bSidebandSyst1 && (fabs(jets.at(0).DeltaPhi(mhtvec)) > 0.5);
	if (jets.size() >= 2) bSidebandSyst1 = bSidebandSyst1 && (fabs(jets.at(1).DeltaPhi(mhtvec)) > 0.5);
	if (jets.size() >= 3) bSidebandSyst1 = bSidebandSyst1 && (fabs(jets.at(2).DeltaPhi(mhtvec)) > 0.1);

	bool bSidebandSyst2 = true;
	if (jets.size() >= 1) bSidebandSyst2 = bSidebandSyst2 && (fabs(jets.at(0).DeltaPhi(mhtvec)) > 0.4);
	if (jets.size() >= 2) bSidebandSyst2 = bSidebandSyst2 && (fabs(jets.at(1).DeltaPhi(mhtvec)) > 0.4);
	if (jets.size() >= 3) bSidebandSyst2 = bSidebandSyst2 && (fabs(jets.at(2).DeltaPhi(mhtvec)) > 0.2);

	//const double mht = t_PFmht;

	const int njet50eta2p5 = CountJets(jets, 50.0, 2.5); 
	const int njet30eta5p0 = CountJets(jets, 30.0, 5.0); 

	const unsigned i_HtBin  = GetVectorIndex(HtBins_, ht);
	const unsigned i_MhtBin = GetVectorIndex(MhtBins_, mht);

	if (bDEBUG) 
	{
		cout << "[jetcoll = " << jetcoll << "] ht/mht/mymht = " 
				<< ht << "/" << mht 
				//<< "/" << mymht 
				<< endl;
	}

	if (jetcoll == 0) //reco hists
	{ 
		Hist.at(i_HtBin).at(i_MhtBin).hv_RecoEvt.h_Mht->Fill(mht, wgt);
		Hist.at(i_HtBin).at(i_MhtBin).hv_RecoEvt.h_Ht->Fill(ht, wgt);
		Hist.at(i_HtBin).at(i_MhtBin).hv_RecoEvt.h_Njet50eta2p5->Fill(njet50eta2p5, wgt);
		Hist.at(i_HtBin).at(i_MhtBin).hv_RecoEvt.h_Njet30eta5p0->Fill(njet30eta5p0, wgt);
		Hist.at(i_HtBin).at(i_MhtBin).hv_RecoEvt.h_DphiMin->Fill(dPhiMin, wgt);
		Hist.at(i_HtBin).at(i_MhtBin).hv_RecoEvt.h_DphiMinVsMht->Fill(mht, dPhiMin, wgt);
		FillJetHistogram(jets, Hist.at(i_HtBin).at(i_MhtBin).hv_RecoJets, jetcoll, mhtvec, wgt);
		/*const double mht_en = Hist.at(i_HtBin).at(i_MhtBin).hv_RecoEvt.h_Mht->GetEntries();
		const double mht_in = Hist.at(i_HtBin).at(i_MhtBin).hv_RecoEvt.h_Mht->Integral();
		const double j1_en = Hist.at(i_HtBin).at(i_MhtBin).hv_RecoJets.at(0).h_Jet_pt->GetEntries();
		const double j1_in = Hist.at(i_HtBin).at(i_MhtBin).hv_RecoJets.at(0).h_Jet_pt->Integral();
		const string mht_name(Hist.at(i_HtBin).at(i_MhtBin).hv_RecoEvt.h_Mht->GetName());
		const string j1_name(Hist.at(i_HtBin).at(i_MhtBin).hv_RecoJets.at(0).h_Jet_pt->GetName());
		cout << "jet coll[=" << jetcoll << "]= mht/j1 = "
					<< mht_en << "|" << mht_in << " / " << j1_en << "|" << j1_in  
					<< " [" << mht_name << "|" << j1_name
					<< endl;
					*/
	} else if (jetcoll == 1)  //gen hists
	{
		Hist.at(i_HtBin).at(i_MhtBin).hv_GenEvt.h_Mht->Fill(mht, wgt);
		Hist.at(i_HtBin).at(i_MhtBin).hv_GenEvt.h_Ht->Fill(ht, wgt);
		Hist.at(i_HtBin).at(i_MhtBin).hv_GenEvt.h_Njet50eta2p5->Fill(njet50eta2p5, wgt);
		Hist.at(i_HtBin).at(i_MhtBin).hv_GenEvt.h_Njet30eta5p0->Fill(njet30eta5p0, wgt);
		Hist.at(i_HtBin).at(i_MhtBin).hv_GenEvt.h_DphiMin->Fill(dPhiMin, wgt);
		Hist.at(i_HtBin).at(i_MhtBin).hv_GenEvt.h_DphiMinVsMht->Fill(mht, dPhiMin, wgt);
		FillJetHistogram(jets, Hist.at(i_HtBin).at(i_MhtBin).hv_GenJets, jetcoll, mhtvec, wgt);
	} else if (jetcoll == 2)  //smeared hists
	{
		Hist.at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.h_Mht->Fill(mht, wgt);
		Hist.at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.h_Ht->Fill(ht, wgt);
		Hist.at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.h_Njet50eta2p5->Fill(njet50eta2p5, wgt);
		Hist.at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.h_Njet30eta5p0->Fill(njet30eta5p0, wgt);
		Hist.at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.h_DphiMin->Fill(dPhiMin, wgt);
		Hist.at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.h_DphiMinVsMht->Fill(mht, dPhiMin, wgt);
		FillJetHistogram(jets, Hist.at(i_HtBin).at(i_MhtBin).hv_SmearedJets, jetcoll, mhtvec, wgt);


		if (bPassRA2dphiCut) 
		{
			Hist.at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.signal->Fill(mht,wgt);
			Hist.at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.signalFineBin->Fill(mht,wgt);
		}

		for (unsigned i=0; i < vDphiVariations.size(); ++i)
		{
			if (dPhiMin > vDphiVariations.at(i)) 
			{
				Hist.at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.pass[i]->Fill(mht, wgt);
				Hist.at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.passFineBin[i]->Fill(mht, wgt);
			} else 
			{
				Hist.at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.fail[i]->Fill(mht, wgt);
				Hist.at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.failFineBin[i]->Fill(mht, wgt);
			}
		}


		if (! bSidebandSyst1) //we want the denominator (or failed events)
		{
			Hist.at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.sidebandSyst[0]->Fill(mht, wgt);
			Hist.at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.sidebandSystFineBin[0]->Fill(mht, wgt);
		}
		if (! bSidebandSyst2)
		{
			Hist.at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.sidebandSyst[1]->Fill(mht, wgt);
			Hist.at(i_HtBin).at(i_MhtBin).hv_SmearedEvt.sidebandSystFineBin[1]->Fill(mht, wgt);
		}

	} else {
		cout << __FUNCTION__ << ": Invalid jet colelctions!" << endl;
		assert (false);
	}
	
}
void FactorizationBySmearing::FillJetHistogram(const vector<TLorentzVector>& jets, 
				vector<JetHist_t> hist, const int jetcoll, const TLorentzVector& mhtvec, const double& wgt) 
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

		if (pt>50 && fabs(eta)<2.5)
		{
			++njetsfound;
			const int j = njetsfound - 1;
			hist.at(j).h_Jet_pt->Fill(pt, wgt);
			hist.at(j).h_Jet_eta->Fill(eta, wgt);
			hist.at(j).h_Jet_phi->Fill(phi, wgt);
			hist.at(j).h_Jet_dphi->Fill(dphi, wgt);
		}
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
					const TLorentzVector& mhtVec)
{
	std::vector<float> vDelPhi_jetmht;

	for (unsigned i = 0 ; i < jets.size() ; ++i)
	{	//use only three leading jets
		//or in the case >=2 min jets, use only 2 jets
		if ( (i+1) > std::min(3, (int)NJet50_min_) ) break; 
		const float delphi_jetmht = fabs(jets.at(i).DeltaPhi(mhtVec));
		vDelPhi_jetmht.push_back(delphi_jetmht);
	}

	std::sort(vDelPhi_jetmht.begin(), vDelPhi_jetmht.end(), sort_using_less_than);	

	//assert (vDelPhi_jetmht.size() == 3 && "ERROR: More than 3 dPhiMin calculations found!");
	//std::cout << "vDelPhi_jetmht size = " << vDelPhi_jetmht.size() << std::endl;
	
	const float dPhiMin = vDelPhi_jetmht.at(0);
	//cout << __FUNCTION__ << ": delphimin = " <<  dPhiMin << endl;
	return dPhiMin;
}
bool FactorizationBySmearing::PassRA2dphiCut(const vector<TLorentzVector>& jets, 
					const TLorentzVector& mht)
{
	//Check if RA2 dphi cuts are satisfied
	bool bPassRA2dphiCut = true;
	if (jets.size() >=1) bPassRA2dphiCut = bPassRA2dphiCut && (std::fabs(jets.at(0).DeltaPhi(mht)) > 0.5);
	if (jets.size() >=2) bPassRA2dphiCut = bPassRA2dphiCut && (std::fabs(jets.at(1).DeltaPhi(mht)) > 0.5);
	if (jets.size() >=3) bPassRA2dphiCut = bPassRA2dphiCut && (std::fabs(jets.at(2).DeltaPhi(mht)) > 0.3);
	
	return bPassRA2dphiCut;
}

void FactorizationBySmearing::BookJerDebugHists()
{
	vector<double> ptBins = smearFunc_->GetPtBinEdges();
	vector<double> etaBins = smearFunc_->GetEtaBinEdges();
	
	TDirectory *dir = (TDirectory*) oFile->FindObject("Hist");
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
			const float ptmin = ptBins.at(ptbin);
			const float ptmax = ptBins.at(ptbin+1);
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

	cout << __FUNCTION__ << ": FATAL ERROR OCCURED!" << endl;
	cout << "val = " << val << " is out of bin ranges[ ";

	for (unsigned bin=0; bin < binEdges.size(); ++bin)
	{
		cout << binEdges.at(bin) << "/";
	}
	cout << "]" << endl;
	
	cout << "bin size = " << binEdges.size() << endl;
	cout << "Assigning the last bin for this extreme value!" << endl;
	return (binEdges.size() - 2);
}
double FactorizationBySmearing::GetLumiWgt(const string& datasetname, const double& dataLumi)
{
	/*************************************************
	 * X sec of Summer 12 Pythis MC samples
	 * **********************************************/

	//cross section for each sample in pt order
	const float xSec[] = {
		1759.549,  //pt 300-470
		113.8791,  //470-600
		26.9921,  //600-800
		3.550036, //8000-1000
		0.737844, //1000-1400
		0.03352235, //1400-1800
		0.001829005 //1800
	};

	const float nEvents[] = {
		5927300, // 300-470  #numbers from DBS, PREP page numbers are approximate
		3994848, // 470-600
		3992760, // 600-800
		3998563, //800-1000
		1964088, //1000-1400
		2000062, //1400-1800
		977586 //1800
	};

	double lumiWgt = 1;

	if( (std::string(datasetname)).find("qcd1") != string::npos)  
	{
		cout << "Lum for qcd 1 used." << endl; 
		lumiWgt = dataLumi/(nEvents[0] / xSec[0]);
	} else if( (std::string(datasetname)).find("qcd2") != string::npos)
	{
		cout << "Lum for qcd 2 used." << endl; 
		lumiWgt = dataLumi/(nEvents[1] / xSec[1]);
	} else if( (std::string(datasetname)).find("qcd3") != string::npos) 
	{
		cout << "Lum for qcd 3 used." << endl; 
		lumiWgt = dataLumi/(nEvents[2] / xSec[2]);
	} else if( (std::string(datasetname)).find("qcd4") != string::npos)  
	{
		cout << "Lum for qcd 4 used." << endl; 
		lumiWgt = dataLumi/(nEvents[3] / xSec[3]);
	} else if( (std::string(datasetname)).find("qcd5") != string::npos)  
	{
		cout << "Lum for qcd 5 used." << endl; 
		lumiWgt = dataLumi/(nEvents[4] / xSec[4]);
	} else if( (std::string(datasetname)).find("qcd6") != string::npos)  
	{ 
		cout << "Lum for qcd 6 used." << endl; 
		lumiWgt = dataLumi/(nEvents[5] / xSec[5]);
	} else if( (std::string(datasetname)).find("qcd7") != string::npos)  
	{
		cout << "Lum for qcd 7 used." << endl; 
		lumiWgt = dataLumi/(nEvents[6] / xSec[6]);
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
void FactorizationBySmearing::BookHistogram(const char *outFileName)
{
	oFile = new TFile(outFileName, "recreate");
	TDirectory *dir = oFile->mkdir("Hist");

	if (bDEBUG) cout << __FUNCTION__ << " dir = " << dir->GetName() << endl;

	for (unsigned htbin=0; htbin < HtBins_.size() -1; ++htbin)
	{
		vector<Hist_t> h;
		for (unsigned mhtbin=0; mhtbin < MhtBins_.size() -1; ++mhtbin)
		{
			const float htmin = HtBins_.at(htbin);
			const float htmax = HtBins_.at(htbin+1);
			const float mhtmin = MhtBins_.at(mhtbin);
			const float mhtmax = MhtBins_.at(mhtbin+1);
			stringstream folder;
			folder << "HT" << htmin << "to" << htmax
						<<  "MHT" << mhtmin << "to" << mhtmax;
			//oFile->ls();
			//cout << __LINE__ << ":" << folder.str() << endl; 
			//dir->cd();
			TDirectory *subdir = dir->mkdir(folder.str().c_str());
			//TFolder *subdir = dir->AddFolder(folder.str().c_str(), folder.str().c_str());
			//cout << "subdir = " << subdir->GetName() << endl;
			//subdir->ls();
			subdir->cd();
			//gDirectory->pwd();

			Hist_t hist;
			GetHist(subdir, hist, htmin, htmax, mhtmin, mhtmax);
			h.push_back(hist);
		}
		Hist.push_back(h);
	}

}

void FactorizationBySmearing::GetHist(TDirectory *dir, Hist_t& hist, 
					const float htmin, const float htmax,
					const float mhtmin, const float mhtmax
					)
{
	stringstream htmhtrange;
	htmhtrange << 	htmin << "<HT<" << htmax << " GeV, " << mhtmin << "<MHT<" << mhtmax << " GeV";
	
	vector<string> jetcoll;
	jetcoll.push_back("reco");
	jetcoll.push_back("gen");
	jetcoll.push_back("smeared");

	for (unsigned i = 0;i < jetcoll.size(); ++i)
	{
		stringstream njet50eta2p5title, njet30eta5p0title, httitle, mhttitle;
		njet50eta2p5title << htmhtrange.str().c_str() << ";Njets [ET>50 GeV,| #eta |<2.5];Events;";
		njet30eta5p0title << htmhtrange.str().c_str() << ";Njets [ET>30 GeV,| #eta |<5.0];Events;";
		httitle << htmhtrange.str().c_str() << ";HT [3 Jets, PT>50 GeV, | #eta |<2.5];Events;";
		mhttitle << htmhtrange.str().c_str() << ";MHT [PT>30 GeV, | #eta |<5.0];Events;";

		stringstream njet50eta2p5name, njet30eta5p0name, htname, mhtname;
		njet50eta2p5name << jetcoll.at(i) << "_njet50eta2p5";
		njet30eta5p0name << jetcoll.at(i) << "_njet30eta5p0";
		htname << jetcoll.at(i) << "_ht";
		mhtname << jetcoll.at(i) << "_mht";

		stringstream dphiminname, dphiminvsmhtname;
		stringstream dphimintitle, dphiminvsmhttitle;
		dphimintitle << htmhtrange.str().c_str() << ";#Delta #Phi_{min};Events;";
		dphiminname << jetcoll.at(i) << "_dphimin";
		dphiminvsmhttitle << htmhtrange.str().c_str() << ";MHT;#Delta #Phi_{min};";
		dphiminvsmhtname << jetcoll.at(i) << "_dphiminVsMht";
		

		if (i==0)
		{
			hist.hv_RecoEvt.h_Njet50eta2p5 = new TH1D(njet50eta2p5name.str().c_str(), njet50eta2p5title.str().c_str(), 20, 0, 20);
			hist.hv_RecoEvt.h_Njet30eta5p0 = new TH1D(njet30eta5p0name.str().c_str(), njet30eta5p0title.str().c_str(), 20, 0, 20);
			hist.hv_RecoEvt.h_Mht = new TH1D(mhtname.str().c_str(), mhttitle.str().c_str(), 500, 0, 1000);
			hist.hv_RecoEvt.h_Ht = new TH1D(htname.str().c_str(), httitle.str().c_str(), 300, 0, 3000);
			hist.hv_RecoEvt.h_DphiMin = new TH1D(dphiminname.str().c_str(), dphimintitle.str().c_str(), 125, 0, 2.5);
			hist.hv_RecoEvt.h_DphiMinVsMht = new TH2D(dphiminvsmhtname.str().c_str(), dphiminvsmhttitle.str().c_str(), 1400, 0, 350,250,0,2.5);
			hist.hv_RecoEvt.h_Njet50eta2p5->Sumw2();
			hist.hv_RecoEvt.h_Njet30eta5p0->Sumw2();
			hist.hv_RecoEvt.h_Mht->Sumw2();
			hist.hv_RecoEvt.h_Ht->Sumw2();
			hist.hv_RecoEvt.h_DphiMin->Sumw2();
			hist.hv_RecoEvt.h_DphiMinVsMht->Sumw2();

			dir->Add(hist.hv_RecoEvt.h_Njet50eta2p5);
			dir->Add(hist.hv_RecoEvt.h_Njet30eta5p0);
			dir->Add(hist.hv_RecoEvt.h_Mht);
			dir->Add(hist.hv_RecoEvt.h_Ht);
			dir->Add(hist.hv_RecoEvt.h_DphiMin);
			dir->Add(hist.hv_RecoEvt.h_DphiMinVsMht);

		} else if (i==1)
		{
			hist.hv_GenEvt.h_Njet50eta2p5 = new TH1D(njet50eta2p5name.str().c_str(), njet50eta2p5title.str().c_str(), 20, 0, 20);
			hist.hv_GenEvt.h_Njet30eta5p0 = new TH1D(njet30eta5p0name.str().c_str(), njet30eta5p0title.str().c_str(), 20, 0, 20);
			hist.hv_GenEvt.h_Mht = new TH1D(mhtname.str().c_str(), mhttitle.str().c_str(), 500, 0, 1000);
			hist.hv_GenEvt.h_Ht = new TH1D(htname.str().c_str(), httitle.str().c_str(), 300, 0, 3000);
			hist.hv_GenEvt.h_DphiMin = new TH1D(dphiminname.str().c_str(), dphimintitle.str().c_str(), 120, 0, 2.5);
			hist.hv_GenEvt.h_DphiMinVsMht = new TH2D(dphiminvsmhtname.str().c_str(), dphiminvsmhttitle.str().c_str(), 1400, 0, 350,250,0,2.5);
			hist.hv_GenEvt.h_Njet50eta2p5->Sumw2();
			hist.hv_GenEvt.h_Njet30eta5p0->Sumw2();
			hist.hv_GenEvt.h_Mht->Sumw2();
			hist.hv_GenEvt.h_Ht->Sumw2();
			hist.hv_GenEvt.h_DphiMin->Sumw2();
			hist.hv_GenEvt.h_DphiMinVsMht->Sumw2();

			dir->Add(hist.hv_GenEvt.h_Njet50eta2p5);
			dir->Add(hist.hv_GenEvt.h_Njet30eta5p0);
			dir->Add(hist.hv_GenEvt.h_Mht);
			dir->Add(hist.hv_GenEvt.h_Ht);
			dir->Add(hist.hv_GenEvt.h_DphiMin);
			dir->Add(hist.hv_GenEvt.h_DphiMinVsMht);
		} else if (i==2)
		{
			hist.hv_SmearedEvt.h_Njet50eta2p5 = new TH1D(njet50eta2p5name.str().c_str(), njet50eta2p5title.str().c_str(), 20, 0, 20);
			hist.hv_SmearedEvt.h_Njet30eta5p0 = new TH1D(njet30eta5p0name.str().c_str(), njet30eta5p0title.str().c_str(), 20, 0, 20);
			hist.hv_SmearedEvt.h_Mht = new TH1D(mhtname.str().c_str(), mhttitle.str().c_str(), 500, 0, 1000);
			hist.hv_SmearedEvt.h_Ht = new TH1D(htname.str().c_str(), httitle.str().c_str(), 300, 0, 3000);
			hist.hv_SmearedEvt.h_DphiMin = new TH1D(dphiminname.str().c_str(), dphimintitle.str().c_str(), 125, 0, 2.5);
			hist.hv_SmearedEvt.h_DphiMinVsMht = new TH2D(dphiminvsmhtname.str().c_str(), dphiminvsmhttitle.str().c_str(), 1400, 0, 350,250,0,2.5);
			hist.hv_SmearedEvt.h_Njet50eta2p5->Sumw2();
			hist.hv_SmearedEvt.h_Njet30eta5p0->Sumw2();
			hist.hv_SmearedEvt.h_Mht->Sumw2();
			hist.hv_SmearedEvt.h_Ht->Sumw2();
			hist.hv_SmearedEvt.h_DphiMin->Sumw2();
			hist.hv_SmearedEvt.h_DphiMinVsMht->Sumw2();

			dir->Add(hist.hv_SmearedEvt.h_Njet50eta2p5);
			dir->Add(hist.hv_SmearedEvt.h_Njet30eta5p0);
			dir->Add(hist.hv_SmearedEvt.h_Mht);
			dir->Add(hist.hv_SmearedEvt.h_Ht);
			dir->Add(hist.hv_SmearedEvt.h_DphiMin);
			dir->Add(hist.hv_SmearedEvt.h_DphiMinVsMht);
		}
	}

	GetJetHist(hist.hv_RecoJets,"reco", htmhtrange.str());
	GetJetHist(hist.hv_GenJets,"gen", htmhtrange.str());
	GetJetHist(hist.hv_SmearedJets,"smeared", htmhtrange.str());

	//for ratio plot
	const double evt_mht_max = 1500, evt_mht_bins = 750;
	const double evt_ht_max = 4000, evt_ht_bins = 800;
	const float npassFailHistBins = 14;
	const float passFailHistBins[] = {50,60,70,80,90,100,110,120,140,160,200,350,500,800,1000};
	//pass variations
	hist.hv_SmearedEvt.pass[0] = new TH1D("pass0","PASS from #Delta#Phi_{min} cut >0.15", npassFailHistBins, passFailHistBins); 
	hist.hv_SmearedEvt.pass[1] = new TH1D("pass1","PASS from #Delta#Phi_{min} cut >0.20", npassFailHistBins, passFailHistBins); 
	hist.hv_SmearedEvt.pass[2] = new TH1D("pass2","PASS from #Delta#Phi_{min} cut >0.25", npassFailHistBins, passFailHistBins); 
	hist.hv_SmearedEvt.pass[3] = new TH1D("pass3","PASS from #Delta#Phi_{min} cut >0.30", npassFailHistBins, passFailHistBins); 
	hist.hv_SmearedEvt.pass[4] = new TH1D("pass4","PASS from #Delta#Phi_{min} cut >0.35", npassFailHistBins, passFailHistBins); 
	hist.hv_SmearedEvt.pass[5] = new TH1D("pass5","PASS from #Delta#Phi_{min} cut >0.40", npassFailHistBins, passFailHistBins); 

	hist.hv_SmearedEvt.fail[0] = new TH1D("fail0","FAIL from #Delta#Phi_{min} cut <0.15", npassFailHistBins, passFailHistBins);
	hist.hv_SmearedEvt.fail[1] = new TH1D("fail1","FAIL from #Delta#Phi_{min} cut <0.20", npassFailHistBins, passFailHistBins);
	hist.hv_SmearedEvt.fail[2] = new TH1D("fail2","FAIL from #Delta#Phi_{min} cut <0.25", npassFailHistBins, passFailHistBins);
	hist.hv_SmearedEvt.fail[3] = new TH1D("fail3","FAIL from #Delta#Phi_{min} cut <0.30", npassFailHistBins, passFailHistBins);
	hist.hv_SmearedEvt.fail[4] = new TH1D("fail4","FAIL from #Delta#Phi_{min} cut <0.35", npassFailHistBins, passFailHistBins);
	hist.hv_SmearedEvt.fail[5] = new TH1D("fail5","FAIL from #Delta#Phi_{min} cut <0.40", npassFailHistBins, passFailHistBins);

	hist.hv_SmearedEvt.passFineBin[0] = new TH1D("passFineBin0","PASS from #Delta#Phi_{min} cut >0.15", evt_mht_bins, 0, evt_mht_max); 
	hist.hv_SmearedEvt.passFineBin[1] = new TH1D("passFineBin1","PASS from #Delta#Phi_{min} cut >0.20", evt_mht_bins, 0, evt_mht_max); 
	hist.hv_SmearedEvt.passFineBin[2] = new TH1D("passFineBin2","PASS from #Delta#Phi_{min} cut >0.25", evt_mht_bins, 0, evt_mht_max); 
	hist.hv_SmearedEvt.passFineBin[3] = new TH1D("passFineBin3","PASS from #Delta#Phi_{min} cut >0.30", evt_mht_bins, 0, evt_mht_max); 
	hist.hv_SmearedEvt.passFineBin[4] = new TH1D("passFineBin4","PASS from #Delta#Phi_{min} cut >0.35", evt_mht_bins, 0, evt_mht_max); 
	hist.hv_SmearedEvt.passFineBin[5] = new TH1D("passFineBin5","PASS from #Delta#Phi_{min} cut >0.40", evt_mht_bins, 0, evt_mht_max); 

	hist.hv_SmearedEvt.failFineBin[0] = new TH1D("failFineBin0","FAIL from #Delta#Phi_{min} cut <0.15", evt_mht_bins, 0, evt_mht_max);
	hist.hv_SmearedEvt.failFineBin[1] = new TH1D("failFineBin1","FAIL from #Delta#Phi_{min} cut <0.20", evt_mht_bins, 0, evt_mht_max);
	hist.hv_SmearedEvt.failFineBin[2] = new TH1D("failFineBin2","FAIL from #Delta#Phi_{min} cut <0.25", evt_mht_bins, 0, evt_mht_max);
	hist.hv_SmearedEvt.failFineBin[3] = new TH1D("failFineBin3","FAIL from #Delta#Phi_{min} cut <0.30", evt_mht_bins, 0, evt_mht_max);
	hist.hv_SmearedEvt.failFineBin[4] = new TH1D("failFineBin4","FAIL from #Delta#Phi_{min} cut <0.35", evt_mht_bins, 0, evt_mht_max);
	hist.hv_SmearedEvt.failFineBin[5] = new TH1D("failFineBin5","FAIL from #Delta#Phi_{min} cut <0.40", evt_mht_bins, 0, evt_mht_max);

	for (int i=0; i<6; ++i) 
	{ 
		hist.hv_SmearedEvt.pass[i]->Sumw2(); hist.hv_SmearedEvt.fail[i]->Sumw2(); 
		hist.hv_SmearedEvt.passFineBin[i]->Sumw2(); hist.hv_SmearedEvt.failFineBin[i]->Sumw2(); 
	}

	hist.hv_SmearedEvt.sidebandSyst[0] = new TH1D("sidebandSyst1"," Sideband Syst1: FAIL from dphi (j1 & j2 , mht) >0.5 and dphi (j3, mht)>0.1 ", npassFailHistBins, passFailHistBins);
	hist.hv_SmearedEvt.sidebandSyst[1] = new TH1D("sidebandSyst2"," Sideband Syst2: FAIL from dphi (j1 & j2 , mht) >0.5 and dphi (j3, mht)>0.1 ", npassFailHistBins, passFailHistBins);
	hist.hv_SmearedEvt.sidebandSystFineBin[0] = new TH1D("sidebandSyst1_fineBin","Sideband Syst1 (FineBinned) : FAIL from dphi (j1 & j2 , mht) >0.5 and dphi (j3, mht)>0.1 ", 1500, 0, 1500);
	hist.hv_SmearedEvt.sidebandSystFineBin[1] = new TH1D("sidebandSyst2_fineBin","Sideband Syst2 (FineBinned) : FAIL from dphi (j1 & j2 , mht) >0.4 and dphi (j3, mht)>0.2 ", 1500, 0, 1500);

	hist.hv_SmearedEvt.signal = new TH1D("signal" ,"Signal Region", npassFailHistBins, passFailHistBins);
	hist.hv_SmearedEvt.signalFineBin = new TH1D("signalFineBin" ,"Signal Region", 1500, 0, 1500);

	hist.hv_SmearedEvt.sidebandSyst[0]->Sumw2();
	hist.hv_SmearedEvt.sidebandSyst[1]->Sumw2();
	hist.hv_SmearedEvt.sidebandSystFineBin[0]->Sumw2();
	hist.hv_SmearedEvt.sidebandSystFineBin[1]->Sumw2();

	hist.hv_SmearedEvt.signal->Sumw2();
	hist.hv_SmearedEvt.signalFineBin->Sumw2();
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
		ptname << jetcoll << "jet" << index << "_pt";
		etaname << jetcoll << "jet" << index << "_eta";
		phiname << jetcoll << "jet" << index << "_phi";
		dphiname << jetcoll << "jet" << index << "_dphi";
		pttitle << htmhtrange << ";Jet-" << index << " P_{T} [GeV];Events;";
		etatitle << htmhtrange << ";Jet-" << index << " #eta;Events;";
		pttitle << htmhtrange << ";Jet-" << index << " #phi ;Events;";
		pttitle << htmhtrange << ";#delta #phi [Jet-" << index << ", MHT];Events;";

		JetHist_t hist;
		hist.h_Jet_pt   = new TH1D(ptname.str().c_str(), pttitle.str().c_str(),   125, 0.0, 2500.0);  
		hist.h_Jet_eta  = new TH1D(etaname.str().c_str(), etatitle.str().c_str(),  120, -6.0, 6.0);  
		hist.h_Jet_phi  = new TH1D(phiname.str().c_str(), phititle.str().c_str(),  70, -3.5, 3.5);  
		hist.h_Jet_dphi  = new TH1D(dphiname.str().c_str(), dphititle.str().c_str(),  70, 0, 3.5);  
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
	for (int bin=1; bin<=h->GetNbinsX(); ++bin)
	{
		const double v = h->GetBinContent(bin);
		const double e = h->GetBinError(bin);
		const double w = h->GetBinWidth(bin);
		const double nv = v/w;
		const double ne = e/w;
		h->SetBinContent(bin, nv);
	}

}


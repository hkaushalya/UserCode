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

using namespace std;

static bool sort_using_less_than(double u, double v)
{
	   return u < v;
}

int main(int argc, char* argv[])
{

	if (argc < 3) {
		cerr <<"Please give 5 arguments " << "runList " << " " << "outputFileName" << " " << "# of evts  nJet50Min  nJet50Max" << endl;
		cerr <<" Valid configurations are " << std::endl;
		cerr << " ./optimize runlist_ttjets.txt isoplots.root qcd.files 100 2 5" << std::endl;
		return -1;
	}

	const char *inputFileList = argv[1];
	const char *outFileName   = argv[2];
	const char *arg3          = argv[3];
	const char *arg4          = argv[4];
	const char *arg5          = argv[5];
	int evts    = -1;

	if (isdigit(arg3[0]))
	{
		evts = atoi(argv[3]);
	} else 
	{
		cout << "argument 3 is not a number. using default value for evts = " << evts << endl;
	}

	//FactorizationBySmearing optimize(inputFileList, outFileName, data);
	FactorizationBySmearing smear(inputFileList, outFileName);
	//cout << "dataset " << dataset << " " << endl;

	smear.EventLoop(inputFileList, evts);

	return 0;
}

void FactorizationBySmearing::EventLoop(const char *datasetname, const int evts2Process) {
//void FactorizationBySmearing::EventLoop(const char *inputFileList) {
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
	absoluteTailScaling_ = false;
	smearedJetPt_        = 0.;
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


	//difference variation of the dPhiMin selections.
	//make sure the book the correct number of histograms 
	//when these are changed!!!
	vDphiVariations.push_back(0.15);
	vDphiVariations.push_back(0.20);
	vDphiVariations.push_back(0.25);
	vDphiVariations.push_back(0.30);
	vDphiVariations.push_back(0.35);
	vDphiVariations.push_back(0.40);



	Long64_t nentries = fChain->GetEntriesFast();
	//cout << "nentries " << nentries << endl;
	if (evts2Process>0 && evts2Process < nentries) 
	{	
	//	nentries = evts2Process;
		cout << "Requested events to process = " << nentries << endl;
	}

	double evtWgt = 1.0;
	const double datasetnameLumi = 10000; // 10 fb-1
	cout << "datasetname = " << datasetname  << endl;
	if( (std::string(datasetname)).find("qcd1") != string::npos)  
	{
		cout << "Lum for qcd 1 used." << endl; 
		evtWgt = datasetnameLumi/(nEvents[0] / xSec[0]);
	} else if( (std::string(datasetname)).find("qcd2") != string::npos)
	{
		cout << "Lum for qcd 2 used." << endl; 
		evtWgt = datasetnameLumi/(nEvents[1] / xSec[1]);
	} else if( (std::string(datasetname)).find("qcd3") != string::npos) 
	{
		cout << "Lum for qcd 3 used." << endl; 
		evtWgt = datasetnameLumi/(nEvents[2] / xSec[2]);
	} else if( (std::string(datasetname)).find("qcd4") != string::npos)  
	{
		cout << "Lum for qcd 4 used." << endl; 
		evtWgt = datasetnameLumi/(nEvents[3] / xSec[3]);
	} else if( (std::string(datasetname)).find("qcd5") != string::npos)  
	{
		cout << "Lum for qcd 5 used." << endl; 
		evtWgt = datasetnameLumi/(nEvents[4] / xSec[4]);
	} else if( (std::string(datasetname)).find("qcd6") != string::npos)  
	{ 
		cout << "Lum for qcd 6 used." << endl; 
		evtWgt = datasetnameLumi/(nEvents[5] / xSec[5]);
	} else if( (std::string(datasetname)).find("qcd7") != string::npos)  
	{
		cout << "Lum for qcd 7 used." << endl; 
		evtWgt = datasetnameLumi/(nEvents[6] / xSec[6]);
	}
	

	std::cout << red << "Dataset = " << datasetname << ": evtWeight " << evtWgt << clearatt << std::endl;

	Long64_t nbytes = 0, nb = 0;
	int decade = 0;
	
	vector<TLorentzVector> recoJets, genJets, smearedGenJets;
	//debug = true;
	debug = false;
	
   smearFunc_ = new SmearFunction();
	if (debug) cout << __FUNCTION__ << ":" << __LINE__ << endl;

	//nentries = 50;
	unsigned n500 = 0, n400 =0, n300 =0, n200=0, n150=0, n80 =0 , n0 = 0;
	unsigned nProcessed = 0;

	NJet50_min_ = 3;
	NJet50_max_ = 1000;

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

		if (debug) cout << red << "run/ls/evt= " << t_EvtRun << " / " 
				<<  t_EvtLS << " / " <<  t_EvtEvent << clearatt  << endl;

		CreateRecoJetVec(recoJets);
		CreateGenJetVec(genJets);
		const int njet50_reco = CountJets(recoJets, 50, 2.5);
		const int njet50_gen = CountJets(genJets, 50, 2.5);
		if (njet50_reco>= NJet50_min_ && njet50_reco< NJet50_max_) FillHistogram(recoJets,0, evtWgt);
		if (njet50_gen>= NJet50_min_ && njet50_gen< NJet50_max_)  FillHistogram(genJets,1, evtWgt);
		const double recoHT= HT(recoJets);
		const double genHT = HT(genJets);
		const double recoMHT= (MHT(recoJets)).Pt();
		const double genMHT = (MHT(genJets)).Pt();

		//do smart sampling
		unsigned nTries_ = 1000; //default
		/*if (debug) nTries_ = 1;
		else
		{
			if (t_PFmht>500) { nTries_ = 2000; ++n500; }
			else if (t_PFmht>400) { nTries_ = 800; ++n400; }
			else if (t_PFmht>300) { nTries_ = 500; ++n300; }
			else if (t_PFmht>200) { nTries_ = 200; ++n200; }
			else if (t_PFmht>150) { nTries_ = 100; ++n150; }
			else if (t_PFmht>80) { nTries_ = 50; ++n80; }
			else { ++n0; }
		}*/

		const double smearingWgt = 1.0/(double)nTries_;

		for (unsigned n = 0; n < nTries_; ++n)
		{
			SmearingGenJets(genJets, smearedGenJets);
			//DumpJets(genJets);
			//DumpJets(smearedGenJets);

			const double smHT = HT(smearedGenJets);
			const double smMHT = (MHT(smearedGenJets)).Pt();
			if (debug) 
			{
				cout << "entry/recoht/genht/smht=" << jentry 
					<< "/" << recoHT << "/" << genHT << "/" << smHT << endl;
				cout << "entry/recomht/genmht/smmht=" << jentry 
					<< "/" << recoMHT << "/" << genMHT << "/" << smMHT << endl;
			}
			const int njet50_smeared = CountJets(smearedGenJets, 50, 2.5);
			if (njet50_smeared>= NJet50_min_ && njet50_smeared< NJet50_max_) 
			{
				FillHistogram(smearedGenJets,2, smearingWgt * evtWgt);
			}
		}
		
		if (debug) cout << "size reco/gen/smeared = " << recoJets.size() << "/" << genJets.size() 
					<< "/" << smearedGenJets.size() << endl;

	} // event loop

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
				TH1* hrecodphimin = (TH1*) Hist.at(htbin).at(mhtbin).hv_RecoEvt.h_Mht->Clone("reco_dphimin");
				hrecodphimin->SetDirectory(0);
				TH1* hsmeardphimin = (TH1*) Hist.at(htbin).at(mhtbin).hv_SmearedEvt.h_Mht->Clone("smear_dphimin");
				hsmeardphimin->SetDirectory(0);


				hrecomht->SetLineWidth(2);
				hsmearmht->SetLineWidth(2);
				hsmearmht->SetLineColor(kRed);
				hrecomht->SetMarkerStyle(kPlus);
				hsmearmht->SetMarkerStyle(24);
				hsmearmht->SetMarkerColor(kRed);

				stringstream recomhtleg,smearmhtleg;
				recomhtleg << "reco (" << hrecomht->Integral() << ")";
				smearmhtleg << "smeared (" << hsmearmht->Integral() << ")";
				
				TLegend *l2 = new TLegend(0.6,0.6,0.9,0.9);
				l2->AddEntry(hrecomht, recomhtleg.str().c_str());
				l2->AddEntry(hsmearmht, smearmhtleg.str().c_str());

				hrecomht->DrawCopy();
				hsmearmht->DrawCopy("same");
				l2->Draw();
				gPad->Print("samples.eps");


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
				hgenjetpt->DrawCopy("same");
				hsmearjetpt->DrawCopy("same");

				stringstream recoleg, genleg, smearleg;
				recoleg << "reco (" << hrecojetpt->Integral() << ")";
				genleg << "gen (" << hgenjetpt->Integral() << ")";
				smearleg << "smeared (" << hsmearjetpt->Integral() << ")";

				TLegend *l = new TLegend(0.6,0.6,0.9,0.9);
				l->AddEntry(hrecojetpt, recoleg.str().c_str());
				l->AddEntry(hgenjetpt, genleg.str().c_str());
				l->AddEntry(hsmearjetpt, smearleg.str().c_str());
				l->Draw();
				gPad->Print("samples.eps");

				TH1* hsmearjetpt_copy = (TH1*) hsmearjetpt->Clone("smear_jet_copy");
				hsmearjetpt_copy->Divide(hrecojetpt);
				hsmearjetpt_copy->SetTitle("Smeared Gen/Reco");
				hsmearjetpt_copy->Draw();
				gPad->Print("samples.eps");

				delete l;
			}
		}
	}
	gPad->Print("samples.eps]");


	//end job summary
	cout << ">>>>>>> " << __FILE__ << ":" << __FUNCTION__ << ": End Job " << endl;
	cout << "nTries[mht>500/400/300/200/150/80 = " << n500 
					<< "/" << n400 << "/" << n200 << "/" << n150 << "/" << n80 << endl; 
	cout << "Entries found/processed = " << nentries << " / " << nProcessed << endl;

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

		//if (debug) cout << __FUNCTION__ << "::jet [" << i 
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
		//if (debug) cout << __FUNCTION__ << "::jet [" << i 
		//		<< "] pt/eta/phi = " << tl.Pt() << "/" 
		//		<< tl.Eta() << "/" << tl.Phi() <<endl;
	}
	std::sort(vjets.begin(), vjets.end(), PtAComparator);
}

//--------------------------------------------------------------------------
void FactorizationBySmearing::SmearingGenJets(const vector<TLorentzVector>& jets_gen, 
		std::vector<TLorentzVector> &genJets_smeared) {

	if (debug) cout << __FUNCTION__ << ":" << __LINE__ << endl;

	genJets_smeared.clear();
	int i_jet = 0;

	if (debug) cout << __FUNCTION__ << ":" << __LINE__ << endl;
	for (vector<TLorentzVector>::const_iterator it = jets_gen.begin(); it != jets_gen.end(); ++it) {

		if (it->Pt() > smearedJetPt_) {
			double newPt = 0;
			newPt = it->Pt() * JetResolutionHist_Pt_Smear(it->Pt(), it->Eta(), i_jet); //what is this i_jet purpose???
			if (debug) cout << "old/new pt = " << it->Pt() << "/" << newPt << endl; 
			double newEta = it->Eta();
			double newPhi = it->Phi();
			double newM   = it->M();

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
	if (debug) cout << __FUNCTION__ << ":" << __LINE__ << endl;
	std::sort(genJets_smeared.begin(), genJets_smeared.end(), PtAComparator);
	if (debug) cout << __FUNCTION__ << ":" << __LINE__ << endl;

	//Fill HT and MHT prediction histos for i-th iteration of smearing
	//FillPredictionHistos_gen(genJets_smeared, i, weight_);

	return;
}

//--------------------------------------------------------------------------
// pt resolution for smearing
double FactorizationBySmearing::JetResolutionHist_Pt_Smear(const double& pt, const double& eta, const int& i) {
	if (debug) cout << __FUNCTION__ << ":" << __LINE__ << endl;
   int i_jet;
   i < 2 ? i_jet = i : i_jet = 2;
   int i_Pt = GetIndex(pt, &PtBinEdges_);
   int i_eta = GetIndex(eta, &EtaBinEdges_);

	const double res = smearFunc_->getSmearFunc(1,i_jet,i_eta,i_Pt)->GetRandom();
	if (debug) cout << __FUNCTION__ << ":" << __LINE__ << ": res = " << res << endl;
	if (debug) cout << __FUNCTION__ << ":" << __LINE__ << endl;

   return res;
}

//--------------------------------------------------------------------------
int FactorizationBySmearing::GetIndex(const double& x, const std::vector<double>* vec) {
   int index = -1;
   // this is a check
   //int index = 0;
	if (debug) cout << __FUNCTION__ << ":" << __LINE__ << ": x / vec size = " << x << " / " << vec->size() << endl; 
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
	
	const double ht = HT(jets);
	const TLorentzVector mhtvec(MHT(jets)); 
	const double mht = mhtvec.Pt(); 
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

	const int njet50eta2p5 = CountJets(jets, 50, 2.5); 
	const int njet30eta5p0 = CountJets(jets, 30, 5.0); 

	const unsigned i_HtBin  = GetVectorIndex(HtBins_, ht);
	const unsigned i_MhtBin = GetVectorIndex(MhtBins_, mht);

	if (debug) 
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

	for (unsigned i=0; i < njet2loop; ++i)
	{
		if (jets.at(i).Pt()>50 && fabs(jets.at(i).Eta())<2.5)
		{
			hist.at(i).h_Jet_pt->Fill(jets.at(i).Pt(), wgt);
			hist.at(i).h_Jet_eta->Fill(jets.at(i).Eta(), wgt);
			hist.at(i).h_Jet_phi->Fill(jets.at(i).Phi(), wgt);
			hist.at(i).h_Jet_dphi->Fill(fabs(jets.at(i).DeltaPhi(mhtvec)), wgt);
		}
	}
}

int FactorizationBySmearing::CountJets(const std::vector<TLorentzVector>& vjets, 
					const double minPt, const double maxEta)
{
	int njet = 0.0;
	for( unsigned int ijet=0; ijet<vjets.size(); ijet++) {    
		if( vjets[ijet].Pt()>minPt && std::fabs(vjets[ijet].Eta())<maxEta ) 
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

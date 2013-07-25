#include <iostream>
#include <string>
#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include <sstream>
#include "TStyle.h"
#include "assert.h"
#include "TF1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include <iomanip>
#include <algorithm>
#include <THStack.h>

using namespace std;

const static float fDATA_LUMI  = 4650.0; //pb-1
//const static float fDATA_LUMI  = 1150.0; //pb-1
//const static float fEWK_LUMI   = 81352581.0/31300.0;        // 31300 pb for NLO xsec and total number of events is 81352581 (sample1+sample2)
//const static float fTTBAR_LUMI = 3701947.0/165.0;
const static float fZNN_LUMI   = 3067017.0/42.8;              //x-sec from note, Anwar hd used 32.92 pb-1 LO
//const static string sDATA_FILE_NAME = "Data_MyPrescaleWgted.root"; //used in AN 02-08-2011
//PU wgted samples
const static string sDATA_FILE_NAME  = "Data_465_PrescaleWgted_03162012.root"; 
const static string sWJET_HT250TO300_FILE_NAME  = "Wjets_HT250to300_Fall11_PUwgted_03162012.root";  //Fall 11 sample
const static string sWJET_HT300TOINF_FILE_NAME  = "Wjets_HT300toInf_Fall11_PUwgted_03162012.root";  //Fall 11 sample
const static string sTTBAR_FILE_NAME = "TTbar_Fall11_PUwgted_03162012.root";  // Fall 11 sample
const static string sZNN_FILE_NAME   = "Zjets_Summ11_PUwgted_03162012.root";    // Summer 11 sample

const static float fWJETS_HT250TO300_LUMI =  9831277.0 / (248.7 * 0.14); 
const static float fWJETS_HT300TOINF_LUMI =  5363746.0 / (317.0 * 0.153); 
const static float fTTBAR_LUMI = 59590147/165.0;




TH1* GetZnnHist(const string histname)
{
	TFile* znnFile = new TFile(sZNN_FILE_NAME.c_str());
	if (znnFile->IsZombie()) { cout << "TTbar file not found!" << endl; assert (false); }

	TH1* znnHist = dynamic_cast<TH1*> (znnFile->Get(histname.c_str()));
	if (znnHist == NULL) { cout << "TTbar hist " << histname << " not found!" << endl; assert (false); }
	znnHist->Sumw2();
	znnHist->SetLineColor(kRed);

	const float int_b4 = znnHist->Integral(); 
	const float scale  = fDATA_LUMI/fZNN_LUMI;
	znnHist->Scale(scale);
	std::cout << __FUNCTION__ << ": scale = " << scale<< ":: Znn Integral b4/a4= "<< int_b4  << "/" << znnHist->Integral() << std::endl; 

/*	new TCanvas();
	znnHist->SetLineColor(kGreen);
	znnHist->DrawCopy();
	hist->DrawCopy("same");
*/
	return znnHist;
}


TH1* GetTTbarHist(const string histname)
{
	TFile* ttbarFile = new TFile(sTTBAR_FILE_NAME.c_str());
	if (ttbarFile->IsZombie()) { cout << "TTbar file not found!" << endl; assert (false); }

	TH1* ttbarHist = dynamic_cast<TH1*> (ttbarFile->Get(histname.c_str()));
	if (ttbarHist == NULL) { cout << "TTbar hist " << histname << " not found!" << endl; assert (false); }
	ttbarHist->Sumw2();
	ttbarHist->SetLineColor(kRed);

	const float int_b4 = ttbarHist->Integral(); 
	const float scale  = fDATA_LUMI/fTTBAR_LUMI;
	ttbarHist->Scale(scale);
	std::cout << __FUNCTION__ << ": scale = " << scale<< ":: TTbar Integral b4/a4= "<< int_b4  << "/" << ttbarHist->Integral() << std::endl; 

/*	new TCanvas();
	ttbarHist->SetLineColor(kGreen);
	ttbarHist->DrawCopy();
	hist->DrawCopy("same");
*/

	return ttbarHist;
}


TH1* GetEWKHist(const string ewkhistname)
{
	TFile* ewkFile1 = new TFile(sWJET_HT250TO300_FILE_NAME.c_str());
	TFile* ewkFile2 = new TFile(sWJET_HT300TOINF_FILE_NAME.c_str());
	if (ewkFile1->IsZombie()) { cout << "EWK file 1 not found!" << endl; assert (false); }
	if (ewkFile2->IsZombie()) { cout << "EWK file 2 not found!" << endl; assert (false); }

	TH1* ewkHist1 = dynamic_cast<TH1*> (ewkFile1->Get(ewkhistname.c_str()));
	if (ewkHist1 == NULL) { cout << "EWK hist " << ewkhistname << " not found in ewkFile 1!" << endl; assert (false); }
	ewkHist1->Sumw2();
	TH1* ewkHist2 = dynamic_cast<TH1*> (ewkFile2->Get(ewkhistname.c_str()));
	if (ewkHist2 == NULL) { cout << "EWK hist " << ewkhistname << " not found in ewkFile 2!" << endl; assert (false); }
	ewkHist2->Sumw2();
	//new TCanvas();
	//ewkHist->DrawCopy();

	const float int_b4 = ewkHist1->Integral(); 
	const float scale1  = fDATA_LUMI/fWJETS_HT250TO300_LUMI;
	const float scale2  = fDATA_LUMI/fWJETS_HT300TOINF_LUMI;

	ewkHist1->Scale(scale1);
	ewkHist2->Scale(scale2);
	ewkHist1->Add(ewkHist2);
	std::cout << __FUNCTION__ << ": scale1/scale2 = " << scale1 << "/" << scale2 << ":: EWK Integral b4/a4= "<< int_b4  << "/" << ewkHist1->Integral() << std::endl; 

	return ewkHist1;
}

TH1* GetQCDHist(const string histname)
{
	//cross section for each sample in pt order
	const float xSec[] = {
		115100,
		24260,
		1168,
		70.22,
		15.55,
		1.844,
		0.3321,
		0.01087,
		0.0003575 
	};

	const float nEvents[] = {
		6127528, //120 //ok
		6220160,	//170 //ok
		6432669,//300  //ok
		3990085,	//470 //ok
		4245695,	//600 //ok
		4053888,//800  //ok
		2093222,//1000 //ok
		2196200, //1400 //ok
		293139 //1800  //ok
	};
	const Int_t nBins = 9;
	TFile *files[nBins];
	TH1 *hists[nBins] = {0,0,0,0,0,0,0,0,0};

	files[0] = new TFile ("QCD_HTbin1.root");
	files[1] = new TFile ("QCD_HTbin2.root");
	files[2] = new TFile ("QCD_HTbin3.root");
	files[3] = new TFile ("QCD_HTbin4.root");
	files[4] = new TFile ("QCD_HTbin5.root");
	files[5] = new TFile ("QCD_HTbin6.root");
	files[6] = new TFile ("QCD_HTbin7.root");
	files[7] = new TFile ("QCD_HTbin8.root");
	files[8] = new TFile ("QCD_HTbin9.root");

	TH1 *res_hist = 0;

	for (int i=0; i<nBins; ++i)
	{
		if (files[i]->IsZombie())
		{
			cout << files[i]->GetName() << " not found!" <<  endl;
			assert (false);
		}
		hists[i] = dynamic_cast<TH1*> (files[i]->Get(histname.c_str()));
		if (hists[i] == 0)
		{
			cout << "hist_pass " << histname << " not found in " << files[i]->GetName() << "!" << endl;
			assert (false);
		} else 
		{
			hists[i]->Sumw2();
			const float scale = fDATA_LUMI/ ( nEvents[i] / xSec[i] );
			hists[i]->Scale(scale);

			if (i == 0) res_hist = dynamic_cast<TH1*> (hists[i]->Clone("histcopy"));
			else res_hist->Add(hists[i]);
		}
	}

	return res_hist;
}

void DumpHist(const TH1* hist)
{
	assert (hist != NULL && "DumpHist:: hist passed is null!"); 
	cout << std::setw(5) << "bin " << std::setw(15) << "[edges]" << std::setw(10) << "content" << std::setw(10) << "error" << endl;
	double sum = 0;
	for (int bin = 0; bin <= hist->GetNbinsX()+1; ++bin)
	{
		if (hist->GetBinContent(bin)>0)
		{
			cout  << std::setprecision(4) << std::setw(5) << bin << std::setw(15) << "[" 
			<< hist->GetBinLowEdge(bin) << ", " << hist->GetXaxis()->GetBinUpEdge(bin)<< "]" 
			<< std::setw(10) << hist->GetBinContent(bin) 
			<< std::setw(10) << hist->GetBinError(bin) << endl;
			if (hist->GetBinCenter(bin) > 500) sum += hist->GetBinContent(bin);
		}
	}

	cout << "Sum (>500) = " << sum << endl;
}

void overLayBackgrounds(const string histname, const string title, const int rebin=1, const string epsname="", const bool logScale =1) 
{

	TH1 *histDataOrig = 0;

	TFile* dataRootFile = new TFile (sDATA_FILE_NAME.c_str());
	if (dataRootFile->IsZombie())
	{
		cout << "Data root file not found!" <<  endl;
		assert (false);
	}

	histDataOrig = dynamic_cast<TH1*> (dataRootFile->Get(histname.c_str()));
	if (histDataOrig == 0 ) { cout << "Hist " << histname << " not found!" << endl; assert (false); }
	histDataOrig->Sumw2();
	histDataOrig->SetLineWidth(2);

	TH1*	histDataFinal = dynamic_cast<TH1*> (histDataOrig->Clone("histpass_copy"));

	/********************************************/
	//substract background contaminations
	/********************************************/
	TH1* ewkHist;
	TH1* ttbarHist;
	TH1* znnHist;

	ewkHist   = GetEWKHist  (histname);
	ttbarHist = GetTTbarHist(histname);
	znnHist   = GetZnnHist  (histname);
	TH1* qcdHist = GetQCDHist(histname);

	if (rebin>1)
	{
		histDataOrig->Rebin(rebin);
		histDataFinal->Rebin(rebin);
		ewkHist->Rebin(rebin);
		ttbarHist->Rebin(rebin);
		znnHist->Rebin(rebin);
		qcdHist->Rebin(rebin);
	}

	histDataFinal->Add(ewkHist,-1);
	histDataFinal->Add(ttbarHist,-1);
	histDataFinal->Add(znnHist,-1);

	const int lineWidth = 2;
	const int ewkColor = 3;
	const int znnColor = 6;
	const int ttbarColor = 4;
	const int qcdColor = 5;

	ewkHist->SetLineColor(ewkColor);
	ewkHist->SetFillColor(ewkColor);
	ewkHist->SetLineWidth(lineWidth);
	znnHist->SetLineColor(znnColor);
	znnHist->SetFillColor(znnColor);
	znnHist->SetLineWidth(lineWidth);
	ttbarHist->SetLineColor(ttbarColor);
	ttbarHist->SetFillColor(ttbarColor);
	ttbarHist->SetLineWidth(lineWidth);
	qcdHist->SetLineColor(qcdColor);
	qcdHist->SetFillColor(qcdColor);
	qcdHist->SetLineWidth(lineWidth);


	histDataFinal->SetLineColor(2);
	histDataFinal->SetLineWidth(lineWidth);
	histDataFinal->SetTitle(title.c_str());
	histDataFinal->SetMinimum(0.5e-1);


	TLegend *leg1  = new TLegend(0.55,0.6,0.9,0.9);
	leg1->AddEntry(histDataOrig,"DATA (before substraction)");
	leg1->AddEntry(histDataFinal,"DATA (after substraction)");
	leg1->AddEntry(ewkHist,"W+Jets (MC)");
	leg1->AddEntry(znnHist,"Z#rightarrow#nu#bar{#nu} (MC)");
	leg1->AddEntry(ttbarHist,"t#bar{t} (MC)");

/*	new TCanvas();
	gStyle->SetOptStat(0);
	gPad->SetLogy();
	gPad->SetTickx();
	gPad->SetTicky();
	histDataFinal->Draw();
	histDataOrig->Draw("same");
	ewkHist->Draw("same");
	znnHist->Draw("same");
	ttbarHist->Draw("same");
	qcdHist->Draw("same");
	leg1->Draw();
*/	


	
	TH1* totBgHist = (TH1*) (qcdHist->Clone("totalBackground"));
	totBgHist->Add(ewkHist);
	totBgHist->Add(znnHist);
	totBgHist->Add(ttbarHist);
	totBgHist->SetLineColor(kRed);
	totBgHist->SetMarkerColor(kRed);
	totBgHist->SetFillColor(10);

	TLegend *leg2  = new TLegend(0.65,0.65,0.9,0.9);
	leg2->SetTextFont(42);
	leg2->AddEntry(histDataOrig,"DATA");
	leg2->AddEntry(totBgHist,"Background");
	leg2->AddEntry(qcdHist,"QCD");
	leg2->AddEntry(ewkHist,"W+Jets");
	leg2->AddEntry(znnHist,"Z#rightarrow#nu#bar{#nu}");
	leg2->AddEntry(ttbarHist,"t#bar{t}");


	THStack *hs = new THStack ("hs", NULL);
	//hs->GetXaxis()->SetTitleFont(42);
	hs->SetMinimum(0.5e-1);
	hs->SetMaximum(1.2 * histDataOrig->GetBinContent(histDataOrig->GetMaximumBin()));
	hs->Add(znnHist);
	hs->Add(ewkHist);
	hs->Add(ttbarHist);
	hs->Add(qcdHist);
	hs->SetTitle(title.c_str());


	new TCanvas();
	if (logScale) gPad->SetLogy();
	hs->Draw("HIST");
	totBgHist->Draw("same");
	histDataOrig->Draw("same");
	leg2->Draw();


	gPad->Print(epsname.c_str());

}

void overLayBackgrounds()
{

	const bool logScale = 0;
	const int rebin = 4;
	overLayBackgrounds("factnomht/HT500to800/dphiMin", "500<HT<800 GeV, MHT>0 GeV;#delta #phi_{min};Events;", rebin, "HT500to800_nomht_dphimin_linear.eps", logScale);
	overLayBackgrounds("factnomht/HT800to1000/dphiMin", "800<HT<1000 GeV, MHT>0 GeV;#delta #phi_{min};Events;", rebin, "HT800to1000_nomht_dphimin_linear.eps", logScale);
	overLayBackgrounds("factnomht/HT1000to1200/dphiMin", "1000<HT<1200 GeV, MHT>0 GeV;#delta #phi_{min};Events;", rebin, "HT1000to1200_nomht_dphimin_linear.eps", logScale);
	overLayBackgrounds("factnomht/HT1200to1400/dphiMin", "1200<HT<1400 GeV, MHT>0 GeV;#delta #phi_{min};Events;", rebin, "HT1200to1400_nomht_dphimin_linear.eps", logScale);
	overLayBackgrounds("factnomht/HT1400to7000/dphiMin", "HT>1400 GeV, MHT>0 GeV;#delta#phi_{min};Events;", rebin, "HT1400to7000_nomht_dphimin_linear.eps", logScale);

	overLayBackgrounds("factmht200/HT500to800/dphiMin", "500<HT<800 GeV, MHT>200 GeV;#delta #phi_{min};Events;", rebin, "HT500to800_MHT200_dphimin_linear.eps", logScale);
	overLayBackgrounds("factmht200/HT800to1000/dphiMin", "800<HT<1000 GeV, MHT>200 GeV;#delta #phi_{min};Events;", rebin, "HT800to1000_MHT200_dphimin_linear.eps", logScale);
	overLayBackgrounds("factmht200/HT1000to1200/dphiMin", "1000<HT<1200 GeV, MHT>200 GeV;#delta #phi_{min};Events;", rebin, "HT1000to1200_MHT200_dphimin_linear.eps", logScale);
	overLayBackgrounds("factmht200/HT1200to1400/dphiMin", "1200<HT<1400 GeV, MHT>200 GeV;#delta #phi_{min};Events;", rebin, "HT1200to1400_MHT200_dphimin_linear.eps", logScale);
	overLayBackgrounds("factmht200/HT1400to7000/dphiMin", "HT>1400 GeV, MHT>200 GeV;#delta#phi_{min};Events;", rebin, "HT1400to7000_MHT200_dphimin_linear.eps", logScale);

	overLayBackgrounds("factmht350/HT500to800/dphiMin", "500<HT<800 GeV, MHT>350 GeV;#delta #phi_{min};Events;", rebin, "HT500to800_MHT350_dphimin_linear.eps", logScale);
	overLayBackgrounds("factmht350/HT800to1000/dphiMin", "800<HT<1000 GeV, MHT>350 GeV;#delta #phi_{min};Events;", rebin, "HT800to1000_MHT350_dphimin_linear.eps", logScale);
	overLayBackgrounds("factmht350/HT1000to1200/dphiMin", "1000<HT<1200 GeV, MHT>350 GeV;#delta #phi_{min};Events;", rebin, "HT1000to1200_MHT350_dphimin_linear.eps", logScale);
	overLayBackgrounds("factmht350/HT1200to1400/dphiMin", "1200<HT<1400 GeV, MHT>350 GeV;#delta #phi_{min};Events;", rebin, "HT1200to1400_MHT350_dphimin_linear.eps", logScale);
	overLayBackgrounds("factmht350/HT1400to7000/dphiMin", "HT>1400 GeV, MHT>350 GeV;#delta#phi_{min};Events;", rebin, "HT1400to7000_MHT350_dphimin_linear.eps", logScale);

//	overLayBackgrounds("factorization_ht500/Pass_RA2dphi_HT500", "HT>500 GeV: Events passing |#Delta#phi(jet 1-2, MHT)|>0.5 & |#Delta#phi(jet 3, MHT)|>0.3;MHT [GeV];Events;", 1, "ht500_pass.eps");
//	overLayBackgrounds("factorization_ht500/Fail_1", "HT>500 GeV: Events passing |#Delta#phi(jet 1-3, MHT)|<0.2;MHT [GeV];Events;", 1, "ht500_fail.eps");
//	overLayBackgrounds("factorization_ht500/mht", "HT>500 GeV: Events passing RA2 Selection Except #Delta#Phi;MHT [GeV];Events;", 8, "ht500_mht.eps");
//	overLayBackgrounds("factorization_ht500/ht", "HT>500 GeV: Events passing RA2 Selection Except #Delta#Phi; HT [GeV];Events;", 1, "ht500_ht.eps");
//	overLayBackgrounds("factorization_ht500/njet_et50eta24", "HT>500 GeV: Events passing RA2 Selection Except #Delta#Phi;NJet [E_{T}>50 GeV & | #eta |<2.5];Events;", 1, "ht500_njet.eps");
//	overLayBackgrounds("factorization_ht500/meteff", "HT>500 GeV: Events passing RA2 Selection Except #Delta#Phi;MEff;Events;", 1, "ht500_meff.eps");
//	overLayBackgrounds("factorization_ht500/pf30_jet1_pt", "HT>500 GeV: Events passing RA2 Selection Except #Delta#Phi;Jet1 P_{T} [GeV];Events;", 4, "ht500_jet1pt.eps");
//	overLayBackgrounds("factorization_ht500/pf30_jet2_pt", "HT>500 GeV: Events passing RA2 Selection Except #Delta#Phi;Jet2 P_{T} [GeV];Events;", 4, "ht500_jet2pt.eps");
//	overLayBackgrounds("factorization_ht500/pf30_jet3_pt", "HT>500 GeV: Events passing RA2 Selection Except #Delta#Phi;Jet3 P_{T} [GeV];Events;", 4, "ht500_jet3pt.eps");
//	overLayBackgrounds("factorization_ht500/pf30_jet1_eta", "HT>500 GeV: Events passing RA2 Selection Except #Delta#Phi;Jet1 #Eta;Events;", 2, "ht500_jet1eta.eps");
//	overLayBackgrounds("factorization_ht500/pf30_jet2_eta", "HT>500 GeV: Events passing RA2 Selection Except #Delta#Phi;Jet2 #Eta;Events;", 2, "ht500_jet2eta.eps");
//	overLayBackgrounds("factorization_ht500/pf30_jet3_eta", "HT>500 GeV: Events passing RA2 Selection Except #Delta#Phi;Jet3 #Eta;Events;", 2, "ht500_jet3eta.eps");

/*
	overLayBackgrounds("factorization_ht800/mht", "HT>800 GeV: Events passing RA2 Selection Except #Delta#Phi;MHT [GeV];Events;", 8, "ht800_mht.eps");
	overLayBackgrounds("factorization_ht800/ht", "HT>800 GeV: Events passing RA2 Selection Except #Delta#Phi; HT [GeV];Events;", 1, "ht800_ht.eps");
	overLayBackgrounds("factorization_ht800/njet_et50eta24", "HT>800 GeV: Events passing RA2 Selection Except #Delta#Phi;NJet [E_{T}>50 GeV & | #eta |<2.5];Events;", 1, "ht800_njet.eps");
	overLayBackgrounds("factorization_ht800/pf30_jet1_pt", "HT>800 GeV: Events passing RA2 Selection Except #Delta#Phi;Jet1 P_{T} [GeV];Events;", 4, "ht800_jet1pt.eps");
	overLayBackgrounds("factorization_ht800/pf30_jet2_pt", "HT>800 GeV: Events passing RA2 Selection Except #Delta#Phi;Jet2 P_{T} [GeV];Events;", 4, "ht800_jet2pt.eps");
	overLayBackgrounds("factorization_ht800/pf30_jet3_pt", "HT>800 GeV: Events passing RA2 Selection Except #Delta#Phi;Jet3 P_{T} [GeV];Events;", 4, "ht800_jet3pt.eps");
*/

/*	overLayBackgrounds("factorization_ht1000/mht", "HT>1000 GeV: Events passing RA2 Selection Except #Delta#Phi;MHT [GeV];Events;", 8, "ht1000_mht.eps");
	overLayBackgrounds("factorization_ht1000/ht", "HT>1000 GeV: Events passing RA2 Selection Except #Delta#Phi; HT [GeV];Events;", 1, "ht1000_ht.eps");
	overLayBackgrounds("factorization_ht1000/njet_et50eta24", "HT>1000 GeV: Events passing RA2 Selection Except #Delta#Phi;NJet [E_{T}>50 GeV & | #eta |<2.5];Events;", 1, "ht1000_njet.eps");
	overLayBackgrounds("factorization_ht1000/pf30_jet1_pt", "HT>1000 GeV: Events passing RA2 Selection Except #Delta#Phi;Jet1 P_{T} [GeV];Events;", 4, "ht1000_jet1pt.eps");
	overLayBackgrounds("factorization_ht1000/pf30_jet2_pt", "HT>1000 GeV: Events passing RA2 Selection Except #Delta#Phi;Jet2 P_{T} [GeV];Events;", 4, "ht1000_jet2pt.eps");
	overLayBackgrounds("factorization_ht1000/pf30_jet3_pt", "HT>1000 GeV: Events passing RA2 Selection Except #Delta#Phi;Jet3 P_{T} [GeV];Events;", 4, "ht1000_jet3pt.eps");

	overLayBackgrounds("factorization_ht1200/mht", "HT>1200 GeV: Events passing RA2 Selection Except #Delta#Phi;MHT [GeV];Events;", 8, "ht1200_mht.eps");
	overLayBackgrounds("factorization_ht1200/ht", "HT>1200 GeV: Events passing RA2 Selection Except #Delta#Phi; HT [GeV];Events;", 1, "ht1200_ht.eps");
	overLayBackgrounds("factorization_ht1200/njet_et50eta24", "HT>1200 GeV: Events passing RA2 Selection Except #Delta#Phi;NJet [E_{T}>50 GeV & | #eta |<2.5];Events;", 1, "ht1200_njet.eps");
	overLayBackgrounds("factorization_ht1200/pf30_jet1_pt", "HT>1200 GeV: Events passing RA2 Selection Except #Delta#Phi;Jet1 P_{T} [GeV];Events;", 4, "ht1200_jet1pt.eps");
	overLayBackgrounds("factorization_ht1200/pf30_jet2_pt", "HT>1200 GeV: Events passing RA2 Selection Except #Delta#Phi;Jet2 P_{T} [GeV];Events;", 4, "ht1200_jet2pt.eps");
	overLayBackgrounds("factorization_ht1200/pf30_jet3_pt", "HT>1200 GeV: Events passing RA2 Selection Except #Delta#Phi;Jet3 P_{T} [GeV];Events;", 4, "ht1200_jet3pt.eps");
*/
}

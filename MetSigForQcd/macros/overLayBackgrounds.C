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
//const static float fEWK_LUMI   = 81352581.0/31300.0;        // 31300 pb for NLO xsec and total number of events is 81352581 (sample1+sample2)
//const static float fTTBAR_LUMI = 3701947.0/165.0;
const static float fZNN_LUMI   = 3067017.0/42.8;              //x-sec from note, Anwar hd used 32.92 pb-1 LO
//const static string sDATA_FILE_NAME = "Data_MyPrescaleWgted.root"; //used in AN 02-08-2011
//PU wgted samples
const static string sDATA_FILE_NAME  = "Data_465_PrescaleWgted_DphiBugFixed_03012012.root"; 
const static string sWJET_HT250TO300_FILE_NAME  = "Wjets1.root";  //Fall 11 sample
const static string sWJET_HT300TOINF_FILE_NAME  = "Wjets2.root";  //Fall 11 sample
const static string sTTBAR_FILE_NAME = "TTbar.root";  // Fall 11 sample
const static string sZNN_FILE_NAME   = "Znn.root";    // Summer 11 sample
//const static string sLM9_FILE_NAME   = "SUSYLM9.root";    // Summer 11 sample
const static string sLM9_FILE_NAME   = "SUSYLM5.root";    // Summer 11 sample
const static float fWJETS_HT250TO300_LUMI =  9831277.0 / (248.7 * 0.14); 
const static float fWJETS_HT300TOINF_LUMI =  5363746.0 / (317.0 * 0.153); 
const static float fTTBAR_LUMI = 59590147/165.0;
//const static float fSUSY_LUMI = 437030/(7.134 	1.48);  //LM9 (seema's)
const static float fSUSY_LUMI = 433070/(0.473 * 1.34); //LM5 (seema's) NLO

std::vector<float> incMHTBins;
const static int nMHTbins = 5;
const static float arrMHTbins[nMHTbins] = {200,350,500,600,7000};
/*to do:
 * include susy file and lumi
 *
 *
 */

struct Predictions_t
{
		float incl_mean[nMHTbins];
		float incl_statErr[nMHTbins];
		float incl_fitErr[nMHTbins];
		float excl_mean[nMHTbins-1];
		float excl_statErr[nMHTbins-1];
		float excl_fitErr[nMHTbins-1];
};

/*******************
 * Find events in each exclusive bin for a given hist
 *******************/
Predictions_t GetExlcusiveBinContents(const TH1* hist)
{
	assert(hist != NULL && "GetExlcusiveBinContents:: hist not found!");

	double mean[incMHTBins.size()-1];
	double statErr[incMHTBins.size()-1];

	Predictions_t results;

	for (unsigned mhtBin = 0; mhtBin < incMHTBins.size(); ++mhtBin)
	{
		if (mhtBin+1<incMHTBins.size()) 
		{
			mean[mhtBin] = 0;
			statErr[mhtBin] = 0;
		}

		for (int bin = 0; bin <= hist->GetNbinsX()+1; ++bin)
		{
			if (hist->GetBinContent(bin)>0)
			{

				const float binCenter = hist->GetBinCenter(bin);
				const float binVal    = hist->GetBinContent(bin);
				const float binErr    = hist->GetBinError(bin);
				const float statErr2  = pow(binErr, 2);

				if (mhtBin+1<incMHTBins.size())
				{
					if (binCenter > incMHTBins.at(mhtBin) 
							&& binCenter < incMHTBins.at(mhtBin+1))
					{

					/*	if (mhtBin ==1)
						{
							cout << std::setprecision(4) << std::setw(5) << bin << std::setw(3) << "[" 
								<< std::setw(10) << hist->GetBinLowEdge(bin) << ", "  << std::setw(10) 
								<< hist->GetXaxis()->GetBinUpEdge(bin)<< "]" 
								<< std::setw(10) << hist->GetBinContent(bin) 
								<< std::setw(10) << hist->GetBinError(bin) 
								<< endl;
						}
						*/
						mean[mhtBin]    += binVal;
						statErr[mhtBin] += statErr2;
					}
				}
			}
		}

		statErr[mhtBin] = sqrt(statErr[mhtBin]);
		//if (mhtBin ==1) { std::cout << "final stat error = " << statErr[mhtBin] << std::endl; }

		results.incl_mean[mhtBin]    = -1;   //to show they are not used or valid for this
		results.incl_statErr[mhtBin] = -1;
		results.incl_fitErr[mhtBin]  = -1;

		
		if (mhtBin+1<incMHTBins.size())
		{
			std::stringstream excl_pred;
			excl_pred << setprecision(3) 
				<< setw(10) << incMHTBins.at(mhtBin) << "<MHT<" << incMHTBins.at(mhtBin+1)
				<< mean[mhtBin] 
				//<< "&$\\pm$" 
				<< " +/- " 
				<< statErr[mhtBin];
			//std::cout << excl_pred.str() << std::endl;
			results.excl_mean[mhtBin]    = mean[mhtBin];
			results.excl_statErr[mhtBin] = statErr[mhtBin];
			results.excl_fitErr[mhtBin]  = -1;
		}
	}
	
	return results;
}

void PrintExclPredictions(const Predictions_t& res, const bool header=true)
{

	if (header)
	{
		cout << setprecision(3) << setw(15) << " MHT "
			<< setw(10) << " mean " 
			<< setw(10) << "stat"
			<< endl;
	}


	for (int i=0; i< nMHTbins -1 ; ++i)
	{
		//cout << "i = " << i  << " [ " << nMHTbins << "]" << endl;
		const float mean    = res.excl_mean[i]; 
		const float statErr = res.excl_statErr[i];

		cout << setprecision(3) << setw(5) << incMHTBins.at(i) << "<MHT<" << incMHTBins.at(i+1) 
					//<< setw(10) << mean << " $\\pm$ " 
					<< setw(10) << mean << " +/- " 
					<< setw(10) << statErr
					<< endl;
	}
}

void PrintResults(const Predictions_t& sig, const Predictions_t& bg, const bool header=true)
{

	if (header)
	{
		cout << setprecision(2) << setw(18) << " MHT "
			<< setw(18) << " signal +/- stat"  
			<< setw(18) << " bg +/- stat"
			//<< setw(10) << " s/b"
			//<< setw(10) << " s/sqrt(b)"
			<< setw(10) << " s/sqrt(s+b)"
			<< endl;
	}


	for (int i=0; i< nMHTbins -1 ; ++i)
	{
		//cout << "i = " << i  << " [ " << nMHTbins << "]" << endl;
		const float sigmean    = sig.excl_mean[i]; 
		const float sigstatErr = sig.excl_statErr[i];
		const float bgmean     = bg.excl_mean[i]; 
		const float bgstatErr  = bg.excl_statErr[i];

		float soverb = 0;
		float soversqrtb = 0;
		if (bgmean==0 && sigmean>0) 
		{
			soverb = 99999.9;
			soversqrtb = 99999.9;
		} else if (bgmean != 0) 
		{
			soverb = sigmean/bgmean;
			soversqrtb = sigmean/sqrt(bgmean);
		}
		
		float soversqrtsb = 0;
		if ( bgmean+sigmean != 0) soversqrtsb = sigmean/sqrt(sigmean+bgmean); 

		cout << setprecision(3) << setw(5) << incMHTBins.at(i) << "<MHT<" << incMHTBins.at(i+1) 
					//<< setw(10) << mean << " $\\pm$ " 
					<< setw(12) << sigmean << " +/- " << sigstatErr 
					<< setw(12) << bgmean << " +/- " << bgstatErr 
					//<< setw(10) << soverb 
					//<< setw(10) << soversqrtb 
					<< setw(10) << soversqrtsb 
					<< endl;
	}
}


void PrintBGbreakdown(const Predictions_t& qcd, 
					const Predictions_t& wjet, 
					const Predictions_t& zinv, 
					const Predictions_t& ttbar, 
					const bool header=true)
{

	if (header)
	{
		cout << setprecision(2) << setw(18) << " MHT "
			<< setw(10) << "qcd "  
			<< setw(10) << "wjet"
			<< setw(10) << "zinv"
			<< setw(10) << "ttbar"
			<< endl;
	}


	for (int i=0; i< nMHTbins -1 ; ++i)
	{
		//cout << "i = " << i  << " [ " << nMHTbins << "]" << endl;

		cout << setprecision(3) << setw(5) << incMHTBins.at(i) << "<MHT<" << incMHTBins.at(i+1) 
					<< setw(10) << qcd.excl_mean[i] 
					<< setw(10) << wjet.excl_mean[i] 
					<< setw(10) << zinv.excl_mean[i] 
					<< setw(10) << ttbar.excl_mean[i] 
					<< endl;
	}
}


void FindHistsYmax(const vector<TH1*> vHists, float& ymax, float& ymin)
{
	ymax = 0;
	ymin = 999999999.0;
	for (unsigned i=0; i< vHists.size(); ++i)
	{
		const float maxbinVal = vHists.at(i)->GetBinContent(vHists.at(i)->GetMaximumBin());
		ymax = std::max(ymax, maxbinVal);
		
		for (unsigned bin =1; bin < vHists.at(i)->GetNbinsX(); ++bin)
		{
			const float binVal = vHists.at(i)->GetBinContent(bin);
			if (binVal==0) continue;
			ymin = std::min(ymin, binVal);
		}
	}
}

TH1* GetZnnHist(const string histname)
{
	TFile* znnFile = new TFile(sZNN_FILE_NAME.c_str());
	if (znnFile->IsZombie()) { cout << "TTbar file not found!" << endl; assert (false); }

	TH1* hist = dynamic_cast<TH1*> (znnFile->Get(histname.c_str()));
	if (hist == NULL) { cout << "Znn hist " << histname << " not found!" << endl; assert (false); }
	TH1* znnHist = dynamic_cast<TH1*> (hist->Clone());
   znnHist->SetDirectory(0);
	znnHist->Sumw2();
	znnHist->SetLineColor(kRed);

	const float int_b4 = znnHist->Integral(); 
	const float scale  = fDATA_LUMI/fZNN_LUMI;
	znnHist->Scale(scale);
//	std::cout << __FUNCTION__ << ": scale = " << scale<< ":: Znn Integral b4/a4= "<< int_b4  << "/" << znnHist->Integral() << std::endl; 

/*	new TCanvas();
	znnHist->SetLineColor(kGreen);
	znnHist->DrawCopy();
	hist->DrawCopy("same");
*/
	delete znnFile;
	return znnHist;
}

TH1* GetSUSYHist(const string histname)
{
	TFile* susyFile = new TFile(sLM9_FILE_NAME.c_str());
	if (susyFile->IsZombie()) { cout << "SUSY file not found!" << endl; assert (false); }

	TH1* hist = dynamic_cast<TH1*> (susyFile->Get(histname.c_str()));
	if (hist == NULL) { cout << "SUSY hist " << histname << " not found!" << endl; assert (false); }
	TH1* susyHist = dynamic_cast<TH1*> (hist->Clone());
   susyHist->SetDirectory(0);
	susyHist->Sumw2();

	const float int_b4 = susyHist->Integral(); 
	const float scale  = fDATA_LUMI/fSUSY_LUMI;
	susyHist->Scale(scale);
//	std::cout << __FUNCTION__ << ": scale = " << scale<< ":: SUSY Integral b4/a4= "<< int_b4  << "/" << susyHist->Integral() << std::endl; 

/*	new TCanvas();
	susyHist->SetLineColor(kGreen);
	susyHist->DrawCopy();
	hist->DrawCopy("same");
*/
	delete susyFile;
	return susyHist;
}


TH1* GetTTbarHist(const string histname)
{
	TFile* ttbarFile = new TFile(sTTBAR_FILE_NAME.c_str());
	if (ttbarFile->IsZombie()) { cout << "TTbar file not found!" << endl; assert (false); }

	TH1* hist = dynamic_cast<TH1*> (ttbarFile->Get(histname.c_str()));
	if (hist == NULL) { cout << "TTbar hist " << histname << " not found!" << endl; assert (false); }
	TH1* ttbarHist = dynamic_cast<TH1*> (hist->Clone());
	ttbarHist->Sumw2();
   ttbarHist->SetDirectory(0);
	ttbarHist->SetLineColor(kRed);

	const float int_b4 = ttbarHist->Integral(); 
	const float scale  = fDATA_LUMI/fTTBAR_LUMI;
	ttbarHist->Scale(scale);
//	std::cout << __FUNCTION__ << ": scale = " << scale<< ":: TTbar Integral b4/a4= "<< int_b4  << "/" << ttbarHist->Integral() << std::endl; 

/*	new TCanvas();
	ttbarHist->SetLineColor(kGreen);
	ttbarHist->DrawCopy();
	hist->DrawCopy("same");
*/
	delete ttbarFile;
	return ttbarHist;
}


TH1* GetEWKHist(const string ewkhistname)
{
	TFile* ewkFile1 = new TFile(sWJET_HT250TO300_FILE_NAME.c_str());
	TFile* ewkFile2 = new TFile(sWJET_HT300TOINF_FILE_NAME.c_str());
	if (ewkFile1->IsZombie()) { cout << "EWK file 1 not found!" << endl; assert (false); }
	if (ewkFile2->IsZombie()) { cout << "EWK file 2 not found!" << endl; assert (false); }

	TH1* hist1 = dynamic_cast<TH1*> (ewkFile1->Get(ewkhistname.c_str()));
	if (hist1 == NULL) { cout << "EWK hist " << ewkhistname << " not found in ewkFile 1!" << endl; assert (false); }
	TH1* ewkHist1 = dynamic_cast<TH1*> (hist1->Clone());
	ewkHist1->Sumw2();
	TH1* hist2 = dynamic_cast<TH1*> (ewkFile2->Get(ewkhistname.c_str()));
	if (hist2 == NULL) { cout << "EWK hist " << ewkhistname << " not found in ewkFile 2!" << endl; assert (false); }
	TH1* ewkHist2 = dynamic_cast<TH1*> (hist2->Clone());
	ewkHist2->Sumw2();
	//new TCanvas();
	//ewkHist->DrawCopy();

   ewkHist1->SetDirectory(0);
   ewkHist2->SetDirectory(0);
	const float int_b4 = ewkHist1->Integral(); 
	const float scale1  = fDATA_LUMI/fWJETS_HT250TO300_LUMI;
	const float scale2  = fDATA_LUMI/fWJETS_HT300TOINF_LUMI;

	ewkHist1->Scale(scale1);
	ewkHist2->Scale(scale2);
	ewkHist1->Add(ewkHist2);
	//std::cout << __FUNCTION__ << ": scale1/scale2 = " << scale1 << "/" << scale2 << ":: EWK Integral b4/a4= "<< int_b4  << "/" << ewkHist1->Integral() << std::endl; 

	delete ewkFile1;
	delete ewkFile2;
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
		TH1* hist = dynamic_cast<TH1*> (files[i]->Get(histname.c_str()));
		if (hist == NULL)
		{
			cout << histname << " not found in file " << files[i]->GetName() << endl;
			assert(false);
		}
		hists[i] = dynamic_cast<TH1*> (hist->Clone());
   	hists[i]->SetDirectory(0);
		//new TCanvas();
		//hists[i]->Draw();
		if (hists[i] == 0)
		{
			cout << "hist_pass " << histname << " not found in " << files[i]->GetName() << "!" << endl;
			assert (false);
		} else 
		{
			hists[i]->Sumw2();
			const float scale = fDATA_LUMI/ ( nEvents[i] / xSec[i] );
			//cout << "hist " << i << " scale = " << scale << endl;
			if (hists[i]->GetEntries() > 1) hists[i]->Scale(scale);

			if (i == 0) res_hist = dynamic_cast<TH1*> (hists[i]->Clone("histcopy"));
			else res_hist->Add(hists[i]);
		}
	}

//	for (int i=0; i<nBins; ++i) delete files[i];


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

void overLayBackgrounds(const string histname, const string title, const int rebin, 
			const string epsname, const bool logScale, const bool drawStack=1) 
{

	//to suppress ROOT runtime warnings (like Sumw2 already created )
	gErrorIgnoreLevel=kError;

	incMHTBins.clear();
	for (int i=0; i<nMHTbins; ++i) incMHTBins.push_back(arrMHTbins[i]);


//	TH1 *histDataOrig = 0;

/*	TFile* dataRootFile = new TFile (sDATA_FILE_NAME.c_str());
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
*/
	/********************************************/
	//substract background contaminations
	/********************************************/
	TH1* ewkHist   = GetEWKHist  (histname);
	TH1* ttbarHist = GetTTbarHist(histname);
	TH1* znnHist   = GetZnnHist  (histname);
	TH1* qcdHist   = GetQCDHist(histname);
	TH1* susyHist  = GetSUSYHist(histname);

	if (rebin>1)
	{
//		histDataOrig->Rebin(rebin);
//		histDataFinal->Rebin(rebin);
		ewkHist->Rebin(rebin);
		ttbarHist->Rebin(rebin);
		znnHist->Rebin(rebin);
		qcdHist->Rebin(rebin);
		susyHist->Rebin(rebin);
	}


/*	histDataFinal->Add(ewkHist,-1);
	histDataFinal->Add(ttbarHist,-1);
	histDataFinal->Add(znnHist,-1);
*/
	const int lineWidth = 2;
	const int ewkColor = 3;
	const int znnColor = 6;
	const int ttbarColor = 4;
	const int qcdColor = kBlack;
	const int susyColor = kRed;

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
	susyHist->SetLineColor(susyColor);
	susyHist->SetFillColor(10);
	susyHist->SetLineWidth(lineWidth);
	susyHist->SetMarkerColor(susyColor);
	susyHist->SetMarkerStyle(23);


	/*histDataFinal->SetLineColor(2);
	histDataFinal->SetLineWidth(lineWidth);
	histDataFinal->SetTitle(title.c_str());
	histDataFinal->SetMinimum(0.5e-1);
*/

	//TLegend *leg1  = new TLegend(0.55,0.6,0.9,0.9);
/*	TLegend *leg1  = new TLegend(0.8,0.8,0.99,0.99);
//	leg1->AddEntry(histDataOrig,"DATA (before substraction)");
//	leg1->AddEntry(histDataFinal,"DATA (after substraction)");
	leg1->AddEntry(ewkHist,"W+Jets");
	leg1->AddEntry(znnHist,"Z#rightarrow#nu#bar{#nu}");
	leg1->AddEntry(ttbarHist,"t#bar{t}");
	leg1->AddEntry(susyHist,"SUSY LM ");
*/
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
	totBgHist->SetLineColor(30);
	totBgHist->SetMarkerColor(30);
	totBgHist->SetFillColor(10);

	TLegend *leg2  = new TLegend(0.65,0.65,0.9,0.9);
	//TLegend *leg2  = new TLegend(0.8,0.8,0.99,0.99);
	leg2->SetTextFont(42);
//	leg2->AddEntry(histDataOrig,"DATA");
//	leg2->AddEntry(totBgHist,"Background");
	leg2->AddEntry(qcdHist,"QCD");
	leg2->AddEntry(ewkHist,"W+Jets");
	leg2->AddEntry(znnHist,"Z#rightarrow#nu#bar{#nu}");
	leg2->AddEntry(ttbarHist,"t#bar{t}");
	if (! drawStack)  leg2->AddEntry(totBgHist,"Sum BGs");
	leg2->AddEntry(susyHist,"SUSY LM5");

	

	vector<TH1*> vHists;
	vHists.push_back(znnHist);
	vHists.push_back(qcdHist);
	vHists.push_back(ewkHist);
	vHists.push_back(ttbarHist);
	vHists.push_back(susyHist);
	//keeping this to get the y-scale correctly in stacked hist
	vHists.push_back(totBgHist);

	float ymax=0, ymin=0;
	FindHistsYmax(vHists, ymax, ymin);
	//cout << "Y min/max = " << ymax << "/ " << ymin << endl;

	new TCanvas();
	if (logScale) gPad->SetLogy();

	if (drawStack)
	{
		THStack *hs = new THStack ("hs", NULL);
		//hs->GetXaxis()->SetTitleFont(42);
		hs->SetMinimum(0.5 * ymin);
		hs->SetMaximum(1.1 * ymax);
		hs->Add(znnHist);
		hs->Add(ewkHist);
		hs->Add(ttbarHist);
		hs->Add(qcdHist);
		hs->SetTitle(title.c_str());
		hs->Draw("HIST");
	} else
	{
		totBgHist->SetMinimum(0);
		totBgHist->SetMaximum(1.2 * ymax);
		totBgHist->SetTitle(title.c_str());
		totBgHist->SetStats(0);
		totBgHist->Draw();
		znnHist->Draw("Lsame");
		ewkHist->Draw("Lsame");
		qcdHist->Draw("Lsame");
		ttbarHist->Draw("Lsame");
	}

//	histDataOrig->Draw("same");
	susyHist->Draw("sameL");
	leg2->Draw();

	gPad->Print(epsname.c_str());

	//PREDICTION FOR EXLUSIVE BINS
	Predictions_t signal = GetExlcusiveBinContents(susyHist);
	Predictions_t bg     = GetExlcusiveBinContents(totBgHist);
	Predictions_t qcd    = GetExlcusiveBinContents(qcdHist);
	Predictions_t ewk    = GetExlcusiveBinContents(ewkHist);
	Predictions_t znn    = GetExlcusiveBinContents(znnHist);
	Predictions_t ttbar    = GetExlcusiveBinContents(ttbarHist);


//	PrintExclPredictions(signal);
//	PrintExclPredictions(bg);

	cout << "====== Results for " << epsname << endl;
	PrintResults(signal, bg);
	cout << ">>>>>>>>>> QCD Predictions " << epsname << endl;
	PrintExclPredictions(qcd);
	/*cout << ">>>>>>>>>> EWK Predictions " << epsname << endl;
	PrintExclPredictions(ewk);
	cout << ">>>>>>>>>> ZNN Predictions " << epsname << endl;
	PrintExclPredictions(znn);
	cout << ">>>>>>>>>> TT  Predictions " << epsname << endl;
	PrintExclPredictions(ttbar);

	cout << "====== Breakdown of backgrounds!!!" << endl;
	PrintBGbreakdown(qcd, ewk, znn,ttbar); 
*/

}

void overLayBackgrounds(const std::string histname, const std::string title, 
				const int rebin, const bool logScale, const bool stackHists)
{
	for (int htBin=0; htBin<6; htBin++)
	{
		//if ( htBin != 1 ) continue;
		string htBinLabel("");
		string htRange("");
		if (htBin==0) { htBinLabel += "HT350to500"; htRange += "350<HT<500";}
		else if (htBin==1) { htBinLabel += "HT500to800"; htRange += "500<HT<800";}
		else if (htBin==2) { htBinLabel += "HT800to1000"; htRange += "800<HT<1000";}
		else if (htBin==3) { htBinLabel += "HT1000to1200"; htRange += "1000<HT<1200";}
		else if (htBin==4) { htBinLabel += "HT1200to1400"; htRange += "1200<HT<1400";}
		else if (htBin==5) { htBinLabel += "HT1400to7000"; htRange += "HT>1400";}

		for (int mhtBin=0; mhtBin <3; ++mhtBin)
		{
			if ( mhtBin != 0 ) continue;
			string mhtBinLabel("");
			string mhtRange("");
			if (mhtBin ==0) { mhtBinLabel += "nomht"; mhtRange += "MHT>0";}
			else if (mhtBin ==1) { mhtBinLabel += "mht200"; mhtRange += "MHT>200";}
			else if (mhtBin ==2) { mhtBinLabel += "mht350"; mhtRange += "MHT>350";}
			else if (mhtBin ==3) { mhtBinLabel += "mht500"; mhtRange += "MHT>500";}

			stringstream path, epsname, histtitle, additional;
			path << "metsig" << mhtBinLabel << "/" << htBinLabel << "/" << histname;
			epsname << htBinLabel << "_" << mhtBinLabel << "_" << histname << ".eps";
			//histtitle << htRange << ", " << mhtRange << ";" << title;
			//histtitle << "#geq 2 Jets: "<< htRange << ", " << mhtRange << ";" << title;
			//histtitle << "#geq 3 Jets and #Delta#Phi cut: "<< htRange << ", " << mhtRange << ";" << title;
			histtitle << "#geq 3 Jets and #Delta#Phi cut and MetSig>300: "<< htRange << ", " << mhtRange << ";" << title;

			overLayBackgrounds(path.str() , histtitle.str(), rebin, epsname.str(), logScale, stackHists);
		}
	}
}

void overLayBackgrounds(const int i)
{
	//to get predictions
	overLayBackgrounds("mht", "MHT;Events;", 1,1,1);
	return;


	if (i==0) overLayBackgrounds("metsig", "METSig;Events;", 2,0,1);
	//return;
	if (i==1) overLayBackgrounds("metsigprob", "METSig Prob;Events;", 4,0,0);
	if (i==2) overLayBackgrounds("ht", "HT;Events;", 10,1,1);
	if (i==3) overLayBackgrounds("mht", "MHT;Events;", 10,1,1);
	

	if (i==4) overLayBackgrounds("jet1_pt", "Jet-1 Pt;Events;", 4,0,1);
	if (i==5) overLayBackgrounds("jet1_dphi", "Jet-1 #Delta#Phi;Events;", 1,0,1);
	if (i==6) overLayBackgrounds("jet2_pt", "Jet-2 Pt;Events;", 4,1,1);
	if (i==7) overLayBackgrounds("njet50", "NJet50;Events;", 1,0,1);
	if (i==8) overLayBackgrounds("metsigVsMHT", "METSig;MHT (GeV);", 5,0,0);
	if (i==10) overLayBackgrounds("metsigVsHT", "METSig;HT (GeV);", 5,0,0);
}

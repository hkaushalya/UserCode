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

using namespace std;

const static float fDATA_LUMI  = 4650.0; //pb-1
//const static float fTTBAR_LUMI = 3701947.0/165.0;
//const static float fZNN_LUMI   = 3067017.0/42.8;     //x-sec from note, Anwar hd used 32.92 pb-1 LO
																	  //2-23-12, Seema said this should be 42.8 not 42.2
//const static string sDATA_FILE_NAME  = "Data.root";
//const static string sDATA_FILE_NAME  = "Data_SmearedJetsNoPrescaleWgts.root";
//const static string sDATA_FILE_NAME  = "Data_PreScaleWgted.root";
//const static string sDATA_FILE_NAME  = "Data_1stfb_12052011.root";
//const static string sWJET_FILE_NAME  = "EWKmc.root";
//const static string sDATA_FILE_NAME = "SmearedJetsWithPrescaleWgts.root";
//const static string sDATA_FILE_NAME = "Data_MyPrescaleWgted.root"; //used in AN 02-08-2011
//const static string sDATA_FILE_NAME = "Data_SmearedRecoJetsWithPrescaleWgts.root";
//const static string sWJET_FILE_NAME  = "WJetsMC.root";
//const static string sTTBAR_FILE_NAME = "TTbarMC.root";
//const static string sZNN_FILE_NAME   = "ZnnMC.root";
//PU weighted, MHT binned, Fall 11 samples, dPhi calculated using all jets Et>30 GeV && |eta|<5.0
const static float fWJETS_HT250TO300_LUMI =  9831277.0 / (248.7 * 0.14); 
const static float fWJETS_HT300TOINF_LUMI =  5363746.0 / (317.0 * 0.153); 
const static float fTTBAR_LUMI = 59590147/165.0;
const static float fZNN_LUMI   = 3067017.0/42.8;     //x-sec from note, Anwar hd used 32.92 pb-1 LO
const static string sDATA_FILE_NAME  = "Data_465_PrescaleWgted_03162012.root"; 
const static string sWJET_HT250TO300_FILE_NAME  = "Wjets_HT250to300_Fall11_PUwgted_03162012.root";  //Fall 11 sample
const static string sWJET_HT300TOINF_FILE_NAME  = "Wjets_HT300toInf_Fall11_PUwgted_03162012.root";  //Fall 11 sample
const static string sTTBAR_FILE_NAME = "TTbar_Fall11_PUwgted_03162012.root";  // Fall 11 sample
const static string sZNN_FILE_NAME   = "Zjets_Summ11_PUwgted_03162012.root";    // Summer 11 sample


const static float passFailRatio_fitrange_xmin = 50.0;
const static float passFailRatio_fitrange_xmax = 150.0;
std::vector<float> incMHTBins;
const static int nMHTbins = 5;
const static float arrMHTbins[nMHTbins] = {200,350,500,600,7000};
float CONST_C      = 0.0217;
//float CONST_C_UP   = 0.07;  //2010 values
//float CONST_C_DOWN = 0.012; //2010 values
//float CONST_C_UP   = 0.0304145;
//oat CONST_C_DOWN = CONST_C_UP - CONST_C;

//ASSIGN 100% ERROR see TWIKI answer we provided (Jan 26th, 2012)
const float CONST_C_UP 	 = CONST_C + CONST_C;
const float CONST_C_DOWN = CONST_C - CONST_C;


struct Predictions_t
{
		float incl_mean[nMHTbins];
		float incl_statErr[nMHTbins];
		float incl_fitErr[nMHTbins];
		float excl_mean[nMHTbins-1];
		float excl_statErr[nMHTbins-1];
		float excl_fitErr[nMHTbins-1];
};

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
			const float mean    = res.excl_mean[i]; 
			const float statErr = res.excl_statErr[i];

		cout << setprecision(3) << setw(15) << incMHTBins.at(i) << "<MHT<" << incMHTBins.at(i+1) 
					//<< setw(10) << mean << " $\\pm$ " 
					<< setw(10) << mean << " +/- " 
					<< setw(10) << statErr
					<< endl;
	}
}

void PrintBackgrounds(
					const Predictions_t& data_b4, 
					const Predictions_t& data_a4, 
					const Predictions_t& ewk, 
					const Predictions_t& ttbar,
					const Predictions_t& znn
					)
{
	cout << ">>>>>> Backgrounds in MHT bins <<<<<<<<" << std::endl;
			cout 
			<< setw(5) << " MHT "
			<< setw(20) << "data+/-stat (b4)" 
			<< setw(20) << "wjet+/-stat" 
			<< setw(20) << "ttbar/-stat" 
			<< setw(20) << "znn +/-stat" 
			<< setw(20) << "data+/-stat (a4)" 
			<< endl;
	//const string seperator("$\\pm$");
	const string seperator(" +/- ");
	//const string elogOption("");
	const string elogOption(" | ");
	for (int i=0; i< nMHTbins -1 ; ++i)
	{
		cout << setprecision(0) 
			<< setw(5) << incMHTBins.at(i) << "<MHT<" << incMHTBins.at(i+1) 
			<< setprecision(1) 
			<< setw(20) << data_b4.excl_mean[i] << seperator << data_b4.excl_statErr[i] << elogOption
			<< setw(15) << ewk.excl_mean[i]     << seperator << ewk.excl_statErr[i] << elogOption
			<< setw(15) << ttbar.excl_mean[i]   << seperator << ttbar.excl_statErr[i] << elogOption
			<< setw(15) << znn.excl_mean[i]     << seperator << znn.excl_statErr[i] << elogOption
			<< setw(20) << data_a4.excl_mean[i] << seperator << data_a4.excl_statErr[i] << elogOption
			<< endl;
	}

}



void PrintResults(const Predictions_t& res, 
	//				const Predictions_t& htSystRes,
					const Predictions_t& qcdSystRes,
					const Predictions_t& C_upSyst,
					const Predictions_t& C_downSyst
					)
{

	cout.precision(2);
	cout << setw(10) << " MHT "
					<< setw(20) << " mean " 
					<< setw(15) << ":statAndFitErr"
	//				<< setw(10) << " htSyst"
					<< setw(10) << " qcdSyst1"
					<< setw(10) << " const C"
					//<< setw(15) << " const C up"
					//<< setw(15) << " const C down"
					<< setw(10) << "total"
					//<< setw(10) << "total up"
					//<< setw(10) << "total down"
					<< endl;
	for (int i=0; i< nMHTbins -1 ; ++i)
	{
			const float mean          = res.excl_mean[i]; 
			const float statAndFitErr = sqrt(pow(res.excl_statErr[i],2)+pow(res.excl_fitErr[i],2));
	//		const float htSys         = fabs(mean - htSystRes.excl_mean[i]);
			const float qcdSys1       = fabs(mean - qcdSystRes.excl_mean[i]);
			const float C_Sys         = max(fabs(mean - C_upSyst.excl_mean[i]), fabs(mean - C_downSyst.excl_mean[i]));
			const float C_Sys_up      = fabs(mean - C_upSyst.excl_mean[i]);
			const float C_Sys_down    = fabs(mean - C_downSyst.excl_mean[i]);

			const float total         = sqrt(  pow(statAndFitErr,2)
	//		                                 + pow(htSys,2)
														+ pow(qcdSys1,2) 
														+ pow(C_Sys,2) 
														);

			const float total_up      = sqrt(  pow(statAndFitErr,2)
	//		                                 + pow(htSys,2)
														+ pow(qcdSys1,2) 
														+ pow(C_Sys_up,2) 
														);

			const float total_down      = sqrt(  pow(statAndFitErr,2)
	//		                                 + pow(htSys,2)
														+ pow(qcdSys1,2) 
														+ pow(C_Sys_down,2) 
														);
		cout << fixed << setw(5) << incMHTBins.at(i) << "<MHT<" << setw(7) << incMHTBins.at(i+1)
					<< " & "
					<< setw(10) << mean 
					<< setw(15) << " & $\\pm$ " 
					<< statAndFitErr
	//				<< setw(15) << " & $\\pm$ " << htSys
					<< setw(15) << " & $\\pm$ " << qcdSys1
					//<< setw(10) << " & " << "$^{+" << C_Sys_up << "}_{-" << C_Sys_down << "}$"
					<< setw(10) << " & $\\pm$ " << C_Sys_up
					//<< setw(10) << " & " << "$^{+" << total_up << "}_{-" << total_down << "}$" << "   \\\\"
					<< setw(10) << " & $\\pm$ " << total_up << "   \\\\"
					<< endl;
	}

	for (int i=0; i< nMHTbins ; ++i)
	{
			const float mean          = res.incl_mean[i]; 
			const float statAndFitErr = sqrt(pow(res.incl_statErr[i],2)+pow(res.incl_fitErr[i],2));
	//		const float htSys         = fabs(mean - htSystRes.incl_mean[i]);
			const float qcdSys1       = fabs(mean - qcdSystRes.incl_mean[i]);
			const float C_Sys         = max(fabs(mean - C_upSyst.incl_mean[i]), fabs(mean - C_downSyst.incl_mean[i]));
			const float C_Sys_up      = fabs(mean - C_upSyst.incl_mean[i]);
			const float C_Sys_down    = fabs(mean - C_downSyst.incl_mean[i]);

			const float total         = sqrt(  pow(statAndFitErr,2)
	//		                                 + pow(htSys,2)
														+ pow(qcdSys1,2) 
														+ pow(C_Sys,2) 
														);

			const float total_up      = sqrt(  pow(statAndFitErr,2)
	//		                                 + pow(htSys,2)
														+ pow(qcdSys1,2) 
														+ pow(C_Sys_up,2) 
														);

			const float total_down    = sqrt(  pow(statAndFitErr,2)
	//		                                 + pow(htSys,2)
														+ pow(qcdSys1,2) 
														+ pow(C_Sys_down,2) 
														);

		cout << fixed << setw(5) << "MHT> " << setw(7) << incMHTBins.at(i) 
					<< setw(10) << " & " << mean 
					<< setw(15) << " & $\\pm$ " 
					<< setw(2) << statAndFitErr
	//				<< setw(10) << " & $\\pm$ " << htSys
					<< setw(10) << " & $\\pm$ " << qcdSys1
					//<< setw(10) << " & " << "$^{+" << C_Sys_up << "}_{-" << C_Sys_down << "}$"
					<< setw(10) << " & $\\pm$ " << C_Sys_up
					//<< setw(10) << " & " << "$^{+" << total_up << "}_{-" << total_down << "}$" << "   \\\\"
					<< setw(10) << " & $\\pm$ " << total_up << "   \\\\"
					<< endl;
	}

}


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

void ScaleTTbarByLumi(TH1* hist)
{
	hist->Scale(fDATA_LUMI/fTTBAR_LUMI);	
}
void ScaleZnnByLumi(TH1* hist)
{
	hist->Scale(fDATA_LUMI/fZNN_LUMI);	
}

//this is to get backgrounds in signal reagions to compare with others
/*void BackgroundsInSignalRegions(const std::string hist, Predictions_t& wjetsPred,
                              Predictions_t& ttbarPred,
                              Predictions_t& znnPred)
{
	TH1* ht500to800_wjets_signal    = dynamic_cast<TH1*> (wjetsRootFile->Get(ht500SignalRegionHist.c_str()));
	assert(ht500to800_wjets_signal != NULL && "sWJET_FILE_NAME:: Fail_lt_point2_500HT800 not found!");
	TH1* ht500to800_ttbar_signal    = dynamic_cast<TH1*> (ttbarRootFile->Get(ht500SignalRegionHist.c_str()));
	assert(ht500to800_ttbar_signal != NULL && "sTTBAR_FILE_NAME:: Fail_lt_point2_500HT800 not found!");
	TH1* ht500to800_znn_signal    = dynamic_cast<TH1*> (znnRootFile->Get(ht500SignalRegionHist.c_str()));
	assert(ht500to800_znn_signal != NULL && "sZNN_FILE_NAME:: Fail_lt_point2_500HT800 not found!");
	ht500to800_wjets_signal->Scale(fDATA_LUMI/fEWK_LUMI);
	ht500to800_ttbar_signal->Scale(fDATA_LUMI/fTTBAR_LUMI);
	ht500to800_znn_signal->Scale(fDATA_LUMI/fZNN_LUMI);

	Predictions_t ht500to800_wjetsInSignalRegions = GetExlcusiveBinContents(ht500to800_wjets_signal);
	Predictions_t ht500to800_ttbarInSignalRegions = GetExlcusiveBinContents(ht500to800_ttbar_signal);
	Predictions_t ht500to800_znnInSignalRegions = GetExlcusiveBinContents(ht500to800_znn_signal);

}
*/

Predictions_t SubstractZnn(TH1* hist, const string histname)
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
	Predictions_t znnBinCont =  GetExlcusiveBinContents(znnHist);
  // PrintExclPredictions(znnBinCont);

	hist->Add(znnHist, -1);

	return znnBinCont;
}


Predictions_t SubstractTTbar(TH1* hist, const string histname)
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
	Predictions_t ttbarBinCont =  GetExlcusiveBinContents(ttbarHist);
//	PrintExclPredictions(ttbarBinCont);

	hist->Add(ttbarHist, -1);

	return ttbarBinCont;
}


Predictions_t SubstractEWK(TH1* hist, const string ewkhistname)
{
	//TFile* ewkFile = new TFile(sWJET_FILE_NAME.c_str());
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
	ewkHist1->Add(ewkHist2,1);
	std::cout << __FUNCTION__ << ": scale1/scale2 = " << scale1 << "/" << scale2 << ":: EWK Integral b4/a4= "<< int_b4  << "/" << ewkHist1->Integral() << std::endl; 

	Predictions_t ewkBinCont =  GetExlcusiveBinContents(ewkHist1);
	PrintExclPredictions(ewkBinCont);

	hist->Add(ewkHist1, -1);

	return ewkBinCont;
}

double GetFitFunctionError(TF1* f1, double x) {

	double err = 0.;
	Double_t val[1] = { x };

	for(Int_t j=0; j<f1->GetNpar(); ++j) 
	{
		err += f1->GradientPar(j,val) * f1->GetParError(j);
	}

	return err;

}


void DumpHist(const TH1* hist)
{
	assert (hist != NULL && "DumpHist:: hist is null!"); 
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
//			if (hist->GetBinCenter(bin) > 500) sum += hist->GetBinContent(bin);
		}
	}

//	cout << "Sum (>500) = " << sum << endl;

}

double expFitFunc(double *x, double *par)
{
	double fitval=0.0;
	if(x[0]>0.0) fitval=par[0] * exp(par[1] * x[0]) + CONST_C;
	else fitval=-1.0E6;
	return fitval;
}
double expFitFunc_Cup(double *x, double *par)
{
	double fitval=0.0;
	if(x[0]>0.0) fitval=par[0] * exp(par[1] * x[0]) + CONST_C_UP;
	else fitval=-1.0E6;
	return fitval;
}
double expFitFunc_Cdown(double *x, double *par)
{
	double fitval=0.0;
	if(x[0]>0.0) fitval=par[0] * exp(par[1] * x[0]) + CONST_C_DOWN;
	else fitval=-1.0E6;
	return fitval;
}



double gausFitFunc(double *x, double *par)
{
	double arg = par[0] * exp(par[1] * x[0]);
	//double fitval = 1.0 / TMath::Erf (arg) - par[2];
	//double fitval = 1.0 / TMath::Erf (arg) - 1 + 0.03;
	//double fitval = 1.0 / TMath::Erf (arg) - 1 + 0.02;
	double fitval = 1.0 / TMath::Erf (arg) - 1 + CONST_C;
	return fitval;
}
double gausFitFunc_Cup(double *x, double *par)
{
	double arg = par[0] * exp(par[1] * x[0]);
	double fitval = 1.0 / TMath::Erf (arg) - 1 + CONST_C_UP;
	return fitval;
}
double gausFitFunc_Cdown(double *x, double *par)
{
	double arg = par[0] * exp(par[1] * x[0]);
	double fitval = 1.0 / TMath::Erf (arg) - 1 + CONST_C_DOWN;
	return fitval;
}


Predictions_t GetPredictions(const TH1* hist, TF1* f1)
{
	assert(hist != NULL && "GetPredictions:: hist not found!");
	assert(f1 != NULL && "GetPredictions:: func not found!");
	//new TCanvas(); gPad->SetLogy(); hist->SetStats(1); hist->DrawCopy(); //gPad->Print("control.eps");

//	std::cout << "INITIAL INFO FOR HIST:"; hist->Print();
	int bin1 = 0;
	for (int bin =0; bin<hist->GetNbinsX(); ++bin) 
	{ 
		if (hist->GetBinLowEdge(bin)<incMHTBins.at(0)) continue;
		else { bin1 = bin; break; }
	}
	int bin2 = hist->GetNbinsX()+1;
	double err =0;
	double integral = hist->IntegralAndError(bin1, bin2, err);
	std::cout << "Intergral for MHT> "<< incMHTBins.at(0) << ":" << integral << "+/-" << err << std::endl;

	double sumGaus[incMHTBins.size()];
	double Gaus_StatErr[incMHTBins.size()];
	double Gaus_FitErr[incMHTBins.size()];
	double sumGaus_excl[incMHTBins.size()-1];
	double Gaus_StatErr_excl[incMHTBins.size()-1];
	double Gaus_FitErr_excl[incMHTBins.size()-1];

	Predictions_t results;

	for (unsigned mhtBin = 0; mhtBin < incMHTBins.size(); ++mhtBin)
	{
		sumGaus[mhtBin] = 0;
		Gaus_StatErr[mhtBin] = 0;
		if (mhtBin+1<incMHTBins.size()) 
		{
			sumGaus_excl[mhtBin] = 0;
			Gaus_StatErr_excl[mhtBin] = 0;
		}

		for (int bin = 0; bin <= hist->GetNbinsX()+1; ++bin)
		{
			if (hist->GetBinContent(bin)>0)
			{
				/*			 cout << std::setprecision(4) << std::setw(5) << bin << std::setw(3) << "[" 
							 << std::setw(10) << hist->GetBinLowEdge(bin) << ", "  << std::setw(10) 
							 << hist->GetXaxis()->GetBinUpEdge(bin)<< "]" 
							 << std::setw(10) << hist->GetBinContent(bin) 
							 << std::setw(10) << hist->GetBinContent(bin) * gausFit2->Eval(hist->GetBinCenter(bin))
							 << endl;
							 */				//inclusive bin stuff

				const float binCenter = hist->GetBinCenter(bin);
				const float binVal    = hist->GetBinContent(bin);
				const float binErr    = hist->GetBinError(bin);
				const float funcVal   = f1->Eval(binCenter);
				const float res       = (binVal * funcVal);
				const float statErr2  = pow(binErr * funcVal, 2);
			 	const float fitErr    = binVal * GetFitFunctionError(f1, binCenter);

				if (hist->GetBinCenter(bin) > incMHTBins.at(mhtBin))
				{
					sumGaus[mhtBin]      += res;
					Gaus_StatErr[mhtBin] += statErr2;
					Gaus_FitErr[mhtBin]  += fitErr;
				}

				//exclsuive bin stuff
				if (mhtBin+1<incMHTBins.size())
				{
					if (hist->GetBinCenter(bin) > incMHTBins.at(mhtBin) 
							&& hist->GetBinCenter(bin) < incMHTBins.at(mhtBin+1))
					{
						sumGaus_excl[mhtBin]      += res;
						Gaus_StatErr_excl[mhtBin] += statErr2;
						Gaus_FitErr_excl[mhtBin]  += fitErr;
					}
				}
			}
		}

		Gaus_StatErr[mhtBin] = sqrt(Gaus_StatErr[mhtBin]);
		Gaus_StatErr_excl[mhtBin] = sqrt(Gaus_StatErr_excl[mhtBin]);

		results.incl_mean[mhtBin] = sumGaus[mhtBin];
		results.incl_statErr[mhtBin] = Gaus_StatErr[mhtBin];
		results.incl_fitErr[mhtBin] = Gaus_FitErr[mhtBin];

		std::stringstream incl_pred;
		incl_pred << setprecision(3) << setw(15) << "MHT > " << incMHTBins.at(mhtBin) 
			<< setw(20) << sumGaus[mhtBin]  << "&$\\pm$" << Gaus_StatErr[mhtBin] << " &$\\pm$ " << Gaus_FitErr[mhtBin];
		//std::cout << incl_pred.str() << std::endl;
		
		if (mhtBin+1<incMHTBins.size())
		{
			std::stringstream excl_pred;
			excl_pred << setprecision(3) << setw(10) << incMHTBins.at(mhtBin) << "<MHT<" << incMHTBins.at(mhtBin+1)
				<< setw(20) << sumGaus_excl[mhtBin]  << "&$\\pm$" << Gaus_StatErr_excl[mhtBin] << " &$\\pm$ " << Gaus_FitErr_excl[mhtBin];
			//std::cout << excl_pred.str() << std::endl;
			results.excl_mean[mhtBin]    = sumGaus_excl[mhtBin];
			results.excl_statErr[mhtBin] = Gaus_StatErr_excl[mhtBin];
			results.excl_fitErr[mhtBin]  = Gaus_FitErr_excl[mhtBin];
		}
	}
	
	return results;
}

TF1* FitPassFailRatio(const int htBin, const int model, TH1* hist, 
				const string title,
				const string printName="factorization")
{
	if (hist == NULL) {cout << __FUNCTION__ << ":: hist is null!" << endl; assert(false);}

	const bool logScale = true;

	//fit range
	stringstream newtitle;
	//newtitle << "Fit Range " << passFailRatio_fitrange_xmin << "--" << passFailRatio_fitrange_xmax << " GeV " << title;
	newtitle << "HT>";
	stringstream epsname;
	epsname << printName;
	//epsname << printName << "_fitrange_" << passFailRatio_fitrange_xmin << "to" << passFailRatio_fitrange_xmax << ".eps";
	if (htBin == 1) { newtitle << 500; epsname << "_HT500" << ".eps"; }
	else if (htBin == 2) { newtitle << 800; epsname << "_HT800" << ".eps"; }
	else if (htBin == 3) { newtitle << 1000; epsname << "_HT1000" << ".eps"; }
	else if (htBin == 4) { newtitle << 1200; epsname << "_HT1200" << ".eps"; }
	else if (htBin == 5) { newtitle << 1400; epsname << "_HT1400" << ".eps"; }

	newtitle << " GeV;#slash{H}_{T} [GeV];Ratio (r);";
	//newtitle << title;
	cout << __FUNCTION__ << ":: Fitting Range = " << passFailRatio_fitrange_xmin << " - " << passFailRatio_fitrange_xmax << endl;

	gStyle->SetOptStat(0);
	//debug stuff
	new TCanvas();
	//gPad->SetGridy();
	gPad->SetTickx();
	if (logScale) gPad->SetLogy();
	hist->SetLineColor(9);
	hist->SetTitle(newtitle.str().c_str());
	hist->SetLineWidth(2);

	gStyle->SetOptFit(1);
	//hist->SetStats(0);
	hist->GetYaxis()->SetRangeUser(1e-2,10);
	hist->GetXaxis()->SetRangeUser(50,490);
	hist->Draw();

	TF1 *expFit=new TF1("fit_1",expFitFunc,passFailRatio_fitrange_xmin,passFailRatio_fitrange_xmax,2);
	expFit->SetParameter(0,0.09);
	expFit->SetParameter(1,-0.0002);
	hist->Fit(expFit,"E0","",passFailRatio_fitrange_xmin, passFailRatio_fitrange_xmax);
	gPad->Modified();
	gPad->Update();
	TPaveStats *e_stats = (TPaveStats*) hist->FindObject("stats");
	e_stats->SetTextColor(kRed);
	TPaveStats *exp_stats = (TPaveStats*) e_stats->Clone("exp_stats");


	TF1 *gausFit = new TF1("gausFit",gausFitFunc, passFailRatio_fitrange_xmin, passFailRatio_fitrange_xmax,2);
	gausFit->SetParameter(0,0.09); 
	gausFit->SetParameter(1,-0.0002);
	gausFit->SetParameter(2,-1.0);
	gausFit->SetLineColor(kGreen);
	hist->Fit(gausFit,"E0","", passFailRatio_fitrange_xmin, passFailRatio_fitrange_xmax);
	gPad->Modified();
	gPad->Update();
	TPaveStats *gaus_stats = (TPaveStats*) hist->FindObject("stats");
	gaus_stats->SetTextColor(kGreen);

	TF1 *expFitResult=new TF1("fit_2",expFitFunc,50,1000.0,2);
	expFitResult->SetParameter(0,expFit->GetParameter(0));
	expFitResult->SetParameter(1,expFit->GetParameter(1));
	expFitResult->SetLineColor(kRed);
	expFitResult->SetLineWidth(2);
	expFitResult->SetParameters( expFit->GetParameters()); 
	expFitResult->SetParErrors ( expFit->GetParErrors() );
	expFitResult->SetChisquare ( expFit->GetChisquare() );
	expFitResult->SetNDF       ( expFit->GetNDF()       );
	expFitResult->SetLineColor ( kRed                    );
	expFitResult->SetLineWidth ( 2                       );
	expFitResult->DrawCopy     ( "same"                  );



	TF1 *gausFitResult = new TF1("GausFitResult",gausFitFunc,50,1000.0,2);
	gausFitResult->SetParameters( gausFit->GetParameters()); 
	gausFitResult->SetParErrors ( gausFit->GetParErrors() );
	gausFitResult->SetChisquare ( gausFit->GetChisquare() );
	gausFitResult->SetNDF       ( gausFit->GetNDF()       );
	gausFitResult->SetLineColor ( kGreen                  );
	gausFitResult->SetLineWidth ( 2                       );
	gausFitResult->DrawCopy     ( "same"                  );

	TLegend *leg  = new TLegend(0.68,0.7,0.9,0.9);
	leg->AddEntry(gausFitResult,"Gaussian Model");
	leg->AddEntry(expFitResult,"Exponential Model");
	leg->Draw();

	const float xmin=0.2, xmax=0.45, ymin=0.7, ymax=0.9;
	gaus_stats->SetX1NDC(xmin);
	gaus_stats->SetX2NDC(xmax);
	gaus_stats->SetY1NDC(ymin);
	gaus_stats->SetY2NDC(ymax);
	gaus_stats->Draw("same");
	exp_stats->SetX1NDC(xmax);
	exp_stats->SetX2NDC(xmax+0.22);
	exp_stats->SetY1NDC(ymin);
	exp_stats->SetY2NDC(ymax);
	exp_stats->Draw("same");
	//overlay the exponential model as well


	gPad->Print(epsname.str().c_str());

	if (model == 1) return gausFitResult;
	else return expFitResult;
}

void makeExpModelPredFromData(const int htBin=1, const int model = 1)
//      const float fitrange_xmin = 50, const float fitrange_xmax = 150) 
{

	incMHTBins.clear();
	for (int i=0; i<nMHTbins; ++i) incMHTBins.push_back(arrMHTbins[i]);

	TH1 *histPassOrig = 0;
	TH1 *histFailOrig = 0;

	TFile* dataRootFile  = new TFile (sDATA_FILE_NAME.c_str());
//	TFile* wjetsRootFile = new TFile (sWJET_FILE_NAME.c_str());
	TFile* ttbarRootFile = new TFile (sTTBAR_FILE_NAME.c_str());
	TFile* znnRootFile   = new TFile (sZNN_FILE_NAME.c_str());
	if (dataRootFile->IsZombie()) { cout << "Data root file not found!" <<  endl; assert (false); }
//	if (wjetsRootFile->IsZombie()) { cout << "wjets root file not found!" <<  endl; assert (false); }
	if (ttbarRootFile->IsZombie()) { cout << "ttbar root file not found!" <<  endl; assert (false); }
	if (znnRootFile->IsZombie()) { cout << "znn root file not found!" <<  endl; assert (false); }

	std::string numerHistName(""), denomHistName("");
	std::string QCDsyst1_numerHistName(""), QCDsyst1_denomHistName("");
	std::string QCDsyst2_numerHistName(""), QCDsyst2_denomHistName("");
	std::string qcdsyst1_sidebandHistName(""), qcdsyst2_sidebandHistName("");
	std::string sideband_hist_name("");
	std::string HTlabel("");
	std::string dirName("");

	if (htBin == 1) { //ht>500 GeV
		HTlabel += "500";
		dirName = "factorization/HT500to800/";
	} else if (htBin == 2) { //ht>800
		HTlabel += "800";
		dirName = "factorization/HT800to1000/";

	} else if (htBin == 3) { //ht>1000
		HTlabel += 1000;
		dirName = "factorization/HT1000to1200/";

	} else if (htBin == 4) { //ht>1200
		HTlabel += 1200;
		dirName = "factorization/HT1200to1400/";

	} else if (htBin == 5) { //ht>1400
		HTlabel += 1400;
		dirName = "factorization/HT1400to7000/";

	} else 
	{
		assert (false && "This HT bin is not defined!!!");
	}

	if (dirName.length()<1) assert(false && "dirName is empty!!");

	numerHistName             += dirName; 
	denomHistName             += dirName; 
	QCDsyst1_denomHistName    += dirName; 
	qcdsyst1_sidebandHistName += dirName; 
	QCDsyst2_denomHistName    += dirName; 
	qcdsyst2_sidebandHistName += dirName; 
	sideband_hist_name        += dirName;   


	numerHistName             += "signal";
	denomHistName             += "fail1";
	QCDsyst1_numerHistName     = numerHistName;
	QCDsyst1_denomHistName    += "sidebandSyst1";
	qcdsyst1_sidebandHistName += "sidebandSyst1_fineBin";
	QCDsyst2_numerHistName     = numerHistName;
	QCDsyst2_denomHistName    += "sidebandSyst2";
	qcdsyst2_sidebandHistName += "sidebandSyst2_fineBin";
	sideband_hist_name        += "failFineBin1";



	cout << "Using numerHistName = " << numerHistName << endl;
	cout << "Using denomHistName = " << denomHistName << endl;

	histPassOrig = dynamic_cast<TH1*> (dataRootFile->Get(numerHistName.c_str()));
	histFailOrig = dynamic_cast<TH1*> (dataRootFile->Get(denomHistName.c_str()));
	if (histPassOrig == 0 ) { cout << "Hist " << numerHistName << " not found!" << endl; assert (false); }
	if (histFailOrig == 0 ) { cout << "Hist " << denomHistName << " not found!" << endl; assert (false); }
	histPassOrig->Sumw2();
	histFailOrig->Sumw2();


	histPassOrig->SetLineWidth(2);
	histPassOrig->SetLineColor(kBlue);
	histFailOrig->SetLineWidth(2);
	histFailOrig->SetLineColor(kBlack);

	TH1*	HIST_numer = dynamic_cast<TH1*> (histPassOrig->Clone("histpass_copy"));
	TH1*	HIST_denom = dynamic_cast<TH1*> (histFailOrig->Clone("histfail_copy"));


/*	cout << __LINE__ << endl;
	TCanvas* cx = new TCanvas();
	gPad->SetLogy();
	gStyle->SetOptStat(0);
	HIST_denom->SetLineColor(kBlue);
	HIST_numer->SetLineColor(kRed);
	HIST_denom->SetTitle(";MHT;");
	HIST_denom->DrawCopy();
	HIST_numer->DrawCopy("Psame");

	TLegend *leg1  = new TLegend(0.7,0.8,0.9,0.9);
	leg1->AddEntry(HIST_numer,"PASS");
	leg1->AddEntry(HIST_denom,"FAIL");
	leg1->Draw();
*/	

	/********************************************/
	//substract EWK contamination
	/********************************************/
	Predictions_t dataInNumer_b4 = GetExlcusiveBinContents(HIST_numer);
	Predictions_t dataInDenom_b4 = GetExlcusiveBinContents(HIST_denom);

	Predictions_t ewkInNumerator   = SubstractEWK  (HIST_numer, numerHistName);
	Predictions_t ttbarInNumerator = SubstractTTbar(HIST_numer, numerHistName);
	Predictions_t znnInNumerator   = SubstractZnn  (HIST_numer, numerHistName);
	Predictions_t ewkInDenom   = SubstractEWK  (HIST_denom, denomHistName);
	Predictions_t ttbarInDenom = SubstractTTbar(HIST_denom, denomHistName);
	Predictions_t znnInDenom   = SubstractZnn  (HIST_denom, denomHistName);

	Predictions_t dataInNumer_a4 = GetExlcusiveBinContents(HIST_numer);
	Predictions_t dataInDenom_a4 = GetExlcusiveBinContents(HIST_denom);

	HIST_numer->Divide(HIST_denom);


	Predictions_t dataRatio = GetExlcusiveBinContents(HIST_numer);

	HIST_numer->GetXaxis()->SetRangeUser(50,500);


	stringstream tmpTitle;
	tmpTitle << "HT>" << HTlabel << " GeV;#slash{H}_{T} [GeV];Ratio (r);";
	const string title = tmpTitle.str();

	TF1 *fitFunc = 0;
	TF1 *fitFunc_Cup = 0;
	TF1 *fitFunc_Cdown = 0;


	if (model == 1) //Gaussian Model
	{
		/***************************
		 * this is all gaussian
		 * ************************/
		fitFunc =  FitPassFailRatio(htBin, model, HIST_numer, title, "Factorization");
		if (fitFunc == NULL) { cout << __LINE__ << "::functions null! " << endl; assert (false); }

		//These two are for const 'c' systematic derivations only!
		cout << "gausFit_C MEAN is using c = " << CONST_C << endl;
		cout << "fitFunc_Cup is using    c+= " << CONST_C_UP << endl;
		fitFunc_Cup = new TF1("fitFunc_Cup", gausFitFunc_Cup, passFailRatio_fitrange_xmin, 1500,2);
		fitFunc_Cup->SetParameters( fitFunc->GetParameters()); 
		cout << "fitFunc_Cdown is using  c-= " << CONST_C_DOWN << endl;
		fitFunc_Cdown = new TF1("gausFit_Cdown", gausFitFunc_Cdown, passFailRatio_fitrange_xmin, 1500,2);
		fitFunc_Cdown->SetParameters( fitFunc->GetParameters()); 	

	} else {
		/***************************
		 * this is all exponetial
		 * ************************/
		fitFunc =  FitPassFailRatio(htBin, model, HIST_numer, title, "Factorization");
		if (fitFunc == NULL) { cout << __LINE__ << "::functions null! " << endl; assert (false); }

		//These two are for const 'c' systematic derivations only!
		cout << "expFit_C MEAN is using c = " << CONST_C << endl;
		cout << "expFit_Cup is using    c+= " << CONST_C_UP << endl;
		fitFunc_Cup = new TF1("fitFunc_Cup",expFitFunc_Cup,passFailRatio_fitrange_xmin, 1500,2);
		fitFunc_Cup->SetParameters( fitFunc->GetParameters()); 
		cout << "expFit_Cdown is using  c-= " << CONST_C_DOWN << endl;
		fitFunc_Cdown = new TF1("gausFit_Cdown",expFitFunc_Cdown, passFailRatio_fitrange_xmin, 1500,2);
		fitFunc_Cdown->SetParameters( fitFunc->GetParameters()); 	
	}



	/**********************************************
	 * Now Get QCD selection dependence syst-type1 fit
	 * *******************************************/

	TH1* QCDsyst1_numerHist = dynamic_cast<TH1*> ((dataRootFile->Get(QCDsyst1_numerHistName.c_str()))->Clone("1"));
	TH1* QCDsyst1_denomHist = dynamic_cast<TH1*> (dataRootFile->Get(QCDsyst1_denomHistName.c_str()));
	assert(QCDsyst1_numerHist != NULL && "QCDsyst1_numerHistName not found!");
	assert(QCDsyst1_denomHist != NULL && "QCDsyst1_denomHistName not found!");
	SubstractEWK  (QCDsyst1_numerHist, QCDsyst1_numerHistName);
	SubstractTTbar(QCDsyst1_numerHist, QCDsyst1_numerHistName);
	SubstractZnn  (QCDsyst1_numerHist, QCDsyst1_numerHistName);
	SubstractEWK  (QCDsyst1_denomHist, QCDsyst1_denomHistName);
	SubstractTTbar(QCDsyst1_denomHist, QCDsyst1_denomHistName);
	SubstractZnn  (QCDsyst1_denomHist, QCDsyst1_denomHistName);
	QCDsyst1_numerHist->Divide(QCDsyst1_denomHist);
	const string QCDsyst1_title = "QCD SYST1: : HT>500 GeV;#slash{H}_{T};Ratio (r) = Pass(RA2 dPhi cuts) / Fail(syst1);";
	TF1 *QCDsyst1_gausFitFunc =  FitPassFailRatio(htBin, model, QCDsyst1_numerHist, QCDsyst1_title, "QCDSyst1");
	if (QCDsyst1_gausFitFunc == NULL) { cout << __LINE__ << "::QCD syst1 functions null! " << endl; assert (false); }

	/**********************************************
	 * Now Get QCD selection dependence syst-type2 fit
	 * *******************************************/

/*H1* QCDsyst2_numerHist = dynamic_cast<TH1*> (dataRootFile->Get(QCDsyst2_numerHistName.c_str()));
//	TH1* QCDsyst2_denomHist = dynamic_cast<TH1*> (dataRootFile->Get(QCDsyst2_denomHistName.c_str()));
	assert(QCDsyst2_numerHist != NULL && "QCDsyst2_numerHistName not found!");
//	assert(QCDsyst2_denomHist != NULL && "QCDsyst2_denomHistName not found!");
	QCDsyst2_numerHist->Divide(QCDsyst2_denomHist);
	SubstractEWK  (QCDsyst2_numerHist, QCDsyst2_numerHistName);
	SubstractTTbar(QCDsyst2_numerHist, QCDsyst2_numerHistName);
	SubstractZnn  (QCDsyst2_numerHist, QCDsyst2_numerHistName);
	SubstractEWK  (QCDsyst2_denomHist, QCDsyst2_denomHistName);
	SubstractTTbar(QCDsyst2_denomHist, QCDsyst2_denomHistName);
	SubstractZnn  (QCDsyst2_denomHist, QCDsyst2_denomHistName);
//	const string QCDsyst2_title = "QCD SYST2: : HT>500 GeV;#slash{H}_{T};Ratio (r) = Pass(RA2 dPhi cuts) / Fail(syst2);";
//	TF1 *QCDsyst2_gausFit =  FitPassFailRatio(model, QCDsyst2_numerHist, QCDsyst2_title, "QCDSyst2_HT500");
//	if (QCDsyst2_gausFit == NULL) { cout << __LINE__ << "::QCD syst2 functions null! " << endl; assert (false); }
*/



	TH1* sideband_hist    = dynamic_cast<TH1*> (dataRootFile->Get(sideband_hist_name.c_str()));
	assert(sideband_hist != NULL && "Fail_lt_point2_500HT800 not found!");
	SubstractEWK  (sideband_hist, sideband_hist_name);
	SubstractTTbar(sideband_hist, sideband_hist_name);
	SubstractZnn  (sideband_hist, sideband_hist_name);
	Predictions_t pred_HT500to800        = GetPredictions(sideband_hist, fitFunc);

	TH1* qcdsyst1_sideband_hist    = dynamic_cast<TH1*> (dataRootFile->Get(qcdsyst1_sidebandHistName.c_str()));
	assert(qcdsyst1_sideband_hist != NULL && "Syst1_Fail_500HT800 not found!");
	SubstractEWK  (qcdsyst1_sideband_hist, qcdsyst1_sidebandHistName);
	SubstractTTbar(qcdsyst1_sideband_hist, qcdsyst1_sidebandHistName);
	SubstractZnn  (qcdsyst1_sideband_hist, qcdsyst1_sidebandHistName);
	Predictions_t pred_HT500to800_qcdSyst1 = GetPredictions(qcdsyst1_sideband_hist, QCDsyst1_gausFitFunc);
	
//	TH1* qcdsyst2_sideband_hist    = dynamic_cast<TH1*> (dataRootFile->Get(qcdsyst2_sideband_hist.c_str()));
//	assert(qcdsyst2_sideband_hist != NULL && "Syst2_Fail_500HT800 not found!");
//	Predictions_t pred_HT500to800_qcdSyst2 = GetPredictions(qcdsyst2_sideband_hist, QCDsyst2_gausFit);

////	cout << "QCD SYST 1 AND SYST2 " << std::endl;
	//PrintExclPredictions(pred_HT500to800_qcdSyst1);
	//PrintExclPredictions(pred_HT500to800_qcdSyst2);

	Predictions_t pred_HT500to800_Cup_syst   = GetPredictions(sideband_hist, fitFunc_Cup);
	Predictions_t pred_HT500to800_Cdown_syst = GetPredictions(sideband_hist, fitFunc_Cdown);


	if (model == 1) cout << " ================ GAUS MODEL ==================";
	else cout << " ================ EXPO MODEL ==================";

	if (htBin == 1) cout << "\n\n >>>>>> 500<HT<800 PREDICTIONS " << endl; 
	else if (htBin == 2) cout << "\n\n >>>>>> 800<HT<1000 PREDICTIONS " << endl; 
	else if (htBin == 3) cout << "\n\n >>>>>> 1000<HT<1200 PREDICTIONS " << endl; 
	else if (htBin == 4) cout << "\n\n >>>>>> 1200<HT<1400 PREDICTIONS " << endl; 
	else if (htBin == 5) cout << "\n\n >>>>>> HT>1400 PREDICTIONS " << endl; 

	PrintResults(pred_HT500to800, pred_HT500to800_qcdSyst1, pred_HT500to800_Cup_syst, pred_HT500to800_Cdown_syst);
	PrintBackgrounds(dataInNumer_b4, dataInNumer_a4, ewkInNumerator, ttbarInNumerator, znnInNumerator);
	PrintBackgrounds(dataInDenom_b4, dataInDenom_a4, ewkInDenom, ttbarInDenom, znnInDenom);
	cout << "====== Final Ratio Plot Values =======" << endl;
	PrintExclPredictions(dataRatio);


}


void makeAllPlots()
{
for (int bin=1; bin <=5; ++bin)
	 makeExpModelPredFromData(bin, 1);
}

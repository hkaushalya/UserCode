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
//const static float fDATA_LUMI  = 1150.0; //pb-1
const static float fEWK_LUMI   = 81352581.0/31300.0;        // 31300 pb for NLO xsec and total number of events is 81352581 (sample1+sample2)
const static float fTTBAR_LUMI = 3701947.0/165.0;
const static float fZNN_LUMI   = 3067017.0/42.2;              //x-sec from note, Anwar hd used 32.92 pb-1 LO
//const static string sDATA_FILE_NAME  = "Data.root";
//const static string sDATA_FILE_NAME  = "Data_SmearedJetsNoPrescaleWgts.root";
//const static string sDATA_FILE_NAME  = "Data_PreScaleWgted.root";
//const static string sDATA_FILE_NAME  = "Data_1stfb_12052011.root";
//const static string sWJET_FILE_NAME  = "EWKmc.root";
//const static string sDATA_FILE_NAME = "SmearedJetsWithPrescaleWgts.root";
const static string sDATA_FILE_NAME = "Data_MyPrescaleWgted.root";
const static string sWJET_FILE_NAME  = "WJetsMC.root";
const static string sTTBAR_FILE_NAME = "TTbarMC.root";
const static string sZNN_FILE_NAME   = "ZnnMC.root";
const static float passFailRatio_fitrange_xmin = 50.0;
const static float passFailRatio_fitrange_xmax = 150.0;
std::vector<float> incMHTBins;
const static int nMHTbins = 5;
const static float arrMHTbins[nMHTbins] = {200,350,500,600,7000};
float CONST_C      = 0.0217;
//float CONST_C      = 0.03;
float CONST_C_UP   = 0.07;
float CONST_C_DOWN = 0.012;
//float CONST_C = 0.03;
struct Predictions_t
{
		float incl_mean[nMHTbins];
		float incl_statErr[nMHTbins];
		float incl_fitErr[nMHTbins];
		float excl_mean[nMHTbins-1];
		float excl_statErr[nMHTbins-1];
		float excl_fitErr[nMHTbins-1];
};

void PrintExclPredictions(const Predictions_t& res)
{

	cout << setprecision(3) << setw(15) << " MHT "
					<< setw(10) << " mean " 
					<< setw(10) << "stat"
					<< endl;
	for (int i=0; i< nMHTbins -1 ; ++i)
	{
			const float mean    = res.excl_mean[i]; 
			const float statErr = res.excl_statErr[i];

		cout << setprecision(3) << setw(15) << incMHTBins.at(i) << "<MHT<" << incMHTBins.at(i+1) 
					<< setw(10) << mean << " $&\\pm$ " 
					<< setw(10) << statErr
					<< endl;
	}
}



void PrintResults(const Predictions_t& res, 
					const Predictions_t& htSystRes,
					const Predictions_t& qcdSystRes,
					const Predictions_t& C_upSyst,
					const Predictions_t& C_downSyst
					)
{

	cout.precision(2);
	cout << setw(10) << " MHT "
					<< setw(20) << " mean " 
					<< setw(15) << ":statAndFitErr"
					<< setw(10) << " htSyst"
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
			const float htSys         = fabs(mean - htSystRes.excl_mean[i]);
			const float qcdSys1       = fabs(mean - qcdSystRes.excl_mean[i]);
			const float C_Sys         = max(fabs(mean - C_upSyst.excl_mean[i]), fabs(mean - C_downSyst.excl_mean[i]));
			const float C_Sys_up      = fabs(mean - C_upSyst.excl_mean[i]);
			const float C_Sys_down    = fabs(mean - C_downSyst.excl_mean[i]);

			const float total         = sqrt(  pow(statAndFitErr,2)
			                                 + pow(htSys,2)
														+ pow(qcdSys1,2) 
														+ pow(C_Sys,2) 
														);

			const float total_up      = sqrt(  pow(statAndFitErr,2)
			                                 + pow(htSys,2)
														+ pow(qcdSys1,2) 
														+ pow(C_Sys_up,2) 
														);

			const float total_down      = sqrt(  pow(statAndFitErr,2)
			                                 + pow(htSys,2)
														+ pow(qcdSys1,2) 
														+ pow(C_Sys_down,2) 
														);
		cout << fixed << setw(5) << incMHTBins.at(i) << "<MHT<" << setw(7) << incMHTBins.at(i+1)
					<< " & "
					<< setw(10) << mean 
					<< setw(15) << " & $\\pm$ " 
					<< statAndFitErr
					<< setw(15) << " & $\\pm$ " << htSys
					<< setw(15) << " & $\\pm$ " << qcdSys1
					<< setw(10) << " & " << "$^{+" << C_Sys_up << "}_{-" << C_Sys_down << "}$"
					<< setw(10) << " & " << "$^{+" << total_up << "}_{-" << total_down << "}$" << "   \\\\"
					<< endl;
	}

	for (int i=0; i< nMHTbins ; ++i)
	{
			const float mean          = res.incl_mean[i]; 
			const float statAndFitErr = sqrt(pow(res.incl_statErr[i],2)+pow(res.incl_fitErr[i],2));
			const float htSys         = fabs(mean - htSystRes.incl_mean[i]);
			const float qcdSys1       = fabs(mean - qcdSystRes.incl_mean[i]);
			const float C_Sys         = max(fabs(mean - C_upSyst.incl_mean[i]), fabs(mean - C_downSyst.incl_mean[i]));
			const float C_Sys_up      = fabs(mean - C_upSyst.incl_mean[i]);
			const float C_Sys_down    = fabs(mean - C_downSyst.incl_mean[i]);

			const float total         = sqrt(  pow(statAndFitErr,2)
			                                 + pow(htSys,2)
														+ pow(qcdSys1,2) 
														+ pow(C_Sys,2) 
														);

			const float total_up      = sqrt(  pow(statAndFitErr,2)
			                                 + pow(htSys,2)
														+ pow(qcdSys1,2) 
														+ pow(C_Sys_up,2) 
														);

			const float total_down    = sqrt(  pow(statAndFitErr,2)
			                                 + pow(htSys,2)
														+ pow(qcdSys1,2) 
														+ pow(C_Sys_down,2) 
														);

		cout << fixed << setw(5) << "MHT> " << setw(7) << incMHTBins.at(i) 
					<< setw(10) << " & " << mean 
					<< setw(15) << " & $\\pm$ " 
					<< setw(2) << statAndFitErr
					<< setw(10) << " & $\\pm$ " << htSys
					<< setw(10) << " & $\\pm$ " << qcdSys1
					<< setw(10) << " & " << "$^{+" << C_Sys_up << "}_{-" << C_Sys_down << "}$"
					<< setw(10) << " & " << "$^{+" << total_up << "}_{-" << total_down << "}$" << "   \\\\"
					<< endl;
	}

}


/*******************
 * Find events in each exclusive bin for a given hist
 *******************/
Predictions_t GetExlcusiveBinContents(const TH1* hist)
{
	assert(hist != NULL && "GetPredictions:: hist not found!");

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
			excl_pred << setprecision(3) << setw(10) << incMHTBins.at(mhtBin) << "<MHT<" << incMHTBins.at(mhtBin+1)
				<< mean[mhtBin]<< "&$\\pm$" << statErr[mhtBin];
			//std::cout << excl_pred.str() << std::endl;
			results.excl_mean[mhtBin]    = mean[mhtBin];
			results.excl_statErr[mhtBin] = statErr[mhtBin];
			results.excl_fitErr[mhtBin]  = -1;
		}
	}
	
	return results;
}





Predictions_t SubstractZnn(TH1* hist, const string histname)
{
	TFile* znnFile = new TFile(sTTBAR_FILE_NAME.c_str());
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
	//PrintExclPredictions(ewkBinCont);

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
	//PrintExclPredictions(ewkBinCont);

	hist->Add(ttbarHist, -1);

	return ttbarBinCont;
}


Predictions_t SubstractEWK(TH1* hist, const string ewkhistname)
{
	TFile* ewkFile = new TFile(sWJET_FILE_NAME.c_str());
	if (ewkFile->IsZombie()) { cout << "EWK file not found!" << endl; assert (false); }

	TH1* ewkHist = dynamic_cast<TH1*> (ewkFile->Get(ewkhistname.c_str()));
	if (ewkHist == NULL) { cout << "EWK hist " << ewkhistname << " not found!" << endl; assert (false); }
	ewkHist->Sumw2();
	//new TCanvas();
	//ewkHist->DrawCopy();

	const float int_b4 = ewkHist->Integral(); 
	const float scale  = fDATA_LUMI/fEWK_LUMI;

	ewkHist->Scale(scale);
	std::cout << __FUNCTION__ << ": scale = " << scale<< ":: EWK Integral b4/a4= "<< int_b4  << "/" << ewkHist->Integral() << std::endl; 

	Predictions_t ewkBinCont =  GetExlcusiveBinContents(ewkHist);
	//PrintExclPredictions(ewkBinCont);

	hist->Add(ewkHist, -1);

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

double expFitFunc(double *x, double *par)
{
	double fitval=0.0;
	//if(x[0]>0.0) fitval=par[0] * exp(par[1] * x[0]) + 0.03;
	//if(x[0]>0.0) fitval=par[0] * exp(par[1] * x[0]) + 0.0217;
	if(x[0]>0.0) fitval=par[0] * exp(par[1] * x[0]) + CONST_C;
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

	std::cout << "INITIAL INFO FOR HIST:"; hist->Print();
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

TF1* FitPassFailRatio(TH1* hist, 
				const string title,
				const string printName="factorization")
{
	if (hist == NULL) {cout << __FUNCTION__ << ":: hist is null!" << endl; assert(false);}

	const bool logScale = true;

	//fit range
	stringstream newtitle;
	newtitle << "Fit Range " << passFailRatio_fitrange_xmin << "--" << passFailRatio_fitrange_xmax << " GeV " << title;
	//newtitle << ";MHT;Ratio (r)";
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
	hist->Draw();

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

	TF1 *gausFitResult = new TF1("GausFitResult",gausFitFunc,50,1000.0,2);
	gausFitResult->SetParameters( gausFit->GetParameters()); 
	gausFitResult->SetParErrors ( gausFit->GetParErrors() );
	gausFitResult->SetChisquare ( gausFit->GetChisquare() );
	gausFitResult->SetNDF       ( gausFit->GetNDF()       );
	gausFitResult->SetLineColor ( kGreen                  );
	gausFitResult->SetLineWidth ( 1                       );
	gausFitResult->DrawCopy     ( "same"                  );

	TLegend *leg  = new TLegend(0.7,0.8,0.9,0.9);
	leg->AddEntry(gausFitResult,"Gaussian Model");
	leg->Draw();

	const float xmin=0.2, xmax=0.45, ymin=0.7, ymax=0.9;
	gaus_stats->SetX1NDC(xmin);
	gaus_stats->SetX2NDC(xmax);
	gaus_stats->SetY1NDC(ymin);
	gaus_stats->SetY2NDC(ymax);
	gaus_stats->Draw("same");

	stringstream epsname;
	epsname << printName << "_fitrange_" << passFailRatio_fitrange_xmin << "to" << passFailRatio_fitrange_xmax << ".eps";
	gPad->Print(epsname.str().c_str());

	return gausFitResult;
}

void makePassFail_data(const float fitrange_xmin = 50, const float fitrange_xmax = 150) 
{

	incMHTBins.clear();
	for (int i=0; i<nMHTbins; ++i) incMHTBins.push_back(arrMHTbins[i]);

	TH1 *histPassOrig = 0;
	TH1 *histFailOrig = 0;

	TFile* dataRootFile = new TFile (sDATA_FILE_NAME.c_str());
	if (dataRootFile->IsZombie())
	{
		cout << "Data root file not found!" <<  endl;
		assert (false);
	}

	const std::string numerHistName("factorization_ht500/Pass_RA2dphi_HT500");
	const std::string denomHistName("factorization_ht500/Fail_1");

	histPassOrig = dynamic_cast<TH1*> (dataRootFile->Get(numerHistName.c_str()));
	histFailOrig = dynamic_cast<TH1*> (dataRootFile->Get(denomHistName.c_str()));
	if (histPassOrig == 0 ) { cout << "Hist " << numerHistName << " not found!" << endl; assert (false); }
	if (histFailOrig == 0 ) { cout << "Hist " << denomHistName << " not found!" << endl; assert (false); }
	histPassOrig->Sumw2();
	histFailOrig->Sumw2();

//	new TCanvas();
//	histPassOrig->DrawCopy();
//	new TCanvas();
//	histFailOrig->DrawCopy();

	//DumpHist(HIST_numer);
	//DumpHist(HIST_denom);
	//std::cout << "Entries pass/ fail " << histPassOrig->GetEntries() << " / " << histFailOrig->GetEntries() << std::endl;
	//new TCanvas();
	//histPassOrig->DrawCopy();
	//histFailOrig->DrawCopy("same");
	//gPad->SetEditable(0);

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

	SubstractEWK  (HIST_numer, numerHistName);
	SubstractTTbar(HIST_numer, numerHistName);
	SubstractZnn  (HIST_numer, numerHistName);
	SubstractEWK  (HIST_denom, denomHistName);
	SubstractTTbar(HIST_denom, denomHistName);
	SubstractZnn  (HIST_denom, denomHistName);
	HIST_numer->Divide(HIST_denom);
	
	const string title = " : HT>500 GeV;#slash{H}_{T};Ratio (r) = Pass(RA2 dPhi cuts) / Fail(#Delta #phi_{min}< 0.2);";
	TF1 *gausFit2 =  FitPassFailRatio(HIST_numer, title, "FactorizationHT500");
	if (gausFit2 == NULL) { cout << __LINE__ << "::functions null! " << endl; assert (false); }

	//These two are for const 'c' systematic derivations only!
	TF1 *gausFit_Cup = new TF1("gausFit_Cup",gausFitFunc_Cup,passFailRatio_fitrange_xmin, 1500,2);
	gausFit_Cup->SetParameters( gausFit2->GetParameters()); 
	TF1 *gausFit_Cdown = new TF1("gausFit_Cup",gausFitFunc_Cdown, passFailRatio_fitrange_xmin, 1500,2);
	gausFit_Cdown->SetParameters( gausFit2->GetParameters()); 	

	/*gausFit_Cup->SetLineColor(kRed);
	gausFit_Cdown->SetLineColor(kRed);
	
	new TCanvas();
	gPad->SetLogy();
	TLegend *tleg = new TLegend(0.7,0.8,0.9,0.9);
	tleg->AddEntry(gausFit2,"Gaussian Fit");
	tleg->AddEntry(gausFit_Cup,"Uncerainty of const 'c'");
	gausFit2->Draw();
	gausFit_Cup->Draw("same");
	gausFit_Cdown->Draw("same");
	tleg->Draw();
	gPad->Print("fit_c_error.eps");
	return;
*/


	/**********************************************
	 * Now Get HT dependence syst fit
	 * *******************************************/
	const std::string HTsyst_numerHistName("factorization_ht600/Pass_RA2dphi");
	const std::string HTsyst_denomHistName("factorization_ht600/Fail_1");
	TH1* HTsyst_numerHist = dynamic_cast<TH1*> (dataRootFile->Get(HTsyst_numerHistName.c_str()));
	TH1* HTsyst_denomHist = dynamic_cast<TH1*> (dataRootFile->Get(HTsyst_denomHistName.c_str()));
	assert(HTsyst_numerHist != NULL && "HTsyst_numerHistName not found!");
	assert(HTsyst_denomHist != NULL && "HTsyst_denomHistName not found!");

	SubstractEWK  (HTsyst_numerHist, HTsyst_numerHistName);
	SubstractTTbar(HTsyst_numerHist, HTsyst_numerHistName);
	SubstractZnn  (HTsyst_numerHist, HTsyst_numerHistName);
	SubstractEWK  (HTsyst_denomHist, HTsyst_denomHistName);
	SubstractTTbar(HTsyst_denomHist, HTsyst_denomHistName);
	SubstractZnn  (HTsyst_denomHist, HTsyst_denomHistName);
	HTsyst_numerHist->Divide(HTsyst_denomHist);

	const string HTsyst_title = " : HT>600 GeV;#slash{H}_{T};Ratio (r) = Pass(RA2 dPhi cuts) / Fail(#Delta #phi_{min}< 0.2);";
	TF1 *HTsyst_gausFit =  FitPassFailRatio(HTsyst_numerHist, HTsyst_title, "FactorizationHT600");
	if (HTsyst_gausFit == NULL) { cout << __LINE__ << "::HT syst functions null! " << endl; assert (false); }


	/**********************************************
	 * Now Get QCD selection dependence syst-type1 fit
	 * *******************************************/
	const std::string QCDsyst1_numerHistName("factorization_ht500/Pass_RA2dphi");
	const std::string QCDsyst1_denomHistName("factorization_ht500/Syst1_Fail_ht500");
	TH1* QCDsyst1_numerHist = dynamic_cast<TH1*> ((dataRootFile->Get(QCDsyst1_numerHistName.c_str()))->Clone("1"));
	TH1* QCDsyst1_denomHist = dynamic_cast<TH1*> (dataRootFile->Get(QCDsyst1_denomHistName.c_str()));
	assert(QCDsyst1_numerHist != NULL && "QCDsyst1_numerHistName not found!");
	assert(QCDsyst1_denomHist != NULL && "QCDsyst1_denomHistName not found!");
	QCDsyst1_numerHist->Divide(QCDsyst1_denomHist);
	//new TCanvas();
	//QCDsyst1_numerHist->DrawCopy();
	const string QCDsyst1_title = "QCD SYST1: : HT>500 GeV;#slash{H}_{T};Ratio (r) = Pass(RA2 dPhi cuts) / Fail(syst1);";
	TF1 *QCDsyst1_gausFit =  FitPassFailRatio(QCDsyst1_numerHist, QCDsyst1_title, "QCDSyst1_HT500");
	if (QCDsyst1_gausFit == NULL) { cout << __LINE__ << "::QCD syst1 functions null! " << endl; assert (false); }


	const std::string QCDsyst2_numerHistName("factorization_ht500/Pass_RA2dphi");
	const std::string QCDsyst2_denomHistName("factorization_ht500/Syst2_Fail_ht500");
	TH1* QCDsyst2_numerHist = dynamic_cast<TH1*> (dataRootFile->Get(QCDsyst2_numerHistName.c_str()));
	TH1* QCDsyst2_denomHist = dynamic_cast<TH1*> (dataRootFile->Get(QCDsyst2_denomHistName.c_str()));
	assert(QCDsyst2_numerHist != NULL && "QCDsyst2_numerHistName not found!");
	assert(QCDsyst2_denomHist != NULL && "QCDsyst2_denomHistName not found!");
	QCDsyst2_numerHist->Divide(QCDsyst2_denomHist);
//	new TCanvas();
//	QCDsyst2_numerHist->DrawCopy();
//	QCDsyst2_denomHist->DrawCopy("same");
	const string QCDsyst2_title = "QCD SYST2: : HT>500 GeV;#slash{H}_{T};Ratio (r) = Pass(RA2 dPhi cuts) / Fail(syst2);";
	TF1 *QCDsyst2_gausFit =  FitPassFailRatio(QCDsyst2_numerHist, QCDsyst2_title, "QCDSyst2_HT500");
	if (QCDsyst2_gausFit == NULL) { cout << __LINE__ << "::QCD syst2 functions null! " << endl; assert (false); }


//	new baseline:
//	(baseline) HT>500, MHT>200

//	for comparison with SUSY11 result:
//		(high HT) HT>800, MHT>200
//		(high MHT) HT>800, MHT>500
//		(medium) HT>500, MHT>350

//		new signal regions with highest expected sensitivity:
//		(high m12, low m0) HT>1000, MHT>600
//		(high m0) HT>1200, MHT>400



	 //make a predicion
/*
	cout << "\n\n >>>>>> HT>500 PREDICTIONS " << endl; 
	TH1* ht500_fail = dynamic_cast<TH1*> (dataRootFile->Get("factorization_ht500/Fail_lt_point2_HT500"));
	assert(ht500_fail != NULL && "Fail_lt_point2_HT500 not found!");
	SubstractEWK(ht500_fail, "factorization_ht500/Fail_lt_point2_HT500");
	GetPredictions(ht500_fail, gausFit2, vHT500Predictions);
	GetPredictions(ht500_fail, HTsyst_gausFit, vHT500Predictions_HTsyst);

	cout << "\n\n >>>>>> HT>800 PREDICTIONS " << endl; 
	TH1* ht800_fail = dynamic_cast<TH1*> (dataRootFile->Get("factorization_ht500/Fail_lt_point2_HT800"));
	assert(ht800_fail != NULL && "Fail_lt_point2_HT800 not found!");
	SubstractEWK(ht800_fail, "factorization_ht500/Fail_lt_point2_HT800");
	GetPredictions(ht800_fail, gausFit2, vHT800Predictions);
	GetPredictions(ht800_fail, HTsyst_gausFit, vHT800Predictions_HTsyst);

	cout << "\n\n >>>>>> HT>1000 PREDICTIONS " << endl; 
	TH1* ht1000_fail = dynamic_cast<TH1*> (dataRootFile->Get("factorization_ht500/Fail_lt_point2_HT1000"));
	assert(ht1000_fail != NULL && "Fail_lt_point2_HT1000 not found!");
	SubstractEWK(ht1000_fail, "factorization_ht500/Fail_lt_point2_HT1000");
	GetPredictions(ht1000_fail, gausFit2, vHT1000Predictions);
	GetPredictions(ht1000_fail, HTsyst_gausFit, vHT1000Predictions_HTsyst);

	cout << "\n\n >>>>>> HT>1200 PREDICTIONS " << endl; 
	TH1* ht1200_fail = dynamic_cast<TH1*> (dataRootFile->Get("factorization_ht500/Fail_lt_point2_HT1200"));
	assert(ht1200_fail != NULL && "Fail_lt_point2_HT1200 not found!");
	SubstractEWK(ht1200_fail, "factorization_ht500/Fail_lt_point2_HT1200");
	GetPredictions(ht1200_fail, gausFit2, vHT1200Predictions);
	GetPredictions(ht1200_fail, HTsyst_gausFit, vHT1200Predictions_HTsyst);

*/	cout << "\n\n >>>>>> HT>1400 PREDICTIONS " << endl; 
	const string ht1400_fail_name("factorization_ht500/Fail_lt_point2_HT1400");
	TH1* ht1400_fail    = dynamic_cast<TH1*> (dataRootFile->Get(ht1400_fail_name.c_str()));
	assert(ht1400_fail != NULL && "Fail_lt_point2_HT1400 not found!");
	SubstractEWK  (ht1400_fail, ht1400_fail_name.c_str());
	SubstractTTbar(ht1400_fail, ht1400_fail_name.c_str());
	SubstractZnn  (ht1400_fail, ht1400_fail_name.c_str());
	Predictions_t pred_HT1400        = GetPredictions(ht1400_fail, gausFit2);
	Predictions_t pred_HT1400_htSyst = GetPredictions(ht1400_fail, HTsyst_gausFit);

	TH1* qcdsyst1_ht1400_fail = dynamic_cast<TH1*> (dataRootFile->Get("factorization_ht500/Syst1_Fail_HT1400_fineBin"));
	assert(qcdsyst1_ht1400_fail != NULL && "Syst1_Fail_HT1400_fineBin not found!");
	Predictions_t pred_HT1400_qcdSyst1 = GetPredictions(qcdsyst1_ht1400_fail, QCDsyst1_gausFit);

	Predictions_t pred_HT1400_Cup_syst   = GetPredictions(ht1400_fail, gausFit_Cup);
	Predictions_t pred_HT1400_Cdown_syst = GetPredictions(ht1400_fail, gausFit_Cdown);

	//exclsuive HT bins
	/***   500 < HT < 800 ******/
	cout << "\n\n >>>>>> 500<HT<800 PREDICTIONS " << endl; 
	const string ht500to800_fail_name("factorization_ht500/Fail_lt_point2_500HT800");
	TH1* ht500to800_fail    = dynamic_cast<TH1*> (dataRootFile->Get(ht500to800_fail_name.c_str()));
	assert(ht500to800_fail != NULL && "Fail_lt_point2_500HT800 not found!");
	SubstractEWK  (ht500to800_fail, ht500to800_fail_name);
	SubstractTTbar(ht500to800_fail, ht500to800_fail_name);
	SubstractZnn  (ht500to800_fail, ht500to800_fail_name);
	Predictions_t pred_HT500to800        = GetPredictions(ht500to800_fail, gausFit2);
	Predictions_t pred_HT500to800_htSyst = GetPredictions(ht500to800_fail, HTsyst_gausFit);

	TH1* qcdsyst1_ht500to800_fail    = dynamic_cast<TH1*> (dataRootFile->Get("factorization_ht500/Syst1_Fail_500HT800_fineBin"));
	assert(qcdsyst1_ht500to800_fail != NULL && "Syst1_Fail_500HT800 not found!");
	Predictions_t pred_HT500to800_qcdSyst1 = GetPredictions(qcdsyst1_ht500to800_fail, QCDsyst1_gausFit);
	
	TH1* qcdsyst2_ht500to800_fail    = dynamic_cast<TH1*> (dataRootFile->Get("factorization_ht500/Syst2_Fail_500HT800_fineBin"));
	assert(qcdsyst2_ht500to800_fail != NULL && "Syst2_Fail_500HT800 not found!");
	Predictions_t pred_HT500to800_qcdSyst2 = GetPredictions(qcdsyst2_ht500to800_fail, QCDsyst2_gausFit);

	cout << "QCD SYST 1 AND SYST2 " << std::endl;
	//PrintExclPredictions(pred_HT500to800_qcdSyst1);
	//PrintExclPredictions(pred_HT500to800_qcdSyst2);


	Predictions_t pred_HT500to800_Cup_syst   = GetPredictions(ht500to800_fail, gausFit_Cup);
	Predictions_t pred_HT500to800_Cdown_syst = GetPredictions(ht500to800_fail, gausFit_Cdown);


	cout << "\n\n >>>>>> 800<HT<1000 PREDICTIONS " << endl; 
	const string ht800to1000_fail_name("factorization_ht500/Fail_lt_point2_800HT1000"); 
	TH1* ht800to1000_fail = dynamic_cast<TH1*> (dataRootFile->Get(ht800to1000_fail_name.c_str()));
	assert(ht800to1000_fail != NULL && "Fail_lt_point2_800HT1000 not found!");
	SubstractEWK  (ht800to1000_fail, ht800to1000_fail_name);
	SubstractTTbar(ht800to1000_fail, ht800to1000_fail_name);
	SubstractZnn  (ht800to1000_fail, ht800to1000_fail_name);
	Predictions_t pred_HT800to1000        = GetPredictions(ht800to1000_fail, gausFit2);
	Predictions_t pred_HT800to1000_htSyst = GetPredictions(ht800to1000_fail, HTsyst_gausFit);

	TH1* qcdsyst1_ht800to1000_fail = dynamic_cast<TH1*> (dataRootFile->Get("factorization_ht800/Syst1_Fail_800HT1000_fineBin"));
	assert(qcdsyst1_ht800to1000_fail       != NULL && "Syst1_Fail_800HT1000 not found!");
	Predictions_t pred_HT800to1000_qcdSyst1 = GetPredictions(qcdsyst1_ht800to1000_fail, QCDsyst1_gausFit);

	Predictions_t pred_HT800to1000_Cup_syst   = GetPredictions(ht800to1000_fail, gausFit_Cup);
	Predictions_t pred_HT800to1000_Cdown_syst = GetPredictions(ht800to1000_fail, gausFit_Cdown);


	cout << "\n\n >>>>>> 1000<HT<1200 PREDICTIONS " << endl; 
	const string ht1000to1200_fail_name("factorization_ht500/Fail_lt_point2_1000HT1200");
	TH1* ht1000to1200_fail = dynamic_cast<TH1*> (dataRootFile->Get(ht1000to1200_fail_name.c_str()));
	assert(ht1000to1200_fail != NULL && "Fail_lt_point2_1000HT1200 not found!");
	SubstractEWK  (ht1000to1200_fail, ht1000to1200_fail_name);
	SubstractTTbar(ht1000to1200_fail, ht1000to1200_fail_name);
	SubstractZnn  (ht1000to1200_fail, ht1000to1200_fail_name);
	Predictions_t pred_HT1000to1200        = GetPredictions(ht1000to1200_fail, gausFit2);
	Predictions_t pred_HT1000to1200_htSyst = GetPredictions(ht1000to1200_fail, HTsyst_gausFit);

	TH1* qcdsyst1_ht1000to1200_fail = dynamic_cast<TH1*> (dataRootFile->Get("factorization_ht1000/Syst1_Fail_1000HT1200_fineBin"));
	assert(qcdsyst1_ht1000to1200_fail != NULL && "Syst1_Fail_1000HT1200 not found!");
	Predictions_t pred_HT1000to1200_qcdSyst1 = GetPredictions(qcdsyst1_ht1000to1200_fail, QCDsyst1_gausFit);

	Predictions_t pred_HT1000to1200_Cup_syst = GetPredictions(ht1000to1200_fail, gausFit_Cup);
	Predictions_t pred_HT1000to1200_Cdown_syst = GetPredictions(ht1000to1200_fail, gausFit_Cdown);



	cout << "\n\n >>>>>> 1200<HT<1400 PREDICTIONS " << endl; 
	const string ht1200to1400_fail_name("factorization_ht500/Fail_lt_point2_1200HT1400");
	TH1* ht1200to1400_fail = dynamic_cast<TH1*> (dataRootFile->Get(ht1200to1400_fail_name.c_str()));
	assert(ht1200to1400_fail != NULL && "Fail_lt_point2_1200HT1400 not found!");
	SubstractEWK  (ht1200to1400_fail, ht1200to1400_fail_name);
	SubstractTTbar(ht1200to1400_fail, ht1200to1400_fail_name);
	SubstractZnn  (ht1200to1400_fail, ht1200to1400_fail_name);
	Predictions_t pred_HT1200to1400        = GetPredictions(ht1200to1400_fail, gausFit2);
	Predictions_t pred_HT1200to1400_htSyst = GetPredictions(ht1200to1400_fail, HTsyst_gausFit);

	TH1* qcdsyst1_ht1200to1400_fail = dynamic_cast<TH1*> (dataRootFile->Get("factorization_ht1200/Syst1_Fail_1200HT1400_fineBin"));
	assert(qcdsyst1_ht1200to1400_fail != NULL && "Syst1_Fail_1200HT1400 not found!");
	Predictions_t pred_HT1200to1400_qcdSyst1 = GetPredictions(qcdsyst1_ht1200to1400_fail, QCDsyst1_gausFit);

	Predictions_t pred_HT1200to1400_Cup_syst = GetPredictions(ht1200to1400_fail, gausFit_Cup);
	Predictions_t pred_HT1200to1400_Cdown_syst = GetPredictions(ht1200to1400_fail, gausFit_Cdown);


	cout << "\n\n >>>>>> 500<HT<800 PREDICTIONS " << endl; 
	PrintResults(pred_HT500to800, pred_HT500to800_htSyst, pred_HT500to800_qcdSyst1, pred_HT500to800_Cup_syst, pred_HT500to800_Cdown_syst);
	cout << "\n >>>>>> 800<HT<1000 PREDICTIONS " << endl; 
	PrintResults(pred_HT800to1000, pred_HT800to1000_htSyst, pred_HT800to1000_qcdSyst1, pred_HT800to1000_Cup_syst, pred_HT800to1000_Cdown_syst);
	cout << "\n >>>>>> 1000<HT<1200 PREDICTIONS " << endl; 
	PrintResults(pred_HT1000to1200, pred_HT1000to1200_htSyst, pred_HT1000to1200_qcdSyst1, pred_HT1000to1200_Cup_syst, pred_HT1000to1200_Cdown_syst);
	cout << "\n >>>>>> 1200<HT<1400 PREDICTIONS " << endl; 
	PrintResults(pred_HT1200to1400, pred_HT1200to1400_htSyst, pred_HT1200to1400_qcdSyst1, pred_HT1200to1400_Cup_syst, pred_HT1200to1400_Cdown_syst);
	cout << "\n >>>>>> HT>1400 PREDICTIONS " << endl; 
	PrintResults(pred_HT1400, pred_HT1400_htSyst, pred_HT1400_qcdSyst1, pred_HT1400_Cup_syst, pred_HT1400_Cdown_syst);


}


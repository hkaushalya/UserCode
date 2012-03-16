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
const static float fZNN_LUMI   = 3067017.0/42.8;     //x-sec from note, Anwar hd used 32.92 pb-1 LO
																	  //2-23-12, Seema said this should be 42.8 not 42.2
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
//float CONST_C_UP   = 0.07;  //2010 values
//float CONST_C_DOWN = 0.012; //2010 values
float CONST_C_UP   = 0.0304145;
float CONST_C_DOWN = 0.0304145;


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
	//newtitle << "Fit Range " << passFailRatio_fitrange_xmin << "--" << passFailRatio_fitrange_xmax << " GeV " << title;
	newtitle << title;
	//newtitle << ";MHT;Ratio (r)";
	cout << __FUNCTION__ << ":: Fitting Range = " << passFailRatio_fitrange_xmin << " - " << passFailRatio_fitrange_xmax << endl;

	gStyle->SetOptStat(0);
	//debug stuff
	new TCanvas();
	//gPad->SetGridy();
	gPad->SetTickx();
	gPad->SetTicky();
	if (logScale) gPad->SetLogy();
	//hist->SetLineColor(9);
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
	gausFitResult->SetLineWidth ( 2                       );
	gausFitResult->DrawCopy     ( "same"                  );

	TLegend *leg  = new TLegend(0.65,0.75,0.9,0.9);
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

void makePassFailHtDep_data(const int htBin = 1, const float fitrange_xmin = 50, const float fitrange_xmax = 150) 
{

	TH1 *histPassOrig = 0;
	TH1 *histFailOrig = 0;

	TFile* dataRootFile = new TFile (sDATA_FILE_NAME.c_str());
	if (dataRootFile->IsZombie())
	{
		cout << "Data root file not found!" <<  endl;
		assert (false);
	}

	/**********************************************
	 * Now Get HT dependence syst fit
	 * *******************************************/
	std::string HTsyst_numerHistName;
	std::string HTsyst_denomHistName;

	if (htBin == 1)
	{
		std::string HTsyst_numerHistName("factorization_ht600/Pass_RA2dphi");
		std::string HTsyst_denomHistName("factorization_ht600/Fail_1");
	} else if (htBin == 2)
	{
		std::string HTsyst_numerHistName("factorization_ht600/Pass_RA2dphi");
		std::string HTsyst_denomHistName("factorization_ht600/Fail_1");
	} else if (htBin == 2)
	{
		std::string HTsyst_numerHistName("factorization_ht600/Pass_RA2dphi");
		std::string HTsyst_denomHistName("factorization_ht600/Fail_1");
	} else if (htBin == 2)
	{
		std::string HTsyst_numerHistName("factorization_ht600/Pass_RA2dphi");
		std::string HTsyst_denomHistName("factorization_ht600/Fail_1");
	} else if (htBin == 2)
	{
		std::string HTsyst_numerHistName("factorization_ht600/Pass_RA2dphi");
		std::string HTsyst_denomHistName("factorization_ht600/Fail_1");
	}



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

	const string HTsyst_title = "HT>600 GeV;MHT [GeV];Ratio (r);";
	TF1 *HTsyst_gausFit =  FitPassFailRatio(HTsyst_numerHist, HTsyst_title, "FactorizationHT600");
	if (HTsyst_gausFit == NULL) { cout << __LINE__ << "::HT syst functions null! " << endl; assert (false); }

}


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

using namespace std;

const static float fDATA_LUMI = 10000; //pb-1
const static float passFailRatio_fitrange_xmin = 50.0;
const static float passFailRatio_fitrange_xmax = 150.0;
std::vector<float> incMHTBins;
//const static int nMHTbins = 5;
const static int nMHTbins = 5;
//const static float arrHTbins[nHTbins] = {0,500,800,1000,1200,1400,7000};
const static float arrMHTbins[nMHTbins] = {200,350,500,600,7000};
float CONST_C = 0.0217;
const static bool bFLATSAMPLE = 0;
const static string sQCD_FLAT_FILE_NAME("QCD_Flat_8TeV.root");
const static float fQCD_FLAT_LUMI = 9998154/2.99913994E10; //pb-1
const Int_t nBins = 7;
TFile *files[nBins];

struct Predictions_t
{
		float incl_mean[nMHTbins];
		float incl_statErr[nMHTbins];
		float incl_fitErr[nMHTbins];
		float excl_mean[nMHTbins-1];
		float excl_statErr[nMHTbins-1];
		float excl_fitErr[nMHTbins-1];
		float incl_signal_mean[nMHTbins];
		float incl_signal_statErr[nMHTbins];
		float excl_signal_mean[nMHTbins-1];
		float excl_signal_statErr[nMHTbins-1];
};

void PrintExclResults(const Predictions_t& gaus, const Predictions_t& exp,
							const Predictions_t& gaus_cplus, Predictions_t& exp_cplus)
{
/*	cout << setprecision(3) 
					<< setw(10) << " HT "
					<< setw(15) << "MHT"
					<< setw(15) << " mean " 
					<< setw(20) << "statFitErr"
					<< setw(30) << "Signal+/-stat"
					<< endl;
*/	for (int i=0; i< nMHTbins -1 ; ++i)
	{
		cout << setprecision(4) 
					//<< setw(5) << HTmin << "-" << HTmax << " & "
					<< setw(10) << incMHTBins.at(i) << "-" << incMHTBins.at(i+1) << " & "
					<< setw(15) << gaus.excl_mean[i]   << "& $\\pm$ " << gaus.excl_statErr[i] << " & ("
					<< setw(15) << gaus_cplus.excl_mean[i] << ") & "
					<< setw(15) << exp.excl_mean[i]   << "& $\\pm$ " << exp.excl_statErr[i] << " & ("
					<< setw(15) << exp_cplus.excl_mean[i] << ") & "
					<< setw(15) << gaus.excl_signal_mean[i] << "& $\\pm$ " << gaus.excl_signal_statErr[i]
					<< endl;
	}
}
void PrintExclResults(const Predictions_t& res, const float HTmin=0, const float HTmax=0)
{
	cout << setprecision(3) 
					<< setw(10) << " HT "
					<< setw(15) << "MHT"
					<< setw(15) << " mean " 
					<< setw(20) << "statFitErr"
					<< setw(30) << "Signal+/-stat"
					<< endl;
	for (int i=0; i< nMHTbins -1 ; ++i)
	{
			const float mean          = res.excl_mean[i]; 
			const float statAndFitErr = sqrt(pow(res.excl_statErr[i],2)+pow(res.excl_fitErr[i],2));
			const float signal        = res.excl_signal_mean[i];
			const float signalErr     = res.excl_signal_statErr[i];

		cout << setprecision(4) 
					<< setw(5) << HTmin << "-" << HTmax << " & "
					<< setw(10) << incMHTBins.at(i) << "-" << incMHTBins.at(i+1) << " & "
					<< setw(15) << mean   << "& $\\pm$ " << statAndFitErr << " & "
					<< setw(15) << signal << "& $\\pm$ " << signalErr
					<< endl;
	}
}
void PrintInclResults(const Predictions_t& res, const float HTmin=0)
{
	cout << setprecision(3) 
					<< setw(10) << " HT "
					<< setw(15) << "MHT"
					<< setw(15) << " mean " 
					<< setw(20) << "statFitErr"
					<< setw(30) << "Signal+/-stat"
					<< endl;
	for (int i=0; i< nMHTbins ; ++i)
	{
			const float mean          = res.incl_mean[i]; 
			const float statAndFitErr = sqrt(pow(res.incl_statErr[i],2)+pow(res.incl_fitErr[i],2));
			const float signal        = res.incl_signal_mean[i];
			const float signalErr     = res.incl_signal_statErr[i];

		cout << setprecision(4) 
					<< setw(10) << HTmin << " & "
					<< setw(10) << incMHTBins.at(i) << " & "
					<< setw(15) << mean   << " & $\\pm$ " << statAndFitErr << " & "
					<< setw(15) << signal << " & $\\pm$ " << signalErr
					<< endl;
	}
	cout << "out of " << __FUNCTION__ << endl;
}


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

double GetFitFunctionError(TF1* f1, double x) {

	double err = 0.;
	Double_t val[1] = { x };

	for(Int_t j=0; j<f1->GetNpar(); ++j) 
	{
		err += f1->GradientPar(j,val) * f1->GetParError(j);
	}

	return err;

}

void OpenFiles()
{
		files[0] = new TFile ("QCD_HT300to470.root");
		files[1] = new TFile ("QCD_HT470to600.root");
		files[2] = new TFile ("QCD_HT600to800.root");
		files[3] = new TFile ("QCD_HT800to1000.root");
		files[4] = new TFile ("QCD_Ht1000to1400.root");
		files[5] = new TFile ("QCD_HT1400to1800.root");
		files[6] = new TFile ("QCD_HT1800.root");

		for (int i=0; i<nBins; ++i)
		{
			if (files[i]->IsZombie())
			{
				cout << "QCD file # " << i << " not found!" <<  endl;
				assert (false);
			}
		}
}

TH1* GetHist(const std::string histname, const float scaleTo=1.0)
{
	TH1 *res_hist = 0;

	if (! bFLATSAMPLE)
	{
		TH1 *hists[nBins] = {0,0,0,0,0,0,0};

		for (int i=0; i<nBins; ++i)
		{
			if (files[i]->IsZombie())
			{
				cout << files[i]->GetName() << " not found!" <<  endl;
				assert (false);
			} else
			{
				hists[i] = dynamic_cast<TH1*> (files[i]->Get(histname.c_str()));
				if (hists[i] == 0 )
				{
					cout << "hist_pass " << histname << " not found in " << files[i]->GetName() << "!" << endl;
					assert (false);
				} else 
				{
					hists[i]->Sumw2();

					const float scale = scaleTo/ ( nEvents[i] / xSec[i] );
					hists[i]->Scale(scale);

					if (i == 0) 
					{ 
						res_hist = dynamic_cast<TH1*> (hists[i]->Clone("histcopy")); 
						res_hist->SetDirectory(0);
					} else { res_hist->Add(hists[i]); }
				}
			}
		}

	} else
	{
		TFile f(sQCD_FLAT_FILE_NAME.c_str());
		if (f.IsZombie())
		{
			cout << __FUNCTION__ << ":" << __LINE__ << ":file " << sQCD_FLAT_FILE_NAME << " not found!"<< endl;
			assert(false);
		}
		TH1* hist = dynamic_cast<TH1*> (f.Get(histname.c_str()));
		if (hist == NULL)
		{
			cout << __FUNCTION__ << ":" << __LINE__ << ": " << histname << " not found!"<< endl;
			assert(false);
		}
		res_hist = dynamic_cast<TH1*> (hist->Clone());
		res_hist->SetDirectory(0);
	}

	return res_hist;
}

double expFitFunc_Cup(double *x, double *par)
{
	double fitval=0.0;
	//if(x[0]>0.0) fitval=par[0] * exp(par[1] * x[0]) + 0.03;
	//if(x[0]>0.0) fitval=par[0] * exp(par[1] * x[0]) + 0.0217;
	if(x[0]>0.0) fitval=par[0] * exp(par[1] * x[0]) + (CONST_C * 2);  //100% error for C
	else fitval=-1.0E6;
	return fitval;
}

double gausFitFunc_Cup(double *x, double *par)
{
	double arg = par[0] * exp(par[1] * x[0]);
	double fitval = 1.0 / TMath::Erf (arg) - 1 + (CONST_C * 2); //100% ERROR FOR C
	return fitval;
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
	//double fitval = 1.0 / TMath::Erf (arg) - 1 + 0.0217;
	double fitval = 1.0 / TMath::Erf (arg) - 1 + CONST_C;
	return fitval;
}

Predictions_t GetPredictions(const TH1* hist, TF1* f1, const TH1* signalHist)
{
	assert(hist != NULL && "GetPredictions:: hist not found!");
	assert(f1 != NULL && "GetPredictions:: func not found!");
	assert(signalHist != NULL && "GetPredictions:: signalHist not found!");
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
/*	double sig_sum = 0;
	double sig_err2 =0;
	for (int bin =0; bin<signalHist->GetNbinsX()+1; ++bin) 
	{
		if (signalHist->GetBinContent(bin)>0)
		{
			sig_sum += signalHist->GetBinContent(bin);
			sig_err2 += pow(signalHist->GetBinError(bin),2);
		}
	}
	std::cout << "Signal hist  := " << sig_sum << " +/- " << sqrt(sig_err2) << std::endl;
*/
	double sumGaus[incMHTBins.size()];
	double Gaus_StatErr[incMHTBins.size()];
	double Gaus_FitErr[incMHTBins.size()];
	double sumGaus_excl[incMHTBins.size()-1];
	double Gaus_StatErr_excl[incMHTBins.size()-1];
	double Gaus_FitErr_excl[incMHTBins.size()-1];
	double signal_mean_excl[incMHTBins.size()-1];
	double signal_statErr_excl[incMHTBins.size()-1];
	double signal_mean_incl[incMHTBins.size()];
	double signal_statErr_incl[incMHTBins.size()];

	Predictions_t results;

	for (unsigned mhtBin = 0; mhtBin < incMHTBins.size(); ++mhtBin)
	{
		sumGaus[mhtBin]      = 0;
		Gaus_StatErr[mhtBin] = 0;
		Gaus_FitErr[mhtBin]  = 0;
		signal_mean_incl[mhtBin]  = 0;
		signal_statErr_incl[mhtBin] = 0;
		if (mhtBin+1<incMHTBins.size()) 
		{
			sumGaus_excl[mhtBin]      = 0;
			Gaus_StatErr_excl[mhtBin] = 0;
			Gaus_FitErr_excl[mhtBin]  = 0;
			signal_mean_excl[mhtBin]  = 0;
			signal_statErr_excl[mhtBin] = 0;
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
				const float binSig    = signalHist->GetBinContent(bin); 
				const float binSigStatErr2 = pow(signalHist->GetBinError(bin),2); 

			//cout << setprecision(3) << setw(15) << "MHT > " << incMHTBins.at(mhtBin) 
			//<< setw(20) << binSig  << "&$\\pm$" << binSigStatErr2 << std::endl;

				if (hist->GetBinCenter(bin) > incMHTBins.at(mhtBin))
				{
					sumGaus[mhtBin]      += res;
					Gaus_StatErr[mhtBin] += statErr2;
					Gaus_FitErr[mhtBin]  += fitErr;
					signal_mean_incl[mhtBin]  += binSig;
					signal_statErr_incl[mhtBin]  += binSigStatErr2;
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
						signal_mean_excl[mhtBin]  += binSig;
						signal_statErr_excl[mhtBin]  += binSigStatErr2;
					}
				}
			}
		}

		Gaus_StatErr[mhtBin] = sqrt(Gaus_StatErr[mhtBin]);
		Gaus_StatErr_excl[mhtBin] = sqrt(Gaus_StatErr_excl[mhtBin]);
		signal_statErr_incl[mhtBin] = sqrt(signal_statErr_incl[mhtBin]);
		signal_statErr_excl[mhtBin] = sqrt(signal_statErr_excl[mhtBin]);

		results.incl_mean[mhtBin] = sumGaus[mhtBin];
		results.incl_statErr[mhtBin] = Gaus_StatErr[mhtBin];
		results.incl_fitErr[mhtBin] = Gaus_FitErr[mhtBin];
		results.incl_signal_mean[mhtBin] = signal_mean_incl[mhtBin];
		results.incl_signal_statErr[mhtBin] = signal_statErr_incl[mhtBin];

		
		if (mhtBin+1<incMHTBins.size())
		{
			std::stringstream excl_pred;
			excl_pred << setprecision(3) << setw(10) << incMHTBins.at(mhtBin) << "<MHT<" << incMHTBins.at(mhtBin+1)
				<< setw(20) << sumGaus_excl[mhtBin]  << "&$\\pm$" << Gaus_StatErr_excl[mhtBin] << " &$\\pm$ " << Gaus_FitErr_excl[mhtBin];
			results.excl_mean[mhtBin]    = sumGaus_excl[mhtBin];
			results.excl_statErr[mhtBin] = Gaus_StatErr_excl[mhtBin];
			results.excl_fitErr[mhtBin]  = Gaus_FitErr_excl[mhtBin];
			results.excl_signal_mean[mhtBin] = signal_mean_excl[mhtBin];
			results.excl_signal_statErr[mhtBin] = signal_statErr_excl[mhtBin];
		}
	}
	
	return results;
}

void makePassFail_QCDMC(const string numeHistName, const string denoHistName, 
			const string title, const string HTbinlabel, 
			const string signalHistName, const string controlHistName,
			const float fitrange_xmin = 50, const float fitrange_xmax = 150) 
{
	incMHTBins.clear();
	for (int i=0; i<nMHTbins; ++i) incMHTBins.push_back(arrMHTbins[i]);
	const bool logScale = true;
	const float scaleTo = fDATA_LUMI; // pb-1

	TH1 *Hist_pass = GetHist(numeHistName, scaleTo);
	TH1 *Hist_fail = GetHist(denoHistName, scaleTo);

	//new TCanvas();
	//Hist_pass->Draw();
	//Hist_fail->Draw("same");
	cout << __LINE__ << ": Integrals= " << Hist_fail->Integral() << "/" << Hist_pass->Integral() << endl;
	if (Hist_fail->GetEntries()<1 || Hist_pass->GetEntries()<1) 
	{
		cout << __LINE__ << ": not enough intries to make the plots!!! " << Hist_fail->Integral() << "/" << Hist_pass->Integral() << endl;
		return;
	}
/*
	new TCanvas();
	gPad->SetLogy();
	gStyle->SetOptStat(0);
	Hist_fail->SetLineColor(kBlue);
	Hist_pass->SetLineColor(kRed);
	Hist_fail->SetTitle(";MHT;");
	Hist_fail->Draw();
	Hist_pass->Draw("Psame");

	TLegend *leg1  = new TLegend(0.7,0.8,0.9,0.9);
	leg1->AddEntry(Hist_pass,"PASS");
	leg1->AddEntry(Hist_fail,"FAIL");
	leg1->Draw();
	return;
*/
	//gPad->SetEditable(0);
	//new TCanvas();
//	new TCanvas(); gPad->SetLogy(); Hist_pass->DrawCopy();
//	new TCanvas(); gPad->SetLogy(); Hist_fail->DrawCopy();
//	new TCanvas(); gPad->SetLogy(); Hist_signal_mht200->DrawCopy();
//	return;

	Hist_pass->Divide(Hist_fail);
	const int maxbin = Hist_pass->GetMaximumBin();
	const double max = Hist_pass->GetBinContent(maxbin);
	//Hist_pass->GetYaxis()->SetRangeUser(-0.05, max+0.05);


	//fit range
	//const float fitrange_xmin = 50, fitrange_xmax = 120;
	stringstream newtitle;
	//newtitle << "Fit Range " << fitrange_xmin << "--" << fitrange_xmax << title;
	newtitle << title << " (fit range = " <<  fitrange_xmin << "--" << fitrange_xmax << "), c = " << CONST_C << ";MHT [GeV];Ratio (r);";


	gStyle->SetOptStat(0);	
	new TCanvas();
	//gPad->SetGridy();
	gPad->SetTickx();
	if (logScale) gPad->SetLogy();
	Hist_pass->SetLineColor(9);
	Hist_pass->SetTitle(newtitle.str().c_str());
	Hist_pass->SetLineWidth(2);

	gStyle->SetOptFit(1);
	//Hist_pass->SetStats(0);
	Hist_pass->SetMinimum(0.5E-2);
	Hist_pass->Draw();
	//return;

	//do fittings exp and gaus

	TF1 *expFit=new TF1("fit_1",expFitFunc,fitrange_xmin,fitrange_xmax,2);
	expFit->SetParameter(0,0.09);
	expFit->SetParameter(1,-0.0002);
	//expFit->SetParameter(2,-1.0);
	Hist_pass->Fit(expFit,"E0","",fitrange_xmin, fitrange_xmax);
	gPad->Modified();
	gPad->Update();
	TPaveStats *e_stats = (TPaveStats*) Hist_pass->FindObject("stats");
	e_stats->SetTextColor(kRed);
	TPaveStats *exp_stats = (TPaveStats*) e_stats->Clone("exp_stats");

	TF1 *expFit2=new TF1("fit_2",expFitFunc,50,1000.0,2);
	expFit2->SetParameter(0,expFit->GetParameter(0));
	expFit2->SetParameter(1,expFit->GetParameter(1));
	//expFit2->SetParameter(2,expFit->GetParameter(2));
	expFit2->SetLineColor(kRed+1);
	expFit2->SetLineWidth(2);

	//for unceratinty on c
	TF1 *expFit_sigma = new TF1("Exp_sigma",expFitFunc_Cup,50,1000.0,2);
	expFit_sigma->SetParameter(0,expFit->GetParameter(0));
	expFit_sigma->SetParameter(1,expFit->GetParameter(1));
	Hist_pass->Fit(expFit_sigma,"E0","",fitrange_xmin, fitrange_xmax);


	TF1 *gausFit=new TF1("fit_3",gausFitFunc,fitrange_xmin, fitrange_xmax,2);
	gausFit->SetParameter(0,0.09); 
	gausFit->SetParameter(1,-0.0002);
	//gausFit->SetParameter(2,-1.0);
	//gausFit->SetParLimits(2,0.01,0.02);
	gausFit->SetLineColor(kGreen-2);
	Hist_pass->Fit(gausFit,"E0","", fitrange_xmin, fitrange_xmax);
	gPad->Modified();
	gPad->Update();
	TPaveStats *gaus_stats = (TPaveStats*) Hist_pass->FindObject("stats");
	gaus_stats->SetTextColor(kGreen);
	//TPaveStats *gaus_stats = (TPaveStats*) g_stats->Clone("gaus_stats");

	TF1 *gausFit2=new TF1("fit_4",gausFitFunc,50,1000.0,2);
	gausFit2->SetParameter(0,gausFit->GetParameter(0)); 
	gausFit2->SetParameter(1,gausFit->GetParameter(1));
	//gausFit2->SetParameter(2,gausFit->GetParameter(2));
	gausFit2->SetLineColor(kGreen);
	gausFit2->SetLineWidth(2);

	//for unceratinty on c
	TF1 *gausFit_sigma=new TF1("Gaus_sigma",gausFitFunc_Cup,50,1000.0,2);
	gausFit_sigma->SetParameter(0,gausFit->GetParameter(0)); 
	gausFit_sigma->SetParameter(1,gausFit->GetParameter(1));
	Hist_pass->Fit(gausFit_sigma,"E0","", fitrange_xmin, fitrange_xmax);




	gausFit2->Draw("same");
	expFit2->Draw("same");

	//TLegend *leg  = new TLegend(0.7,0.8,0.9,0.9);
	TLegend *leg  = new TLegend(0.7,0.7,0.9,0.9);
	leg->AddEntry(gausFit2,"Gaussian");
	leg->AddEntry(expFit2,"Exponential");
	leg->SetTextFont(42);
	leg->Draw();

	const float xmin=0.2, xmax=0.45, ymin=0.7, ymax=0.9;
	gaus_stats->SetX1NDC(xmin);
	gaus_stats->SetX2NDC(xmax);
	gaus_stats->SetY1NDC(ymin);
	gaus_stats->SetY2NDC(ymax);
	gaus_stats->Draw("same");
	exp_stats->SetX1NDC(xmax);
	exp_stats->SetX2NDC(xmax+0.25);
	exp_stats->SetY1NDC(ymin);
	exp_stats->SetY2NDC(ymax);
	exp_stats->Draw("same");

	//fit quality
	stringstream gaus_fit_res, exp_fit_res;
	const float gausFit_goodness =  gausFit->GetChisquare()/gausFit->GetNDF();
	gaus_fit_res << "#chi^{2}/ndof = " << gausFit->GetChisquare()
			 << "/"<< gausFit->GetNDF() << " = " << gausFit_goodness;
	const float expFit_goodness =  expFit->GetChisquare()/expFit->GetNDF();
	exp_fit_res << "#chi^{2}/ndof = " << expFit->GetChisquare()
			 << "/"<< expFit->GetNDF() << " = " << expFit_goodness;

	cout << gaus_fit_res.str() << endl;
	cout << exp_fit_res.str() << endl;

	TPaveText *pt1 = new TPaveText(0.5,0.8,0.7,0.9);
	pt1->AddText(gaus_fit_res.str().c_str());
	pt1->AddText(exp_fit_res.str().c_str());
//	pt1->Draw("same");

	//stringstream epsname;
	//epsname << "factnomht/HTbin" << HTbin << "_fitrange_" << fitrange_xmin << "to" << fitrange_xmax << ".eps";
	//epsname << "factnomht_" << HTbinlabel << ".eps";
	//gPad->Print(epsname.str().c_str());
	gPad->Print("ratios.eps");

	 //make a predicion
	TH1 *sigHist      = GetHist(signalHistName, scaleTo);
	TH1 *controlHist  = GetHist(controlHistName, scaleTo);
	Predictions_t pred_gaus       = GetPredictions(controlHist, gausFit2, sigHist); 
	Predictions_t pred_gaus_sigma = GetPredictions(controlHist, gausFit_sigma, sigHist); 
	Predictions_t pred_exp        = GetPredictions(controlHist, expFit2, sigHist); 
	Predictions_t pred_exp_sigma  = GetPredictions(controlHist, expFit_sigma, sigHist); 

	//PrintExclResults(pred_gaus);
	//PrintExclResults(pred_exp);
	PrintExclResults(pred_gaus, pred_exp, pred_gaus_sigma, pred_exp_sigma);

}

void makePassFail_QCDMC() 
{
	vector<string> folders, htbinlabels;
/*	folders.push_back("factnomht");
	folders.push_back("factnomhtnjet4");
	folders.push_back("factnomhtnjet5");
	folders.push_back("factnomhtnjet6");
	folders.push_back("factnomhtnjet7");
	folders.push_back("factnomhtnjet8");
*/
	folders.push_back("factnomhtnjet2to3");
	folders.push_back("factnomhtnjet4to5");
	folders.push_back("factnomhtnjet6to7");
	folders.push_back("factnomhtnjet8");

	htbinlabels.push_back("HT500to800");
	htbinlabels.push_back("HT800to1000");
	htbinlabels.push_back("HT1000to1200");
	htbinlabels.push_back("HT1200to1400");
	htbinlabels.push_back("HT1400to7000");


	OpenFiles();
	TCanvas *c = new TCanvas("print");	
	c->Draw();
	c->Print("ratios.eps[");
	
	for (unsigned i = 0; i < folders.size(); ++i)
	{	
		string njet("");
		if (i==0) njet += "[2-3]";
		else if (i==1) njet += "[4-5]";
		else if (i==2) njet += "[6-7]";
		else if (i==3) njet += "#geq 8";
		//else if (i==4) njet += "7";
		//else if (i==5) njet += "8";
		
		cout << " >>>>> " << njet << " <<<<< " << endl;

		for (unsigned j = 0; j < htbinlabels.size(); ++j)
		{
			string htrange("");
			if (j==0) htrange += "500<HT<800 GeV";
			else if (j==1) htrange += "800<HT<1000 GeV";
			else if (j==2) htrange += "1000<HT<1200 GeV";
			else if (j==3) htrange += "1200<HT<1400 GeV";
			else if (j==4) htrange += "HT>1400 GeV";

			cout << " >>>>>>>>>>> " << htrange << endl; 
			stringstream title, numeHistName, denoHistName, signalHistName, controlHistName;
			//title << htrange << ";#slash{H}_{T};Ratio (r) = Pass(RA2 dPhi cuts) / Fail(#Delta #phi_{min}< 0.2);";
			title << "Njet " << njet << ", " << htrange;
			numeHistName << folders.at(i) << "/" << htbinlabels.at(j) << "/signal";
			denoHistName << folders.at(i) << "/" << htbinlabels.at(j) << "/fail1";
			signalHistName << folders.at(i) << "/" << htbinlabels.at(j) << "/signalFineBin";
			controlHistName << folders.at(i) << "/" << htbinlabels.at(j) << "/failFineBin1";
			makePassFail_QCDMC(numeHistName.str(), denoHistName.str(), title.str(), htbinlabels.at(j), signalHistName.str(), controlHistName.str()); 
		}
	}
	c->Print("ratios.eps]");
}

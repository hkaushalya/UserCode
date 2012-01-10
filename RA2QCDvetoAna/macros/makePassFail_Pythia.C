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

const static float fDATA_LUMI = 4650; //pb-1
const static float passFailRatio_fitrange_xmin = 50.0;
const static float passFailRatio_fitrange_xmax = 150.0;
std::vector<float> incMHTBins;
const static int nMHTbins = 5;
const static float arrMHTbins[nMHTbins] = {200,350,500,600,7000};
float CONST_C = 0.0217;

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

double GetFitFunctionError(TF1* f1, double x) {

	double err = 0.;
	Double_t val[1] = { x };

	for(Int_t j=0; j<f1->GetNpar(); ++j) 
	{
		err += f1->GradientPar(j,val) * f1->GetParError(j);
	}

	return err;

}


TH1* GetHist(const std::string histname, const float scaleTo=1.0)
{
	const Int_t nBins = 9;
	TFile *files[nBins];
	TH1 *hists[nBins] = {0,0,0,0,0,0,0,0,0};

	files[0] = new TFile ("qcd1/Merged.root");
	files[1] = new TFile ("qcd2/Merged.root");
	files[2] = new TFile ("qcd3/Merged.root");
	files[3] = new TFile ("qcd4/Merged.root");
	files[4] = new TFile ("qcd5/Merged.root");
	files[5] = new TFile ("qcd6/Merged.root");
	files[6] = new TFile ("qcd7/Merged.root");
	files[7] = new TFile ("qcd8/Merged.root");
	files[8] = new TFile ("qcd9/Merged.root");

	TH1 *res_hist = 0;

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

				if (i == 0) { res_hist = dynamic_cast<TH1*> (hists[i]->Clone("histcopy")); } 
				else { res_hist->Add(hists[i]); }
			}
		}
	}
	
	return res_hist;

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

void makePassFail_Pythia(const int HTbin = 1, const float fitrange_xmin = 50, const float fitrange_xmax = 150) 
{

	incMHTBins.clear();
	for (int i=0; i<nMHTbins; ++i) incMHTBins.push_back(arrMHTbins[i]);
	const bool logScale = true;
	const float scaleTo = 4650; // pb

	string title("");
	//const std::string histname1("factorization_ht500/Pass_RA2dphi_HT500");
	std::string numeHistName("");
	std::string denoHistName("");

	if (HTbin == 1)   //500gev
	{
		title +=" : HT>500 GeV;#slash{H}_{T};Ratio (r) = Pass(RA2 dPhi cuts) / Fail(#Delta #phi_{min}< 0.2);";
		numeHistName += "factorization_ht500/Pass_RA2dphi";
		denoHistName += "factorization_ht500/Fail_1";
	} else if (HTbin == 2)  //600gev
	{
		title +=" : HT>600 GeV;#slash{H}_{T};Ratio (r) = Pass(RA2 dPhi cuts) / Fail(#Delta #phi_{min}< 0.2);";
		numeHistName += "factorization_ht600/Pass_RA2dphi";
		denoHistName += "factorization_ht600/Fail_1";
	} else if (HTbin == 3)  //800gev
	{
		title +=" : HT>800 GeV;#slash{H}_{T};Ratio (r) = Pass(RA2 dPhi cuts) / Fail(#Delta #phi_{min}< 0.2);";
		numeHistName += "factorization_ht800/Pass_RA2dphi";
		denoHistName += "factorization_ht800/Fail_1";
	} else if (HTbin == 4)  //1000gev
	{
		title +=" : HT>1000 GeV;#slash{H}_{T};Ratio (r) = Pass(RA2 dPhi cuts) / Fail(#Delta #phi_{min}< 0.2);";
		numeHistName += "factorization_ht1000/Pass_RA2dphi";
		denoHistName += "factorization_ht1000/Fail_1";
	} else if (HTbin == 5)  //1200gev
	{
		title +=" : HT>1200 GeV;#slash{H}_{T};Ratio (r) = Pass(RA2 dPhi cuts) / Fail(#Delta #phi_{min}< 0.2);";
		numeHistName += "factorization_ht1200/Pass_RA2dphi";
		denoHistName += "factorization_ht1200/Fail_1";
	}

	TH1 *Hist_pass = GetHist(numeHistName, scaleTo);
	TH1 *Hist_fail = GetHist(denoHistName, scaleTo);
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
	newtitle << "MC Smeared Jets (Fit Range " << fitrange_xmin << "--" << fitrange_xmax << " GeV)";
	//<< title;
	newtitle << ";MHT [GeV];Ratio (r);";


	gStyle->SetOptStat(0);
	//debug stff
	new TCanvas();
	//gPad->SetGridy();
	gPad->SetTickx();
	gPad->SetTicky();
	if (logScale) gPad->SetLogy();
	Hist_pass->SetLineColor(9);
	Hist_pass->SetTitle(newtitle.str().c_str());
	Hist_pass->SetLineWidth(2);

	gStyle->SetOptFit(1);
	//Hist_pass->SetStats(0);
	Hist_pass->Draw();

	//do fittings exp and gaus

	cout << __LINE__ << endl;
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
	cout << __LINE__ << endl;

	TF1 *expFit2=new TF1("fit_2",expFitFunc,50,1000.0,2);
	expFit2->SetParameter(0,expFit->GetParameter(0));
	expFit2->SetParameter(1,expFit->GetParameter(1));
	//expFit2->SetParameter(2,expFit->GetParameter(2));
	expFit2->SetLineColor(kRed);
	expFit2->SetLineWidth(1);
	cout << __LINE__ << endl;


	TF1 *gausFit=new TF1("fit_3",gausFitFunc,fitrange_xmin, fitrange_xmax,2);
	gausFit->SetParameter(0,0.09); 
	gausFit->SetParameter(1,-0.0002);
	gausFit->SetParameter(2,-1.0);
	//gausFit->SetParLimits(2,0.01,0.02);
	gausFit->SetLineColor(kGreen);
	Hist_pass->Fit(gausFit,"E0","", fitrange_xmin, fitrange_xmax);
	gPad->Modified();
	gPad->Update();
	TPaveStats *gaus_stats = (TPaveStats*) Hist_pass->FindObject("stats");
	gaus_stats->SetTextColor(kGreen);
	//TPaveStats *gaus_stats = (TPaveStats*) g_stats->Clone("gaus_stats");

	TF1 *gausFit2=new TF1("fit_4",gausFitFunc,50,1000.0,2);
	gausFit2->SetParameter(0,gausFit->GetParameter(0)); 
	gausFit2->SetParameter(1,gausFit->GetParameter(1));
	gausFit2->SetParameter(2,gausFit->GetParameter(2));
	gausFit2->SetLineColor(kGreen);
	gausFit2->SetLineWidth(1);
	gausFit2->Draw("same");
	expFit2->Draw("same");

	//TLegend *leg  = new TLegend(0.7,0.8,0.9,0.9);
	TLegend *leg  = new TLegend(0.65,0.8,0.9,0.9);
	leg->AddEntry(gausFit2,"Gaus Model");
	leg->AddEntry(expFit2,"Exp Model");
	leg->Draw();

	const float xmin=0.3, xmax=0.65, ymin=0.7, ymax=0.9;
	gaus_stats->SetX1NDC(xmin);
	gaus_stats->SetX2NDC(xmax);
	gaus_stats->SetY1NDC(ymin);
	gaus_stats->SetY2NDC(ymax);
	gaus_stats->Draw("same");
	//exp_stats->SetX1NDC(xmax);
	//exp_stats->SetX2NDC(xmax+0.25);
	//exp_stats->SetY1NDC(ymin);
	//exp_stats->SetY2NDC(ymax);
	//exp_stats->Draw("same");


	//fit quality
	stringstream gaus_fit_res, exp_fit_res;
	const float gausFit_goodness =  gausFit->GetChisquare()/gausFit->GetNDF();
	gaus_fit_res << "#chi^{2}/ndof = " << gausFit->GetChisquare()
			 << "/"<< gausFit->GetNDF() << " = " << gausFit_goodness;
	const float expFit_goodness =  expFit->GetChisquare()/expFit->GetNDF();
	exp_fit_res << "#chi^{2}/ndof = " << expFit->GetChisquare()
			 << "/"<< expFit->GetNDF() << " = " << expFit_goodness;

	TPaveText *pt1 = new TPaveText(0.5,0.8,0.7,0.9);
	pt1->AddText(gaus_fit_res.str().c_str());
	pt1->AddText(exp_fit_res.str().c_str());
//	pt1->Draw("same");

	stringstream epsname;
	epsname << "factorization_HTbin" << HTbin << "_fitrange_" << fitrange_xmin << "to" << fitrange_xmax << ".eps";
	gPad->Print(epsname.str().c_str());
	/*return;

	 //make a predicion
	cout << "\n\n >>>>>> HT>500 PREDICTIONS " << endl; 
	TH1 *sigHist_ht500       = GetHist("factorization_ht500/Signal_HT500MHT200", scaleTo);
	TH1 *sidebandHist_ht500  = GetHist("factorization_ht500/Fail_lt_point2_HT500", scaleTo);
	Predictions_t pred_HT500 = GetPredictions(sidebandHist_ht500, gausFit2, sigHist_ht500); 

	cout << "\n\n >>>>>> HT>800 PREDICTIONS " << endl; 
	TH1 *sigHist_ht800       = GetHist("factorization_ht500/Signal_HT800MHT200", scaleTo);
	TH1 *sidebandHist_ht800  = GetHist("factorization_ht500/Fail_lt_point2_HT800", scaleTo);
	Predictions_t pred_HT800 = GetPredictions(sidebandHist_ht800, gausFit2, sigHist_ht800); 

	cout << "\n\n >>>>>> HT>1000 PREDICTIONS " << endl; 
	TH1 *sigHist_ht1000       = GetHist("factorization_ht500/Signal_HT1000MHT200", scaleTo);
	TH1 *sidebandHist_ht1000  = GetHist("factorization_ht500/Fail_lt_point2_HT1000", scaleTo);
	Predictions_t pred_HT1000 = GetPredictions(sidebandHist_ht1000, gausFit2, sigHist_ht1000); 

	cout << "\n\n >>>>>> HT>1200 PREDICTIONS " << endl; 
	TH1 *sigHist_ht1200       = GetHist("factorization_ht500/Signal_HT1200MHT200", scaleTo);
	TH1 *sidebandHist_ht1200  = GetHist("factorization_ht500/Fail_lt_point2_HT1200", scaleTo);
	Predictions_t pred_HT1200 = GetPredictions(sidebandHist_ht1200, gausFit2, sigHist_ht1200); 



	cout << "\n\n >>>>>> 500<HT<800 PREDICTIONS " << endl; 
	TH1 *sigHist_ht500to800       = GetHist("factorization_ht500/Signal_HT500to800MHT200", scaleTo);
	TH1 *sidebandHist_ht500to800  = GetHist("factorization_ht500/Fail_lt_point2_500HT800", scaleTo);
	Predictions_t pred_HT500to800 = GetPredictions(sidebandHist_ht500to800, gausFit2, sigHist_ht500to800); 

	cout << "\n\n >>>>>> 800<HT<1000 PREDICTIONS " << endl; 
	TH1 *sigHist_ht800to1000       = GetHist("factorization_ht500/Signal_HT800to1000MHT200", scaleTo);
	TH1 *sidebandHist_ht800to1000  = GetHist("factorization_ht500/Fail_lt_point2_800HT1000", scaleTo);
	Predictions_t pred_HT800to1000 = GetPredictions(sidebandHist_ht800to1000, gausFit2, sigHist_ht800to1000); 

	cout << "\n\n >>>>>> 1000<HT<1200 PREDICTIONS " << endl; 
	TH1 *sigHist_ht1000to1200       = GetHist("factorization_ht500/Signal_HT1000to1200MHT200", scaleTo);
	TH1 *sidebandHist_ht1000to1200  = GetHist("factorization_ht500/Fail_lt_point2_1000HT1200", scaleTo);
	Predictions_t pred_HT1000to1200 = GetPredictions(sidebandHist_ht1000to1200, gausFit2, sigHist_ht1000to1200); 

	cout << "\n\n >>>>>> 1200<HT<1400 PREDICTIONS " << endl; 
	TH1 *sigHist_ht1200to1400       = GetHist("factorization_ht500/Signal_HT1200to1400MHT200", scaleTo);
	TH1 *sidebandHist_ht1200to1400  = GetHist("factorization_ht500/Fail_lt_point2_1200HT1400", scaleTo);
	Predictions_t pred_HT1200to1400 = GetPredictions(sidebandHist_ht1200to1400, gausFit2, sigHist_ht1200to1400); 

	cout << "\n\n >>>>>> HT>500 PREDICTIONS " << endl; 
	PrintInclResults(pred_HT500, 500);
	cout << "\n\n >>>>>> HT>800 PREDICTIONS " << endl; 
	PrintInclResults(pred_HT800, 800);
	cout << "\n\n >>>>>> HT>1000 PREDICTIONS " << endl; 
	PrintInclResults(pred_HT1000, 1000);
	cout << "\n\n >>>>>> HT>1200 PREDICTIONS " << endl; 
	PrintInclResults(pred_HT1200, 1200);


	cout << "\n\n >>>>>> 500<HT<800 PREDICTIONS " << endl; 
	PrintExclResults(pred_HT500to800, 500, 800);
	cout << "\n\n >>>>>> 800<HT<1000 PREDICTIONS " << endl; 
	PrintExclResults(pred_HT800to1000, 800, 1000);
	cout << "\n\n >>>>>> 1000<HT<1200 PREDICTIONS " << endl; 
	PrintExclResults(pred_HT1000to1200, 1000, 1200);
	cout << "\n\n >>>>>> 1200<HT<1400 PREDICTIONS " << endl; 
	PrintExclResults(pred_HT1200to1400, 1200, 1400);

*/
	//save results to a hist
	TFile f("Results.root","RECREATE");
	Hist_pass->Write();
	expFit2->Write();
	gausFit2->Write();
	leg->Write();
	exp_stats->Write();
	gaus_stats->Write();
	f.Close();


}



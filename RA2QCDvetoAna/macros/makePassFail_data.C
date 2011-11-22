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
	if(x[0]>0.0) fitval=par[0] * exp(par[1] * x[0]) + 0.02;
	else fitval=-1.0E6;
	return fitval;
}

double gausFitFunc(double *x, double *par)
{
	double arg = par[0] * exp(par[1] * x[0]);
	//double fitval = 1.0 / TMath::Erf (arg) - par[2];
	//double fitval = 1.0 / TMath::Erf (arg) - 1 + 0.03;
	double fitval = 1.0 / TMath::Erf (arg) - 1 + 0.02;
	return fitval;
}



void makePassFail_data(const float fitrange_xmin = 50, const float fitrange_xmax = 150) 
{
	const bool logScale = true;
	const string title = " : HT>350 GeV;#slash{H}_{T};Ratio (r) = Pass(RA2 dPhi cuts) / Fail(#Delta #phi_{min}< 0.2);";

	TH1 *hist_pass = 0;
	TH1 *hist_fail = 0;
	TH1 *hist_control = 0;

	TFile* rootFile = new TFile ("Data.root");
	if (rootFile->IsZombie())
	{
		cout << "Data root file not found!" <<  endl;
		assert (false);
	}


	cout << __LINE__ << endl;
	const std::string histname1("factorization_ht500/Pass_RA2dphi_HT500");
	const std::string histname2("factorization_ht500/Fail_1");
	const std::string control1("factorization_ht500/Fail_lt_point2_HT500");

	cout << __LINE__ << endl;
	hist_pass = dynamic_cast<TH1*> (rootFile->Get(histname1.c_str()));
	hist_fail = dynamic_cast<TH1*> (rootFile->Get(histname2.c_str()));
	hist_control = dynamic_cast<TH1*> (rootFile->Get(control1.c_str()));
	if (hist_pass == 0 ) { cout << "Hist " << histname1 << " not found!" << endl; assert (false); }
	if (hist_fail == 0 ) { cout << "Hist " << histname2 << " not found!" << endl; assert (false); }
	if (hist_control == 0 ) { cout << "Hist " << control1 << " not found!" << endl; assert (false); }
	hist_pass->Sumw2();
	hist_fail->Sumw2();
	hist_control->Sumw2();
	cout << __LINE__ << endl;

	//DumpHist(Hist_pass);
	//DumpHist(Hist_fail);
	//std::cout << "Entries pass/ fail " << hist_pass->GetEntries() << " / " << hist_fail->GetEntries() << std::endl;
	//new TCanvas();
	//hist_pass->DrawCopy();
	//hist_fail->DrawCopy("same");
	//gPad->SetEditable(0);

	hist_pass->SetLineWidth(2);
	hist_pass->SetLineColor(kBlue);
	hist_fail->SetLineWidth(2);
	hist_fail->SetLineColor(kBlack);

	TH1*	Hist_pass = dynamic_cast<TH1*> (hist_pass->Clone("histpass_copy"));
	TH1*	Hist_fail = dynamic_cast<TH1*> (hist_fail->Clone("histfail_copy"));

	cout << __LINE__ << endl;
/*	new TCanvas();
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
	Hist_pass->Divide(Hist_fail);
	const int maxbin = Hist_pass->GetMaximumBin();
	const double max = Hist_pass->GetBinContent(maxbin);
	//Hist_pass->GetYaxis()->SetRangeUser(-0.05, max+0.05);

	//fit range
	//const float fitrange_xmin = 50, fitrange_xmax = 120;
	stringstream newtitle;
	newtitle << "Fit Range " << fitrange_xmin << "--" << fitrange_xmax << title;

	gStyle->SetOptStat(0);
	//debug stuff
	new TCanvas();
	//gPad->SetGridy();
	gPad->SetTickx();
	if (logScale) gPad->SetLogy();
	Hist_pass->SetLineColor(9);
	Hist_pass->SetTitle(newtitle.str().c_str());
	Hist_pass->SetLineWidth(2);

	gStyle->SetOptFit(1);
	//Hist_pass->SetStats(0);
	Hist_pass->GetYaxis()->SetRangeUser(1e-2,10);
	Hist_pass->Draw();

	cout << __LINE__ << endl;
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
	//expFit2->SetParameter(0,3.725);
	//expFit2->SetParameter(1,-0.0116);
	expFit2->SetLineColor(kRed);
	expFit2->SetLineWidth(1);

	cout << __LINE__ << endl;

	TF1 *gausFit=new TF1("fit_3",gausFitFunc,fitrange_xmin, fitrange_xmax,2);
	gausFit->SetParameter(0,0.09); 
	gausFit->SetParameter(1,-0.0002);
	gausFit->SetParameter(2,-1.0);
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
	//gausFit2->SetParameter(0,0.2016); 
	//gausFit2->SetParameter(1,0.00778);
	gausFit2->SetLineColor(kGreen);
	gausFit2->SetLineWidth(1);
	gausFit2->Draw("same");
	expFit2->Draw("same");

	cout << __LINE__ << endl;
	TLegend *leg  = new TLegend(0.7,0.8,0.9,0.9);
	leg->AddEntry(gausFit2,"Gaus");
	leg->AddEntry(expFit2,"Exp");
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


	cout << __LINE__ << endl;
	//fit quality
/*	stringstream gaus_fit_res, exp_fit_res;
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
*/
	stringstream epsname;
	epsname << "factorization_fitrange_" << fitrange_xmin << "to" << fitrange_xmax << ".eps";
	gPad->Print(epsname.str().c_str());


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
	double sumGaus = 0, sumExp = 0;
	double sumGaus_excl = 0, sumExp_excl = 0;
	double Gaus_StatErr = 0, Exp_Stat_Err = 0;
	double Gaus_StatErr_excl = 0, Exp_Stat_Err_excl = 0;
	std::vector<float> incMHTBins;
	incMHTBins.push_back(200);
	incMHTBins.push_back(350);
	incMHTBins.push_back(500);

	
	//HT>500 predictions
	cout << "\n\n >>>>>> HT>500 PREDICTIONS " << endl; 
	incMHTBins.clear();
	incMHTBins.push_back(200);
	incMHTBins.push_back(350);

	//TH1* ht500_fail = dynamic_cast<TH1*> (rootFile->Get("factorization_ht500/Fail_lt_point2_HT500"));
	TH1* ht500_fail = dynamic_cast<TH1*> (rootFile->Get("factorization_ht500/Fail_1"));
	//TH1* ht500_fail = dynamic_cast<TH1*> (rootFile->Get("Fail_lt_point2_HT800"));
	assert(ht500_fail != NULL && "Fail_lt_point2_HT500 not found!");
	new TCanvas(); gPad->SetLogy(); ht500_fail->DrawCopy();

	for (unsigned mhtBin = 0; mhtBin < incMHTBins.size(); ++mhtBin)
	{
		sumGaus = 0, sumExp = 0;
		sumGaus_excl = 0, sumExp_excl = 0;
		Gaus_StatErr = 0, Exp_Stat_Err = 0;
		Gaus_StatErr_excl = 0, Exp_Stat_Err_excl = 0;
		for (int bin = 0; bin <= ht500_fail->GetNbinsX()+1; ++bin)
		{
			if (ht500_fail->GetBinContent(bin)>0)
			{
				/*			 cout << std::setprecision(4) << std::setw(5) << bin << std::setw(3) << "[" 
							 << std::setw(10) << ht500_fail->GetBinLowEdge(bin) << ", "  << std::setw(10) 
							 << ht500_fail->GetXaxis()->GetBinUpEdge(bin)<< "]" 
							 << std::setw(10) << ht500_fail->GetBinContent(bin) 
							 << std::setw(10) << ht500_fail->GetBinContent(bin) * gausFit2->Eval(ht500_fail->GetBinCenter(bin))
							 << std::setw(10) << ht500_fail->GetBinContent(bin) * expFit2->Eval(ht500_fail->GetBinCenter(bin)) << endl;
							 */				//inclusive bin stuff
				//if (ht500_fail->GetBinCenter(bin) > 500.)

				if (ht500_fail->GetBinCenter(bin) > incMHTBins.at(mhtBin))
				{
					sumGaus      += ht500_fail->GetBinContent(bin) * gausFit2->Eval(ht500_fail->GetBinCenter(bin));
					sumExp       += ht500_fail->GetBinContent(bin) * expFit2->Eval(ht500_fail->GetBinCenter(bin));
					Gaus_StatErr += ht500_fail->GetBinError(bin)   * gausFit2->Eval(ht500_fail->GetBinCenter(bin));
					Exp_Stat_Err += ht500_fail->GetBinError(bin)   * expFit2->Eval(ht500_fail->GetBinCenter(bin));
				}

				//exclsuive bin stuff
				if (mhtBin+1<incMHTBins.size())
				{
					if (ht500_fail->GetBinCenter(bin) > incMHTBins.at(mhtBin) 
							&& ht500_fail->GetBinCenter(bin) < incMHTBins.at(mhtBin+1))
					{
						sumGaus_excl      += ht500_fail->GetBinContent(bin) * gausFit2->Eval(ht500_fail->GetBinCenter(bin));
						sumExp_excl       += ht500_fail->GetBinContent(bin) * expFit2->Eval(ht500_fail->GetBinCenter(bin));
						Gaus_StatErr_excl += ht500_fail->GetBinError(bin)   * gausFit2->Eval(ht500_fail->GetBinCenter(bin));
						Exp_Stat_Err_excl += ht500_fail->GetBinError(bin)   * expFit2->Eval(ht500_fail->GetBinCenter(bin));
					}

				}

			}
		}
		cout<< setprecision(3) << setw(15) << "MHT > " << incMHTBins.at(mhtBin) 
			<< setw(20) << sumGaus  << "+/-" << Gaus_StatErr 
			<< setw(20) << sumExp  << " +/-" << Exp_Stat_Err << endl;
		if (mhtBin+1<incMHTBins.size())
		{
			cout<< setprecision(3) << setw(10) << incMHTBins.at(mhtBin) << "<MHT<" << incMHTBins.at(mhtBin+1)
				<< setw(20) << sumGaus_excl  << "+/-" << Gaus_StatErr_excl 
				<< setw(20) << sumExp_excl  << " +/-" << Exp_Stat_Err_excl << endl;
		}

	}

	//HT>800 predictions
	cout << "\n\n >>>>>> HT>800 PREDICTIONS " << endl; 
	incMHTBins.clear();
	incMHTBins.push_back(200);
	incMHTBins.push_back(500);

	//TH1* ht800_fail = dynamic_cast<TH1*> (rootFile->Get("factorization_ht500/Fail_lt_point2_HT800"));
	TH1* ht800_fail = dynamic_cast<TH1*> (rootFile->Get("factorization_ht800/Fail_1"));
	assert(ht800_fail != NULL && "Fail_lt_point2_HT800 not found!");
	new TCanvas(); gPad->SetLogy(); ht800_fail->DrawCopy();

	for (unsigned mhtBin = 0; mhtBin < incMHTBins.size(); ++mhtBin)
	{
		sumGaus = 0, sumExp = 0;
		sumGaus_excl = 0, sumExp_excl = 0;
		Gaus_StatErr = 0, Exp_Stat_Err = 0;
		Gaus_StatErr_excl = 0, Exp_Stat_Err_excl = 0;
		for (int bin = 0; bin <= ht800_fail->GetNbinsX()+1; ++bin)
		{
			if (ht800_fail->GetBinContent(bin)>0)
			{
				/*			 cout << std::setprecision(4) << std::setw(5) << bin << std::setw(3) << "[" 
							 << std::setw(10) << ht800_fail->GetBinLowEdge(bin) << ", "  << std::setw(10) 
							 << ht800_fail->GetXaxis()->GetBinUpEdge(bin)<< "]" 
							 << std::setw(10) << ht800_fail->GetBinContent(bin) 
							 << std::setw(10) << ht800_fail->GetBinContent(bin) * gausFit2->Eval(ht800_fail->GetBinCenter(bin))
							 << std::setw(10) << ht800_fail->GetBinContent(bin) * expFit2->Eval(ht800_fail->GetBinCenter(bin)) << endl;
							 */				//inclusive bin stuff

				if (ht800_fail->GetBinCenter(bin) > incMHTBins.at(mhtBin))
				{
					sumGaus      += ht800_fail->GetBinContent(bin) * gausFit2->Eval(ht800_fail->GetBinCenter(bin));
					sumExp       += ht800_fail->GetBinContent(bin) * expFit2->Eval(ht800_fail->GetBinCenter(bin));
					Gaus_StatErr += ht800_fail->GetBinError(bin)   * gausFit2->Eval(ht800_fail->GetBinCenter(bin));
					Exp_Stat_Err += ht800_fail->GetBinError(bin)   * expFit2->Eval(ht800_fail->GetBinCenter(bin));
				}

				//exclsuive bin stuff
				if (mhtBin+1<incMHTBins.size())
				{
					if (ht800_fail->GetBinCenter(bin) > incMHTBins.at(mhtBin) 
							&& ht800_fail->GetBinCenter(bin) < incMHTBins.at(mhtBin+1))
					{
						sumGaus_excl      += ht800_fail->GetBinContent(bin) * gausFit2->Eval(ht800_fail->GetBinCenter(bin));
						sumExp_excl       += ht800_fail->GetBinContent(bin) * expFit2->Eval(ht800_fail->GetBinCenter(bin));
						Gaus_StatErr_excl += ht800_fail->GetBinError(bin)   * gausFit2->Eval(ht800_fail->GetBinCenter(bin));
						Exp_Stat_Err_excl += ht800_fail->GetBinError(bin)   * expFit2->Eval(ht800_fail->GetBinCenter(bin));
					}

				}

			}
		}
		cout << setprecision(3) << setw(15) << "MHT > " << incMHTBins.at(mhtBin) 
			<< setw(20) << sumGaus  << "+/-" << Gaus_StatErr 
			<< setw(20) << sumExp  << " +/-" << Exp_Stat_Err << endl;
		if (mhtBin+1<incMHTBins.size())
		{
			cout<< setprecision(3) << setw(10) << incMHTBins.at(mhtBin) << "<MHT<" << incMHTBins.at(mhtBin+1)
				<< setw(20) << sumGaus_excl  << "+/-" << Gaus_StatErr_excl 
				<< setw(20) << sumExp_excl  << " +/-" << Exp_Stat_Err_excl << endl;
		}

	}



}


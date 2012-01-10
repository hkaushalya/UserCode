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

double expFitFunc(double *x, double *par)
{
	double fitval=0.0;
	//if(x[0]>0.0) fitval=par[0] * exp(par[1] * x[0]) + 0.03;
	if(x[0]>0.0) fitval=par[0] * exp(par[1] * x[0]) + 0.0217;
	else fitval=-1.0E6;
	return fitval;
}

double gausFitFunc(double *x, double *par)
{
	double arg = par[0] * exp(par[1] * x[0]);
	//double fitval = 1.0 / TMath::Erf (arg) - par[2];
	//double fitval = 1.0 / TMath::Erf (arg) - 1 + 0.03;
	double fitval = 1.0 / TMath::Erf (arg) - 1 + 0.0217;
	return fitval;
}

void DumpHist(const TH1* h)
{
	h->Print();
	cout << "bin" << "val +/- err" << endl;

	for (int bin =1; bin <= h->GetNbinsX(); ++bin)
	{
	 cout << "[" << h->GetBinLowEdge(bin) << "] = " << h->GetBinContent(bin) << " / "
	 		<< h->GetBinError(bin) << endl;
	}

}


//void makePassFail(const string histname1, 
//			string histname2, 
//			const string title, const bool logScale)
void makePassFail_MG(const float fitrange_xmin = 50, const float fitrange_xmax = 150) 
{

	const bool logScale = true;
	const string title = " : HT>500 GeV;#slash{H}_{T};Ratio (r) = Pass(RA2 dPhi cuts) / Fail(#Delta #phi_{min}< 0.2);";
	const float scaleTo = 4650; // pb

	//cross section for each sample in pt order
/*	const float cs[] = {	
								2.78e+11 ,
								833057,
								17132,
								409.6
								};
								*/
	const float cs[] = {	2.7800E+08 * 0,15134,
								833057.0 * 0.238 ,
								17132 * 0.34,
								409.6 * 0.29932
								};


	const float nevts[] = {23739003,
								  20674219,
								  14437469,
								  6294851
	};
	
	
	

	const Int_t nBins = 4;
	TFile *files[nBins];
	TH1 *hists_pass[nBins] = {0,0,0,0};
	TH1 *hists_fail[nBins] = {0,0,0,0};
	TH1 *hists_signal_mht200[nBins] = {0,0,0,0};
	TH1 *hists_control[nBins] = {0,0,0,0,};

	files[0] = new TFile ("qcd1/Merged.root");
	files[1] = new TFile ("qcd2/Merged.root");
	files[2] = new TFile ("qcd3/Merged.root");
	files[3] = new TFile ("qcd4/Merged.root");

	TH1 *Hist_pass = 0;
	TH1 *Hist_fail = 0;
	TH1 *Hist_signal_mht200 = 0;
	TH1 *Hist_control = 0;
	const std::string histname1("factorization_ht500/Pass_RA2dphi_HT500");
	const std::string histname2("factorization_ht500/Fail_1");
	//const std::string histname1("factorization_ht350/Pass_RA2dphi");
	//const std::string histname2("factorization_ht350/Fail_1");
	//const std::string signal1_mht200("factorization_ht500/Signal_0");
	//const std::string control1("factorization_ht500/Fail_lt_point2_HT500");
	const std::string signal1_mht200("factorization_ht500/Signal_HT500MHT200");
	const std::string control1("factorization_ht500/Fail_lt_point2_HT500");

	for (int i=0; i<nBins; ++i)
	{
		if (files[i]->IsZombie())
		{
			cout << files[i]->GetName() << " not found!" <<  endl;
			assert (false);
		} else
		{
			//cout << files[i]->GetName() << " opened." <<  endl;

			hists_pass[i] = dynamic_cast<TH1*> (files[i]->Get(histname1.c_str()));
			hists_fail[i] = dynamic_cast<TH1*> (files[i]->Get(histname2.c_str()));
			hists_signal_mht200[i] = dynamic_cast<TH1*> (files[i]->Get(signal1_mht200.c_str()));
			hists_control[i] = dynamic_cast<TH1*> (files[i]->Get(control1.c_str()));
			if (hists_pass[i] == 0)
			{
				cout << histname1 << " not found in " << files[i]->GetName() << "!" << endl;
				assert (false);
			}
			if (hists_fail[i] == 0)
			{
				cout << histname1 << " not found in " << files[i]->GetName() << "!" << endl;
				assert (false);
			}
			if (hists_signal_mht200[i] == 0)
			{
				cout << signal1_mht200 << " not found in " << files[i]->GetName() << "!" << endl;
				assert (false);
			}
			if (hists_control[i] == 0)
			{
				cout << control1 << " not found in " << files[i]->GetName() << "!" << endl;
				assert (false);
			}
		


			{
				//std::cout << "Entries pass/ fail " << hists_pass[i]->GetEntries() << " / " << hists_fail[i]->GetEntries() << std::endl;
				//hists_pass[i]->Print();
				//new TCanvas();
				//hists_pass[i]->Draw();
				//gPad->SetEditable(0);
				hists_pass[i]->Sumw2();
				hists_fail[i]->Sumw2();
				hists_signal_mht200[i]->Sumw2();
				hists_control[i]->Sumw2();
						
        		//print last bin events
				/*const int lastbin = hists_pass[i]->GetNbinsX();
				std::cout << "pass: last bin events [" << hists_pass[i]->GetBinLowEdge(lastbin) << "] " <<
					hists_pass[i]->GetBinContent(lastbin)
					<< " ][overflow = " << hists_pass[i]->GetBinContent(lastbin+1)
					<< std::endl;
				std::cout << "fail: last bin events [" << hists_fail[i]->GetBinLowEdge(lastbin) << "] " <<
					hists_fail[i]->GetBinContent(lastbin)
					<< " ][overflow = " << hists_fail[i]->GetBinContent(lastbin+1)
					<< std::endl;
				*/
				DumpHist(hists_pass[i]);
				DumpHist(hists_fail[i]);


				const float scale = scaleTo/ ( nevts[i] / cs[i] );
				hists_pass[i]->Scale(scale);
				hists_pass[i]->SetLineWidth(2);
				hists_pass[i]->SetLineColor(i);
				hists_fail[i]->Scale(scale);
				hists_fail[i]->SetLineWidth(2);
				hists_fail[i]->SetLineColor(i+11);
				
				hists_signal_mht200[i]->Scale(scale);
				hists_control[i]->Scale(scale);

				//if (i != 0) hists_pass[0]->Add(hists_pass[i]);
				if (i == 0) 
				{
					Hist_pass = dynamic_cast<TH1*> (hists_pass[i]->Clone("histpass_copy"));
					Hist_fail = dynamic_cast<TH1*> (hists_fail[i]->Clone("histfail_copy"));
					Hist_signal_mht200 = dynamic_cast<TH1*> (hists_signal_mht200[i]->Clone("histsignal_mth200_copy"));
					Hist_control = dynamic_cast<TH1*> (hists_control[i]->Clone("histscontrol_copy"));
				} else
				{
					Hist_pass->Add(hists_pass[i]);
					Hist_fail->Add(hists_fail[i]);
					Hist_signal_mht200->Add(hists_signal_mht200[i]);
					Hist_control->Add(hists_control[i]);
				}

			}
		}
	}

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
	/*new TCanvas(); gPad->SetLogy(); Hist_pass->DrawCopy();
	new TCanvas(); gPad->SetLogy(); Hist_fail->DrawCopy();
	new TCanvas(); gPad->SetLogy(); Hist_signal_mht200->DrawCopy();
*/

	Hist_pass->Divide(Hist_fail);
	const int maxbin = Hist_pass->GetMaximumBin();
	const double max = Hist_pass->GetBinContent(maxbin);
	//Hist_pass->GetYaxis()->SetRangeUser(-0.05, max+0.05);



	//fit range
	//const float fitrange_xmin = 50, fitrange_xmax = 120;
	stringstream newtitle;
	newtitle << "Fit Range " << fitrange_xmin << "--" << fitrange_xmax << title;


	gStyle->SetOptStat(0);
	//debug stff
	new TCanvas();
	//gPad->SetGridy();
	gPad->SetTickx();
	if (logScale) gPad->SetLogy();
	Hist_pass->SetLineColor(9);
	Hist_pass->SetTitle(newtitle.str().c_str());
	Hist_pass->SetLineWidth(2);

	gStyle->SetOptFit(1);
	//Hist_pass->SetStats(0);
	Hist_pass->Draw();

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
	expFit2->SetLineColor(kRed);
	expFit2->SetLineWidth(1);


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
	epsname << "factorization_fitrange_" << fitrange_xmin << "to" << fitrange_xmax << ".eps";
	gPad->Print(epsname.str().c_str());

	 //make a predicion
	double sumGaus = 0, sumExp = 0;
	double sumGaus_excl = 0, sumExp_excl = 0;
	double Gaus_StatErr = 0, Exp_Stat_Err = 0;
	double Gaus_StatErr_excl = 0, Exp_Stat_Err_excl = 0;
	std::vector<float> incMHTBins;
	incMHTBins.push_back(200);
	incMHTBins.push_back(350);
	incMHTBins.push_back(500);

	
	cout << "\n\n >>>>>> HT>500 PREDICTIONS " << endl; 
	cout << setw(15) << "MHT bin "<< setw(20) << " Val" << setw(20) << " Gaus " << setw(20) << " Exp "  << setw(20) << "  Signal "<< endl;

	for (unsigned mhtBin = 0; mhtBin < incMHTBins.size(); ++mhtBin)
	{
		sumGaus = 0, sumExp = 0;
		sumGaus_excl = 0, sumExp_excl = 0;
		Gaus_StatErr = 0, Exp_Stat_Err = 0;
		Gaus_StatErr_excl = 0, Exp_Stat_Err_excl = 0;
		double signal = 0, signal_excl = 0;


		for (int bin = 0; bin <= Hist_control->GetNbinsX()+1; ++bin)
		{
							 				//inclusive bin stuff
				if (Hist_control->GetBinCenter(bin) > incMHTBins.at(mhtBin))
				{
							/* cout << std::setprecision(4) << std::setw(5) << bin << std::setw(3) << "[" 
							 << std::setw(10) << Hist_control->GetBinLowEdge(bin) << ", "  << std::setw(10) 
							 << Hist_control->GetXaxis()->GetBinUpEdge(bin)<< "]" 
							 << std::setw(15) << Hist_control->GetBinContent(bin) 
							 << std::setw(15) << Hist_control->GetBinContent(bin) * gausFit2->Eval(Hist_control->GetBinCenter(bin))
							 << std::setw(15) << Hist_control->GetBinContent(bin) * expFit2->Eval(Hist_control->GetBinCenter(bin))
							 << std::setw(15) << Hist_signal_mht200->GetBinContent(bin) 
							 << endl;
							 */

					sumGaus      += Hist_control->GetBinContent(bin) * gausFit2->Eval(Hist_control->GetBinCenter(bin));
					sumExp       += Hist_control->GetBinContent(bin) * expFit2->Eval(Hist_control->GetBinCenter(bin));
					Gaus_StatErr += Hist_control->GetBinError(bin)   * gausFit2->Eval(Hist_control->GetBinCenter(bin));
					Exp_Stat_Err += Hist_control->GetBinError(bin)   * expFit2->Eval(Hist_control->GetBinCenter(bin));
					signal       += Hist_signal_mht200->GetBinContent(bin);
				}

									//exclsuive bin stuff
				if (mhtBin+1<incMHTBins.size())
				{
					if (Hist_control->GetBinCenter(bin) > incMHTBins.at(mhtBin) 
							&& Hist_control->GetBinCenter(bin) < incMHTBins.at(mhtBin+1))
					{
						sumGaus_excl      += Hist_control->GetBinContent(bin) * gausFit2->Eval(Hist_control->GetBinCenter(bin));
						sumExp_excl       += Hist_control->GetBinContent(bin) * expFit2->Eval(Hist_control->GetBinCenter(bin));
						Gaus_StatErr_excl += Hist_control->GetBinError(bin)   * gausFit2->Eval(Hist_control->GetBinCenter(bin));
						Exp_Stat_Err_excl += Hist_control->GetBinError(bin)   * expFit2->Eval(Hist_control->GetBinCenter(bin));
						signal_excl       += Hist_signal_mht200->GetBinContent(bin);
					}

				}

		}

		cout<< setprecision(3) << setw(15) << "MHT > " << incMHTBins.at(mhtBin) 
			<< setw(20) << sumGaus  << "+/-" << Gaus_StatErr 
			<< setw(20) << sumExp  << " +/-" << Exp_Stat_Err  << setw(30) << signal << endl;

		if (mhtBin+1<incMHTBins.size())
		{
			cout<< setprecision(3) << setw(10) << incMHTBins.at(mhtBin) << "<MHT<" << incMHTBins.at(mhtBin+1)
				<< setw(20) << sumGaus_excl  << "+/-" << Gaus_StatErr_excl 
				<< setw(20) << sumExp_excl  << " +/-" << Exp_Stat_Err_excl  << setw(30) << signal_excl << endl;
		}

	}

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



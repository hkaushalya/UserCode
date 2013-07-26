#include<iostream>
#include<sstream>
#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include "TPad.h"
#include <vector>
#include "TStyle.h"
#include "TLegend.h"
#include <string>
#include "assert.h"
#include "TPaveText.h"
#include <iomanip>

using namespace std;

const static double LargeNegNum = -99999.99;


double StatErr(const TH1* h)
{
	double toterr2 = 0;
	for (int bin = 1; bin <= h->GetNbinsX(); ++bin)
	{
		const double binerr = h->GetBinError(bin);
		toterr2 += (binerr*binerr);
	}
	const double toterr = sqrt(toterr2);
	return toterr;
}

TCanvas* GetCanvas(TPad *p1, TPad *p2)
{
   //TCanvas *c1 = new TCanvas("c1", "c1",15,60,550,600);
   TCanvas *c1 = new TCanvas("c1");
   c1->Range(0,0,1,1);
   c1->SetBorderSize(2);
   c1->SetFrameFillColor(0);
  
// ------------>Primitives in pad: c1_1
   TPad *c1_1 = new TPad("c1_1", "c1_1",0.01,0.30,0.99,0.99);
   c1_1->Draw();
   c1_1->cd();
   c1_1->SetBorderSize(2);
   c1_1->SetTickx(1);
   c1_1->SetTicky(1);
   c1_1->SetTopMargin(0.1);
   c1_1->SetBottomMargin(0.0);
   //c1_1->SetFrameFillColor(3);
	//c1_1->SetLogy();
  
// ------------>Primitives in pad: c1_2
   TPad *c1_2 = new TPad("c1_2", "c1_2",0.01,0.01,0.99,0.30);
   c1_2->Draw();
   c1_2->cd();
   c1_2->SetBorderSize(2);
   c1_2->SetTickx(1);
   c1_2->SetTicky(1);
   c1_2->SetTopMargin(0.0);
   c1_2->SetBottomMargin(0.24);
   c1_2->SetFrameFillColor(0);
	c1_2->SetGridx();
	c1_2->SetGridy();

	p1 = c1_1;
	p2 = c1_2;
	return c1;
};

class Hist
{
	public:
	Hist(const string name, const string title="", const int rebin=1, 
			const double xmin=LargeNegNum, const double xmax=LargeNegNum)
	{
		name_  = name;
		title_ = title;
		rebin_ = rebin;
		xmin_  = xmin;
		xmax_  = xmax;
	};
	
	string Name() const { return name_; }
	string Title() const { return title_; }
	int Rebin() const { return rebin_; }
	double Xmin() const { return xmin_; }
	double Xmax() const { return xmax_; }

	private:
		string name_, title_;
		int rebin_;
		double xmin_, xmax_;
};


//void makeMetRatio(const int sample=0)
void makeMetRatio(const string title="")
{
	vector<Hist> hist2print;
	TPaveText *tx = new TPaveText(.05,.1,.95,.8);
//	tx->AddText("Using Deafult JERs for all jets");
//	tx->AddText("Using b-Jet JERs");

	string mtitle("");

	//if (sample==1) mtitle += "NJet(70/50/30>=2/4/5), #slash{E}_{T}>175, Triplet>1, 80<TopMass<270, TOP+0.5*BJET>500, MT2>300, #Delta#Phi(.5,.5,.3), BJets>=1";
/*	if (sample==1)      mtitle += "All Stop cuts applied (use default JERs for all jets)";
	else if (sample==2) mtitle += "All Stop cuts applied + Inverted #Delta#Phi (use default JERs for all jets)";
	else if (sample==3) mtitle += "All Stop cuts applied (use b-Jet JERs)";
	else if (sample==4) mtitle += "All Stop cuts applied + Inverted #Delta#Phi (use b-Jet JERs)";
*/


	vector<pair<unsigned, unsigned> > vBitMasks;
			
	//unsigned bitMaskArray[] = {0,1,2,3,129,130,131,195,257,258,259,323};
	vBitMasks.push_back(make_pair(128,256));
	vBitMasks.push_back(make_pair(129,257));
	vBitMasks.push_back(make_pair(130,258));
	vBitMasks.push_back(make_pair(131,259));
	vBitMasks.push_back(make_pair(195,323));


	//hist2print.push_back(Hist("met","Smeared Gen-Jets (QCD MG)",2,150.0, 400.0));
	hist2print.push_back(Hist("met",title,4,0.0, 500.0));
	hist2print.push_back(Hist("met1",title,4,0.0, 500.0));
	hist2print.push_back(Hist("met2",title,4,0.0, 500.0));
	hist2print.push_back(Hist("met3",title,4,0.0, 500.0));
	hist2print.push_back(Hist("met4",title,4,0.0, 500.0));
	hist2print.push_back(Hist("dphiMin_met",title));
	hist2print.push_back(Hist("dphiMin_mht",title));
	hist2print.push_back(Hist("mht",title,4,0.0, 500.0));
	hist2print.push_back(Hist("unclmet",title,4,0.0, 500.0));
	hist2print.push_back(Hist("jet1_pt",title,4,0.0, 500.0));
	hist2print.push_back(Hist("jet1_eta",title,2));
	hist2print.push_back(Hist("jet1_phi",title,2));
	hist2print.push_back(Hist("jet1_dphimet",title,2));
	hist2print.push_back(Hist("jet1_dphimht",title,2));
	hist2print.push_back(Hist("jet2_pt",title,4,0.0, 500.0));
	hist2print.push_back(Hist("jet2_eta",title,2));
	hist2print.push_back(Hist("jet2_phi",title,2));
	hist2print.push_back(Hist("jet2_dphimet",title,2));
	hist2print.push_back(Hist("jet2_dphimht",title,2));
	hist2print.push_back(Hist("jet3_pt",title,4,0.0, 500.0));
	hist2print.push_back(Hist("jet3_eta",title,2));
	hist2print.push_back(Hist("jet3_phi",title,2));
	hist2print.push_back(Hist("jet3_dphimet",title,2));
	hist2print.push_back(Hist("jet3_dphimht",title,2));

	TFile *outRootFile_num = new TFile("Merged.root");
//	TFile *outRootFile_den = new TFile("StopInvrtDphiCut.root");

   TCanvas *c = new TCanvas("c1");
   c->Range(0,0,1,1);
   c->SetBorderSize(2);
   c->SetFrameFillColor(0);
  
// ------------>Primitives in pad: c1_1
   TPad *c1_1 = new TPad("c1_1", "c1_1",0.01,0.30,0.99,0.99);
   c1_1->Draw();
   c1_1->cd();
   c1_1->SetBorderSize(2);
   c1_1->SetTickx(1);
   c1_1->SetTicky(1);
   c1_1->SetTopMargin(0.1);
   c1_1->SetBottomMargin(0.0);
   //c1_1->SetFrameFillColor(3);
	//c1_1->SetLogy();
  
  c->cd();
// ------------>Primitives in pad: c1_2
   TPad *c1_2 = new TPad("c1_2", "c1_2",0.01,0.01,0.99,0.30);
   c1_2->Draw();
   c1_2->cd();
   c1_2->SetBorderSize(2);
   c1_2->SetTickx(1);
   c1_2->SetTicky(1);
   c1_2->SetTopMargin(0.0);
   c1_2->SetBottomMargin(0.24);
   c1_2->SetFrameFillColor(0);
	c1_2->SetGridx();
	c1_2->SetGridy();
	c->cd();
	gStyle->SetOptStat(0);
	gPad->Print("met_ratio.eps[");

	
	for (unsigned j=0;j<hist2print.size(); ++j)
	{
		for (unsigned i=0;i<vBitMasks.size(); ++i)
		{
			const unsigned mask1 = vBitMasks.at(i).first;
			const unsigned mask2 = vBitMasks.at(i).second;

			stringstream path, reco_num_hist_name;
			stringstream smear_num_hist, reco_num_hist;
			stringstream smear_den_hist, reco_den_hist;
			stringstream num_folder, den_folder;
			num_folder << "Hist/Mask"<< mask1 <<"HT0to8000MHT0to8000/";
			den_folder << "Hist/Mask"<< mask2 <<"HT0to8000MHT0to8000/";
			//cout << "folder = " << folder.str() << endl;

			smear_num_hist << num_folder.str() << "smeared_" << hist2print.at(j).Name();
			reco_num_hist  << num_folder.str() << "reco_"    << hist2print.at(j).Name();
			smear_den_hist << den_folder.str() << "smeared_" << hist2print.at(j).Name();
			reco_den_hist  << den_folder.str() << "reco_"    << hist2print.at(j).Name();

			TH1* hsmear_num = (TH1*) (outRootFile_num->Get(smear_num_hist.str().c_str()));
			if (hsmear_num == NULL) { cout << "hsmear_num = " << smear_num_hist.str() << " was not found!" << endl; assert(false); } 
			hsmear_num->SetDirectory(0);
			hsmear_num->Sumw2();

			TH1* hsmear_den = (TH1*) (outRootFile_num->Get(smear_den_hist.str().c_str()));
			if (hsmear_den == NULL) { cout << "hsmear_den = " << smear_den_hist.str() << " was not found!" << endl; assert(false); } 
			hsmear_den->SetDirectory(0);
			hsmear_den->Sumw2();


			TH1* hreco_num = (TH1*) (outRootFile_num->Get(reco_num_hist.str().c_str()));
			if (hreco_num == NULL) { cout << "hreco_num = " << reco_num_hist.str() << " was not found!" << endl; assert(false); } 
			hreco_num->SetDirectory(0);
			hreco_num->Sumw2();

			TH1* hreco_den = (TH1*) (outRootFile_num->Get(reco_den_hist.str().c_str()));
			if (hreco_den == NULL) { cout << "hreco_den = " << reco_den_hist.str() << " was not found!" << endl; assert(false); } 
			hreco_den->SetDirectory(0);
			hreco_den->Sumw2();



			const int rebin = hist2print.at(j).Rebin();
			//const string title = hist2print.at(j).Title();
			const string title("");

			const double xmin = hist2print.at(j).Xmin();
			const double xmax = hist2print.at(j).Xmax();

			if (rebin>1)
			{
				hsmear_num->Rebin(rebin);
				hsmear_den->Rebin(rebin);
				hreco_num->Rebin(rebin);
				hreco_den->Rebin(rebin);
			}
			if (title.length()>0)
			{
				hsmear_num->SetTitle(title.c_str());
				hsmear_den->SetTitle(title.c_str());
				hreco_num->SetTitle(title.c_str());
				hreco_den->SetTitle(title.c_str());
			}
			if (xmin != LargeNegNum || xmax != LargeNegNum)
			{
				hsmear_num->GetXaxis()->SetRangeUser(xmin,xmax);
				hsmear_den->GetXaxis()->SetRangeUser(xmin,xmax);
				hreco_num->GetXaxis()->SetRangeUser(xmin,xmax);
				hreco_den->GetXaxis()->SetRangeUser(xmin,xmax);
			}

			hsmear_num->SetLineColor(kRed);
			hsmear_num->SetMarkerColor(kRed);
			hsmear_num->SetMarkerStyle(26);
			hsmear_num->SetLineWidth(2);
			//hsmear_num->GetXaxis()->SetRangeUser(0,300);

			hsmear_den->SetLineColor(kRed);
			hsmear_den->SetMarkerColor(kRed);
			hsmear_den->SetMarkerStyle(20);
			hsmear_den->SetLineWidth(2);
			//hsmear_den->GetXaxis()->SetRangeUser(0,300);


			hreco_num->SetLineColor(kBlack);
			hreco_num->SetMarkerColor(kBlack);
			hreco_num->SetMarkerStyle(26);
			hreco_num->SetLineWidth(2);
			//hreco_num->GetXaxis()->SetRangeUser(0,300);

			hreco_den->SetLineColor(kBlack);
			hreco_den->SetMarkerColor(kBlack);
			hreco_den->SetMarkerStyle(20);
			hreco_den->SetLineWidth(2);
			//hreco_den->GetXaxis()->SetRangeUser(0,300);


			hsmear_num->GetYaxis()->CenterTitle(1);
			hsmear_num->SetLabelFont(42,"XYZ");
			hsmear_num->SetTitleFont(42,"XYZ");
			hsmear_num->GetYaxis()->SetTitleOffset(0.8);
			hsmear_num->SetLabelSize(0.05,"XYZ");
			hsmear_num->SetTitleSize(0.06,"XYZ");


			TH1 *hratio = (TH1*) (hsmear_num->Clone("hnum_copy"));
			hratio->Divide(hsmear_den);
			hratio->SetTitle("");
			hratio->GetYaxis()->SetTitle("#frac{Stop cuts}{Stop cuts+Invt. #Delta#Phi}");
			hratio->GetYaxis()->SetRangeUser(0,2);

			TH1 *hratio_reco = (TH1*) (hreco_num->Clone("hnum_copy"));
			hratio_reco->Divide(hreco_den);
			hratio_reco->SetTitle("");
			hratio_reco->GetYaxis()->SetRangeUser(0,0.2);


			hratio->GetYaxis()->SetTitleOffset(0.4);
			hratio->GetXaxis()->SetTitleOffset(0.9);
			hratio->GetYaxis()->CenterTitle(1);
			hratio->GetXaxis()->CenterTitle(1);
			hratio->SetLabelSize(0.125,"XYZ");
			hratio->SetTitleSize(0.1,"XYZ");
			//	hratio->SetLabelFont(labelfont,"XYZ");
			//	hratio->SetTitleFont(titlefont,"XYZ");
			hratio->GetXaxis()->SetTickLength(0.08);

			stringstream numleg, denleg, reco_numleg, reco_denleg;
			const double sum_num = hsmear_num->Integral(1, hsmear_num->GetNbinsX()+1);
			const double sum_den = hsmear_den->Integral(1, hsmear_den->GetNbinsX()+1);
			const double sum_reco_num = hreco_num->Integral(1, hreco_num->GetNbinsX()+1);
			const double sum_reco_den = hreco_den->Integral(1, hreco_den->GetNbinsX()+1);

			const double staterr_num = StatErr(hsmear_num);
			const double staterr_den = StatErr(hsmear_den);
			const double staterr_reco_num = StatErr(hreco_num);
			const double staterr_reco_den = StatErr(hreco_den);

			numleg << "Smear[#Delta#Phi] (" << setprecision(1) << fixed << sum_num << "#pm" << staterr_num << ")";
			denleg << "Smear[Inv.#Delta#Phi] ("<< setprecision(1) << fixed  << sum_den  << "#pm" << staterr_den << ")";
			reco_numleg << "Reco[#Delta#Phi] ("<< setprecision(1) << fixed  << sum_reco_num << "#pm" << staterr_reco_num  << ")";
			reco_denleg << "Reco[Inv.#Delta#Phi] ("<< setprecision(1) << fixed  << sum_reco_den  << "#pm" << staterr_reco_den << ")";

			TLegend *l2 = new TLegend(0.6,0.6,0.9,0.9);
			l2->AddEntry(hsmear_num, numleg.str().c_str());
			l2->AddEntry(hsmear_den, denleg.str().c_str());
			l2->AddEntry(hreco_num, reco_numleg.str().c_str());
			l2->AddEntry(hreco_den, reco_denleg.str().c_str());

			c1_1->cd();
			gPad->SetLogy();
			hsmear_num->DrawCopy();
			hsmear_den->DrawCopy("same");
			hreco_num->DrawCopy("same");
			hreco_den->DrawCopy("same");
			l2->Draw();
			//tx->Draw();
			c1_2->cd();
			hratio->DrawCopy();
			hratio_reco->DrawCopy("same");
			c->cd();
			gPad->Print("met_ratio.eps");

		}
	}

	gPad->Print("met_ratio.eps]");
}

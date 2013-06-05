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

using namespace std;

const static double LargeNegNum = -99999.99;

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
			const double xmin=LargeNegNum, const double xmax=LargeNegNum, const bool logyScale=1)
	{
		name_  = name;
		title_ = title;
		rebin_ = rebin;
		xmin_  = xmin;
		xmax_  = xmax;
		logy_   = logyScale;
	};
	
	string Name() const { return name_; }
	string Title() const { return title_; }
	int Rebin() const { return rebin_; }
	double Xmin() const { return xmin_; }
	double Xmax() const { return xmax_; }
	bool LogY() const { return logy_; }

	private:
		string name_, title_;
		int rebin_;
		double xmin_, xmax_;
		bool logy_;
};


//void makeHist(const int sample, const int dataset=1)
void makeHist(const string title="")
{
	vector<Hist> hist2print;
	TPaveText *tx = new TPaveText(.05,.1,.95,.8);
//	tx->AddText("Using Deafult JERs for all jets");
//	tx->AddText("Using b-Jet JERs");

/*	string title("QCD MG:");

	//if (sample==1) title += "NJet(70/50/30>=2/4/5), #slash{E}_{T}>175, Triplet>1, 80<TopMass<270, TOP+0.5*BJET>500, MT2>300, #Delta#Phi(.5,.5,.3), BJets>=1";
	if (sample==1)      title += "All Stop cuts applied (use default JERs for all jets)";
	else if (sample==2) title += "All Stop cuts applied + Inverted #Delta#Phi (use default JERs for all jets)";
	else if (sample==3) title += "All Stop cuts applied (use b-Jet JERs)";
	else if (sample==4) title += "All Stop cuts applied + Inverted #Delta#Phi (use b-Jet JERs)";
	else if (sample==5) title += "No cuts applied";
*/
	hist2print.push_back(Hist("met",title,2,100.0, 400.0,1));
	hist2print.push_back(Hist("unclmet",title,2,100.0, 400.0,1));
	hist2print.push_back(Hist("mht",title,2,100.0, 400.0,1));
	hist2print.push_back(Hist("ht",title,2,0,2000,1));
	hist2print.push_back(Hist("njet30eta5p0",title,1,0,15,1));
	hist2print.push_back(Hist("nbjets",title,1,0,10,1));
//	hist2print.push_back(Hist("bjetPt",title,2));
	hist2print.push_back(Hist("M123",title,2));
//	hist2print.push_back(Hist("M23overM123",title));
	hist2print.push_back(Hist("MT2",title,2));
	hist2print.push_back(Hist("MTb",title));
	hist2print.push_back(Hist("MTt",title));
	hist2print.push_back(Hist("MTb_p_MTt",title,2,400,1000,1));
	//hist2print.push_back(Hist("jet1_pt",title,2));
	//hist2print.push_back("bjetPt");
//	hist2print.push_back(Hist("bjetMass",title,2,0,200));
//	hist2print.push_back(Hist("dphimin",title,4));


	TFile *outRootFile = new TFile("Merged.root");

	/*TPad *c1=0, *c2=0;
	TCanvas *c = GetCanvas(c1, c2);
	if (c ==NULL|| c1 == 0 ||c2 == 0)
	{
		cout << " A drawing pad is null !"<< endl;
		cout << "c = " << c << endl;
		cout << "c1 = " << c1 << endl;
		cout << "c2 = " << c2 << endl;
		assert(false);
	}*/
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
	gPad->Print("samples.eps[");

	for (unsigned ihist=0; ihist < hist2print.size(); ++ihist)
	{
		stringstream path, reco_hist_name, gen_hist_name, smear_hist_name;
		stringstream reco_hist, gen_hist, smear_hist;
		stringstream folder;
		folder << "Hist/Njet0to1000HT0to8000MHT0to8000/";
		//cout << "folder = " << folder.str() << endl;

	/*	if ((hist2print.at(ihist).Name()).find("Jet"))
		{
			reco_hist_name << folder.str() << "reco" << hist2print.at(ihist).Name() << "_copy";
			reco_hist << folder.str() << "reco" << hist2print.at(ihist).Name();
			smear_hist_name << folder.str() << "smeared" << hist2print.at(ihist).Name() << "_copy";
			smear_hist << folder.str() << "smeared" << hist2print.at(ihist).Name();
			gen_hist_name << folder.str() << "gen" << hist2print.at(ihist).Name() << "_copy";
			gen_hist << folder.str() << "gen" << hist2print.at(ihist).Name();
		} else
	*/	{
			reco_hist_name << folder.str() << "reco_" << hist2print.at(ihist).Name() << "_copy";
			reco_hist << folder.str() << "reco_" << hist2print.at(ihist).Name();
			smear_hist_name << folder.str() << "smeared_" << hist2print.at(ihist).Name() << "_copy";
			smear_hist << folder.str() << "smeared_" << hist2print.at(ihist).Name();
			gen_hist_name << folder.str() << "gen_" << hist2print.at(ihist).Name() << "_copy";
			gen_hist << folder.str() << "gen_" << hist2print.at(ihist).Name();
		}

		TH1* hreco = (TH1*) (outRootFile->Get(reco_hist.str().c_str()));
		if (hreco == NULL) { cout << "hreco = " << reco_hist.str() << " was not found!" << endl; assert(false); } 
		hreco->SetDirectory(0);
		TH1* hsmear = (TH1*) (outRootFile->Get(smear_hist.str().c_str()));
		if (hsmear == NULL) { cout << "hsmear = " << smear_hist.str() << " was not found!" << endl; assert(false); } 
		hsmear->SetDirectory(0);
		TH1* hgen = (TH1*) (outRootFile->Get(gen_hist.str().c_str()));
		//->Clone(gen_hist_name.str().c_str()));
		if (hgen == NULL) { cout << "hgen = " << gen_hist.str() << " was not found!" << endl; assert(false); } 
		hgen->SetDirectory(0);

		hreco->Sumw2();
		hsmear->Sumw2();
		hgen->Sumw2();

		const int rebin = hist2print.at(ihist).Rebin();
		const string title = hist2print.at(ihist).Title();
		const double xmin = hist2print.at(ihist).Xmin();
		const double xmax = hist2print.at(ihist).Xmax();

		if (rebin>1)
		{
			hreco->Rebin(rebin);
			hsmear->Rebin(rebin);
			hgen->Rebin(rebin);
		}
		//if (title.length()>0)
		{
			hreco->SetTitle(title.c_str());
			hsmear->SetTitle(title.c_str());
			hgen->SetTitle(title.c_str());
		}
		if (xmin != LargeNegNum || xmax != LargeNegNum)
		{
			hreco->GetXaxis()->SetRangeUser(xmin,xmax);
			hsmear->GetXaxis()->SetRangeUser(xmin,xmax);
			hgen->GetXaxis()->SetRangeUser(xmin,xmax);
		}

		hgen->SetLineColor(kBlue);
		hgen->SetMarkerColor(kBlue);
		hgen->SetMarkerStyle(24);
		hgen->SetLineWidth(2);
		hsmear->SetLineColor(kRed);
		hsmear->SetMarkerColor(kRed);
		hsmear->SetMarkerStyle(24);
		hsmear->SetLineWidth(2);
		hreco->SetLineWidth(2);
		hreco->SetMarkerStyle(kDot);
		hreco->SetLineColor(kBlack);
		hreco->SetMarkerColor(kBlack);
		//hreco->GetXaxis()->SetRangeUser(0,300);
		//hsmear->GetXaxis()->SetRangeUser(0,300);


		hreco->GetYaxis()->CenterTitle(1);
		hreco->SetLabelFont(42,"XYZ");
		hreco->SetTitleFont(42,"XYZ");
		hreco->GetYaxis()->SetTitleOffset(0.8);
		hreco->SetLabelSize(0.05,"XYZ");
		hreco->SetTitleSize(0.06,"XYZ");



		TH1 *hsmeartoreco_ratio = (TH1*) (hsmear->Clone("hsmear_copy"));
		hsmeartoreco_ratio->Divide(hreco);
		hsmeartoreco_ratio->SetTitle("");
		hsmeartoreco_ratio->GetYaxis()->SetTitle("Smear/Reco");
		hsmeartoreco_ratio->GetYaxis()->SetRangeUser(0,2.);

		hsmeartoreco_ratio->GetYaxis()->SetTitleOffset(0.4);
		hsmeartoreco_ratio->GetXaxis()->SetTitleOffset(0.9);
		hsmeartoreco_ratio->GetYaxis()->CenterTitle(1);
		hsmeartoreco_ratio->GetXaxis()->CenterTitle(1);
		hsmeartoreco_ratio->SetLabelSize(0.125,"XYZ");
		hsmeartoreco_ratio->SetTitleSize(0.125,"XYZ");
		//	hsmeartoreco_ratio->SetLabelFont(labelfont,"XYZ");
		//	hsmeartoreco_ratio->SetTitleFont(titlefont,"XYZ");
		hsmeartoreco_ratio->GetXaxis()->SetTickLength(0.07);

		stringstream recoleg,smearleg, genleg;
		const double sum_reco  = hreco->Integral(1, hreco->GetNbinsX()+1);
		const double sum_smear = hsmear->Integral(1, hsmear->GetNbinsX()+1);
		const double sum_gen   = hgen->Integral(1, hgen->GetNbinsX()+1);
		recoleg << "reco (" << sum_reco << ")";
		smearleg << "smear (" << sum_smear << ")";
		genleg << "gen (" << sum_gen << ")";

		TLegend *l2 = new TLegend(0.6,0.6,0.9,0.9);
		l2->AddEntry(hreco, recoleg.str().c_str());
		//l2->AddEntry(hgen, genleg.str().c_str());
		l2->AddEntry(hsmear, smearleg.str().c_str());

		c1_1->cd();
		gPad->SetLogy(hist2print.at(ihist).LogY());
	
		hreco->DrawCopy();
		//hgen->DrawCopy("same");
		hsmear->DrawCopy("same");
		l2->Draw();
		//tx->Draw();
		c1_2->cd();
		hsmeartoreco_ratio->DrawCopy();
		c->cd();
		gPad->Print("samples.eps");

	}

	gPad->Print("samples.eps]");
}

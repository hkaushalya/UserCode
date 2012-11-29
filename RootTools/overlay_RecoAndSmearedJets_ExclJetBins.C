#include <iostream>
#include <string>
#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include <sstream>
#include "TStyle.h"
#include "TLegend.h"
#include "assert.h"

/*************************************
 * overlays RECO and SMEARED spectrums
 * and plot the ratios of all exclusive
 * njet bins.
 * **********************************/

using namespace std;

const static float fDATA_LUMI = 10000; //pb-1
const Int_t nBins = 4;
TFile *files[nBins];
struct hist_t {
	string name;
	string title;
	int rebin;
	int normalizeByBinWidth;
};

void OpenFiles()
{
	/*files[0] = new TFile ("11012012_Excl2J_MhtFromEta5p0/qcd_all.root");
	files[1] = new TFile ("11012012_Excl3to5J_MhtFromEta5p0/qcd_all.root");
	files[2] = new TFile ("11012012_Excl6to7J_MhtFromEta5p0/qcd_all.root");
	files[3] = new TFile ("11012012_Excl8toInfJ_MhtFromEta5p0/qcd_all.root");
	*/

	/*files[0] = new TFile ("11012012_Excl2J_MhtFromEta2p5/qcd_all.root");
	files[1] = new TFile ("11012012_Excl3to5J_MhtFromEta2p5/qcd_all.root");
	files[2] = new TFile ("11012012_Excl6to7J_MhtFromEta2p5/qcd_all.root");
	files[3] = new TFile ("11012012_Excl8toInfJ_MhtFromEta2p5/qcd_all.root");
	*/

	files[0] = new TFile ("Excl2Jets/qcd_all.root");
	files[1] = new TFile ("Excl3to5Jets/qcd_all.root");
	files[2] = new TFile ("Excl6to7Jets/qcd_all.root");
	files[3] = new TFile ("Excl8toInfJets/qcd_all.root");
	
	for (int i = 0; i < nBins; ++i)
	{
		if (files[i]->IsZombie())
		{
			cout << "QCD file # " << i << " not found!" <<  endl;
			assert (false);
		}
	}
}
void NormByBinWidth(TH1* h)
{
	//assuming sumw2() already called upon
	//no error propgation will be done
	for (int bin=1; bin<= h->GetNbinsX(); ++bin)
	{
		const float binw = h->GetBinWidth(bin);
		const float val = h->GetBinContent(bin); 
		const float err = h->GetBinError(bin);
		const float nval = val/binw;
		const float nerr = err/binw;
		h->SetBinContent(bin, nval);
		h->SetBinError(bin, nerr);
	}
}

vector<TH1*> GetHist(const string histname)
{
	//hists are already scaled to 10fb-1
	vector<TH1*> hists;
	for (int i=0; i< nBins; ++i)
	{
		TH1* h = dynamic_cast<TH1*> (files[i]->Get(histname.c_str()));
		if (h == NULL)
		{
			cout << "hist " << histname << " not found in file # " << i <<  " !" << endl;
			assert (false);
		}
		TH1* hist = dynamic_cast<TH1*> (h->Clone());
		hist->Sumw2();
		hist->SetLineWidth(2);
		hists.push_back(hist);
	}

	return hists;
}
void GetYRange(vector<vector<TH1* > > hists, double& ymin, double& ymax)
{
	ymin=99999999.9; ymax=0;
	for (int i=0; i < hists.size(); ++i)
	{
		for (int j=0; j< hists.at(0).size(); ++j)
		{
			double min = hists.at(i).at(j)->GetBinContent(hists.at(i).at(j)->GetMinimumBin());
			double max = hists.at(i).at(j)->GetBinContent(hists.at(i).at(j)->GetMaximumBin());
			if (max>ymax) ymax = max;
			if (min<ymin) ymin = min;
		}
	}

}


TCanvas* overlay_RecoSmeared_ExclJetBins(const string& folder, 
								const string& htrange, const hist_t& h
								)
{

	TLegend *leg  = new TLegend(0.6,0.65,0.9,0.9);
	leg->SetTextFont(42);
	TCanvas* c1 = new TCanvas("c1");
	gPad->SetLogy();
	gPad->SetTickx();
	gPad->SetTicky();

	vector<string> jetcoll;
	jetcoll.push_back("reco");
	jetcoll.push_back("smeared");

	stringstream title;
	title << htrange << ", L = 10 fb^{-1}" << ";" << h.title ;
	
	vector< vector<TH1* > > HISTS;
	for (unsigned j=0; j< jetcoll.size(); ++j)
	{
		stringstream histname;
		histname << folder << "/" << jetcoll.at(j) << h.name;
		cout << __LINE__ << ": Looking for hist: " << histname.str().c_str() << endl;
		vector<TH1*> hists = GetHist(histname.str());
		HISTS.push_back(hists);
	}

	double ymin=0, ymax = 0;
	GetYRange(HISTS,ymin, ymax); 
	if (ymin==0) ymin = 0.0005;

	for (int j=0; j<jetcoll.size(); ++j)
	{
		stringstream jetc;
		if (j==0) jetc << "RECO"; 
		else if (j==1) jetc << "SMEAR"; 

		for (int bin=0; bin<nBins; ++bin)
		{
			HISTS.at(j).at(bin)->Rebin(h.rebin);
			HISTS.at(j).at(bin)->SetTitle(title.str().c_str());
			HISTS.at(j).at(bin)->SetMarkerStyle(20+j);
			HISTS.at(j).at(bin)->SetStats(0);
			

			if (j==1) //all smeared plots 
			{
				HISTS.at(j).at(bin)->SetLineColor(kRed);
				HISTS.at(j).at(bin)->SetMarkerColor(kRed);
				
			}

			int mstyle = 21+bin;
			HISTS.at(j).at(bin)->SetMarkerStyle(mstyle);
			stringstream legname;
			if (bin==0) legname << "==2 Jets";
			else if (bin==1) legname << "3-5 Jets";
			else if (bin==2) legname << "6-7 Jets";
			else if (bin==3) legname << "#geq 8 Jets";

			if (j==0 && bin==0) 
			{
				HISTS.at(j).at(bin)->GetYaxis()->SetRangeUser(ymin*0.5, ymax*10);
				HISTS.at(j).at(bin)->Draw();
			} else  HISTS.at(j).at(bin)->Draw("same");

			const double sum = HISTS.at(j).at(bin)->Integral(1, HISTS.at(j).at(bin)->GetNbinsX()+1); 
			legname << "[" << jetc.str() << "] (" << sum << ")";
			leg->AddEntry(HISTS.at(j).at(bin), legname.str().c_str());
		}
	}
	leg->Draw();

	
	/*TH1* ratio = dynamic_cast<TH1*> (hists.at(2)->Clone("ratio"));
	ratio->GetYaxis()->SetTitle("Smeared/Reco");
	ratio->SetTitle("");
	ratio->Divide(hists.at(0));
	//ratio->GetYaxis()->SetRangeUser(-0.01,2.01);
	ratio->GetYaxis()->SetRangeUser(0.49,1.51);
	//ratio->SetTickLength (+0.01,"Y");

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
	c1_1->SetLogy();
   
	hists.at(0)->GetYaxis()->CenterTitle(1);
	hists.at(0)->SetLabelFont(42,"XYZ");
	hists.at(0)->SetTitleFont(42,"XYZ");
	hists.at(0)->GetYaxis()->SetTitleOffset(0.8);
	hists.at(0)->SetLabelSize(0.05,"XYZ");
	hists.at(0)->SetTitleSize(0.06,"XYZ");
   hists.at(0)->Draw("P");
   hists.at(2)->Draw("same P");
	leg->Draw();

   c1_1->Modified();
   c1->cd();
  
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
   
	ratio->GetYaxis()->SetTitleOffset(0.4);
	ratio->GetXaxis()->SetTitleOffset(0.9);
	ratio->GetYaxis()->CenterTitle(1);
	ratio->GetXaxis()->CenterTitle(1);
	ratio->SetLabelSize(0.125,"XYZ");
	ratio->SetTitleSize(0.125,"XYZ");
//	ratio->SetLabelFont(labelfont,"XYZ");
//	ratio->SetTitleFont(titlefont,"XYZ");
   ratio->GetXaxis()->SetTickLength(0.07);
   ratio->Draw("");
   
   c1_2->Modified();
   c1->cd();
   //c1->Modified();
   //c1->cd();
   //c1->SetSelected(c1);
	*/
	return c1;
	//return gPad;

}

void overlay_RecoSmeared_ExclJetBins()
{
	vector<string> folders, htbinlabels, hists;
	//folders.push_back("Hist/HT0to7000MHT0to7000");
	//folders.push_back("Hist/HT500to7000MHT0to7000");
	folders.push_back("Hist/HT500to750MHT0to8000");
	folders.push_back("Hist/HT750to1000MHT0to8000");
	folders.push_back("Hist/HT1000to1250MHT0to8000");
	folders.push_back("Hist/HT1250to1500MHT0to8000");
	folders.push_back("Hist/HT1500to8000MHT0to8000");
	/*folders.push_back("Hist/HT0to500MHT0to7000");
	folders.push_back("Hist/HT500to900MHT0to7000");
	folders.push_back("Hist/HT900to1300MHT0to7000");
	folders.push_back("Hist/HT1300to7000MHT0to7000");
	*/

	vector<string> jetcoll;
	jetcoll.push_back("reco");
	jetcoll.push_back("smeared");

	vector<hist_t> histlist;
	for (int i=0; ;++i)
	{
		hist_t h;
		h.name = "";
		h.normalizeByBinWidth = 0;
		if (i==0) {
			h.name = "_mht";
			h.title = "MHT;Events";
			h.rebin = 10;
		} else if (i==1) {
			h.name = "_ht";
			h.title = "HT;Events";
			h.rebin = 2;
		} else if (i==2) {
			h.name = "jet1_pt";
			h.title = "P_{T} [Jet1^{E_{T}>30, |#eta |<5.0}];Events;";
			h.rebin = 1;
		} else if (i==3) {
			h.name = "jet2_pt";
			h.title = "P_{T} [Jet2^{E_{T}>30, |#eta |<5.0}];Events;";
			h.rebin = 1;
		} else if (i==4) {
			h.name = "jet3_pt";
			h.title = "P_{T} [Jet3^{E_{T}>30, |#eta |<5.0}];Events;";
			h.rebin = 1;
		} else if (i==5) {
			h.name = "jet1_eta";
			h.title = "#eta [Jet1^{E_{T}>30, |#eta |<5.0}];Events;";
			h.rebin = 2;
		} else if (i==6) {
			h.name = "jet2_eta";
			h.title = "#eta [Jet2^{E_{T}>30, |#eta |<5.0}];Events;";
			h.rebin = 2;
		} else if (i==7) {
			h.name = "jet3_eta";
			h.title = "#eta [Jet3^{E_{T}>30, |#eta |<5.0}];Events;";
			h.rebin = 2;
		} else if (i==8) {
			h.name = "_njet50eta2p5";
			h.title = "N Jets [E_{T}>50 GeV && |#eta |<2.5];Events";
			h.rebin = 1;
		} else if (i==9) {
			h.name = "_dphimin";
			h.title = "#Delta #Phi_{min};Events";
			h.rebin = 5;
			h.normalizeByBinWidth = 1;
		} else if (i==10) {
			h.name = "_njet30eta5p0";
			h.title = "N Jets [E_{T}>30 GeV && |#eta |<5.0];Events";
			h.rebin = 1;
		} else if (i==11) {
			h.name = "jet1_dphi";
			h.title = "#Delta #Phi (Jet 1^{ pt>50,GeV | #eta | < 2.5}, MHT);Events";
			h.rebin = 2;
		} else if (i==12) {
			h.name = "jet2_dphi";
			h.title = "#Delta #Phi (Jet 2^{ pt>50 GeV, | #eta | < 2.5}, MHT);Events";
			h.rebin = 2;
		} else if (i==13) {
			h.name = "jet3_dphi";
			h.title = "#Delta #Phi (Jet 3^{ pt>50 GeV, |#eta | < 2.5}, MHT);Events";
			h.rebin = 2;

	/*	} else if (i==10) {
			h.name = "pass0";
			h.title = "MHT;Events with #Delta #Phi_{min}>0.15";
			h.rebin = 1;
			h.normalizeByBinWidth = 1;
		} else if (i==11) {
			h.name = "pass1";
			h.title = "MHT;Events with #Delta #Phi_{min}>0.20";
			h.rebin = 1;
			h.normalizeByBinWidth = 1;
		} else if (i==12) {
			h.name = "pass3";
			h.title = "MHT;Events with #Delta #Phi_{min}>0.30";
			h.rebin = 1;
			h.normalizeByBinWidth = 1;
		} else if (i==13) {
			h.name = "fail0";
			h.title = "MHT;Events with #Delta #Phi_{min}<0.15";
			h.rebin = 1;
			h.normalizeByBinWidth = 1;
		} else if (i==14) {
			h.name = "fail1";
			h.title = "MHT;Events with #Delta #Phi_{min}<0.20";
			h.rebin = 1;
			h.normalizeByBinWidth = 1;
		} else if (i==15) {
			h.name = "fail3";
			h.title = "MHT;Events with #Delta #Phi_{min}<0.30";
			h.rebin = 1;
			h.normalizeByBinWidth = 1;
		} else if (i==16) {
			h.name = "signal";
			h.title = "MHT;Events Passing RA2 #Delta #Phi cut";
			h.rebin = 1;
			h.normalizeByBinWidth = 1;
		} else if (i==17) {
			h.name = "_njet30eta5p0";
			h.title = "N Jets [E_{T}>30 GeV && |#eta |<5.0];Events";
			h.rebin = 1;
*/		} else { break; }
		histlist.push_back(h);
	}

	OpenFiles();
	TCanvas *c = new TCanvas("print");	
	c->Draw();
	c->Print("reco_smear_compare_ExclJetBins.eps[");

	for (unsigned h = 0; h < histlist.size(); ++h)
	{
		//if (h!=0) continue;
		cout << histlist.at(h).name << endl;

		for (unsigned j = 0; j < folders.size(); ++j)
		{
			string htrange("");
			//if (j==0) htrange += "0<HT<500 GeV";
			if (j==0) htrange += "500<HT<750 GeV, MHT>0 GeV";
			else if (j==1) htrange += "750<HT<1000 GeV";
			else if (j==2) htrange += "1000<HT<1250 GeV";
			else if (j==3) htrange += "1250<HT<1500 GeV";
			else if (j==4) htrange += "1500<HT<8000 GeV";


			//stringstream histName;
			//histName << folders.at(j) << "/" << histlist.at(h).name;;
			//histName << histlist.at(h).name;;
			//cout << "hist name = " << histName.str() << endl;
			TCanvas* c1 = overlay_RecoSmeared_ExclJetBins(folders.at(j), htrange, histlist.at(h)); 
			c1->Print("reco_smear_compare_ExclJetBins.eps");
		}
	}
	c->Print("reco_smear_compare_ExclJetBins.eps]");
}


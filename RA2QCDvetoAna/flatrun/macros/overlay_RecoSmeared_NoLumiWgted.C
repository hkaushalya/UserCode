#include <iostream>
#include <string>
#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include <sstream>
#include "TStyle.h"
#include "TLegend.h"
#include "assert.h"

const static float fDATA_LUMI = 10000; //pb-1
const Int_t nBins = 7;
TFile *files[nBins];
struct hist_t {
	string name;
	string title;
	int rebin;
	int normalizeByBinWidth;
	int log;
};

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


using namespace std;


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
TH1* GetHist(const string histname)
{
	const float scaleTo = fDATA_LUMI; // pb

	TH1 *hists[nBins] = {0,0,0,0,0,0,0};
	TH1 *Hist = 0;

	for (int i=0; i<nBins; ++i)
	{
		hists[i] = dynamic_cast<TH1*> (files[i]->Get(histname.c_str()));
		if (hists[i] != 0)
		{
			hists[i]->Sumw2();
			const float scale = scaleTo/ ( nEvents[i] / xSec[i] );
			hists[i]->Scale(scale);
			hists[i]->SetLineWidth(2);
			hists[i]->SetLineColor(i);

			if (i == 0) 
			{
				Hist = dynamic_cast<TH1*> (hists[i]->Clone());
				Hist->SetDirectory(0);
			} else
			{
				Hist->Add(hists[i]);
			}

		} else {
			cout << "hist " << histname << " not found in " << files[i]->GetName() << "!" << endl;
			assert (false);
		}
	}

	return Hist;
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

TCanvas* overlay_RecoSmeared_NoLumiWgted(const string folder, const string htrange,
								const hist_t histinfo)
{
	TLegend *leg  = new TLegend(0.6,0.65,0.9,0.9);
	leg->SetTextFont(42);
	vector<TH1*> hists;
	
	stringstream title;
	title << htrange << " [" << histinfo.name << "];" << histinfo.title;
	cout << title.str() << endl;
	stringstream recohistname, smearhistname;
	recohistname << folder << "/reco" << histinfo.name; 
	smearhistname << folder << "/smeared" << histinfo.name; 
	TH1* recoHist = GetHist(recohistname.str());
	TH1* smearHist = GetHist(smearhistname.str());

	cout << recoHist->GetName() << ": int  = " << recoHist->Integral() << endl;
	cout << smearHist->GetName() << ": int  = " << smearHist->Integral() << endl;

	stringstream epsname;
	epsname << folder << "_" << histinfo.name;

	//Hist->SetTitle(title.str().c_str());
	recoHist->SetMarkerStyle(20);
	smearHist->SetMarkerStyle(22);
	recoHist->SetLineColor(kRed);
	recoHist->SetMarkerColor(kRed);
	recoHist->SetStats(0);
	recoHist->Rebin(histinfo.rebin);




	TCanvas *c1 = new TCanvas("c1");
	if (histinfo.log) gPad->SetLogy();
	gPad->SetTickx();
	gPad->SetTicky();

	recoHist->Draw("P");
	smearHist->Draw("same P");
	
	return c1;
}

void overlay_RecoSmeared_NoLumiWgted()
{
	vector<string> folders, htbinlabels, hists;
	folders.push_back("Hist/HT0to7000MHT0to7000");
	/*folders.push_back("Hist/HT0to500MHT0to7000");
	folders.push_back("Hist/HT500to900MHT0to7000");
	folders.push_back("Hist/HT900to1300MHT0to7000");
	folders.push_back("Hist/HT1300to7000MHT0to7000");
	*/

	vector<string> jetcoll;
	jetcoll.push_back("reco");
	jetcoll.push_back("gen");
	jetcoll.push_back("smeared");

	vector<hist_t> histlist;
		for (int i=0; ;++i)
		{
			hist_t h;
			h.name = "";
			h.normalizeByBinWidth = 0;
			h.log = 0;
			if (i==0) {
				h.name += "_mht";
				h.title = "MHT;Events";
				h.rebin = 10;
			} else if (i==1) {
				h.name += "_ht";
				h.title = "HT;Events";
				h.rebin = 2;
			} else if (i==2) {
				h.name += "jet1_pt";
				h.title = "P_{T} [Jet1^{E_{T}>30, |#eta |<5.0}];Events;";
				h.rebin = 1;
			} else if (i==3) {
				h.name += "jet2_pt";
				h.title = "P_{T} [Jet2^{E_{T}>30, |#eta |<5.0}];Events;";
				h.rebin = 1;
			} else if (i==4) {
				h.name += "jet3_pt";
				h.title = "P_{T} [Jet3^{E_{T}>30, |#eta |<5.0}];Events;";
				h.rebin = 1;
			} else if (i==5) {
				h.name += "jet1_eta";
				h.title = "#eta [Jet1^{E_{T}>30, |#eta |<5.0}];Events;";
				h.rebin = 2;
			} else if (i==6) {
				h.name += "jet2_eta";
				h.title = "#eta [Jet2^{E_{T}>30, |#eta |<5.0}];Events;";
				h.rebin = 2;
			} else if (i==7) {
				h.name += "jet3_eta";
				h.title = "#eta [Jet3^{E_{T}>30, |#eta |<5.0}];Events;";
				h.rebin = 2;
			} else if (i==8) {
				h.name += "_njet50eta2p5";
				h.title = "N Jets [E_{T}>50 GeV && |#eta |<2.5];Events";
				h.rebin = 1;
				h.log = 1;
			} else if (i==9) {
				h.name += "_dphimin";
				h.title = "#Delta #Phi_{min};Events";
				h.rebin = 5;
				h.normalizeByBinWidth = 1;
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
				h.name += "_njet30eta5p0";
				h.title = "N Jets [E_{T}>30 GeV && |#eta |<5.0];Events";
				h.rebin = 1;
		*/	} else { break; }
			histlist.push_back(h);
		}

	OpenFiles();
	TCanvas *c = new TCanvas("print");	
	c->Draw();
	c->Print("reco_smear_compare_fromNoLumiWgtd.eps[");

	for (unsigned h = 0; h < histlist.size(); ++h)
	{
	//	if (h!=0) continue;
		cout << histlist.at(h).name << endl;

		for (unsigned j = 0; j < folders.size(); ++j)
		{
			string htrange("");
			if (j==0) htrange += "0<HT<500 GeV";
			else if (j==1) htrange += "500<HT<900 GeV";
			else if (j==2) htrange += "900<HT<1300 GeV";
			else if (j==3) htrange += "1300<HT<7000 GeV";


			//stringstream histName;
			//histName << folders.at(j) << "/" << histlist.at(h).name;;
			//histName << histlist.at(h).name;;
			//cout << "hist name = " << histName.str() << endl;
			TCanvas* c1 = overlay_RecoSmeared_NoLumiWgted(folders.at(j), htrange, histlist.at(h)); 
			c1->Print("reco_smear_compare_fromNoLumiWgtd.eps");
		}
	}
	c->Print("reco_smear_compare_fromNoLumiWgtd.eps]");
}


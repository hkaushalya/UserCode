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
const Int_t nBins = 1;
TFile *files[nBins];
struct hist_t {
	string name;
	string title;
	int rebin;
	int normalizeByBinWidth;
};

void OpenFiles()
{

	files[0] = new TFile ("qcd_all.root");
	
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


TCanvas* overlay_RecoAndSmearedJets_ExclJetBins(const vector<pair<unsigned, unsigned> >& jetBins,
								const pair<float, float> htBin, const pair<float, float> mhtBin,
								
								const hist_t& h
								)
{

	TLegend *leg  = new TLegend(0.6,0.65,0.9,0.9);
	leg->SetTextFont(42);

//   TCanvas *c1 = new TCanvas("c1", "c1",15,60,550,600);
   TCanvas *c1 = new TCanvas("c1");
   c1->Range(0,0,1,1);
   c1->SetBorderSize(2);
	//TCanvas* c1 = new TCanvas("c1");
	//gPad->SetLogy();
	gPad->SetTickx();
	gPad->SetTicky();

	
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
	//c1_1->SetFillColor(7);

   c1_1->Modified();
   c1->cd();
// ------------>Primitives in pad: c1_2
	c1->cd();
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
	//c1_2->SetFillColor(9);


	stringstream htrange;
	htrange << htBin.first << "<HT<" << htBin.second << " GeV, "
		<< mhtBin.first << "<MHT<" << mhtBin.second << " GeV";

	stringstream title;
	title << htrange.str() << ", L = 10 fb ^{ -1}" << ";" << h.title ;

	vector<vector<TH1*> > HISTS;
	vector<TH1*> reco_hists, smeared_hists, ratio_hists;

	for (int jetbin=0; jetbin< jetBins.size(); ++jetbin)
	{
		stringstream folder;
		folder << "Hist/Njet" << jetBins.at(jetbin).first << "to" << jetBins.at(jetbin).second 
			<< "HT"   << htBin.first << "to" << htBin.second 
			<< "MHT"  << mhtBin.first << "to" << mhtBin.second;

		stringstream reco_histname, smear_histname;
		reco_histname << folder.str() << "/reco" << h.name;
		smear_histname << folder.str() << "/smeared" << h.name;
		cout << __LINE__ << ": Looking for hist: " << reco_histname.str().c_str() << endl;
		cout << __LINE__ << ": Looking for hist: " << smear_histname.str().c_str() << endl;

		TH1 *temp_recohist=0, *temp_smearhist=0, *temp_ratiohist=0;

		temp_recohist = dynamic_cast<TH1*> (files[nBins-1]->Get(reco_histname.str().c_str()));
		if (temp_recohist != 0)
		{
			temp_recohist->SetDirectory(0);
			reco_hists.push_back(temp_recohist);
		} else  { cout << __LINE__ << ": " << reco_histname.str() << " not found!" << endl; assert(false); }

		temp_smearhist = dynamic_cast<TH1*> (files[nBins-1]->Get(smear_histname.str().c_str()));
		if (temp_smearhist != 0) {
			temp_smearhist->SetDirectory(0);
			smeared_hists.push_back(temp_smearhist);
		} else  { cout << __LINE__ << ": " << smear_histname.str() << " not found!" << endl; assert(false); }

		//this is a copy of smearhist, so no need for another check
		cout << __LINE__ << endl;
		temp_ratiohist = dynamic_cast<TH1*> (files[nBins-1]->Get(smear_histname.str().c_str()));
		temp_ratiohist->SetDirectory(0);

		ratio_hists.push_back(temp_ratiohist);
		cout << __LINE__ << endl;
	}
	cout << __LINE__ << endl;
	HISTS.push_back(reco_hists);
	HISTS.push_back(smeared_hists);

	double ymin=0, ymax = 0;
	GetYRange(HISTS,ymin, ymax); 
	if (ymin==0) ymin = 0.0005;

	cout << __LINE__ << endl;
	for (int jetbin=0; jetbin<jetBins.size(); ++jetbin)
	{
		reco_hists.at(jetbin)->Rebin(h.rebin);
		reco_hists.at(jetbin)->SetTitle(title.str().c_str());
		reco_hists.at(jetbin)->SetMarkerStyle(24+jetbin);
		reco_hists.at(jetbin)->SetStats(0);

		smeared_hists.at(jetbin)->Rebin(h.rebin);
		smeared_hists.at(jetbin)->SetTitle(title.str().c_str());
		smeared_hists.at(jetbin)->SetMarkerStyle(24+jetbin);
		smeared_hists.at(jetbin)->SetStats(0);
		smeared_hists.at(jetbin)->SetLineColor(kRed);
		smeared_hists.at(jetbin)->SetMarkerColor(kRed);

		ratio_hists.at(jetbin)->Rebin(h.rebin);
		ratio_hists.at(jetbin)->SetMarkerStyle(24+jetbin);
		ratio_hists.at(jetbin)->SetStats(0);
		ratio_hists.at(jetbin)->SetLineColor(kRed);
		ratio_hists.at(jetbin)->SetMarkerColor(kRed);
		stringstream ratio_title;
		ratio_title << ";" << reco_hists.at(jetbin)->GetXaxis()->GetTitle() << "; Smeared/Reco;";
		ratio_hists.at(jetbin)->SetTitle(ratio_title.str().c_str());
		ratio_hists.at(jetbin)->SetStats(0);

		ratio_hists.at(jetbin)->Divide(reco_hists.at(jetbin));

		stringstream legname;
		legname << jetBins.at(jetbin).first << "-" << jetBins.at(jetbin).second << " Jets";

		stringstream reco_legname, smeared_legname;
		const double reco_sum = reco_hists.at(jetbin)->Integral(1, reco_hists.at(jetbin)->GetNbinsX()+1); 
		const double smeared_sum = smeared_hists.at(jetbin)->Integral(1, smeared_hists.at(jetbin)->GetNbinsX()+1); 
		reco_legname << legname.str() << "[RECO] (" << reco_sum << ")";
		smeared_legname << legname.str() << "[SMEAR] (" << smeared_sum << ")";
		leg->AddEntry(reco_hists.at(jetbin), reco_legname.str().c_str());
		leg->AddEntry(smeared_hists.at(jetbin), smeared_legname.str().c_str());
	}


	cout << __LINE__ << endl;
	c1_1->cd();
	for (int bin=0; bin<jetBins.size(); ++bin)
	{
		if (bin==0) 
		{
			reco_hists.at(bin)->GetYaxis()->SetRangeUser(ymin*0.5, ymax*10);
			reco_hists.at(bin)->Draw();
			smeared_hists.at(bin)->Draw("same");
		} else  
		{
			reco_hists.at(bin)->Draw("same");
			smeared_hists.at(bin)->Draw("same");
		}
	}

	cout << __LINE__ << endl;
	leg->Draw();
	c1->cd();
	c1_2->cd();
	for (int bin=0; bin<jetBins.size(); ++bin)
	{
		if (bin==0) 
		{
			ratio_hists.at(bin)->GetYaxis()->SetRangeUser(0.01, 1.99);
			ratio_hists.at(bin)->GetYaxis()->CenterTitle(1);
			ratio_hists.at(bin)->GetXaxis()->CenterTitle(1);
			ratio_hists.at(bin)->SetLabelFont(42,"XYZ");
			ratio_hists.at(bin)->SetTitleFont(42,"XYZ");
			ratio_hists.at(bin)->GetYaxis()->SetTitleOffset(0.4);
			ratio_hists.at(bin)->GetXaxis()->SetLabelOffset(0.02);
			ratio_hists.at(bin)->SetLabelSize(0.09,"XYZ");
			ratio_hists.at(bin)->SetTitleSize(0.10,"XYZ");
			ratio_hists.at(bin)->Draw();
		} else  ratio_hists.at(bin)->Draw("same");
	}

	cout << __LINE__ << endl;


	
	/*TH1* ratio = dynamic_cast<TH1*> (hists.at(2)->Clone("ratio"));
	ratio->GetYaxis()->SetTitle("Smeared/Reco");
	ratio->SetTitle("");
	ratio->Divide(hists.at(0));
	//ratio->GetYaxis()->SetRangeUser(-0.01,2.01);
	ratio->GetYaxis()->SetRangeUser(0.49,1.51);
	//ratio->SetTickLength (+0.01,"Y");

   
	hists.at(0)->GetYaxis()->CenterTitle(1);
	hists.at(0)->SetLabelFont(42,"XYZ");
	hists.at(0)->SetTitleFont(42,"XYZ");
	hists.at(0)->GetYaxis()->SetTitleOffset(0.8);
	hists.at(0)->SetLabelSize(0.05,"XYZ");
	hists.at(0)->SetTitleSize(0.06,"XYZ");
   hists.at(0)->Draw("P");
   hists.at(2)->Draw("same P");
	leg->Draw();

  
   
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

void overlay_RecoAndSmearedJets_ExclJetBins()
{
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
			h.rebin = 4;
		} else if (i==2) {
			h.name = "jet1_pt";
			h.title = "P_{T} [Jet1^{E_{T}>30, |#eta |<5.0}];Events;";
			h.rebin = 2;
		} else if (i==3) {
			h.name = "jet2_pt";
			h.title = "P_{T} [Jet2^{E_{T}>30, |#eta |<5.0}];Events;";
			h.rebin = 2;
		} else if (i==4) {
			h.name = "jet3_pt";
			h.title = "P_{T} [Jet3^{E_{T}>30, |#eta |<5.0}];Events;";
			h.rebin = 2;
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
			h.title = "#Delta #Phi [Jet 1^{ pt>50,GeV | #eta | < 2.5}, MHT];Events";
			h.rebin = 2;
		} else if (i==12) {
			h.name = "jet2_dphi";
			h.title = "#Delta #Phi [Jet 2^{ pt>50 GeV, | #eta | < 2.5}, MHT];Events";
			h.rebin = 2;
		} else if (i==13) {
			h.name = "jet3_dphi";
			h.title = "#Delta #Phi [Jet 3^{ pt>50 GeV, |#eta | < 2.5}, MHT];Events";
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

	//jetbins
	//
	vector<pair<unsigned, unsigned> > jetBins;
	vector<pair<float, float> > htBins, mhtBins;
	pair<float, float> jetbin1(2,2);	
	pair<float, float> jetbin2(3,5);	
	pair<float, float> jetbin3(6,7);	
	pair<float, float> jetbin4(8,1000);	

	pair<float, float> htbin1(500,750);	
	pair<float, float> htbin2(750,1000);	
	pair<float, float> htbin3(1000,1250);	
	pair<float, float> htbin4(1500,8000);	


	pair<float, float> mhtbin1(0,8000);	
	
	jetBins.push_back(jetbin1);
	jetBins.push_back(jetbin2);
	jetBins.push_back(jetbin3);
	jetBins.push_back(jetbin4);
	
	htBins.push_back(htbin1);
	htBins.push_back(htbin2);
	htBins.push_back(htbin3);
	htBins.push_back(htbin4);

	mhtBins.push_back(mhtbin1);

	for (unsigned h = 0; h < histlist.size(); ++h)
	{
			for (unsigned htbin =0; htbin < htBins.size(); ++htbin)
			{
				for (unsigned mhtbin =0; mhtbin < mhtBins.size(); ++mhtbin)
				{
					//if (h!=0) continue;
					cout << histlist.at(h).name << endl;
					TCanvas* c1 = overlay_RecoAndSmearedJets_ExclJetBins(jetBins, htBins.at(htbin), mhtBins.at(mhtbin), histlist.at(h)); 
					c1->Print("reco_smear_compare_ExclJetBins.eps");

				}
			}
	}
	c->Print("reco_smear_compare_ExclJetBins.eps]");
}


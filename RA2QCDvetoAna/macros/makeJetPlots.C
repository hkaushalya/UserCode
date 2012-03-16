#include <iostream>
#include <sstream>
#include <string>
#include "TFile.h"
#include "TH1.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
/*
 * To compare  data with smeared jets from data.
 *
 */

void makeJetPlots(const std::string histname)
{

	TFile *f1 = new TFile("Data.root");
	TH1 *h1 = (TH1*) f1->Get(histname.c_str());

	TFile *f2 = new TFile("Data_SmearedJetsNoPrescaleWgts.root");
	TH1 *h2 = (TH1*) f2->Get(histname.c_str());

	h2->SetLineColor(kRed);
	h2->SetLineWidth(2);
	h1->SetLineWidth(2);
	TLegend *leg = new TLegend(0.65,0.7,0.9,0.9);
	leg->AddEntry(h1,"DATA");
	leg->AddEntry(h2,"DATA-Smeared");

	h1->Rebin(4);
	h2->Rebin(4);

/*	new TCanvas();
	gPad->SetLogy();
	h1->Fit("gaus","","",200,350);
	h1->DrawCopy();
	
	new TCanvas();
	gPad->SetLogy();
	h2->Fit("gaus","","",200,350);
	h2->DrawCopy();
*/
	//return;
	gStyle->SetOptStat(0);
	new TCanvas();
	gPad->SetLogy();
	h1->Draw();
	h2->Draw("same");
	leg->Draw();
	
	stringstream epsname;
	epsname << h1->GetName() << ".eps";
	gPad->Print(epsname.str().c_str());

}

void makeJetPlots()
{

	makeJetPlots("factorization_ht500/pf30_jet1_pt");
	makeJetPlots("factorization_ht500/pf30_jet2_pt");
	makeJetPlots("factorization_ht500/pf30_jet3_pt");
	return;
	makeJetPlots("factorization_ht500/pf30_jet1_eta");
	makeJetPlots("factorization_ht500/pf30_jet2_eta");
	makeJetPlots("factorization_ht500/pf30_jet3_eta");
}

#include <iostream>
#include <sstream>
#include <vector>
#include "TCanvas.h"
#include "TPad.h"
#include "TFile.h"
#include "TH1.h"
#include "TStyle.h"

void printJERs() {

	vector<double> ptBins;
	vector<double> etaBins;
	ptBins.push_back(0);
	ptBins.push_back(20);
	ptBins.push_back(30);
	ptBins.push_back(50);
	ptBins.push_back(80);
	ptBins.push_back(120);
	ptBins.push_back(170);
	ptBins.push_back(230);
	ptBins.push_back(300);
	ptBins.push_back(380);
	ptBins.push_back(470);
	ptBins.push_back(570);
	ptBins.push_back(680);
	ptBins.push_back(800);
	ptBins.push_back(1000);
	ptBins.push_back(1300);
	ptBins.push_back(1700);
	ptBins.push_back(2200);
	ptBins.push_back(2800);
	ptBins.push_back(3500);

	etaBins.push_back(0); 
	etaBins.push_back(0.3); 
	etaBins.push_back(0.5); 
	etaBins.push_back(0.8); 
	etaBins.push_back(1.1); 
	etaBins.push_back(1.4); 
	etaBins.push_back(1.7); 
	etaBins.push_back(2.0); 
	etaBins.push_back(2.3); 
	etaBins.push_back(2.8); 
	etaBins.push_back(3.2); 
	etaBins.push_back(4.1); 
	etaBins.push_back(5.0);


	TFile f("qcd_all.root");

	TCanvas *c1 = new TCanvas("c1");
	gStyle->SetOptStat(11111);
	gPad->Print("JERs.eps[");

	for (unsigned ptbin=0; ptbin < ptBins.size() -1; ++ptbin)
	{
		for (unsigned etabin=0; etabin < etaBins.size() -1; ++etabin)
		{
			stringstream name;
			name << "Hist/JERS/jer_pt"<< ptBins.at(ptbin) << "to" << ptBins.at(ptbin+1) 
					<< "_eta" << etaBins.at(etabin) << "to" << etaBins.at(etabin+1);
			TH1* hist = (TH1*) f.Get(name.str().c_str());
			hist->Draw();
			gPad->Print("JERs.eps");
		}
	}
	gPad->Print("JERs.eps]");
}

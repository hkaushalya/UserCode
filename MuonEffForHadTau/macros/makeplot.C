#include<iostream>
#include<sstream>
#include "TFile.h"
#include "TCanvas.h"
#include <vector>
#include <string>
#include <utility>
#include "TH1.h"
void makeplot(const string filename)
{
	
	TFile f(filename.c_str());
	if (f.IsZombie()) return;

	vector <pair<double, double> > vJetBins_, vHtBins_;
	vector<string> histNames;

	vJetBins_.push_back(make_pair(2,2));
	vJetBins_.push_back(make_pair(3,5));
	vJetBins_.push_back(make_pair(6,1000));
	vHtBins_.push_back(make_pair(500,8000));

	histNames.push_back("nvtx");
	histNames.push_back("ht");
	histNames.push_back("mht");
	histNames.push_back("njet");
	histNames.push_back("muEff_num");
	histNames.push_back("muEff_den");
	histNames.push_back("matched_MuPtRatio");
	histNames.push_back("matched_MuERatio");
	histNames.push_back("matched_MuEDiff");
	histNames.push_back("matched_MuDelR");

	stringstream epsname;
	TCanvas *c = new TCanvas();
	c->Print("epsname.eps[");
	for (unsigned ihist=0; ihist<histNames.size(); ++ihist)
	{
		for (unsigned jetbin=0; jetbin<vJetBins_.size(); ++jetbin)
		{
			unsigned minjet= vJetBins_.at(jetbin).first;
			unsigned maxjet= vJetBins_.at(jetbin).second;
			for (unsigned htbin=0; htbin<vHtBins_.size(); ++htbin)
			{
				unsigned minht= vHtBins_.at(htbin).first;
				unsigned maxht= vHtBins_.at(htbin).second;
				stringstream folder, histname;
				folder << "Njet" << minjet << "to" << maxjet << "Ht" << minht << "to" << maxht;
				histname << folder.str() << histNames.at(ihist);
				TH1* hist = (TH1*) f.Get(histname.str().c_str());
				hist->Draw();
				c->Print("epsname.eps");
			}
		}
	}
	c->Print("epsname.eps]");
}

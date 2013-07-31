//#include<commonheader.h>

#include<iostream>
#include<string>
#include<vector>
#include<utility>
#include<iomanip>
#include<sstream>
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TStyle.h"
using namespace std;

void makeMuEff(const string filename="") 
{
	if (filename.length()<1)
	{
		cout << "Need input root file!" << endl;
		return;
	}

	vector< pair<double,double> >  vJetBins, vHtBins;
	vJetBins.push_back(make_pair(2,2));
	vJetBins.push_back(make_pair(3,5));
	vJetBins.push_back(make_pair(6,7));
	vJetBins.push_back(make_pair(8,1000));
	vHtBins.push_back(make_pair(500,8000));

	TFile f(filename.c_str());
	if (f.IsZombie())
	{
		cout << "File with name " << filename << " cannot be opened or does not exist!" << endl;
		return;
	}

	const string numname("muEff_num");
	const string denname("muEff_den");

	TCanvas *c = new TCanvas();
	size_t dotpos = filename.find(".");
	string substr = filename.substr(0,dotpos);
	stringstream epsname_o, epsname, epsname_c;
	epsname_o << substr << ".eps[";
	epsname << substr << ".eps";
	epsname_c << substr << ".eps]";
	c->Print(epsname_o.str().c_str());
	gStyle->SetOptStat(0);
	gStyle->SetPaintTextFormat("1.1f");
	for (unsigned jetbin=0; jetbin < vJetBins.size(); ++jetbin)
	{
		double minjet= vJetBins.at(jetbin).first;
		double maxjet= vJetBins.at(jetbin).second;
		for (unsigned htbin=0; htbin < vHtBins.size(); ++htbin)
		{
			double minht= vHtBins.at(htbin).first;
			double maxht= vHtBins.at(htbin).second;
			stringstream folder, numhistname, denhistname;
			folder << "muoneff/Hist/Njet" << minjet << "to" << maxjet << "Ht" << minht << "to" << maxht << "/";
			numhistname << folder.str() << numname;
			denhistname << folder.str() << denname;
			TH1* numhist = dynamic_cast<TH1*>( f.Get(numhistname.str().c_str()));
			TH1* denhist = dynamic_cast<TH1*> (f.Get(denhistname.str().c_str()));

			if (numhist == NULL) 
			{
				cout << "Numerator hist " << numhistname.str() << " is not found in file " << f.GetName() << endl;
				return;
			}

			if (denhist == NULL) 
			{
				cout << "Denominator hist " << denhistname.str() << " is not found in file " << f.GetName() << endl;
				return;
			}
			numhist->Draw("colz");
			c->Print(epsname.str().c_str());
			denhist->Draw("colz");
			c->Print(epsname.str().c_str());
			numhist->Sumw2();
			numhist->Divide(denhist);
			stringstream title;
			title << "Jets " << minjet << "-" << maxjet << ", HT>" << minht <<  ": Muon Reco+ID efficiency from " << substr << endl;
			numhist->SetTitle(title.str().c_str());
			//numhist->GetXaxis()->SetRangeUser(0,500);
			numhist->SetAxisRange(0,500,"X");
			numhist->SetAxisRange(0,3.5,"Y");
			//numhist->GetYaxis()->SetRangeUser(0,3.5);
			numhist->Draw("colzTEXT90E");
			c->Print(epsname.str().c_str());
		}
	}
	c->Print(epsname_c.str().c_str());

}

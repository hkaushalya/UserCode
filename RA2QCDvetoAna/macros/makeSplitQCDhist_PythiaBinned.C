#include <iostream>
#include <string>
#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include <sstream>
#include "TStyle.h"
#include "assert.h"

using namespace std;


void makeSplitQCDhist_PythiaBinned(const string histname, const string title, const bool logScale)
{

	const float scaleTo = 1; // pb

	//cross section for each sample in pt order
	const float cs[] = {	115100,
								24260,
								1168,
								70.22,
								15.55,
								1.844,
								0.3321,
								0.01087,
								0.0003575 
								};

	const float nevts[] = {6127528, //120 //ok
									6220160,	//170 //ok
									6432669,//300  //ok
									3990085,	//470 //ok
									4245695,	//600 //ok
									4053888,//800  //ok
									2093222,//1000 //ok
									2196200, //1400 //ok
									293139 //1800  //ok
	};
	
	
	

	const Int_t nBins = 9;
	TFile *files[nBins];
	TH1 *hists[nBins] = {0,0,0,0,0,0,0,0,0};

	files[0] = new TFile ("qcd1/Merged.root");
	files[1] = new TFile ("qcd2/Merged.root");
	files[2] = new TFile ("qcd3/Merged.root");
	files[3] = new TFile ("qcd4/Merged.root");
	files[4] = new TFile ("qcd5/Merged.root");
	files[5] = new TFile ("qcd6/Merged.root");
	files[6] = new TFile ("qcd7/Merged.root");
	files[7] = new TFile ("qcd8/Merged.root");
	files[8] = new TFile ("qcd9/Merged.root");

	TH1 *Hist = 0;
	int entries = 0;

	for (int i=0; i<nBins; ++i)
	{
		if (files[i]->IsZombie())
		{
			cout << files[i]->GetName() << " not found!" <<  endl;
			assert (false);
		} else
		{
			//cout << files[i]->GetName() << " opened." <<  endl;

			hists[i] = dynamic_cast<TH1*> (files[i]->Get(histname.c_str()));
			if (hists[i] != 0)
			{
				//hists[i]->Print();
				//new TCanvas();
				//hists[i]->Draw();
				//gPad->SetEditable(0);
				hists[i]->Sumw2();
				entries += hists[i]->GetEntries();
						
				const float scale = scaleTo/ ( nevts[i] / cs[i] );
				hists[i]->Scale(scale);
				hists[i]->SetLineWidth(2);
				hists[i]->SetLineColor(i);
				
				//if (i != 0) hists[0]->Add(hists[i]);
				if (i == 0) 
				{
					Hist = dynamic_cast<TH1*> (hists[i]->Clone("hist0_copy"));
				} else
				{
					Hist->Add(hists[i]);
				}

			} else {
				cout << "hist " << histname << " not found in " << files[i]->GetName() << "!" << endl;
				assert (false);
			}
		}
	}
	cout << "Entries found in " << histname << " =  " << entries << endl;

	gStyle->SetOptStat(0);
	//debug stff
	new TCanvas();
	gPad->SetGridy();
	if (logScale) gPad->SetLogy();
	Hist->SetLineColor(9);
	Hist->SetTitle(title.c_str());
	Hist->SetLineWidth(2);
	
	if (histname.find("phi") != std::string::npos)
	{
		Hist->GetXaxis()->SetRangeUser(0,4);
	}


	Hist->Draw();
	//debug
/*	for (int i=0; i<nBins; ++i)
	{
		hists[i]->SetLineWidth(2);
		hists[i]->SetLineColor(i);
		hists[i]->Draw("same");
	}
*/


	//debug
/*	hists[0]->SetTitle(title.c_str());
	hists[0]->SetLineWidth(2);
	new TCanvas();
	if (logScale) gPad->SetLogy();
	hists[0]->Draw();
*/
/*	Hist->SetTitle(title.c_str());
	Hist->SetLineWidth(2);
	new TCanvas();
	if (logScale) gPad->SetLogy();
	Hist->Draw();

*/
	stringstream epsname;
	epsname << hists[0]->GetName() << ".eps";
	gPad->Print(epsname.str().c_str());

	Hist->SetName(hists[0]->GetName());

	TFile f("Results.root","UPDATE");
	Hist->Write();
	f.Close();
	


}

void makeSplitQCDhist_PythiaBinned()
{
	makeSplitQCDhist_PythiaBinned("factorization_ht500/pf30_jet1_pt","HT>500 GeV;Jet1-P_{T};Events;", 1);
	makeSplitQCDhist_PythiaBinned("factorization_ht500/pf30_jet2_pt","HT>500 GeV;Jet2-P_{T};Events;", 1);
	makeSplitQCDhist_PythiaBinned("factorization_ht500/pf30_jet3_pt","HT>500 GeV;Jet3-P_{T};Events;", 1);
	makeSplitQCDhist_PythiaBinned("factorization_ht500/ht" , "HT>500 GeV;HT;Events" , 1);
	makeSplitQCDhist_PythiaBinned("factorization_ht500/mht", "HT>500 GeV;MHT;Events", 1);
	makeSplitQCDhist_PythiaBinned("factorization_ht500/meteff","HT>500 GeV;MEff;Events", 1);
	makeSplitQCDhist_PythiaBinned("factorization_ht500/njet_et50eta24"  , "HT>500 GeV;N Jets [E_{T}>50 GeV && |#eta|<2.5];Events", 1);
	makeSplitQCDhist_PythiaBinned("factorization_ht500/pf30_jet1_delphi", "HT>500 GeV: Jet-1 ;#Delta#Phi(jet-1, MHT);Events;", 0);
	makeSplitQCDhist_PythiaBinned("factorization_ht500/pf30_jet2_delphi", "HT>500 GeV: Jet-2 ;#Delta#Phi(jet-2, MHT);Events;", 0);
	makeSplitQCDhist_PythiaBinned("factorization_ht500/pf30_jet3_delphi", "HT>500 GeV: Jet-3 ;#Delta#Phi(jet-3, MHT);Events;", 0);
	makeSplitQCDhist_PythiaBinned("factorization_ht500/delPhiMin","HT>500 GeV;#Delta#Phi_{min};Events", 1);

}

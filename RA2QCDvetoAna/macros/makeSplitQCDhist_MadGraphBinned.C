#include <iostream>
#include <string>
#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include <sstream>
#include "TStyle.h"
#include "assert.h"

using namespace std;


void makeSplitQCDhist_MadGraphBinned(const string histname, const string title, const bool logScale)
{

	const float scaleTo = 1; // pb

	//cross section for each sample in pt order
	const float cs[] = {	2.78e8 * 0.15143 ,
								833057 * 0.23835,
								17132 * 0.34085,
								409.6 * 0.29932
								};

	/*const float nevts[] = {23739003,
								  20674219,
								  14437469,
								  6294851
	};*/
	const float nevts[] = {2.12e7,
								  1.99e7,
								  1.47e7,
								  7.33e6
	};
	
	
	

	const Int_t nBins = 4;
	TFile *files[nBins];
	TH1 *hists[nBins] = {0,0,0,0};

	files[0] = new TFile ("qcd1/Merged.root");
	files[1] = new TFile ("qcd2/Merged.root");
	files[2] = new TFile ("qcd3/Merged.root");
	files[3] = new TFile ("qcd4/Merged.root");

	TH1 *Hist = 0;

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
				//hists[i]->DrawCopy();
				//gPad->SetEditable(0);
				hists[i]->Sumw2();
						
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

	//gStyle->SetOptStat(0);
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

void makeSplitQCDhist_MadGraphBinned()
{
	//RA2b plots

	makeSplitQCDhist_MadGraphBinned("factorization_ht500/pf30_jet1_pt","HT>500 GeV;Jet1-P_{T};Events;",1);
	makeSplitQCDhist_MadGraphBinned("factorization_ht500/ht","HT>500 GeV;HT;Events",1);
	makeSplitQCDhist_MadGraphBinned("factorization_ht500/mht","HT>500 GeV;MHT;Events",1);
	makeSplitQCDhist_MadGraphBinned("factorization_ht500/meteff","HT>500 GeV;MEff;Events",1);
	makeSplitQCDhist_MadGraphBinned("factorization_ht500/pf30_jet2_pt","RA2;Jet2-P_{T};Events;",1);
	makeSplitQCDhist_MadGraphBinned("factorization_ht500/pf30_jet3_pt","RA2;Jet3-P_{T};Events;",1);
	makeSplitQCDhist_MadGraphBinned("factorization_ht500/njet_et50eta24","HT>500 GeV;N Jets [E_{T}>50 GeV && |#eta|<2.5];Events",1);

//   makeSplitQCDhist_MadGraphBinned("QCDvetoAna/metsig_delphinorm_fail","#slash{E}_{T}-Significance of Events Failing #Delta#Phi_{min}^{norm} cut;#slash{E}_{T}-Sig;Events;",0);
//	makeSplitQCDhist_MadGraphBinned("QCDvetoAna/metsig_delphinorm_pass","#slash{E}_{T}-Significance of Events Passing #Delta#Phi_{min}^{norm} cut;#slash{E}_{T}-Sig;Events;",0);
	//makeSplitQCDhist_MadGraphBinned("QCDvetoAna/pat_njet_et50eta25",";Njet (Et50Eta2.5);Events;",0);


}

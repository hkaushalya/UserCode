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
	makeSplitQCDhist_PythiaBinned("factorization_ht500/pf30_jet1_pt","HT>500 GeV;Jet1-P_{T};Events;",1);
	makeSplitQCDhist_PythiaBinned("factorization_ht500/ht","HT>500 GeV;HT;Events",1);
	makeSplitQCDhist_PythiaBinned("factorization_ht500/mht","HT>500 GeV;MHT;Events",1);
	makeSplitQCDhist_PythiaBinned("factorization_ht500/meteff","HT>500 GeV;MEff;Events",1);
	makeSplitQCDhist_PythiaBinned("factorization_ht500/pf30_jet2_pt","RA2;Jet2-P_{T};Events;",1);
	makeSplitQCDhist_PythiaBinned("factorization_ht500/pf30_jet3_pt","RA2;Jet3-P_{T};Events;",1);
	makeSplitQCDhist_PythiaBinned("factorization_ht500/njet_et50eta24","HT>500 GeV;N Jets [E_{T}>50 GeV && |#eta|<2.5];Events",1);



//	makeSplitQCDhist_PythiaBinned("factorization_ht350/pf30_jet1_pt","RA2;Jet1-P_{T};Events;",1);
//	makeSplitQCDhist_PythiaBinned("factorization_ht350/pf30_jet2_pt","RA2;Jet2-P_{T};Events;",1);
//	makeSplitQCDhist_PythiaBinned("factorization_ht350/pf30_jet3_pt","RA2;Jet3-P_{T};Events;",1);
/*	makeSplitQCDhist_PythiaBinned("factorization_ht500/ra2b_pf30_jet1_eta","RA2b;Jet1-#{Eta};Events;",0);
	makeSplitQCDhist_PythiaBinned("factorization_ht500/ra2b_pf30_jet2_eta","RA2b;Jet2-#{Eta};Events;",0);
	makeSplitQCDhist_PythiaBinned("factorization_ht500/ra2b_pf30_jet3_eta","RA2b;Jet3-#{Eta};Events;",0);
*/
/*	makeSplitQCDhist_PythiaBinned("factorization_ht500/ra2b_pf30_jet1_delTJetPt","RA2b;#Delta T;Events;",1);
	makeSplitQCDhist_PythiaBinned("factorization_ht500/ra2b_pf30_jet2_delTJetPt","RA2b;#Delta T;Events;",1);
	makeSplitQCDhist_PythiaBinned("factorization_ht500/ra2b_pf30_jet3_delTJetPt","RA2b;#Delta T;Events;",1);
*/
//	makeSplitQCDhist_PythiaBinned("factorization_ht500/ra2b_pf30_jet1_delTvsJetPt","RA2b;P_{T}^{Jet1};#Delta T;",1);
//	makeSplitQCDhist_PythiaBinned("factorization_ht500/ra2b_pf30_jet2_delTvsJetPt","RA2b;P_{T}^{Jet2};#Delta T;",1);
//	makeSplitQCDhist_PythiaBinned("factorization_ht500/ra2b_pf30_jet3_delTvsJetPt","RA2b;P_{T}^{Jet3};#Delta T;",1);


/*	makeSplitQCDhist_PythiaBinned("factorization_ht500/ra2b_pf30_jet1_delTDevidedByJetPtvsJetPt","RA2b;P_{T}^{Jet1};#Delta T;",1);
	makeSplitQCDhist_PythiaBinned("factorization_ht500/ra2b_pf30_jet2_delTDevidedByJetPtvsJetPt","RA2b;P_{T}^{Jet2};#Delta T;",1);
	makeSplitQCDhist_PythiaBinned("factorization_ht500/ra2b_pf30_jet3_delTDevidedByJetPtvsJetPt","RA2b;P_{T}^{Jet3};#Delta T;",1);
*/
//	makeSplitQCDhist_PythiaBinned("factorization_ht500/ra2b_delPhiMin","RA2b;#Delta#Phi_{min};Events;",0);
//	makeSplitQCDhist_PythiaBinned("factorization_ht500/ra2b_delPhiMinNorm","RA2b;#Delta#Phi_{min}^{norm};Events;",0);
//	makeSplitQCDhist_PythiaBinned("factorization_ht500/ra2b_delPhiMinVsMET","RA2b;MET; Avg. #Delta#Phi_{min};",1);
//	makeSplitQCDhist_PythiaBinned("factorization_ht500/ra2b_delPhiMinNormVsMET","RA2b;MET; Avg. #Delta#Phi_{min}^{norm};",1);

/*	makeSplitQCDhist_PythiaBinned("factorization_ht500/ra2b_delPhiMinVsMET","RA2b;#slash{E}_{T};#Delta#Phi_{min};",0);
	makeSplitQCDhist_PythiaBinned("factorization_ht500/ra2b_delPhiMinNormVsMET","RA2b;#slash{E}_{T};#Delta#Phi_{min}^{norm};",0);

	makeSplitQCDhist_PythiaBinned("factorization_ht500/ra2b_PassFail","RA2b: Using #Delta#Phi_{min} ;#slash{E}_{T};PASS/FAIL Ratio;",0);
	makeSplitQCDhist_PythiaBinned("factorization_ht500/ra2b_PassFail_Norm","RA2b: Using #Delta#Phi_{min}^{norm} ;#slash{E}_{T};PASS/FAIL Ratio;",0);
*/
	//ra2 plots

/*	makeSplitQCDhist_PythiaBinned("factorization_ht500/pat_mht",";MHT;Events;",1);
	makeSplitQCDhist_PythiaBinned("factorization_ht500/pat_ht",";HT;Events",1);
	makeSplitQCDhist_PythiaBinned("factorization_ht500/pat_meff",";MEff;Events;",1);
*/
//	makeSplitQCDhist_PythiaBinned("factorization_ht500/metsig_delphi_fail","#slash{E}_{T}-Significance of Events Failing #Delta#Phi_{min} cut;#slash{E}_{T}-Sig;Events;",0);
//	makeSplitQCDhist_PythiaBinned("factorization_ht500/metsig_delphi_pass","#slash{E}_{T}-Significance of Events Passing #Delta#Phi_{min} cut;#slash{E}_{T}-Sig;Events;",0);

//   makeSplitQCDhist_PythiaBinned("factorization_ht500/metsig_delphinorm_fail","#slash{E}_{T}-Significance of Events Failing #Delta#Phi_{min}^{norm} cut;#slash{E}_{T}-Sig;Events;",0);
//	makeSplitQCDhist_PythiaBinned("factorization_ht500/metsig_delphinorm_pass","#slash{E}_{T}-Significance of Events Passing #Delta#Phi_{min}^{norm} cut;#slash{E}_{T}-Sig;Events;",0);
	//makeSplitQCDhist_PythiaBinned("factorization_ht500/pat_njet_et50eta25",";Njet (Et50Eta2.5);Events;",0);

/*	makeSplitQCDhist_PythiaBinned("factorization_ht500/delPhiMin_mht",";#Delta#Phi_{min};Events;",0);
	makeSplitQCDhist_PythiaBinned("factorization_ht500/delPhiMinNorm_mht",";#Delta#Phi_{min}^{norm};Events;",0);
*/
//	makeSplitQCDhist_PythiaBinned("factorization_ht500/hPassFail_mht","Using #Delta#Phi_{min} ;MHT;PASS/FAIL Ratio;",0);
//	makeSplitQCDhist_PythiaBinned("factorization_ht500/hPassFail_Norm_mht","Using #Delta#Phi_{min}^{norm};MHT;PASS/FAIL Ratio;",0);


//	makeSplitQCDhist_PythiaBinned("factorization_ht500/delPhiMin_met",";#Delta#Phi_{min};Events;",1);
//	makeSplitQCDhist_PythiaBinned("factorization_ht500/delPhiMinNorm_met",";#Delta#Phi_{min}^{norm};Events;",1);

//	makeSplitQCDhist_PythiaBinned("factorization_ht500/hPassFail_met","Using #Delta#Phi_{min} ;MET;PASS/FAIL Ratio;",0);
//	makeSplitQCDhist_PythiaBinned("factorization_ht500/hPassFail_Norm_met","Using #Delta#Phi_{min}^{norm};MET;PASS/FAIL Ratio;",0);


/*
	makeSplitQCDhist_PythiaBinned("factorization_ht500/pat_pf50eta25_jet1_pt","Jet-1 (PF50Eta2.5);E_{T};Events;",1);
	makeSplitQCDhist_PythiaBinned("factorization_ht500/pat_pf50eta25_jet2_pt","Jet-2 (PF50Eta2.5);E_{T};Events;",1);
	makeSplitQCDhist_PythiaBinned("factorization_ht500/pat_pf50eta25_jet3_pt","Jet-3 (PF50Eta2.5);E_{T};Events;",1);
	makeSplitQCDhist_PythiaBinned("factorization_ht500/pat_pf50eta25_jet1_eta","Jet-1 (PF50Eta2.5);#Eta;Events;",0);
	makeSplitQCDhist_PythiaBinned("factorization_ht500/pat_pf50eta25_jet2_eta","Jet-2 (PF50Eta2.5);#Eta;Events;",0);
	makeSplitQCDhist_PythiaBinned("factorization_ht500/pat_pf50eta25_jet3_eta","Jet-3 (PF50Eta2.5);#Eta;Events;",0);

	makeSplitQCDhist_PythiaBinned("factorization_ht500/pat_pf50eta25_jet1_delphi", "Jet-1 (PF50Eta2.5);#Delta#Phi(jet-1, MHT);Events;",0);
	makeSplitQCDhist_PythiaBinned("factorization_ht500/pat_pf50eta25_jet2_delphi", "Jet-2 (PF50Eta2.5);#Delta#Phi(jet-2, MHT);Events;",0);
	makeSplitQCDhist_PythiaBinned("factorization_ht500/pat_pf50eta25_jet3_delphi", "Jet-3 (PF50Eta2.5);#Delta#Phi(jet-3, MHT);Events;",0);
*/

}

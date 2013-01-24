#include "abcdMethod.hh"
#include <iostream>
#include "TPad.h"
#include <assert.h>
#include "TStyle.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TMath.h"
#include <sstream>
#include <ctype.h>
#include <iomanip>
#include <algorithm>
#include "exception"
#include <cmath>
#include <memory>
#include "TMatrixDSym.h"
#include "TFitResult.h"

string SystName(const int i)
{
	switch (i)
	{
		case 0: return "UpperTailScaling x 2.0";
		case 1: return "UpperTailScaling x 0.5";
		case 2: return "LowerTailScaling x 2.0";
		case 3: return "LowerTailScaling x 0.5";
		case 4: return "AdditionalSmearing +10%";
		case 5: return "AdditionalSmearing -10%";
		default: return "Unspecified systematic type!!!";
	}

}

vector<float> 
GetVarBinVector (float down, float up5, float up10, float up20, float up50, 
					float width1, float width2, float width3, float width4)
{
	vector<float> result;
	float point = down;
	const unsigned nregion = 4;
	const float up [] = {up5, up10, up20, up50};
	const float step [] = {width1, width2, width3, width4};		//1j pet

	result.push_back (point);
	for (unsigned region = 0; region != nregion; ++ region)
	{
		while (point + step[region] < up[region] + 1)
		{
			point += step[region];
			result.push_back (point);
	 	};
	};
	
	if (point + 1 < up[nregion-1])
		result.push_back (up[nregion-1]);
	
	return result;
};

void NormalizeBinWidth(TH1* hist)
{
	assert (hist != NULL && "requirement failed"); //spec
	assert (hist->GetDimension() == 1 && "CommonTools::NormalizeBinWidth:: hist is not 1-D");

	for (unsigned bin = 1; bin <= unsigned (hist->GetNbinsX()); ++ bin)
	{
		const float width = hist->GetXaxis()->GetBinWidth (bin);
		hist->SetBinContent (bin, hist->GetBinContent (bin) / width);
		hist->SetBinError (bin, hist->GetBinError (bin) / width);
	};

}

auto_ptr<TH1> FillVarBinHist (const vector<float>& bins, TH1 *input)
{
	assert (input != NULL && "requirement failed"); //spec
	assert (input->GetDimension() == 1 && "requirement failed"); //spec
	const unsigned nbin = unsigned (input->GetNbinsX ());

	auto_ptr<TH1> result (new TH1F (input->GetName(), input->GetTitle(),
						 bins.size() - 1, &bins[0]));

	result->SetBinContent (0, input->GetBinContent (0));
	result->SetBinError (0, input->GetBinError (0));
	result->SetBinContent (bins.size(), input->GetBinContent (nbin));
	result->SetBinError (bins.size(), input->GetBinError (nbin));
	for (unsigned bin = 1; bin != nbin; ++ bin)
	{
		const float low = input->GetXaxis()->GetBinLowEdge (bin);
		const float high = input->GetXaxis()->GetBinUpEdge (bin);
		const unsigned bin1 = result->FindBin (0.99 * low + 0.01 * high);
		const unsigned bin2 = result->FindBin (0.01 * low + 0.99 * high);
		assert (bin1 == bin2 && "requirement failed: histogram bin edges don't match"); //spec

		const float va = input->GetBinContent (bin);
		const float ea = input->GetBinError (bin);
		const float vb = result->GetBinContent (bin2);
		const float eb = result->GetBinError (bin2);
		const float vc = va + vb;
		const float ec = sqrt (ea * ea + eb * eb);
		result->SetBinContent (bin2, vc);
		result->SetBinError (bin2, ec);
	};

	assert (result.get() != NULL && "postcondition failed"); //spec
	return result;
}

TH1* MakeVariableBinHist (TH1 *hist, const float xmin, const float xpoint1,
							const float xpoint2, const float xpoint3, const float xpoint4,
				 			const float width1, const float width2, const float width3, 
							const float width4)
{
	assert (hist != NULL && "CommonTools::MakeVariableBinHist:: hist is NULL!");
	assert (hist->GetDimension() == 1 && "CommonTools::MakeVariableBinHist:: hist is not 1-D");
	
  	auto_ptr<TH1> result = FillVarBinHist ( 
											GetVarBinVector(xmin, xpoint1, xpoint2, 
															xpoint3, xpoint4, width1, 
															width2, width3, width4)
											, hist);
  	NormalizeBinWidth(result.get());
  	return result.release ();
};


//tokenizer
vector<string> FindAllOf(const string str, const string sub=":")
{
	const size_t len = str.length();
	vector<size_t> brac_pos;

	brac_pos.push_back(0);
	size_t pos = str.find(sub, 0);
	while(pos != string::npos)
	{
		brac_pos.push_back(pos);
		pos = str.find(sub,pos+1);
	}
	//check if last position overlaps with end of the string
	unsigned pos_size = brac_pos.size();
	if (brac_pos.at(pos_size-1) + 1 < len)
	{
		brac_pos.push_back(len);
	}
	
	vector<string> substrs;

	for (unsigned i=0; i < brac_pos.size()-1; ++i)
	{
		string ss;
		if (i==0) { 
			ss = str.substr(0, brac_pos.at(i+1));
		} else {
			ss = str.substr(brac_pos.at(i)+1, brac_pos.at(i+1)- brac_pos.at(i) - 1);
		}
		substrs.push_back(ss);
	}

	return substrs;
}

//one all strings are broken to x-y type use this to
//get numerical values
vector<pair<float, float> > GetVals(vector<string> strs, const string sep="-")
{
	vector <pair<float, float> > bins;
	for (unsigned i=0; i < strs.size(); ++i)
	{
		const string str = strs.at(i);
		const size_t len = str.length();
		const size_t pos = str.find_first_of(sep);
		
		const string str_1 = str.substr(0,pos);
		const string str_2 = str.substr(pos+1, len - pos);
		
		const float min = atof(str_1.c_str());
		const float max = atof(str_2.c_str());
		//cout << __FUNCTION__ << ": range for " << str <<  ": " << min << "->" << max <<endl;
		bins.push_back(make_pair(min, max));
	}

	return bins;
}

double GetFitFunctionError(TF1* f1, TFitResultPtr fitResPtr, double val) {
	TMatrixDSym covMat_ratio1pol2 = fitResPtr->GetCorrelationMatrix();
//	covMat_ratio1pol2.Print();

/*	for(int ir=0; ir<2; ir++){
		for(int ic=0; ic<2; ic++){
			std::cout<<"ir : "<<ir<<"  ic : "<<ic<<"  -->  "<<covMat_ratio1pol2(ir, ic)<<std::endl;
		}
	}
	std::cout<<std::endl;
*/
//	std::cout<<"\n\nprint out the errors (taking into account of correlation)...??? NOT complete next line!!!!!"<<std::endl;
	double parErrp0 = fitResPtr->ParError(0), parErrp1 = fitResPtr->ParError(1);
//	printf("parErrp0 : %9.5e  parErrp1 : %9.5e\n", parErrp0, parErrp1);

	const double err = sqrt( val*val*parErrp1*parErrp1*covMat_ratio1pol2(1, 1)*covMat_ratio1pol2(1, 1) + parErrp0*parErrp0*covMat_ratio1pol2(0, 0)*covMat_ratio1pol2(0, 0) + 2*val*parErrp0*parErrp1*covMat_ratio1pol2(0, 1) );
	double centralVal = f1->Eval(val);

//	cout << __FUNCTION__ << ": val/err = " << centralVal << "/" << err << endl;

	return err;
}

double expFitFunc(double *x, double *par)
{
	double fitval=0.0;
	if(x[0]>0.0) fitval=par[0] * exp(par[1] * x[0]) + par[2];
	else fitval=-1.0E6;
	return fitval;
}

double gausFitFunc(double *x, double *par)
{
	double arg = par[0] * exp(par[1] * x[0]);
	double fitval = 1.0 / TMath::Erf (arg) - 1 + par[2];
	return fitval;
}

double GetAvgVal(TH1* hist, const float& xmin)
{
	assert(hist != NULL && "GetAvgVal:: hist is null!");
	assert(hist->GetDimension() == 1 && "GetAvgVal:: hist is not 1D!");
	
	double sum = 0, N = 0;

	for (int bin=1; bin<= hist->GetNbinsX(); ++bin)
	{
		if (hist->GetBinCenter(bin)< xmin) continue;
		
		//cout << "Adding bin " << bin << " = " << hist->GetBinContent(bin) <<endl; 
		sum += hist->GetBinContent(bin); 
		++N;
	}

	//cout << "avg = " << (sum/N) << endl;
	return (sum/N);

}


abcdMethod::abcdMethod()
{

	Nsyst = 6;
	b_doSystematics = true;
	//numHistName = "smeared_failFineBin1";
	//denHistName = "smeared_signalFineBin";
	numHistName = "smeared_signal";
	denHistName = "smeared_fail1";
	controlHistName = "smeared_failFineBin1";
	signalHistName = "smeared_signalFineBin";
	controlsyst1HistName = "smeared_sidebandSyst1_fineBin";
	//controlsyst2HistName = "smeared_sidebandSyst2_fineBin";
	
	canvas_for_ratios = new TCanvas("canvas_for_ratios","Canvas for Ratio Plots");
	gPad->SetLogy();
	gPad->SetTickx();
	//open eps files for printing all ratio hists
	epsfile_for_ratios = "ratioplots.eps";
	stringstream epsopen;
	epsopen << epsfile_for_ratios << "[";
	canvas_for_ratios->Print(epsopen.str().c_str());

	canvas_for_debug = new TCanvas("canvas_for_debug","Canvas for debug Plots");
	gPad->SetLogy();
	gPad->SetTickx();
	//open eps files for printing all debug hists
	epsfile_for_debug = "debugplots.eps";
	stringstream epsopenfordebug;
	epsopenfordebug << epsfile_for_debug << "[";
	canvas_for_debug->Print(epsopenfordebug.str().c_str());

}

abcdMethod::~abcdMethod()
{
	vg_Inputs.clear();
	stringstream epsclose, epsclosefordebug;
	epsclose << epsfile_for_ratios << "]";
	epsclosefordebug << epsfile_for_debug << "]";
	canvas_for_ratios->Clear();
	canvas_for_ratios->Print(epsclose.str().c_str());
	canvas_for_debug->Print(epsclosefordebug.str().c_str());
	delete canvas_for_ratios;
	delete canvas_for_debug;
}


TH1* abcdMethod::GetRatioHist(TFile* f, const string path_to_hist, const bool debug)
{

	stringstream numHistPath, denHistPath;
	denHistPath << path_to_hist << "/" << denHistName;

	int rebin = 1;

	string filepath(f->GetName());
	if (filepath.find("Mean") == string::npos )
	{
		cout << __FUNCTION__ << ": mean with diff name found " << endl;
		numHistPath << path_to_hist << "/" << numHistName;
	} else {
		numHistPath << path_to_hist << "/smeared_signal";
		rebin = 2;
	}
	TH1* h_num = dynamic_cast<TH1*> (f->Get(numHistPath.str().c_str()));
	TH1* h_den = dynamic_cast<TH1*> (f->Get(denHistPath.str().c_str()));

	if (h_num == NULL) { 
		cout << "numerator hist not found in file " << f->GetName() << " in path " << numHistPath.str() << endl; 
		assert (false);
	}
	if (h_den == NULL) { 
		cout << "denominator hist not found in file " << f->GetName() << " in path " << denHistPath.str() << endl; 
		assert (false);
	}

	if (h_num->GetNbinsX() != h_den->GetNbinsX()) {
		cout << __FUNCTION__ << ": Numerator hist (" << h_num->GetName() << ") and Denominator hist (" << h_den->GetName() 
				<< " found in file (" << f->GetName() << " have different number of bins!!! "
				<< h_num->GetNbinsX() << "/" << h_den->GetNbinsX() << endl;
		assert(false);
	}


	//const float xmin = 50, xpoint1 = 120, xpoint2 = 200, xpoint3 = 500, xpoint4 = 1100, width1 = 10, width2 = 20, width3 = 150, width4 = 300;
	//const float xmin = 50, xpoint1 = 160, xpoint2 = 300, xpoint3 = 400, xpoint4 = 1100, width1 = 10, width2 = 20, width3 = 100, width4 = 350;
	const float xmin = 50, xpoint1 = 160, xpoint2 = 300, xpoint3 = 400, xpoint4 = 1100, width1 = 10, width2 = 20, width3 = 50, width4 = 50;
	h_num = MakeVariableBinHist (h_num, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4);
	h_den = MakeVariableBinHist (h_den, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4);



	TH1* ratio2 = dynamic_cast<TH1*> (h_num->Clone("TEMP"));
	ratio2->SetDirectory(0);
	ratio2->Divide(h_den);

	//const float xmin = 50, xpoint1 = 150, xpoint2 = 250, xpoint3 = 400, xpoint4 = 1200, width1 = 10, width2 = 10, width3 = 50, width4 = 200;
	//const float xmin = 50, xpoint1 = 150, xpoint2 = 250, xpoint3 = 400, xpoint4 = 1200, width1 = 20, width2 = 10, width3 = 20, width4 = 250;
	//const float xmin = 50, xpoint1 = 150, xpoint2 = 250, xpoint3 = 400, xpoint4 = 1200, width1 = 20, width2 = 10, width3 = 20, width4 = 100;
//	const float xmin = 50, xpoint1 = 160, xpoint2 = 200, xpoint3 = 500, xpoint4 = 1100, width1 = 10, width2 = 40, width3 = 150, width4 = 300;
	//TH1* ratio2 = MakeVariableBinHist (ratio2, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4);
	TH1* ratio = dynamic_cast<TH1*> (ratio2->Clone("final"));
	ratio->SetDirectory(0);

	return ratio;
}

void abcdMethod::AddInput(const string rootFileName,const string path_to_hist, const string legendText,
						const string searchbins, const bool debug)
{
	
	input_t input;
		
	input.file = new TFile (rootFileName.c_str());
	TFile f(rootFileName.c_str());
	if (f.IsZombie()) {
		cout << "File with name " << rootFileName << " not found!" << endl;
		assert(false);
	}

	input.legText = legendText;
	input.hist_ratio = GetRatioHist(input.file, path_to_hist, debug);
	input.gaus_func = 0;
	input.exp_func = 0;
	input.gaus_fitResPtr = 0;
	input.exp_fitResPtr = 0;

	input.hist_ratio->Print();

	//control and signal hists are selected base on the search bins specified.
	SetSearchBinInfo(input.njetrange, input.htrange, input.mhtrange , input.results, searchbins);
	//GetControlHist(input.file, input.hist_control_finebin, input.njetrange, input.htrange, input.mhtrange, path_to_hist, debug);  
	TFile *fsigcon = new TFile("Mean/qcd_all.root");
	if (fsigcon->IsZombie()) {
		cout << __FUNCTION__ << ": control/signal root file Mean.qcd_all.root is not found!! " << endl;
		assert(false);
	}
	GetControlHist(fsigcon, input.hist_control_finebin, input.njetrange, input.htrange, input.mhtrange, path_to_hist, debug);  
	//GetSignalHist(input.file, input.hist_signal_finebin, input.njetrange, input.htrange, input.mhtrange, path_to_hist, debug);  
	GetSignalHist(fsigcon, input.hist_signal_finebin, input.njetrange, input.htrange, input.mhtrange, path_to_hist, debug);  
	if (b_doSystematics) InitSysts(input.Systs, input.njetrange, input.htrange, input.mhtrange, path_to_hist, debug);  

	vg_Inputs.push_back(input);
}

void abcdMethod::Run()
{

	for (unsigned i=0; i < vg_Inputs.size(); ++i)
	{
		//one fit per jet bin.
		stringstream str_njetrange, legText;
		str_njetrange <<  vg_Inputs.at(i).njetrange.at(0).first << "-" << vg_Inputs.at(i).njetrange.at(0).second << " Jets "; 
		legText <<  str_njetrange.str() <<  ": Default conditions for 'C'";
		//GetFits(vg_Inputs.at(i).hist_ratio, vg_Inputs.at(i).gaus_func,  vg_Inputs.at(i).exp_func, vg_Inputs.at(i).legText);
		GetFits(vg_Inputs.at(i).hist_ratio, 
					vg_Inputs.at(i).gaus_func,  
					vg_Inputs.at(i).gaus_fitResPtr,  
					vg_Inputs.at(i).exp_func, 
					vg_Inputs.at(i).exp_fitResPtr, 
					legText.str());

		if (b_doSystematics)
		{
			for (int s=0; s <Nsyst; ++s)
			{
				stringstream systtype;
				systtype << str_njetrange.str() << ": " << SystName(s) << " (Syst " << (s+1) << ")";
				GetFits(vg_Inputs.at(i).Systs.at(s).hist_ratio, 
						vg_Inputs.at(i).Systs.at(s).gaus_func, 
						vg_Inputs.at(i).Systs.at(s).gaus_fitResPtr, 
						vg_Inputs.at(i).Systs.at(s).exp_func,
						vg_Inputs.at(i).Systs.at(s).exp_fitResPtr,
						systtype.str() );
			}
		}



		//GET RESULTS for each of the search bins
		for (unsigned jetbin = 0; jetbin < vg_Inputs.at(i).njetrange.size(); ++jetbin)
		{
			for (unsigned htbin = 0; htbin < vg_Inputs.at(i).htrange.size(); ++htbin)
			{
				for (unsigned mhtbin = 0; mhtbin < vg_Inputs.at(i).mhtrange.size(); ++mhtbin)
				{

					const float mht_min = vg_Inputs.at(i).mhtrange.at(mhtbin).first; 
					const float mht_max = vg_Inputs.at(i).mhtrange.at(mhtbin).second; 

					//gaus results
					GetPredictions(
							mht_min, mht_max, 
							vg_Inputs.at(i).gaus_func,	
							vg_Inputs.at(i).gaus_fitResPtr,	
							vg_Inputs.at(i).hist_signal_finebin.at(jetbin).at(htbin),
							vg_Inputs.at(i).hist_control_finebin.at(jetbin).at(htbin),
							vg_Inputs.at(i).results.at(jetbin).at(htbin).at(mhtbin).first );
					//exponential results
					GetPredictions(
							mht_min, mht_max, 
							vg_Inputs.at(i).exp_func,	
							vg_Inputs.at(i).exp_fitResPtr,	
							vg_Inputs.at(i).hist_signal_finebin.at(jetbin).at(htbin),
							vg_Inputs.at(i).hist_control_finebin.at(jetbin).at(htbin),
							vg_Inputs.at(i).results.at(jetbin).at(htbin).at(mhtbin).second );

					const double cen_pred_gaus = vg_Inputs.at(i).results.at(jetbin).at(htbin).at(mhtbin).first.pred_mean;
					const double cen_pred_expo = vg_Inputs.at(i).results.at(jetbin).at(htbin).at(mhtbin).second.pred_mean;

					if (b_doSystematics)
					{

						//total systematic for this njet/ht/mht bin
						double c_systTot_gaus2 = 0, c_systTot_expo2 = 0;
						for (int s=0; s <Nsyst; ++s)
						{
							stringstream systtype;
							systtype << "Syst " << (s+1);

							GetPredictions(
									mht_min, mht_max, 
									//vg_Inputs.at(i).Systs.at(jetbin).at(htbin).at(mhtbin).at(s).gaus_func,	
									vg_Inputs.at(i).Systs.at(s).gaus_func,	
									vg_Inputs.at(i).Systs.at(s).gaus_fitResPtr,	
									vg_Inputs.at(i).hist_signal_finebin.at(jetbin).at(htbin),
									vg_Inputs.at(i).hist_control_finebin.at(jetbin).at(htbin),
									//vg_Inputs.at(i).Systs.at(jetbin).at(htbin).at(mhtbin).at(s).res_gaus
									vg_Inputs.at(i).Systs.at(s).res_gaus
									);

							GetPredictions(
									mht_min, mht_max, 
									//vg_Inputs.at(i).Systs.at(jetbin).at(htbin).at(mhtbin).at(s).exp_func,	
									vg_Inputs.at(i).Systs.at(s).exp_func,	
									vg_Inputs.at(i).Systs.at(s).exp_fitResPtr,	
									vg_Inputs.at(i).hist_signal_finebin.at(jetbin).at(htbin),
									vg_Inputs.at(i).hist_control_finebin.at(jetbin).at(htbin),
									//vg_Inputs.at(i).Systs.at(jetbin).at(htbin).at(mhtbin).at(s).res_expo
									vg_Inputs.at(i).Systs.at(s).res_expo
									);

							//const double alt_pred_gaus =  vg_Inputs.at(i).Systs.at(jetbin).at(htbin).at(mhtbin).at(s).res_gaus.pred_mean;
							//const double alt_pred_expo =  vg_Inputs.at(i).Systs.at(jetbin).at(htbin).at(mhtbin).at(s).res_expo.pred_mean;
							const double alt_pred_gaus =  vg_Inputs.at(i).Systs.at(s).res_gaus.pred_mean;
							const double alt_pred_expo =  vg_Inputs.at(i).Systs.at(s).res_expo.pred_mean;


							c_systTot_gaus2 += pow( fabs(cen_pred_gaus - alt_pred_gaus), 2); 
							c_systTot_expo2 += pow( fabs(cen_pred_expo - alt_pred_expo), 2); 
						}

						vg_Inputs.at(i).results.at(jetbin).at(htbin).at(mhtbin).first.pred_sysErr = sqrt(c_systTot_gaus2);
						vg_Inputs.at(i).results.at(jetbin).at(htbin).at(mhtbin).second.pred_sysErr = sqrt(c_systTot_expo2);

					} // doSystematics


					const double pred_statErr_gaus =  vg_Inputs.at(i).results.at(jetbin).at(htbin).at(mhtbin).first.pred_statErr;
					const double pred_fitErr_gaus  =  vg_Inputs.at(i).results.at(jetbin).at(htbin).at(mhtbin).first.pred_fitErr;
					const double pred_sysErr_gaus  =  vg_Inputs.at(i).results.at(jetbin).at(htbin).at(mhtbin).first.pred_sysErr;

					const double pred_statErr_expo =  vg_Inputs.at(i).results.at(jetbin).at(htbin).at(mhtbin).second.pred_statErr;
					const double pred_fitErr_expo  =  vg_Inputs.at(i).results.at(jetbin).at(htbin).at(mhtbin).second.pred_fitErr;
					const double pred_sysErr_expo  =  vg_Inputs.at(i).results.at(jetbin).at(htbin).at(mhtbin).second.pred_sysErr;

					const double tot_err_gaus  = sqrt( pow(pred_statErr_gaus,2) + pow(pred_fitErr_gaus, 2) + pow(pred_sysErr_gaus, 2) ); 
					const double tot_err_expo  = sqrt( pow(pred_statErr_expo,2) + pow(pred_fitErr_expo, 2) + pow(pred_sysErr_expo, 2) ); 

					vg_Inputs.at(i).results.at(jetbin).at(htbin).at(mhtbin).first.pred_totsysErr = tot_err_gaus;
					vg_Inputs.at(i).results.at(jetbin).at(htbin).at(mhtbin).second.pred_totsysErr = tot_err_expo;

					//do the sideband selection systematic




				} // mhtbin loop

			} //htbin loop
		} //jetbin loop
	} // input loop


	//print results

	for (unsigned i=0; i < vg_Inputs.size(); ++i)
	{
		PrintResults(vg_Inputs.at(i).njetrange, vg_Inputs.at(i).htrange, vg_Inputs.at(i).mhtrange, vg_Inputs.at(i).results);
	}



}

void abcdMethod::GetFits(TH1* ratio_hist, 
									TF1*& gausFit2, TFitResultPtr& gausFitResPtr,
									TF1*& expFit2, TFitResultPtr& expFitResPtr,
									const string optTitleText)
{

	if (ratio_hist == NULL)
	{
		cout << __FUNCTION__ << ": ratio his is null!" << endl;
		assert(false);
	}

	stringstream newtitle;
	if (optTitleText.length()>0) newtitle << optTitleText << ";MHT [GeV];Ratio (r);";
	newtitle << ";MHT [GeV];Ratio (r);";

	canvas_for_ratios->cd();
	canvas_for_ratios->Clear();

	gStyle->SetOptStat(0);	
	gPad->SetTickx();
	gPad->SetTicky();
	ratio_hist->SetLineColor(9);
	ratio_hist->SetTitle(newtitle.str().c_str());
	ratio_hist->SetLineWidth(2);

	gStyle->SetOptFit(1);
	//ratio_hist->SetStats(0);
	ratio_hist->Draw();
	//return;


	//use mean value of the last several bin as the 'C'
	const double CONST_C = GetAvgVal(ratio_hist, 400);
	cout << "FOR " << optTitleText << ":CONST_C = " << CONST_C  << endl;

	//do fittings exp and gaus
	//const float C_UPLIMIT = CONST_C+0.00001; //this only to get this values on the stat box
	const float C_LOLIMIT = CONST_C-0.00001;

	ratio_hist->SetMinimum(C_LOLIMIT/10);

	const float fitrange_xmin = 50.0;
	const float fitrange_xmax = 150.0;

	TF1 *expFit=new TF1("exp_func",expFitFunc,fitrange_xmin,fitrange_xmax,3);
	expFit->SetParameter(0,0.09);
	expFit->SetParameter(1,-0.0002);
	//expFit->SetParameter(2,CONST_C);
	//expFit->SetParLimits(2,C_LOLIMIT, C_UPLIMIT);
	expFit->FixParameter(2,CONST_C);
	expFitResPtr = ratio_hist->Fit(expFit,"E0S","",fitrange_xmin, fitrange_xmax);
	gPad->Modified();
	gPad->Update();
	TPaveStats *e_stats = (TPaveStats*) ratio_hist->FindObject("stats");
	e_stats->SetTextColor(kRed);
	TPaveStats *exp_stats = (TPaveStats*) e_stats->Clone("exp_stats");

	//TF1 *expFit2=new TF1("fit_2",expFitFunc,50,1000.0,3);
	//expFit2 = new TF1((*expFit));
	expFit2 = ratio_hist->GetFunction("exp_func");
	expFit2->SetLineColor(kRed+1);
	expFit2->SetLineWidth(2);
	expFit2->SetRange(50,1000);

	TF1 *gausFit=new TF1("gaus_func",gausFitFunc,fitrange_xmin, fitrange_xmax,3);
	gausFit->SetParameter(0,0.09); 
	gausFit->SetParameter(1,-0.0002);
	gausFit->FixParameter(2,CONST_C);
	//gausFit->SetParameter(2,CONST_C);
	//gausFit->SetParLimits(2,C_LOLIMIT,C_UPLIMIT);
	gausFit->SetLineColor(kGreen-2);
	gausFitResPtr = ratio_hist->Fit(gausFit,"E0S","", fitrange_xmin, fitrange_xmax);
	gPad->Modified();
	gPad->Update();
	TPaveStats *gaus_stats = (TPaveStats*) ratio_hist->FindObject("stats");
	gaus_stats->SetTextColor(kGreen);
	//TPaveStats *gaus_stats = (TPaveStats*) g_stats->Clone("gaus_stats");

	gausFit2 = new TF1(*gausFit);
	gausFit2->SetLineColor(kGreen);
	gausFit2->SetRange(50,1000);
	gausFit2->SetLineWidth(2);

	gausFit2->Draw("same");
	expFit2->Draw("same");
	//TLegend *leg  = new TLegend(0.7,0.8,0.9,0.9);
	TLegend *leg  = new TLegend(0.7,0.7,0.9,0.9);
	leg->AddEntry(gausFit2,"Gaussian");
	leg->AddEntry(expFit2,"Exponential");
	leg->SetTextFont(42);
	leg->Draw();

	const float xmin=0.2, xmax=0.45, ymin=0.7, ymax=0.9;
	gaus_stats->SetX1NDC(xmin);
	gaus_stats->SetX2NDC(xmax);
	gaus_stats->SetY1NDC(ymin);
	gaus_stats->SetY2NDC(ymax);
	gaus_stats->Draw("same");
	exp_stats->SetX1NDC(xmax);
	exp_stats->SetX2NDC(xmax+0.25);
	exp_stats->SetY1NDC(ymin);
	exp_stats->SetY2NDC(ymax);
	exp_stats->Draw("same");

	//fit quality
	stringstream gaus_fit_res, exp_fit_res;
	const float gausFit_goodness =  gausFit->GetChisquare()/gausFit->GetNDF();
	gaus_fit_res << "#chi^{2}/ndof = " << gausFit->GetChisquare()
			 << "/"<< gausFit->GetNDF() << " = " << gausFit_goodness;
	const float expFit_goodness =  expFit->GetChisquare()/expFit->GetNDF();
	exp_fit_res << "#chi^{2}/ndof = " << expFit->GetChisquare()
			 << "/"<< expFit->GetNDF() << " = " << expFit_goodness;

	//cout << gaus_fit_res.str() << endl;
	//cout << exp_fit_res.str() << endl;


//	gPad->SetEditable(false);

	canvas_for_ratios->Modified();
	canvas_for_ratios->Update();
	canvas_for_ratios->SetEditable(false);
	PrintRatioCanvas();
	canvas_for_ratios->Print(epsfile_for_ratios.c_str());
	canvas_for_ratios->SetEditable(true);

	delete expFit;
	delete gausFit;
	//delete exp_stats;
	//delete gaus_stats;
	//delete leg;

}


void abcdMethod::GetPredictions(const float& min_mht, const float& max_mht, 
									const TF1* func, const TFitResultPtr fitResPtr, 
									const TH1* signalHist, const TH1* controlHist, Predictions_t& pred)
{
	assert(controlHist != NULL && "GetPredictions:: controlHist not found!");
	assert(func != NULL && "GetPredictions:: function not found!");
	assert(signalHist != NULL && "GetPredictions:: signalHist not found!");
	//new TCanvas(); gPad->SetLogy(); controlHist->SetStats(1); controlHist->DrawCopy(); //gPad->Print("control.eps");
	
	if (controlHist->GetNbinsX() != signalHist->GetNbinsX())
	{
		cout<< __FUNCTION__<< ":: controlHist and signalHist have different binning of " << controlHist->GetNbinsX() << " / " << signalHist->GetNbinsX() << " !!!" << endl;
		assert(false);
	}


	//need this as this need to be passed for fit error calculation function which 
	//will change this a bit
	TF1 f1(*func);

	Predictions_t res;
	res.obs_mean = 0;
	res.obs_statErr = 0;
	res.pred_mean = 0;
	res.pred_statErr = 0;
	res.pred_fitErr  = 0;
	res.pred_sysErr = 0;
	res.pred_totsysErr = 0;


	for (int bin = 0; bin <= controlHist->GetNbinsX()+1; ++bin)
	{
		if (controlHist->GetBinContent(bin)>0)
		{
			const float binCenter = controlHist->GetBinCenter(bin);
			if (binCenter >= min_mht && binCenter < max_mht)
			{
			 /*cout << setprecision(4) << setw(5) << bin << setw(3) << "[" 
				<< setw(10) << controlHist->GetBinLowEdge(bin) << ", "  << setw(10) 
				<< controlHist->GetXaxis()->GetBinUpEdge(bin)<< "]" 
				<< setw(10) << controlHist->GetBinContent(bin) 
				<< setw(10) << f1.Eval(controlHist->GetBinCenter(bin))
				<< setw(10) << controlHist->GetBinContent(bin) * f1.Eval(controlHist->GetBinCenter(bin))
				<< setw(5) << " = " 
				<< setw(10) << signalHist->GetBinWidth(bin)
				<< setw(10) << signalHist->GetBinContent(bin)
				<< endl;
			*/
				const float binVal    = controlHist->GetBinContent(bin);
				const float binErr    = controlHist->GetBinError(bin);
				const float funcVal   = f1.Eval(binCenter);
				const float expt      = (binVal * funcVal);
				const float expt_statErr2  = pow(binErr * funcVal, 2);
				const float expt_fitErr    = binVal * GetFitFunctionError(&f1, fitResPtr, binCenter);
				const float binSig         = signalHist->GetBinContent(bin); 
				const float binSigStatErr2 = pow(signalHist->GetBinError(bin),2); 

				//cout << setprecision(3) << setw(15) << "MHT > " << incMHTBins.at(mhtBin) 
				//<< setw(20) << binSig  << "&$\\pm$" << binSigStatErr2 << endl;

				res.obs_mean += binSig;
				res.obs_statErr += binSigStatErr2;

				res.pred_mean += expt;
				res.pred_statErr += expt_statErr2;
				res.pred_fitErr  += expt_fitErr;
			}
		}
	}
	res.obs_statErr = sqrt(res.obs_statErr);
	res.pred_statErr = sqrt(res.pred_statErr);

	pred = res;
}


void abcdMethod::PrintRatioCanvas()
{
	cout << "Printing ratio canvas.. " << endl;
	canvas_for_ratios->cd();
	canvas_for_ratios->Print(epsfile_for_ratios.c_str());
}
void abcdMethod::PrintDebugCanvas()
{
	cout << "Printing debug canvas.. " << endl;
	canvas_for_debug->cd();
	canvas_for_debug->Print(epsfile_for_debug.c_str());
}


void abcdMethod::PrintResults(vector< pair<float, float> >& njetrange,
										vector< pair<float, float> >& htrange,
										vector< pair<float, float> >& mhtrange,
const vector< vector<vector< pair<abcdMethod::Predictions_t, abcdMethod::Predictions_t> > > >& results)
{
	cout <<  "\tNJET""\tHT "
					<< setw(15) << "MHT"
					<< setw(15) << " mean " 
					<< setw(20) << "statFitErr"
					<< setw(30) << "Signal+/-stat"
					<< endl;


	for (unsigned jetbin = 0; jetbin < njetrange.size(); ++jetbin)
	{
		const float jetmin = njetrange.at(jetbin).first;
		const float jetmax = njetrange.at(jetbin).second;

		double var_g = 0, var_e = 0;

		for (unsigned htbin = 0; htbin < htrange.size(); ++htbin)
		{

			const float htmin = htrange.at(htbin).first;
			const float htmax = htrange.at(htbin).second;
			for (unsigned mhtbin = 0; mhtbin < mhtrange.size(); ++mhtbin)
			{
				const float mhtmin = mhtrange.at(mhtbin).first;
				const float mhtmax = mhtrange.at(mhtbin).second;

				const float obs_mean             = results.at(jetbin).at(htbin).at(mhtbin).first.obs_mean;
				const float obs_statErr          = results.at(jetbin).at(htbin).at(mhtbin).first.obs_statErr;
				const float pred_mean_gaus       = results.at(jetbin).at(htbin).at(mhtbin).first.pred_mean;
				const float pred_statErr_gaus    = results.at(jetbin).at(htbin).at(mhtbin).first.pred_statErr;
				const float pred_fitErr_gaus     = results.at(jetbin).at(htbin).at(mhtbin).first.pred_fitErr;
				const float pred_sysErr_gaus     = results.at(jetbin).at(htbin).at(mhtbin).first.pred_sysErr;
				const float pred_totsysErr_gaus  = results.at(jetbin).at(htbin).at(mhtbin).first.pred_totsysErr;


				const float pred_mean_expo       = results.at(jetbin).at(htbin).at(mhtbin).second.pred_mean;
				const float pred_statErr_expo    = results.at(jetbin).at(htbin).at(mhtbin).second.pred_statErr;
				const float pred_fitErr_expo     = results.at(jetbin).at(htbin).at(mhtbin).second.pred_fitErr;
				const float pred_sysErr_expo     = results.at(jetbin).at(htbin).at(mhtbin).second.pred_sysErr;
				const float pred_totsysErr_expo  = results.at(jetbin).at(htbin).at(mhtbin).second.pred_totsysErr;
				
				var_g += pow(obs_mean-pred_mean_gaus,2);
				var_e += pow(obs_mean-pred_mean_expo,2);
	
				cout << setprecision(4) 
					<< "\t" << jetmin << "-" << jetmax 
					<< "\t" << htmin << "-" << htmax 
					<< "\t" << mhtmin << "-" << mhtmax
					<< setprecision(3) 
					<< "\t" << pred_mean_gaus   << "[" << pred_statErr_gaus << "/" << pred_fitErr_gaus << "/" << pred_sysErr_gaus << "=" << pred_totsysErr_gaus << "]"
					<< "\t" << pred_mean_expo   << "[" << pred_statErr_expo << "/" << pred_fitErr_expo << "/" << pred_sysErr_expo << "=" << pred_totsysErr_expo << "]"
					<< "\t" << obs_mean   << "[" << obs_statErr << "]"
					<< endl;
			}
		}
		double sig_g = var_g /double (htrange.size()+ mhtrange.size());
		double sig_e = var_e /double (htrange.size()+ mhtrange.size());
		cout << "Variation for " << jetmin << "-" << jetmax << ": gaus/exp = " << sig_g << "/" << sig_e << endl; 

	}
}

void abcdMethod::SetSearchBinInfo( 
			vector< pair<float, float> >& njetrange,
			vector< pair<float, float> >& htrange,
			vector< pair<float, float> >& mhtrange,
			vector< vector<vector< pair<abcdMethod::Predictions_t, abcdMethod::Predictions_t> > > >& results, 
			const string searchbins)
{
	string str(searchbins);

	//search for spaces and remove them
	//str.erase(remove_if(str.begin(), str.end(), isspace), str.end());
	if (str.find(" ") != string::npos)
	{
		cout << __FUNCTION__ << ": Input search bins has a space. please remove it!!!" << endl;
		assert(false);
	}

	//first break
	vector<string> vstrs = FindAllOf(str, ":");

	//expecting njet/ht/mht bin
	if (vstrs.size() != 3)
	{
		cout << __FUNCTION__ << ":Input parse error. please check search bin input format!" << endl;
		assert(false);
	}
	
	string str_njets = vstrs.at(0);
	string str_hts   = vstrs.at(1);
	string str_mhts  = vstrs.at(2);

	vector<pair<float, float> > nJetBins, htBins, mhtBins;
	
	//njet setup: this is assuming there will be only one input njet range per input
	vector<string> njstrs;
	njstrs.push_back(str_njets);
	nJetBins = GetVals(njstrs);

	vector<string> htstrs; 
	if (str_hts.find(",") != string::npos) //has more than one search bin
	{
		htstrs = FindAllOf(str_hts, ",");
	} else {
		htstrs.push_back(str_hts);
	}
	htBins = GetVals(htstrs);

	vector<string> mhtstrs;
	if (str_mhts.find(",") != string::npos) //has more than one search bin
	{
		mhtstrs = FindAllOf(str_mhts, ",");
	} else {
		mhtstrs.push_back(str_mhts);
	}
	mhtBins = GetVals(mhtstrs);
	
	//SearchBins sb;
	njetrange = nJetBins;
	htrange = htBins;
	mhtrange = mhtBins;
	//also use the file name here and setup the control/signal hists
	
	

	//now setup prediction bins

	for (unsigned jetbin = 0; jetbin < njetrange.size(); ++jetbin)
	{
		vector<vector<pair <Predictions_t, Predictions_t> > > vHtBinPred;
		for (unsigned htbin = 0; htbin < htrange.size(); ++htbin)
		{
			vector<pair <Predictions_t, Predictions_t> > vMhtBinPred;
			for (unsigned mhtbin = 0; mhtbin < mhtrange.size(); ++mhtbin)
			{
					Predictions_t gaus_pred, exp_pred;
					gaus_pred.obs_mean = 0;
					gaus_pred.obs_statErr = 0;
					gaus_pred.pred_mean = 0;
					gaus_pred.pred_statErr = 0;
					gaus_pred.pred_fitErr = 0;
					gaus_pred.pred_sysErr = 0;
					gaus_pred.pred_totsysErr = 0;

					exp_pred.obs_mean = 0;
					exp_pred.obs_statErr = 0;
					exp_pred.pred_mean = 0;
					exp_pred.pred_statErr = 0;
					exp_pred.pred_fitErr = 0;
					exp_pred.pred_sysErr = 0;
					exp_pred.pred_totsysErr = 0;
					vMhtBinPred.push_back(make_pair(gaus_pred, exp_pred));
			} //mht loop
			vHtBinPred.push_back(vMhtBinPred);
		} //ht loop
		results.push_back(vHtBinPred);
	} //jet loop
	

}

void abcdMethod::GetSignalHist(TFile *f, vector< vector<TH1*> > & signalHists, 
										vector< pair<float, float> >& njetrange,
										vector< pair<float, float> >& htrange,
										vector< pair<float, float> >& mhtrange,
										const string path_to_hist, const bool debug)
{

	//vector< vector<TH1*> > vHist;
	for (unsigned jetbin = 0; jetbin < njetrange.size(); ++jetbin)
	{
		const float jetmin = njetrange.at(jetbin).first;
		const float jetmax = njetrange.at(jetbin).second;
		vector<TH1*> vHistSignal;

		for (unsigned htbin = 0; htbin < htrange.size(); ++htbin)
		{

			const float htmin = htrange.at(htbin).first;
			const float htmax = htrange.at(htbin).second;

			stringstream commonpath, controlHist, signalHist;

			commonpath << "Hist/Njet" << jetmin << "to" << jetmax 
				<< "HT" << htmin << "to" << htmax
				<< "MHT0to8000";

			//const string temp_signal_hist_name ("smear_signalFineBin");
			const string temp_signal_hist_name ("smeared_signalFineBin");
			cout << __FUNCTION__ << ":" << __LINE__ <<  "::WARN!!! Using TEMP signal name "<< temp_signal_hist_name << endl;
			signalHist << commonpath.str() << "/" << temp_signal_hist_name;

			TH1* s_hist = dynamic_cast<TH1*> (f->Get(signalHist.str().c_str()));

			if (s_hist == NULL)
			{
				cout << __FUNCTION__ << ": Signal hist (" << signalHist.str() << ") not found in  file "
					<< f->GetName() << endl;
				assert (false);
			} else { 
				//temp fix to make euqal number of bins with the control hist
				//s_hist->Rebin(2);
				TH1* hist = dynamic_cast<TH1*> (s_hist->Clone());
				hist->SetDirectory(0); 
				vHistSignal.push_back(hist);
			}
		} //ht loop
		
		//vHist.push_back(vHistSignal);
		signalHists.push_back(vHistSignal);
	}

}
//will need to merge these two methods. for now separate because of naming issues
//in the smear_signal (should be smeared_signal

void abcdMethod::GetControlHist(TFile *f, vector< vector<TH1*> > & controlHists, 
										vector< pair<float, float> >& njetrange,
										vector< pair<float, float> >& htrange,
										vector< pair<float, float> >& mhtrange,
										const string path_to_hist, const bool debug)
{

	//vector< vector<TH1*> > vHist;
	for (unsigned jetbin = 0; jetbin < njetrange.size(); ++jetbin)
	{
		const float jetmin = njetrange.at(jetbin).first;
		const float jetmax = njetrange.at(jetbin).second;
		vector<TH1*> vHist;

		for (unsigned htbin = 0; htbin < htrange.size(); ++htbin)
		{

			const float htmin = htrange.at(htbin).first;
			const float htmax = htrange.at(htbin).second;

			stringstream commonpath, histname;

			commonpath << "Hist/Njet" << jetmin << "to" << jetmax 
				<< "HT" << htmin << "to" << htmax
				<< "MHT0to8000";

			histname << commonpath.str() << "/" << controlHistName;

			TH1* s_hist = dynamic_cast<TH1*> (f->Get(histname.str().c_str()));

			if (s_hist == NULL)
			{
				cout << __FUNCTION__ << ": Signal hist(" << histname.str() << ") not found in  file "
					<< f->GetName() << endl;
				assert (false);
			} else { 
				TH1* hist = dynamic_cast<TH1*> (s_hist->Clone());
				//cout << __FUNCTION__ << " : " << histname.str() << ": bin w = " << hist->GetBinWidth(10) << endl;  
				hist->SetDirectory(0); 
				vHist.push_back(hist);
			}
		} //ht loop
		
		controlHists.push_back(vHist);
	}

}

void abcdMethod::InitSysts(
			//vector< vector< vector< vector< abcdMethod::Syst_t > > > > & systs, 
			vector< abcdMethod::Syst_t > & systs, 
			vector< pair<float, float> >& njetrange,
			vector< pair<float, float> >& htrange,
			vector< pair<float, float> >& mhtrange,
			const string path_to_hist, const bool debug)
{

//	for (unsigned jetbin = 0; jetbin < njetrange.size(); ++jetbin)
//	{
//		vector<vector<vector< Syst_t > > > vSyst_ht;
//		for (unsigned htbin = 0; htbin < htrange.size(); ++htbin)
//		{
//			vector< vector<Syst_t> > vSyst_mht;
//			for (unsigned mhtbin = 0; mhtbin < mhtrange.size(); ++mhtbin)
//			{

				//const string systfilepath("Systematics/Syst");
				const string systfilepath("Syst");
				//const string systfilename("qcd_all.root");
				const string systfilename("qcd_all_HTranged.root");

				vector<Syst_t> vSyst;
				for (int i=0; i <Nsyst; ++i)
				{
					Syst_t syst;
					stringstream filename, legtext;
					const int j = i + 1;
					filename << systfilepath << j << "/" << systfilename;
					legtext << "Syst" << j;
					syst.file = new TFile(filename.str().c_str());
					if (syst.file->IsZombie())
					{
						cout << "File with name " << filename.str() << " not found!" << endl;
						assert(false);
					}

					syst.legText = legtext.str();
					syst.hist_ratio = GetRatioHist(syst.file, path_to_hist, debug);
					syst.gaus_func = 0;
					syst.exp_func  = 0;
					syst.gaus_fitResPtr = 0;
					syst.exp_fitResPtr  = 0;
					syst.gaus_sys  = 0;
					syst.exp_sys   = 0;
					
					const double large_neg_number = -99999.99; //so variable should not be used in the final calcualtions
					syst.res_gaus.obs_mean = large_neg_number;
					syst.res_gaus.obs_statErr = large_neg_number;
					syst.res_gaus.pred_mean = 0;  //only useful variable 
					syst.res_gaus.pred_statErr = large_neg_number;
					syst.res_gaus.pred_fitErr = large_neg_number;
					syst.res_gaus.pred_sysErr = large_neg_number;
					syst.res_gaus.pred_totsysErr = large_neg_number;

					syst.res_expo.obs_mean = large_neg_number;
					syst.res_expo.obs_statErr = large_neg_number;
					syst.res_expo.pred_mean = 0;  //only useful variable 
					syst.res_expo.pred_statErr = large_neg_number;
					syst.res_expo.pred_fitErr = large_neg_number;
					syst.res_expo.pred_sysErr = large_neg_number;
					syst.res_expo.pred_totsysErr = large_neg_number;

					systs.push_back(syst);
//					vSyst.push_back(syst);
				}
//				vSyst_mht.push_back(vSyst);
//			}
//			vSyst_ht.push_back(vSyst_mht);
//		}
//		systs.push_back(vSyst_ht);
//	} //jet loop
	

}

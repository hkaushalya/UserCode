/*******************************************************************
*  This has a bunch of utilities that is used under Stntuple
*  framework and when I make the final hitsts.
*  Author: Sam Hewamanage <samantha@fnal.gov>
********************************************************************
*	$Id: CommonTools.cc,v 1.3 2010/03/18 15:34:45 samantha Exp $
*	$Log: CommonTools.cc,v $
*	Revision 1.3  2010/03/18 15:34:45  samantha
*	ADDED: 1. Two overloaded GetMaxBinContent() method to find the contents of the
*	maximum bin in a hist or a set of hists.
*		2. class StringVector, which can hold a set of strings to which can be
*		printed at anytime during a job. Same info can be written to a canvas
*		and saved to the root file at the end of the job. I created this so  do
*		not have to search for the log file for a job. The root file will have
*		all the info I need about the job.
*	
*	Revision 1.2  2010/02/13 22:11:41  samantha
*	MODIFIED: FindNonZeroXrange() is modified to inly search the visible
*	region of the x-axis (i.e. exclude underflow and overflow bins.)
*	
*	Revision 1.1  2009/12/15 23:23:54  samantha
*	A collection of utitlity functions, color couts, Draw hists, rebin, etc.
*	
*
********************************************************************/

#include "CommonTools.hh"
#include <iostream>
#include "TStyle.h"
#include <iomanip>
#include <cmath>
#include <algorithm>
#include "TLatex.h"


//=================================================================//
//dump a bin info
//=================================================================//
void DumpBin(const TH1* hist, const int iBin,const std::string sMsg)
{
	assert (hist != NULL && "DumpBin::Hist is null!");
	assert (hist->GetDimension() == 1 && "CommonTools::DumpBin:: hist is not 1-D");
	assert ( (iBin>=1 && iBin<= hist->GetNbinsX()) && "DumpBin::Bin requested is out of range");
	
	if (sMsg.length()>0) std::cout << sMsg << std::endl;
	std::cout << "\t\t" << __FUNCTION__ << "::bin[loedge] : [val/err] = "
		<< std::setw(3) << iBin << "[" << hist->GetXaxis()->GetBinLowEdge(iBin) << "]" << std::setw(20)
		<< hist->GetBinContent(iBin) << " / " << hist->GetBinError(iBin) << std::endl;
}

//=================================================================//
void DrawClone(const TH1* h, const bool log)
{
	assert (h != NULL && "hist is null");
	assert (h->GetDimension() == 1 && "CommonTools::DrawClone:: hist is not 1-D");
	
	TH1* hist = dynamic_cast<TH1*> (h->Clone("copy"));
	new TCanvas();
	gStyle->SetOptStat("neou");
	gPad->SetGridx();
	gPad->SetGridy();
	if (log) gPad->SetLogy();
	if (hist->GetLineColor() == 5) hist->SetLineColor(46);
	hist->SetStats(1);
	hist->DrawClone();
}

//=================================================================//
// divide bins by bin width
//=================================================================//
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

//=================================================================//
//dumps hist bins including underflow and overflow
//=================================================================//
void DumpHist(const TH1* hist, const std::string sMsg)
{
	assert (hist != NULL && "requirement failed"); //spec
	assert (hist->GetDimension() == 1 && "CommonTools::DumpHist:: hist is not 1-D");

	if (sMsg.length()>0) std::cout << sMsg << std::endl;
	std::cout << std::setw(4) << "Bin" << std::setw(11) << "Low Edge" 
		<< std::setw(15) << "Bin Val"
		<< std::setw(15) << "Bin Err"
		<< std::endl;

	for (int bin = 0; bin <= hist->GetNbinsX()+1; ++bin)
	{
		if (hist->GetBinContent(bin)>0 || hist->GetBinError(bin)>0)
		{
			std::cout << "[" << std::setw(3) << bin << "]" 
				<< std::setw(9) << hist->GetXaxis()->GetBinLowEdge(bin)
				<< std::setw(15) << hist->GetBinContent(bin) 
				<< std::setw(15) << hist->GetBinError(bin) << std::endl;
		}
	}

}
//=================================================================//
//finds the bin number corresponding to a value
//=================================================================//
int GetBin(const TH1* hist, const double val)
{
	assert (hist != NULL && "GetBin::Hist is null!");
	for (int bin = 1; bin <= hist->GetNbinsX(); ++bin)
	{
		if (val > hist->GetXaxis()->GetBinLowEdge(bin) 
			  && val <= hist->GetXaxis()->GetBinUpEdge(bin)) return bin;
	}
	std::cout << __FUNCTION__ << "::WARNING! Did not find a bin for Value =" << val
		<< ". Hist bounds are [" << hist->GetXaxis()->GetBinLowEdge(1)
		<< ", " << hist->GetXaxis()->GetBinUpEdge(hist->GetNbinsX()) 
		<< "] ! returning -1" << std::endl;
	return -1;
}


//=================================================================//
// find the nozero x range two two 1-D hists
//=================================================================//
void FindCommonNonZeroXrange(const TH1* h1, const TH1* h2, float& xmin, float& xmax)
{
	assert (h1 != NULL && "requirement failed"); //spec
	assert (h2 != NULL && "requirement failed"); //spec
	assert (h1->GetDimension() == 1 && "CommonTools::FindNonZeroXrange:: hist1 is not 1-D");
	assert (h2->GetDimension() == 1 && "CommonTools::FindNonZeroXrange:: hist2 is not 1-D");
	for (int bin=0; bin<=h1->GetNbinsX(); ++bin)
	{
		if (h1->GetBinContent(bin)>0 || h2->GetBinContent(bin)>0)
		{
			xmin = h1->GetXaxis()->GetBinLowEdge(bin);
			break;
		}
	}
	for (int bin=h1->GetNbinsX(); bin>=1; --bin)
	{
		if (h1->GetBinContent(bin)>0 || h2->GetBinContent(bin)>0)
		{
			xmax = h1->GetXaxis()->GetBinUpEdge(bin);
			break;
		}
	}
}
//=================================================================//
// find the nozero x range for a 1-D hist in visible region
//=================================================================//
void FindNonZeroXrange(const TH1* h1, int& xminbin, int& xmaxbin)
{
	assert (h1 != NULL && "CommonTools::FindNonZeroXrange:: hist1 is NULL!");
	assert (h1->GetDimension() == 1 && "CommonTools::FindNonZeroXrange:: hist1 is not 1-D!");
	for (int bin=1; bin<=h1->GetNbinsX(); ++bin)
	{
		if (h1->GetBinContent(bin)>0)
		{
			xminbin = bin;
			break;
		}
	}
	for (int bin=h1->GetNbinsX()+1; bin>=0; --bin)
	{
		if (h1->GetBinContent(bin)>0)
		{
			xmaxbin = bin;
			break;
		}
	}
}

//-----------------------------------------------------------------------------
// effects: make a new binning with variable bin sizes, and the given
//   cutoff points between bin sizes
// returns: the binning (by bin edges)
// guarantee: strong
//// throws: std::bad_alloc
//-----------------------------------------------------------------------------
std::vector<float> 
GetVarBinVector (float down, float up5, float up10, float up20, float up50, 
					float width1, float width2, float width3, float width4)
{
	std::vector<float> result;
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
//-----------------------------------------------------------------------------
// make a variable binned hist with given specification (overloaded method)
//-----------------------------------------------------------------------------
TH1F* MakeVariableBinHist (const std::string name, const std::string title,
							const float xmin, const float xpoint1,
							const float xpoint2, const float xpoint3, const float xpoint4,
				 			const float width1, const float width2, const float width3, 
							const float width4)
{

	std::vector<float> vBins = GetVarBinVector(xmin, xpoint1, xpoint2, 
															xpoint3, xpoint4, width1, 
															width2, width3, width4);
	return new TH1F(name.c_str(), title.c_str(),vBins.size()-1, &vBins[0]);
}

//-----------------------------------------------------------------------------
// effects: turn this histogram into a histogram with variable bin
//   size
// returns: the new histogram
// guarantee: strong
// throws: std::bad_alloc
// requires: input != NULL
// requires: input->GetDimension() == 1
// requires: histogram bin edges must match
// postcondition: result != NULL
//-----------------------------------------------------------------------------
std::auto_ptr<TH1> FillVarBinHist (const std::vector<float>& bins, TH1 *input)
{
	assert (input != NULL && "requirement failed"); //spec
	assert (input->GetDimension() == 1 && "requirement failed"); //spec
	const unsigned nbin = unsigned (input->GetNbinsX ());

	std::auto_ptr<TH1> result (new TH1F (input->GetName(), input->GetTitle(),
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

//-----------------------------------------------------------------------------
// make variable binned hist from a given hist. this method is overloaded
//-----------------------------------------------------------------------------
TH1* MakeVariableBinHist (TH1 *hist, const float xmin, const float xpoint1,
							const float xpoint2, const float xpoint3, const float xpoint4,
				 			const float width1, const float width2, const float width3, 
							const float width4)
{
	assert (hist != NULL && "CommonTools::MakeVariableBinHist:: hist is NULL!");
	assert (hist->GetDimension() == 1 && "CommonTools::MakeVariableBinHist:: hist is not 1-D");
	
  	std::auto_ptr<TH1> result = FillVarBinHist ( 
											GetVarBinVector(xmin, xpoint1, xpoint2, 
															xpoint3, xpoint4, width1, 
															width2, width3, width4)
											, hist);
  	NormalizeBinWidth(result.get());
  	return result.release ();
};


void WilsonInterval(float k, float n, float &center, float &error)
{
	// See https://hep.baylor.edu/elog/mfrank/125 or search
	// "Wilson Score Interval" for a derivation of the below.
	// We pick z = 1 for a 68.3% confidence interval

	float p;

	if (n > 0)
	{
		p = k/n;
		center = ( p + (1/(2*n)) ) / (1 + (1/n));
		error = sqrt( (1/n) * (p*(1-p) + (1/(4*n))) ) / (1+(1/n));
	} else {
		center = 0;
		error = 0;
	}
}



void WilsonInterval(TH1* &nom, TH1* denom)
{
	for (int bin = 1; bin <= nom->GetNbinsX(); ++bin)
	{
		float k = nom->GetBinContent(bin);
		float n = denom->GetBinContent(bin);
		float center, error;

		WilsonInterval(k, n, center, error);

		nom->SetBinContent(bin, center);
		nom->SetBinError(bin, error);
	}
}

Clock::Clock()
{
	start = -1;
	end = -1;
};
Clock::Clock(const std::string n)
{
	name = n;
	Clock();
};
void Clock::Start() 
{
	std::cout << "in start " << std::endl;
	start = time_t(NULL);	//CINT seem to return 0 for this! 
	std::cout << "start time = " << start << std::endl;
}
void Clock::Stop() 
{ 
	std::cout << "in stop " << std::endl;
	end = time_t(NULL);  
	std::cout << "end  time = " << end << std::endl;
}

void Clock::ShowEndTimeHours() const 
{ 
	if (end > 0)
		std::cout << name << ": Total Time " <<  (end - start)/3600. << " (hrs)" << std::endl; 
	else 
		std::cout << "Clock has not stopped!" << start << ", " << end << std::endl;

	std::cout << name << ": Total Time " <<  (end - start)/3600. << " (hrs)" << std::endl; 
}

//find the value of the bin with the largest number of entries
//in a given hist
double GetMaxBinContent(const TH1* hist)
{
	assert (hist != NULL && "CommonTools::GetMaxBinContent()::Given hist is null!");
	assert (hist->GetDimension() == 1 && "CommonTools::GetMaxBinContent()::Hist must be 1-D!");
	return hist->GetBinContent(hist->GetMaximumBin());
}
//find the value of the bin with the largest number of entries
//in a set of given hists
double GetMaxBinContent(const std::vector<TH1*> hist)
{
	double max = 0;
	for (unsigned int i=0; i< hist.size(); ++i)
	{
		max = std::max(max, GetMaxBinContent(hist.at(i)));
	}
	return max;
}

StringVector::StringVector()
{
	c = 0;
};
StringVector::~StringVector()
{
	if (c != 0 || c != NULL) delete c;
};

TCanvas* StringVector::writeToCanvas()
{
	if (strVec.size())
	{
		if (strVec.size()>50) std::cout << __FUNCTION__ 
			<< ": Warning! Too many lines to print. Some may not print correctly. Think of creating a new StringVector..." << std::endl;
		if (c == 0 || c == NULL) c = new TCanvas();
		//const float h = c->GetWindowHeight();
	//	const float w = c->GetWindowWidth();

		float lineheight = 0.02;
		float ypos = 0.95;
		//float ypos = 0.1;
		TLatex l;
		l.SetTextSize(0.02);
		l.SetTextFont(80);
		l.SetTextAlign(12);
		
		for (std::vector<std::string>::const_iterator it = strVec.begin(); it != strVec.end(); ++it)
		{
			l.DrawLatex(0.03,ypos, (*it).c_str());
			ypos -= lineheight;
		}
		return dynamic_cast<TCanvas*> (c->Clone());
	} else 
	{
		std::cout << __FUNCTION__ << ":Vector is empty! Nothing to write.! returning null pointer!" << std::endl;
		return c;
	}

};


void SetSmallNumberErrors(TH1* hist, const int confInteval)
{
	//confInteval =0: 68.3% 
	//confInteval =1: 95% 
	
	assert (hist != NULL && "SetSmallNumberErrors:: hist is null");
	assert (hist->Dimensions() ==1 && "SetSmallNumberErrors:: hist is not 1-D !!!");

	vector<float> errLo, errHi;

	if (confInteval == 0) //68%
	{
		errLo.push_back(0.0);   errHi.push_back(1.15); //n=0
		errLo.push_back(0.268); errHi.push_back(2.5); //n=1
		errLo.push_back(0.864); errHi.push_back(3.86); //n=2
		errLo.push_back(1.55); errHi.push_back(5.15); //n=3
		errLo.push_back(2.29); errHi.push_back(6.4); //n=4
		errLo.push_back(3.06); errHi.push_back(7.63); //n=5
		errLo.push_back(3.85); errHi.push_back(8.84); //n=6
		errLo.push_back(4.65); errHi.push_back(10.0); //n=7
		errLo.push_back(5.47); errHi.push_back(11.2); //n=8
		errLo.push_back(6.30); errHi.push_back(12.4); //n=9
		errLo.push_back(7.14); errHi.push_back(13.5); //n=10
		errLo.push_back(8.84); errHi.push_back(15.8); //n=12
		errLo.push_back(10.6); errHi.push_back(18.1); //n=14
		errLo.push_back(12.3); errHi.push_back(20.4); //n=16
		errLo.push_back(14.1); errHi.push_back(22.6); //n=18
		errLo.push_back(15.8); errHi.push_back(24.8); //n=20

	} else if (confInteval == 1) //98%
	{
		errLo.push_back(0.0);    errHi.push_back(3.0); //n=0
		errLo.push_back(0.042);  errHi.push_back(4.77); //n=1
		errLo.push_back(0.304);  errHi.push_back(6.4); //n=2
		errLo.push_back(0.713);  errHi.push_back(7.95); //n=3
		errLo.push_back(1.21);   errHi.push_back(9.43); //n=4
		errLo.push_back(1.76);   errHi.push_back(10.9); //n=5
		errLo.push_back(2.35);   errHi.push_back(12.3); //n=6
		errLo.push_back(2.97);   errHi.push_back(13.6); //n=7
		errLo.push_back(3.62);   errHi.push_back(15.0); //n=8
		errLo.push_back(4.29);   errHi.push_back(16.3); //n=9
		errLo.push_back(4.98);   errHi.push_back(17.6); //n=10
		errLo.push_back(6.40);   errHi.push_back(20.2); //n=12
		errLo.push_back(7.86);   errHi.push_back(22.7); //n=14
		errLo.push_back(9.36);   errHi.push_back(25.2); //n=16
		errLo.push_back(10.9);   errHi.push_back(27.7); //n=18
		errLo.push_back(12.4);   errHi.push_back(30.1); //n=20

	}

	//set errors including underflow and overflow bins
	for (int bin = 0; bin < hist->GetNbinsX()+1; ++bin)
		
	

}

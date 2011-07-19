#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "TH1.h"
#include "TFile.h"
#include <assert.h>
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPad.h"

using namespace std;

// this is a to overlay a hist from different
// datasets found in different root files.
// Must compiled and loaded.
class HistMaker
{
	public:
		HistMaker();
		~HistMaker();
		void AddHist(const string completeName, const string title="", const bool logScale=false,
						const double yMaxMultiplier=1);
		void AddFile(const string fileName, const string legendText,
				const int marker = 0, const int lineWidth=2, const int lineColor=2); 
		void LegendOptions(const float x1=0.6, const float x2=0.9, const float y1=0.7, const float y2=0.9, const float font=42);
		void DrawAll();
		void Draw(const int);
		//void DrawAll(const int);
		void PrintFormat(const string format);
		struct histProperties_t
		{
			string fileName;
			int lineWidth;
			int lineColor;
			string legendText;
			int marker;
			TH1* hist; // pointer to hist
		};
		struct histNames_t
		{
			string name;
			string title;
			bool logScale;
			double yMaxMultiplier;
		};

		struct legendOptions_t
		{
			float x1;
			float x2;
			float y1;
			float y2;
			float font;
		};

		struct statOptions_t
		{
			bool optStat;
		};
	private:
		vector<histProperties_t> vHistPropeties;
		vector<histNames_t> vHistNames;
		TH1* GetHist(const string, const string);
		legendOptions_t legendOptions;
		string printFormat;  //eps, gif, pdf

};

HistMaker::HistMaker()
{
	legendOptions.x1 = 0.6;
	legendOptions.x2 = 0.9;
	legendOptions.y1 = 0.7;
	legendOptions.y2 = 0.9;
	legendOptions.font = 42;
	
	printFormat="";
}

HistMaker::~HistMaker()
{
}

void HistMaker::PrintFormat(const string format)
{
	if (format == "eps" || format == "pdf" || format == "png"
			|| format == "gif" || format == "jpg")
	{ 
		printFormat = format;
	} else cout << __FUNCTION__ <<":: Invalid print format " << format << " specified!";
}


void HistMaker::AddHist(const string histName, const string title, const bool logScale, const double yMaxMultiplier)
{ 
	histNames_t hName;
	hName.name = histName;
	hName.title = title;
	hName.logScale = logScale;
	hName.yMaxMultiplier = yMaxMultiplier;
	vHistNames.push_back(hName);
}

void HistMaker::AddFile(const string fileName, const string legendText,
				 const int marker, const int lineWidth, const int lineColor)
{
	histProperties_t hist;
	hist.fileName = fileName;
	hist.lineWidth = lineWidth;
	hist.lineColor = lineColor;
	hist.legendText = legendText;
	hist.marker = marker;
	hist.hist = 0;
	vHistPropeties.push_back(hist);
}


TH1* HistMaker::GetHist(const string fileName, const string completePath)
{
	TFile file(fileName.c_str());
	if (file.IsZombie())
	{
		cout << __FUNCTION__ << ":: file " << fileName 
			<< " not found ! returning 0!" << std::endl;
		return 0;
	} else
	{
		TH1* hist = dynamic_cast<TH1*> (file.Get(completePath.c_str()));
		if (hist != NULL)
		{
			TH1* hist_clone = dynamic_cast<TH1*> (hist->Clone());
			hist_clone->SetDirectory(0);
			return hist_clone;
		} else
		{
			cout << __FUNCTION__ << ":: hist " << completePath 
				<< " not found  in file " << fileName 
				<< "! returning 0!" << std::endl;
			return 0;
		}
	}

}

void HistMaker::LegendOptions(const float x1, const float x2, const float y1, const float y2, const float font)
{
	legendOptions.x1 = x1;
	legendOptions.x2 = x2;
	legendOptions.y1 = y1;
	legendOptions.y2 = y2;
	legendOptions.font = font;

}

void HistMaker::Draw(const int iHist)
{
	if (iHist>=0 && iHist< vHistNames.size())
	{
		double yMax = 0;
		double yMin = 0;
		if (vHistNames.at(iHist).logScale) yMin = 0.005;
		TLegend *leg  = new TLegend(legendOptions.x1, legendOptions.y1,legendOptions.x2,legendOptions.y2);
		for (int i=0; i< vHistPropeties.size(); ++i)
		{
			vHistPropeties.at(i).hist = GetHist(vHistPropeties.at(i).fileName, vHistNames.at(iHist).name);
			if (vHistPropeties.at(i).hist == 0)
			{
				cout << __FUNCTION__ << "::Hist not found! exiting!" << endl;
				return;
			}
			const double ymax = vHistPropeties.at(i).hist->GetBinContent(vHistPropeties.at(i).hist->GetMaximumBin());
			if (ymax> yMax) yMax = ymax;
			vHistPropeties.at(i).hist->Print();
			vHistPropeties.at(i).hist->SetLineColor(vHistPropeties.at(i).lineColor);
			vHistPropeties.at(i).hist->SetLineWidth(vHistPropeties.at(i).lineWidth);
			if (vHistNames.at(iHist).title.length()>0)
			{
				vHistPropeties.at(i).hist->SetTitle(vHistNames.at(iHist).title.c_str());
			}

			leg->AddEntry(vHistPropeties.at(i).hist, vHistPropeties.at(i).legendText.c_str());
		}

		new TCanvas();
		if (vHistNames.at(iHist).logScale) gPad->SetLogy();
		//need another loop to set maximum y
		for (int i=0; i< vHistPropeties.size(); ++i)
		{
			vHistPropeties.at(i).hist->SetMaximum(yMax * vHistNames.at(iHist).yMaxMultiplier);
			vHistPropeties.at(i).hist->SetMinimum(yMin);
			vHistPropeties.at(i).hist->SetStats(0);
			if (i==0) vHistPropeties.at(i).hist->Draw();
			else vHistPropeties.at(i).hist->Draw("same");
		}
		leg->Draw();
		//print only if a format is specified
		if (printFormat.length()>0)
		{
			stringstream printFile;
			printFile << vHistPropeties.at(0).hist->GetName() << "." << printFormat; 
			gPad->Print(printFile.str().c_str());
		}

	}

}
void HistMaker::DrawAll()
{
		for (int ihist=0; ihist < vHistNames.size(); ++ihist) Draw(ihist);
}

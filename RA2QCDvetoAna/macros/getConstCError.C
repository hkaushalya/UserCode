#include <algorithm>
#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include  "TCanvas.h"

using namespace std;

void getConstCError() {
	TFile f("ReponseSmearingSyst.root");

	TH1* h[5];
 
	h[0] = (TH1*) (f.Get("passfail_QCDPythia_SmearedDefault"));
	h[1] = (TH1*) (f.Get("passfail_QCDPythia_lowTailMin"));
	h[2] = (TH1*) (f.Get("passfail_QCDPythia_lowTailMax"));
	h[3] = (TH1*) (f.Get("passfail_QCDPythia_upperTailMin"));
	h[4] = (TH1*) (f.Get("passfail_QCDPythia_upperTailMax"));

	float lastBinCont[5];
	new TCanvas();
	for (int i=0; i<5; ++i) 
	{
		if ( i == 0) h[i]->DrawCopy();
		else h[i]->DrawCopy("same");
		lastBinCont[i] = h[i]->GetBinContent(h[i]->GetNbinsX());
		cout << h[i]->GetName() << " = " << lastBinCont[i] << endl;
	}

	float d1 = fabs(lastBinCont[0] - lastBinCont[1]);
	float d2 = fabs(lastBinCont[0] - lastBinCont[2]);
	float d3 = fabs(lastBinCont[0] - lastBinCont[3]);
	float d4 = fabs(lastBinCont[0] - lastBinCont[4]);

	float max1 = std::max(d1,d2);
	float max2 = std::max(d3,d4);

	float max = std::max(max1, max2);

	cout << "max = " << max << endl;


}

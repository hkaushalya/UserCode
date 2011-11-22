#include <TChain.h>
#include "TChainElement.h"
#include "TCollection.h"// for TIter
#include <sstream>
#include <iostream>
#include "TSystem.h"
#include <cmath>
#include <exception>
#include <assert.h>
#include "TFile.h"
#include "Factorization.hh"
#include<string>
#include "TObject.h"
#include "TObjArray.h"

void DoRunMyJob (const int events, const int dataset=0)
{
	TChain *ch = new TChain("AUX");
 	const std::string fileprefix="PhoJets";
  	std::stringstream rootfile;
	
	int MCFLAG = 0;
	if (dataset != 1) MCFLAG = 1;

	if (dataset == 1)
	{
		std::cout << __FILE__ << ":: Processing data" << std::endl;
		ch->Add("/uscms_data/d3/samantha/SUSYRA2work/CMSSW_4_2_5/src/UserCode/RA2QCDvetoAna/test/FactorizationResults/11172011_DataTree/1162011_Data1/StdRA2FlatTree.root");
	} else if (dataset == 2)
	{
		std::cout << __FILE__ << ":: Processing data" << std::endl;
		ch->Add("/uscms_data/d3/samantha/SUSYRA2work/CMSSW_4_2_5/src/UserCode/RA2QCDvetoAna/test/FactorizationResults/FlatTrees/StdRA2FlatTree_LM9.root");
	} else 
	{
		std::cout << __FILE__ << ":: NO dataset found!" << std::endl;
		return;
	}

	TObjArray *fileElements=ch->GetListOfFiles();
	TIter next(fileElements);
	TChainElement *chEl=0;
	while (( chEl=(TChainElement*)next() )) 
	{
		TFile f(chEl->GetTitle());
		assert ( (! f.IsZombie()) && " a file is not found. pl check!");
		std::cout << "Attached data file: ";  f.Print();
	}
	std::cout << __FILE__ << ": Total Entries in chain = " << ch->GetEntries() << std::endl;

	//std::cout << __FILE__ << " ch = " << ch << std::endl;
	TObjArray *brs = ch->GetListOfBranches();
	TIter n(brs);
	TChainElement *b=0;
	while (( b=(TChainElement*)n()))
	{
		//b->Print();
	}
	
//	ch->Draw("vtx_NClass12");

	
/*	const float fMinJetEt = 15.;
	const float fMaxJetEt = 1200.;
	const float fMinJetEta = -3.0;
	const float fMaxJetEta = 3.0;
	const int iMinNjet15 = 1;
	const int iMaxNjet15 = 100;
	const float fMinMet = 0;
	const float fMaxMet = 1500;
	int iUseNvtxWgts = 0;
	if (dataset != 1 && dataset != 50) iUseNvtxWgts = 1;
	std::cout << "iUseNvtxWgts = " << iUseNvtxWgts << std::endl;
*/
	Factorization *fact = new Factorization;
	fact->SetHistFileName("Data.root");
	fact->Run(ch,events);

}

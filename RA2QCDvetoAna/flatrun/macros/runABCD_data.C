#include "abcdMethod_data.hh"
#include <iostream>
#include "TSystem.h"


using namespace std;


int main()
{
	abcdMethod_data *ab; 

return 1;
};

void runABCD_data()
{

	//format of input searhc bins:
	//njet:ht:mht, separate mutiple search bins per njet or ht bin by commas
	//

	//gSystem->CompileMacro("abcdMethod_data.cc","=");
	const string BGfile("TotBG_MCpuWgted_Lumi_12187pb.root");
	abcdMethod_data *ab = new abcdMethod_data(BGfile);

	//const string projfile("01292013_data_nodphi_trigsel_notrigwgt/qcd_all_HTranged.root");
	const string projfile("01292013_data_nodphi_trigsel_trigwgt_noBGsub/qcd_all_HTranged.root");

	const float c_for_3to5jets = 0.0309481;
//	ab->AddInput(projfile.c_str(),"Hist/Njet3to5", "Central Val",
					//"3-5:500-800:400-600,600-5000:", c_for_3to5jets);
//					"3-5:500-800:300-350:", c_for_3to5jets);
	ab->AddInput(projfile.c_str(),"Hist/Njet3to5", "Central Val",
					"3-5:500-800,800-1000,1000-1250:200-300,300-450,450-600,600-5000:", c_for_3to5jets);
	ab->AddInput(projfile.c_str(),"Hist/Njet3to5", "Central Val",
					"3-5:1250-1500:200-300,300-450,450-5000:", c_for_3to5jets);
	ab->AddInput(projfile.c_str(),"Hist/Njet3to5", "Central Val",
					"3-5:1500-8000:200-300,300-5000:", c_for_3to5jets);

	const float c_for_6to7jets = 0.193871;
	ab->AddInput(projfile.c_str(),"Hist/Njet6to7", "Central Val",
					"6-7:500-800,800-1000,1000-1250,1250-1500:200-300,300-450,450-5000:", c_for_6to7jets);
	ab->AddInput(projfile.c_str(),"Hist/Njet6to7", "Central Val",
					"6-7:1500-8000:100-200,200-300,300-5000:", c_for_6to7jets);

	const float c_for_8toInfjets = 0.23578; 

	ab->AddInput(projfile.c_str(),"Hist/Njet8to1000", "Central Val",
					"8-1000:500-800,800-1000,1000-1250,1250-1500,1500-8000:200-5000:", c_for_8toInfjets);

	
	ab->Run();
	
}

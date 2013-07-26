#include "abcdMethod_MC.hh"
#include <iostream>
#include "TSystem.h"


using namespace std;


int main()
{
	abcdMethod_MC *ab; 

return 1;
};

void runABCD_MC()
{
	//gSystem->CompileMacro("abcdMethod_MC.cc","=");
	abcdMethod_MC *ab = new abcdMethod_MC();

	///ab->AddInput("Mean/qcd_all_HTranged.root","Hist/Njet3to5HT500to750MHT0to8000", "Central Val",
//	ab->AddInput("Mean/qcd_all_HTranged.root","Hist/Njet6to7", "Central Val",
//					"6-7:500-750:200-350,350-500:");

//	ab->AddInput("Mean/qcd_all_HTranged.root","Hist/Njet2to2", "Central Val",
//					"2-2:500-750:200-300,300-450,450-600,600-5000:");
/*	ab->AddInput("Mean/qcd_all_HTranged.root","Hist/Njet2to2", "Central Val",
					"2-2:750-1000,1000-1250:100-200,200-300,300-450,450-600,600-5000:");
	ab->AddInput("Mean/qcd_all_HTranged.root","Hist/Njet2to2", "Central Val",
					"2-2:1250-1500:100-200,200-300,300-450,450-5000:");
	ab->AddInput("Mean/qcd_all_HTranged.root","Hist/Njet2to2", "Central Val",
					"2-2:1500-8000:100-200,200-300,300-5000:");
*/


	ab->AddInput("Mean/qcd_all_HTranged.root","Hist/Njet3to5", "Central Val",
					"3-5:500-750,750-1000,1000-1250:200-300,300-450,450-600,600-5000:");
	ab->AddInput("Mean/qcd_all_HTranged.root","Hist/Njet3to5", "Central Val",
					"3-5:1250-1500:200-300,300-450,450-5000:");
	ab->AddInput("Mean/qcd_all_HTranged.root","Hist/Njet3to5", "Central Val",
					"3-5:1500-8000:200-300,300-5000:");


	ab->AddInput("Mean/qcd_all_HTranged.root","Hist/Njet6to7", "Central Val",
					"6-7:500-750,750-1000,1000-1250,1250-1500:200-300,300-450,450-5000:");
	ab->AddInput("Mean/qcd_all_HTranged.root","Hist/Njet6to7", "Central Val",
					"6-7:1500-8000:100-200,200-300,300-5000:");

	ab->AddInput("Mean/qcd_all_HTranged.root","Hist/Njet8to1000", "Central Val",
					"8-1000:500-750,750-1000,1000-1250,1250-1500,1500-8000:200-5000:");

	//format of input searhc bins:
	//njet:ht:mht, separate mutiple search bins per njet or ht bin by commas
	//

	
	ab->Run();
	
}

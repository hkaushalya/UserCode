#include <iostream>
#include <vector>
#include <string>

using namespace std;

int main(int argc, char* argv[])
{
	//cout << __FUNCTION__ << ": number of arguments = " << argc << endl;
	
	if (argc <= 3) {
		cerr <<"Please give 3 arguments " << "runList " << " " << "outputFileName" << " " << "# of evts  nJet50Min  nJet50Max" << endl;
		cerr <<" Valid configurations are " << std::endl;
		cerr << " ./optimize filelist.txt out.root nevts2process[-1 = all] smearingSyst[0=mean, 1-6 systs]" << std::endl;
		cerr << "Eg:  ./optimize filelist.txt out.root 100" << std::endl;
		return -1;
	}

	for (int i=0; )
	{

	}

	const char *inputFileList = argv[1];
	const char *outFileName   = argv[2];
	const char *evt2Proc      = argv[3];
		
	int evts    = -1;
	int systematic_var = 0;
	unsigned cutmask = 128+64+32+16+8+4+2+1;  //enable all cuts

	if (isdigit(evt2Proc[0]))
	{
		int num = atoi(argv[3]);
		if (num>=1) evts = num;
		//cout << __FUNCTION__ << ": evts = " << evts << endl;
	} else {
			cout << "argument 4 is not a number. using default value for evts = " << evts << endl;
	}
	
	if (argc>4)
	{
		//cout << __FUNCTION__ << ": processing arg4 ..." << endl;
		const char *g4 = argv[4];
		if (isdigit(g4[0])) systematic_var = atoi(g4);
		else {
			cout << "argument 4 is not a number. using default value for systematic_var = " << systematic_var << endl;
		}
	}

	if (argc>5)
	{
		const char *g5 = argv[5];
		if (isdigit(g5[0])) cutmask = atoi(g5);
		else {
			cout << "argument 5 is not a number. using default value for cutmask = " << cutmask << endl;
		}
	}


	return 0;
}

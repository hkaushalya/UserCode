void RunMyJob (int events = -1, int dataset=0)
{
	std::cout << "events, dataset = " << events << "," << dataset << std::endl;
	gROOT->Reset();
	gSystem->Load("$ROOTSYS/lib/libPhysics.so");
	std::vector<std::string> vLibs;
	vLibs.push_back("Factorization.cc");
	vLibs.push_back("DoRunMyJob.cc");

	bool bReady = true;
	for (unsigned i=0; i<vLibs.size(); ++i)
	{
		std::cout << "Compiling file " << vLibs.at(i);
		if (gSystem->CompileMacro(vLibs[i].c_str(),"k") == 1) 
		{
			std::cout << " success." << vLibs.at(i);
			continue;
		}
		bReady = false;
		break;
	}

	std::cout << __LINE__ << std::endl;
	if (! bReady) 
	{
		//gSystem->ListLibraries();
		std::cout << __FILE__ << ":: NOT ALL REQUIRED LIBRARIES ARE LOADED PROPERLY!" << std::endl;
		return;
	}
	
	std::cout << __FILE__ << ":" << __LINE__ << std::endl;
	DoRunMyJob (events, dataset);
	std::cout << __FILE__ << ":" << __LINE__ << std::endl;
	//gObjectTable->Print();
}

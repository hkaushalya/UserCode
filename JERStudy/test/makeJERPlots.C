


void makeJETPlots()
{

	TFile *f = new TFile ("JER_FullMCQCDSample.root");
	std::string path("jerana/");

	float eta_bin_step = 0.1;
	const float eta_max = 5.0;
	const float pt_max = 1500;
	
	const float pt_lessthan100_resolution = 0.2;
	const float pt_morethan100_resolution = 0.1;

	for (float f=0.0; f<= eta_max; f+=eta_bin_step)
	{
		pt_bin_step = 20;
		for (float pt=0.0; pt<= pt_max; pt+=pt_bin_step)
		{
			if (iVerbose)std::cout <<  pt << " | "; 
			if (pt>=960) pt_bin_step = 100.0;
			else if (pt>=800) pt_bin_step = 80.0;
			else if (pt>=620) pt_bin_step = 60.0;
			else if (pt>=420) pt_bin_step = 40.0;
			else if (pt>=300) pt_bin_step = 30.0;

			const float eta_min = f; 
			const float eta_max = f + eta_bin_step;
			const float pt_min = pt; 
			const float pt_max = pt+pt_bin_step;

			std::stringstream name;
			name << "eta_" << etpt.eta_max << "<|#Eta|<" << etpt.eta_max;
			pt_range << etpt.pt_min << "<pt<" << etpt.pt_max;
			etpt.str_eta_range = eta_range.str();
			etpt.str_pt_range = pt_range.str();

			vEtaPtBins.push_back(etpt);

			MatchJetHist_t hist;
			vJERHists.push_back(hist);
			BookHistograms(fs, vJERHists.end()-1, etpt);
		}
	}

}

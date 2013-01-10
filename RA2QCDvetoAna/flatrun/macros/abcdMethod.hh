#ifndef ABCDMETHOD_MC_HH
#define ABCDMETHOD_MC_HH

#include <string>
#include "TF1.h"
#include "TH1.h"
#include "TCanvas.h"
#include <vector>
#include <utility>
#include "TFile.h"

class abcdMethod
{
	public:
		abcdMethod();
		~abcdMethod();

	
		struct Predictions_t
		{
			double obs_mean;
			double obs_statErr;
			double pred_mean;
			double pred_statErr;
			double pred_fitErr;
			double pred_sysErr; //total error for 'C'
			double pred_totsysErr; //total error stat+fit+c-err
		};

		struct Syst_t {
			TFile* file;
			std::string legText;
			TH1* hist_ratio;
			TF1* gaus_func;
			TF1* exp_func;
			double gaus_sys;
			double exp_sys;
			Predictions_t res_gaus;
			Predictions_t res_expo;
		};

		struct input_t {
			TFile* file;
			std::string legText;
			TH1* hist_ratio;
			TF1* gaus_func;
			TF1* exp_func;
			std::vector< std::pair<float, float> > njetrange;
			std::vector< std::pair<float, float> > htrange;
			std::vector< std::pair<float, float> > mhtrange;
			std::vector< std::vector<TH1* > > hist_control_finebin; 
			std::vector< std::vector<TH1* > > hist_signal_finebin; 
			std::vector< std::vector<std::vector< std::pair<abcdMethod::Predictions_t, abcdMethod::Predictions_t> > > > results;  //gaus/exp 
			//njet/ht/mht/systs
			std::vector< std::vector< std::vector<std::vector< abcdMethod::Syst_t > > > > Systs; //holds 6 individual systematics. pred_mean is the only used variable
		};


	private:

		int Nsyst;
		bool b_doSystematics;
		std::string numHistName, denHistName, controlHistName, signalHistName;
		std::vector<TH1*> vRatioHist; //0=mean, >0 are all systematics
		TH1* GetRatioHist(TFile* f, const std::string path_to_hist, const bool debug=false);
		TH1* GetControlHist(TFile* f, const std::string path_to_hist, const bool debug=false);
		void GetSignalHist(TFile *f, std::vector< std::vector<TH1*> > & signalHists, 
										std::vector< std::pair<float, float> >& njetrange,
										std::vector< std::pair<float, float> >& htrange,
										std::vector< std::pair<float, float> >& mhtrange,
										const std::string path_to_hist, const bool debug);
		void GetControlHist(TFile *f, std::vector< std::vector<TH1*> > & controlHists, 
										std::vector< std::pair<float, float> >& njetrange,
										std::vector< std::pair<float, float> >& htrange,
										std::vector< std::pair<float, float> >& mhtrange,
										const std::string path_to_hist, const bool debug);
		void GetFits(TH1* ratio_hist, TF1*& gausFunc, TF1*& expFunc, const std::string optTitleText="");
		void SetSearchBinInfo( std::vector< std::pair<float, float> >& njetrange,
										std::vector< std::pair<float, float> >& htrange,
										std::vector< std::pair<float, float> >& mhtrange,
										std::vector< std::vector<std::vector< std::pair<abcdMethod::Predictions_t, abcdMethod::Predictions_t> > > >& results, 
										const std::string searchbins);
		void InitSysts(
				std::vector< std::vector< std::vector< std::vector< abcdMethod::Syst_t > > > > & systs, 
				std::vector< std::pair<float, float> >& njetrange,
				std::vector< std::pair<float, float> >& htrange, 
				std::vector< std::pair<float, float> >& mhtrange, 
				const std::string path_to_hist, const bool debug);
		void PrintRatioCanvas();
		void PrintDebugCanvas();
		std::vector<input_t> vg_Inputs;
		TCanvas *canvas_for_ratios, *canvas_for_debug;
		std::string epsfile_for_ratios, epsfile_for_debug;
		void PrintResults(std::vector< std::pair<float, float> >& njetrange,
										std::vector< std::pair<float, float> >& htrange,
										std::vector< std::pair<float, float> >& mhtrange,
					const std::vector< std::vector<std::vector< std::pair<abcdMethod::Predictions_t, abcdMethod::Predictions_t> > > >& results);
		void  GetPredictions(const float& min_mht, const float& max_mht, const TF1* f1, const TH1* signalHist, const TH1* controlHist, Predictions_t& pred);

	public:
		void AddInput(const std::string rootFileName,const std::string path_to_hist, const std::string legendText, 
							const std::string searchbins, const bool debug=false); 
		void DoSystematics(const bool b) { b_doSystematics = b; }
		void Run();
};

#endif

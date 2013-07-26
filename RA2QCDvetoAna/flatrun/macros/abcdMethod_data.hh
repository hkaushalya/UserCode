#ifndef ABCDMETHOD_DATA_HH
#define ABCDMETHOD_DATA_HH

#include <string>
#include "TF1.h"
#include "TH1.h"
#include "TCanvas.h"
#include <vector>
#include <utility>
#include "TFile.h"
#include "TFitResultPtr.h"

using namespace std;

class abcdMethod_data
{
	public:
		abcdMethod_data(const string bgFileName);
		~abcdMethod_data();

	
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
			string legText;
			TH1* hist_ratio;
			TF1* gaus_func;
			TF1* exp_func;
			float c_value;
			TFitResultPtr gaus_fitResPtr;
			TFitResultPtr exp_fitResPtr;
			double gaus_sys;
			double exp_sys;
			Predictions_t res_gaus;
			Predictions_t res_expo;
		};

		struct input_t {
			TFile* file;
			string legText;
			TH1* hist_ratio;
			TF1* gaus_func;
			TF1* exp_func;
			float c_value;
			TFitResultPtr gaus_fitResPtr;
			TFitResultPtr exp_fitResPtr;
			vector< pair<float, float> > njetrange;
			vector< pair<float, float> > htrange;
			vector< pair<float, float> > mhtrange;
			vector< vector<TH1* > > hist_control_finebin; 
			vector< vector<TH1* > > hist_signal_finebin; 
			vector< vector<TH1* > > hist_controlsyst1_finebin; 
			vector< vector<TH1* > > hist_controlsyst2_finebin; 
		
			vector< vector<vector< pair<abcdMethod_data::Predictions_t, abcdMethod_data::Predictions_t> > > > results;  //gaus/exp 
			//njet/ht/mht/systs
			//vector< vector< vector<vector< abcdMethod_data::Syst_t > > > > Systs; //holds 6 individual systematics. pred_mean is the only used variable
			//one fit per njet bin
			vector< abcdMethod_data::Syst_t > Systs; //holds 6 individual systematics. pred_mean is the only used variable
		};


	private:

		int Nsyst;
		bool b_doSystematics;
		string numHistName, denHistName, controlHistName, signalHistName;
		string controlsyst1HistName, constrolsyst2HistName;
		vector<TH1*> vRatioHist; //0=mean, >0 are all systematics
		TH1* GetRatioHist(TFile* f, const string path_to_hist, const bool debug=false);
		TH1* GetControlHist(TFile* f, const string path_to_hist, const bool debug=false);
		void GetSignalHist(TFile *f, vector< vector<TH1*> > & signalHists, 
										vector< pair<float, float> >& njetrange,
										vector< pair<float, float> >& htrange,
										vector< pair<float, float> >& mhtrange,
										const string path_to_hist, const bool debug);
		void GetControlHist(TFile *f, vector< vector<TH1*> > & controlHists, 
										vector< pair<float, float> >& njetrange,
										vector< pair<float, float> >& htrange,
										vector< pair<float, float> >& mhtrange,
										const string path_to_hist, const bool debug);
		void GetFits(TH1* ratio_hist, 
							TF1*& gausFunc, TFitResultPtr& gausFitResPtr, 
							TF1*& expFunc, TFitResultPtr& expFitResPtr,
							const float c_value,
							const string optTitleText="");
		void SetSearchBinInfo( vector< pair<float, float> >& njetrange,
										vector< pair<float, float> >& htrange,
										vector< pair<float, float> >& mhtrange,
										vector< vector<vector< pair<abcdMethod_data::Predictions_t, abcdMethod_data::Predictions_t> > > >& results, 
										const string searchbins);
		void InitSysts(
				//vector< vector< vector< vector< abcdMethod_data::Syst_t > > > > & systs, 
				vector< abcdMethod_data::Syst_t > & systs, 
				vector< pair<float, float> >& njetrange,
				vector< pair<float, float> >& htrange, 
				vector< pair<float, float> >& mhtrange, 
				const string path_to_hist, const bool debug);
		void PrintRatioCanvas();
		void PrintDebugCanvas();
		vector<input_t> vg_Inputs;
		TCanvas *canvas_for_ratios, *canvas_for_debug;
		string epsfile_for_ratios, epsfile_for_debug;
		void PrintResults(vector< pair<float, float> >& njetrange,
										vector< pair<float, float> >& htrange,
										vector< pair<float, float> >& mhtrange,
					const vector< vector<vector< pair<abcdMethod_data::Predictions_t, abcdMethod_data::Predictions_t> > > >& results);
		void  GetPredictions(const float& min_mht, const float& max_mht, 
									const TF1* f1, const TFitResultPtr fitResPtr, 
									const TH1* signalHist, const TH1* controlHist, Predictions_t& pred);
		double GetCerrValue(const int jetbin_low, const int jetbin_hi, const double x);
		double GetCsyst(TH1* h, const int jetbin_lo, const int jetbin_hi, const double mht_min, const double mht_max);
		void SubtractBackgrounds(TH1* h, const string histname, const bool bUseHtMerged=0);
		TFile *bgRootFile;
		TFile *bgHTmergedRootFile;
		string bgRootFileName, bgHTmergedRootFileName;

	public:
		void AddInput(const string rootFileName,const string path_to_hist, const string legendText, 
							const string searchbins, const float c_value, const bool debug=false); 
		void DoSystematics(const bool b) { b_doSystematics = b; }
		void Run();
};

#endif

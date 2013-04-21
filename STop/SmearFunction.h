#ifndef SMEAR_FUNCTION_H
#define SMEAR_FUNCTION_H

#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TKey.h"
#include "TFile.h"
#include "TMath.h"
#include "TArray.h"

#include <memory>
#include <string>
#include <vector>
#include <cassert>
#include <cmath>
#include <iostream>

using namespace std;

/***************************************************************
 * See cc file for class description.
 * Sam Hewamanage, Florida International University
 **************************************************************/

class SmearFunction {

public:
   SmearFunction();
   ~SmearFunction();

	void Init();
   TF1* getSigmaPtForRebalancing(int i_jet, int i_eta) const;
   TF1* getSigmaPtScaledForRebalancing(int i_jet, int i_eta) const;
   TH1F* getSmearFunc(int i_flav, int i_jet, int i_eta, int i_Pt) const;

private:
   typedef std::vector<std::string>::const_iterator StrIter;

   void FillSigmaHistsForRebalancing();
   void ResizeSmearFunctions();
   void CalculateSmearFunctions();

   double GetAdditionalSmearing(const double&, const double&);
   double GetLowerTailScaling(const double&, const double&);
   double GetUpperTailScaling(const double&, const double&);
   void FoldWithGaussian(const TH1&, TH1&, const double&);
   void StretchHisto(const TH1&, TH1&, const double&);
   int GetIndex(const double&, const std::vector<double>*);
  
   std::vector<double> PtBinEdges_scaling_;
   std::vector<double> EtaBinEdges_scaling_;
   std::vector<double> AdditionalSmearing_;
   std::vector<double> LowerTailScaling_;
   std::vector<double> UpperTailScaling_;
   double AdditionalSmearing_variation_;
   double LowerTailScaling_variation_;
   double UpperTailScaling_variation_;

   std::string inputhist1HF_;
   std::string inputhist2HF_;
   std::string inputhist3pHF_;
   std::string inputhist1NoHF_;
   std::string inputhist2NoHF_;
   std::string inputhist3pNoHF_;
   std::string smearingfile_;
   std::string bprobabilityfile_;
   std::string outputfile_;

   int NRebin_;
   bool absoluteTailScaling_;
   double A0RMS_;
   double A1RMS_;
   double probExtreme_;

   std::string uncertaintyName_;

   std::vector<double> PtBinEdges_;
   std::vector<double> EtaBinEdges_;

   //// vectors of response functions
   std::vector<std::vector<std::vector<std::vector<TH1F*> > > >smearFunc;
   std::vector<std::vector<std::vector<std::vector<TH1F*> > > >smearFunc_Core;
   std::vector<std::vector<std::vector<std::vector<TH1F*> > > >smearFunc_LowerTail;
   std::vector<std::vector<std::vector<std::vector<TH1F*> > > >smearFunc_UpperTail;
   std::vector<std::vector<std::vector<std::vector<TH1F*> > > >smearFunc_scaled;
   std::vector<std::vector<std::vector<TH1F*> > >SigmaPtHist;
   std::vector<std::vector<std::vector<TF1*> > >SigmaPt;
   std::vector<std::vector<std::vector<TH1F*> > >SigmaPtHist_scaled;
   std::vector<std::vector<std::vector<TF1*> > >SigmaPt_scaled;

   std::vector<std::vector<std::vector<TH1F*> > > smearFunc_total;
   std::vector<std::vector<TH1F*> > SigmaPtHist_total;
   std::vector<std::vector<TF1*> > SigmaPt_total;
   std::vector<std::vector<TH1F*> > SigmaPtHist_scaled_total;
   std::vector<std::vector<TF1*> > SigmaPt_scaled_total;

	unsigned uResFuncCollType_; //default=0 (JetAll), jet Ranked=1

public:


	string SmearingFile() const { return smearingfile_; }
	void SetSmearingFile(const string name) { smearingfile_ = name; }
	void SetResFuncColl(const unsigned i);
	unsigned GetResFuncCollType() const { return uResFuncCollType_; }

	//setters (override default settings)
	void SetPtBinEdges(const vector<double> ptBins);
	void SetEtaBinEdges(const vector<double> etaBins);

	void SetAbsoluteTailScaling        (const bool  b) { absoluteTailScaling_          = b; }
   void SetAdditionalSmearingVariation(const float v) { AdditionalSmearing_variation_ = v; }
   void SetLowerTailScalingVariation  (const float v) { LowerTailScaling_variation_   = v; }
   void SetUpperTailScalingVariation  (const float v) { UpperTailScaling_variation_   = v; }

	void SetPtBinEdges_scaling     (const vector<double>& vec) { PtBinEdges_scaling_  = vec; }
	void SetEtaBinEdges_scaling    (const vector<double>& vec) { EtaBinEdges_scaling_ = vec; }
	void SetLowerTail_scaling      (const vector<double>& vec) { LowerTailScaling_    = vec; }
	void SetUpperTail_scaling      (const vector<double>& vec) { UpperTailScaling_    = vec; }
	void SetAdditionalSmear_scaling(const vector<double>& vec) { AdditionalSmearing_  = vec; }

	vector<double>* GetPtBinEdges()  { return &PtBinEdges_;  }
	vector<double>* GetEtaBinEdges() { return &EtaBinEdges_; }
	bool   GetAbsoluteTailScaling() const { return absoluteTailScaling_; }
   double GetAdditionalSmearingVariation() const { return AdditionalSmearing_variation_; }
   double GetLowerTailScalingVariation  () const { return LowerTailScaling_variation_;   }
   double GetUpperTailScalingVariation  () const { return UpperTailScaling_variation_;   }

	vector<double> GetPtBinEdges_scaling     () const { return PtBinEdges_scaling_;  }
	vector<double> GetEtaBinEdges_scaling    () const { return EtaBinEdges_scaling_; }
	vector<double> GetLowerTail_scaling      () const { return LowerTailScaling_;    }
	vector<double> GetUpperTail_scaling      () const { return UpperTailScaling_;    }
	vector<double> GetAdditionalSmear_scaling() const { return AdditionalSmearing_;  }

	void UncertaintyName(const string s) { uncertaintyName_ = s; }
	string GetUncertaintyName() const { return uncertaintyName_; }

};

#endif

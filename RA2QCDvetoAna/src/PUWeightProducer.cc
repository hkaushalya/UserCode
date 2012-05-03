#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include <ctype.h>
#include <stdlib.h>
#include <exception>
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"


using namespace edm;
using namespace std;

class LumiWeightProducer : public edm::EDProducer{

  public:

    explicit LumiWeightProducer(const edm::ParameterSet & iConfig);
    ~LumiWeightProducer();
 private:
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void beginJob();
    virtual void endJob();
    virtual void beginRun(edm::Run&, const edm::EventSetup&);
    virtual void endRun(edm::Run&, const edm::EventSetup&);

    string PUscenario;

 	 double sumLumiWeights; //sum of lumi weights for cross checks
	 //to be added to event
	 int debug;

 	 edm::InputTag vertexSrc_;
    // histograms
	 TH1D* TNPUInTime_;
	 TH1D* TNPUTrue_;
	 TH1D* TNVTX_;

};

LumiWeightProducer::LumiWeightProducer(const edm::ParameterSet & iConfig) {

	PUscenario    = iConfig.getParameter<string> ("puScenario");
	debug          = iConfig.getUntrackedParameter<int> ("debug",0);

  	produces<double> ( "puWeight" );

	sumLumiWeights   = 0;
}

LumiWeightProducer::~LumiWeightProducer() {
}

// ------------ method called on each new Event  ------------
void LumiWeightProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

	//cout << "-------------------------------------------" << endl;

	lumiWeight = LumiWeights_.weight( iEvent );
	edm::Handle<edm::TriggerResults> triggerResults;
	iEvent.getByLabel(triggerResults_,triggerResults);

	if (! triggerResults.isValid()) 
	{
		cout << __FUNCTION__ << ":: No valid triggerResult found!" << endl;
		assert(false); 
	}

	bool accept = false;

   auto_ptr<double> pOut(new double(puWeight));
	iEvent.put(pOut,"puWeight");
}

// ------------ method called once each job just before starting event loop  ------------
void LumiWeightProducer::beginJob() {


 /*
	2011 Samples

	-	The 2011 samples have various designations that correspond to different pileup treatments. Here is a list

	-	PU_S1: 311X MC with the 2010 Pileup conditions
	-	PU_S2: 413X MC where the proper (correlated) out-of-time pileup was introduced. It was a test towards the Summer11 production... Made with 25ns bunch spacing. 'Flat10+tail distribution' (see below)
	-	PU_S3: 42X MC. The first round of Summer11 MC. Has random (no correlation) in-time and out-of-time pileup. (Flat10+tail)
	-	PU_S4: 42X MC. The remainder of the Summer11 MC. Has properly-correlated in-time and out-of-time pileup. (Flat10+tail)
	-	PU_S6: 42X MC. Fall11 MC. Has properly-correlated in-time and out-of-time pileup and Pileup Truth information. (Fall11 Distribution) 
*/



	//2011 Pileup Scenarios https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupInformation#2011_Pileup_Scenarios
	//"Flat to 10 plus tail" scenario shown above, the relative normalization of each bin is: 
	vMCNvtxDist.push_back(0.069286816); //1
	vMCNvtxDist.push_back(0.069286816); //2
	vMCNvtxDist.push_back(0.069286816); //3
	vMCNvtxDist.push_back(0.069286816); //4
	vMCNvtxDist.push_back(0.069286816); //5
	vMCNvtxDist.push_back(0.069286816); //6
	vMCNvtxDist.push_back(0.069286816); //7
	vMCNvtxDist.push_back(0.069286816); //8
	vMCNvtxDist.push_back(0.069286816); //9
	vMCNvtxDist.push_back(0.069286816); //10
	vMCNvtxDist.push_back(0.069286816); //11
	vMCNvtxDist.push_back(0.06518604 ); //12
	vMCNvtxDist.push_back(0.053861878); //13
	vMCNvtxDist.push_back(0.040782032); //14
	vMCNvtxDist.push_back(0.030135062); //15
	vMCNvtxDist.push_back(0.019550796); //16
	vMCNvtxDist.push_back(0.012264707); //17
	vMCNvtxDist.push_back(0.007449117); //18
	vMCNvtxDist.push_back(0.004502075); //19
	vMCNvtxDist.push_back(0.002194605); //20
	vMCNvtxDist.push_back(0.001166276); //21
	vMCNvtxDist.push_back(0.000476543); //22
	vMCNvtxDist.push_back(0.000188109); //23
	vMCNvtxDist.push_back(7.52436E-05); //24
	vMCNvtxDist.push_back(1.25406E-05); //25

	//DATA Nvtx distribution from Pileup_2011_EPS_8_jul.root
	vDATANvtxDist.push_back(1.45417e+07); //1
	vDATANvtxDist.push_back(3.47743e+07); //2
	vDATANvtxDist.push_back(7.89247e+07); //3
	vDATANvtxDist.push_back(1.26467e+08); //4
	vDATANvtxDist.push_back(1.59329e+08); //5
	vDATANvtxDist.push_back(1.67603e+08); //6
	vDATANvtxDist.push_back(1.52684e+08); //7
	vDATANvtxDist.push_back(1.23794e+08); //8
	vDATANvtxDist.push_back(9.09462e+07); //9
	vDATANvtxDist.push_back(6.13973e+07); //10
	vDATANvtxDist.push_back(3.8505e+07); //11
	vDATANvtxDist.push_back(2.2628e+07); //12
	vDATANvtxDist.push_back(1.25503e+07); //13
	vDATANvtxDist.push_back(6.61051e+06); //14
	vDATANvtxDist.push_back(3.32403e+06); //15
	vDATANvtxDist.push_back(1.60286e+06); //16
	vDATANvtxDist.push_back(743920); //17
	vDATANvtxDist.push_back(333477); //18
	vDATANvtxDist.push_back(144861); //19
	vDATANvtxDist.push_back(61112.7); //20
	vDATANvtxDist.push_back(25110.2); //21
	vDATANvtxDist.push_back(10065.1); //22
	vDATANvtxDist.push_back(3943.98); //23
	vDATANvtxDist.push_back(1513.54); //24
	vDATANvtxDist.push_back(896.161); //25

	LumiWeights_ = edm::LumiReWeighting(vMCNvtxDist, vDATANvtxDist);


}

// ------------ method called once each job just after ending the event loop  ------------
void LumiWeightProducer::endJob() {
}

// ------------ method called once each run just before starting event loop  ------------
void LumiWeightProducer::beginRun(edm::Run &run, const edm::EventSetup& iSetup) {


}

// ------------ method called once each run just after starting event loop  ------------
void LumiWeightProducer::endRun(edm::Run &run, const edm::EventSetup& iSetup) {
}


#include "FWCore/Framework/interface/MakerMacros.h"
//define this as a plug-in
DEFINE_FWK_MODULE(LumiWeightProducer);

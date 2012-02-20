#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"

using namespace edm;
using namespace std;

class RunFilter : public edm::EDFilter{

  public:

    explicit RunFilter(const edm::ParameterSet & iConfig);
    ~RunFilter();
    virtual bool filter(edm::Event&, const edm::EventSetup&);
    virtual void beginJob();
    virtual void endJob();
    virtual bool beginRun(edm::Run&, const edm::EventSetup&);
    virtual bool endRun(edm::Run&, const edm::EventSetup&);


    double minRunNumber, maxRunNumber;

};

RunFilter::RunFilter(const edm::ParameterSet & iConfig) {

   minRunNumber = iConfig.getParameter<double>("lowRunNumber");
   maxRunNumber = iConfig.getParameter<double>("hiRunNumber");
 
}

RunFilter::~RunFilter() {
}

// ------------ method called on each new Event  ------------
bool RunFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {

	const double run  = iEvent.id().run();
	bool pass = true;
	//std::cout << "run = " << run << std::endl;
	if (run < minRunNumber || run > maxRunNumber) pass = false;

   return pass;

}

// ------------ method called once each job just before starting event loop  ------------
void RunFilter::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void RunFilter::endJob() {
		std::cout << __FILE__ << ":: Processed Run Range : " << minRunNumber << " to " << maxRunNumber << std::endl;
}

// ------------ method called once each run just before starting event loop  ------------
bool RunFilter::beginRun(edm::Run &run, const edm::EventSetup& iSetup) {
  return true;
}

// ------------ method called once each run just after starting event loop  ------------
bool RunFilter::endRun(edm::Run &run, const edm::EventSetup& iSetup) {
  return true;
}


#include "FWCore/Framework/interface/MakerMacros.h"
//define this as a plug-in
DEFINE_FWK_MODULE(RunFilter);

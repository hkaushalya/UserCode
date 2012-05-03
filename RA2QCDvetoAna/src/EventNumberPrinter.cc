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

class EventNumberPrinter : public edm::EDFilter{

  public:

    explicit EventNumberPrinter(const edm::ParameterSet & iConfig);
    ~EventNumberPrinter();
    virtual bool filter(edm::Event&, const edm::EventSetup&);
    virtual void beginJob();
    virtual void endJob();
    virtual bool beginRun(edm::Run&, const edm::EventSetup&);
    virtual bool endRun(edm::Run&, const edm::EventSetup&);



};

EventNumberPrinter::EventNumberPrinter(const edm::ParameterSet & iConfig) {
 
}

EventNumberPrinter::~EventNumberPrinter() {
}

// ------------ method called on each new Event  ------------
bool EventNumberPrinter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {

	cout << ">>>>>>>>>>> Run:Lumi:Event:: " 
		<< std::setw(10) << iEvent.id().run()
	   << std::setw(6)  << iEvent.id().luminosityBlock()
	   << std::setw(15) << iEvent.id().event()
		<<   "<<<<<<<<<<< "
		<< endl;

		return true;
}

// ------------ method called once each job just before starting event loop  ------------
void EventNumberPrinter::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void EventNumberPrinter::endJob() {
}

// ------------ method called once each run just before starting event loop  ------------
bool EventNumberPrinter::beginRun(edm::Run &run, const edm::EventSetup& iSetup) {
  return true;
}

// ------------ method called once each run just after starting event loop  ------------
bool EventNumberPrinter::endRun(edm::Run &run, const edm::EventSetup& iSetup) {
  return true;
}


#include "FWCore/Framework/interface/MakerMacros.h"
//define this as a plug-in
DEFINE_FWK_MODULE(EventNumberPrinter);

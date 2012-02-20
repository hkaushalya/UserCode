// -*- C++ -*-
//
// Package:    RA2QCDvetoAna
// Class:      GenPtHatFilter
// 
/**\class GenPtHatFilter GenPtHatFilter.cc UserCode/RA2QCDvetoAna/src/GenPtHatFilter.cc

 Description: Studies the DelPhi(jet,MET) for QCD veto.

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  samantha hewamanage
//         Created:  Tue Jul 26 22:15:44 CDT 2011
// $Id$
//
//


// system include files
#include <memory>
#include <vector>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
 
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/VertexReco/interface/print.h"
#include "FWCore/Utilities/interface/Verbosity.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

//ROOT
#include "TH1F.h"
#include "TH2F.h"
#include "assert.h"
#include "TVector3.h"
#include <iomanip>

//
// class declaration
//

class GenPtHatFilter : public edm::EDFilter {
   public:
      explicit GenPtHatFilter(const edm::ParameterSet&);
      ~GenPtHatFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
		

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------
		//jet collections
//		edm::InputTag genInputTag_;

		double minPtHat;
		TH1F *pthat_b4, *pthat_a4;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GenPtHatFilter::GenPtHatFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
//	genInputTag_ = iConfig.getParameter<edm::InputTag>("genInfoInputTag");
	minPtHat = iConfig.getUntrackedParameter<double>("dMinPtHat",0.0);

	edm::Service<TFileService> fs;

	//hists before any cuts
	pthat_b4 = fs->make<TH1F> ("pthat_b4" ,"pt-hat before cuts ;p_{T}-hat [GeV];Events;", 4000, 0, 4000);
	pthat_a4 = fs->make<TH1F> ("pthat_a4" ,"pt-hat after cuts ;p_{T}-hat [GeV];Events;", 4000, 0, 4000);

	//these are general event hist to check the PAT tuples cuts
}


GenPtHatFilter::~GenPtHatFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
GenPtHatFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

	//iEvent.getByLabel(genInputTag_, genHandle);
	edm::Handle<GenEventInfoProduct> genHandle;

	if (iEvent.isRealData()) return true;

	const bool genSource = iEvent.getByLabel("generator", genHandle);
 
	if (! genSource) 
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":genHandle handle not found!" << std::endl;
		assert(false);
	}

	const float ptHat = genHandle->binningValues()[0];
	pthat_b4->Fill(ptHat);
	bool pass = true;
	if (ptHat < minPtHat) pass = false;
	else { pthat_a4->Fill(ptHat); }

	//std::cout << minPtHat << " / " << ptHat  << " (pass = " << pass << ")" << std::endl;
   return pass;
}

// ------------ method called once each job just before starting event loop  ------------
void 
GenPtHatFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenPtHatFilter::endJob() {
	std::cout << __FILE__ << ":"<<  __FUNCTION__ << ":: minPtHat = " << minPtHat << std::endl;
}

// ------------ method called when starting to processes a run  ------------
bool 
GenPtHatFilter::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
GenPtHatFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
GenPtHatFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
GenPtHatFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenPtHatFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenPtHatFilter);

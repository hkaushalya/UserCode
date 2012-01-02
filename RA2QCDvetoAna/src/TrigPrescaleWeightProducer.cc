#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTfilters/interface/HLTHighLevel.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include <ctype.h>
#include <stdlib.h>
#include <exception>

using namespace edm;
using namespace std;

class TrigPrescaleWeightProducer : public edm::EDProducer{

  public:

    explicit TrigPrescaleWeightProducer(const edm::ParameterSet & iConfig);
    ~TrigPrescaleWeightProducer();
 private:
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void beginJob();
    virtual void endJob();
    virtual void beginRun(edm::Run&, const edm::EventSetup&);
    virtual void endRun(edm::Run&, const edm::EventSetup&);

	 void WildCardRemove(vector<string>& newHLTPathsByName, const edm::TriggerNames&);
	 float GetHTthreshold(const string trigName);
    vector<string> TokenizeByAsterix(const string str);
		
	 edm::InputTag triggerResults_;
    std::vector<std::string > HLTPathsByName_;
    std::vector<unsigned int> HLTPathsByIndex_;
    HLTConfigProvider hltConfig_;
    string _processName;

	 double prescaleWeight;
	 string highestPrescaledTriggerName;
	 //to be added to event
	 vector<string> trigsFired;
	 int debug;

};

TrigPrescaleWeightProducer::TrigPrescaleWeightProducer(const edm::ParameterSet & iConfig) {

	triggerResults_ = iConfig.getParameter<edm::InputTag>("triggerResults");
	HLTPathsByName_ = iConfig.getParameter<std::vector<std::string > >("hltPaths");
	_processName    = iConfig.getParameter<string> ("HLTProcess");
	debug           = iConfig.getUntrackedParameter<int> ("debug",0);

  produces<double> ( "prescaleWeight" );
  produces<string> ( "highestPrescaledTriggerName" );

}

TrigPrescaleWeightProducer::~TrigPrescaleWeightProducer() {
}

// ------------ method called on each new Event  ------------
void TrigPrescaleWeightProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

	//cout << "-------------------------------------------" << endl;

	// retreive TriggerResults from the event
	edm::Handle<edm::TriggerResults> triggerResults;
	iEvent.getByLabel(triggerResults_,triggerResults);

	if (! triggerResults.isValid()) 
	{
		cout << __FUNCTION__ << ":: No valid triggerResult found!" << endl;
		assert(false); 
	}

	bool accept = false;

	if (debug) std::cout << "TriggerResults found, number of HLT paths: " << triggerResults->size() << std::endl;
	// get trigger names
	const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerResults);
	if (debug) {
		for (unsigned int i=0; i < triggerNames.size(); i++) {
			if (debug) std::cout << "trigger path= " << triggerNames.triggerName(i) << std::endl;
		}
	}

	//triggerIndex method DOES NOT RECOGNIZE WIDLCARDS!!!!
	//It requires a string. If input name is not found, it returns
	//a values equal to the size of the trigger list.
	//So I need to strip out the wild crads and replace with matching trigger names!

	vector<string> newHLTPathsByName;
	WildCardRemove(newHLTPathsByName, triggerNames);

	HLTPathsByIndex_.resize(newHLTPathsByName.size());

	unsigned int n = newHLTPathsByName.size();

	for (unsigned int i=0; i!=n; i++) {
		HLTPathsByIndex_[i] = triggerNames.triggerIndex(newHLTPathsByName[i]);
	}

	vector<string> firedTrigNames;
	vector<float> firedTrigPrescale;
	unsigned highestHTtrigIndex = 99999;
	int  highestHTprescale = 1.0;
	float highestHT = 0.0;
	// count number of requested HLT paths which have fired
	unsigned int fired=0;
	for (unsigned int i=0; i < HLTPathsByIndex_.size() ; i++) {
		if (HLTPathsByIndex_.at(i) < triggerResults->size()) {
			if (triggerResults->accept(HLTPathsByIndex_.at(i))) {
				fired++;
				//std::cout << "Fired HLT path= " << newHLTPathsByName[i] << std::endl;
				accept = true;
				const int currentPrescale = hltConfig_.prescaleValue(iEvent, iSetup, newHLTPathsByName.at(i));
				if (debug) std::cout << "Fired HLT path= " << newHLTPathsByName[i] << " with prescale = " 
							<< currentPrescale << std::endl;

				const float ht = GetHTthreshold(newHLTPathsByName.at(i));
				if (ht > highestHT)
				{
					highestHT = ht;
					highestHTtrigIndex = i;
					highestHTprescale = currentPrescale;
				}
			}
		}
	}
	if (highestHTtrigIndex < 99999)
	{
		if (debug) cout << "Highest trigger is: " << newHLTPathsByName.at(highestHTtrigIndex) << " with prescale = " << highestHTprescale << endl;
		prescaleWeight = highestHTprescale;
		highestPrescaledTriggerName = newHLTPathsByName.at(highestHTtrigIndex);
	} else {
		cout << __FILE__ << "::" << __FUNCTION__ << ":: WARNING!! NO triggers fired for this event!!!!" << endl;
		prescaleWeight = 1;
		highestPrescaledTriggerName = "";
	}
	
   auto_ptr<double> pOut(new double(prescaleWeight));
   auto_ptr<string> pOut2(new string(highestPrescaledTriggerName));
	iEvent.put(pOut,"prescaleWeight");
	iEvent.put(pOut2,"highestPrescaledTriggerName");

	//if (fired > 1) cout <<  "\t :: accept = " << accept  << " /nFired = " << fired << endl;
}

// ------------ method called once each job just before starting event loop  ------------
void TrigPrescaleWeightProducer::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void TrigPrescaleWeightProducer::endJob() {
}

// ------------ method called once each run just before starting event loop  ------------
void TrigPrescaleWeightProducer::beginRun(edm::Run &run, const edm::EventSetup& iSetup) {

   bool changed(true);
   if (hltConfig_.init(run, iSetup, _processName, changed)) {
      if (changed) {
         //         hltConfig_.dump("Streams");
         //         hltConfig_.dump("Datasets");
         //         hltConfig_.dump("PrescaleTable");
         //         hltConfig_.dump("ProcessPSet");
      }
   } else {
      // if init returns FALSE, initialisation has NOT succeeded, which indicates a problem
      // with the file and/or code and needs to be investigated!
      cout << " HLT config extraction failure with process name " << _processName << endl;
      // In this case, all access methods will return empty values!
   }

}

// ------------ method called once each run just after starting event loop  ------------
void TrigPrescaleWeightProducer::endRun(edm::Run &run, const edm::EventSetup& iSetup) {
}


void TrigPrescaleWeightProducer::WildCardRemove(vector<string>& newHLTPathsByName, const edm::TriggerNames &trigNames)
{
	for (unsigned i=0; i< HLTPathsByName_.size(); ++i)
	{
		if (debug) cout << " HLTPathsByName_["<< i << "] = " << HLTPathsByName_.at(i) << endl;
		
		vector<string> name_parts = TokenizeByAsterix(HLTPathsByName_.at(i));
		
		//now fine the actual trigger path name that matches closley
		for (unsigned k =0; k < trigNames.size(); ++k)
		{
			bool match = true;
			for (unsigned p = 0; p < name_parts.size(); p++)
			{
				if ((trigNames.triggerName(k)).find(name_parts.at(p)) != string::npos)
				{
					match = match & true;
					//cout << ">>Match found for : " << name_parts.at(p) << endl;
				} else
				{
					match = match & false;
				}
			}

		   if (match)
			{
				newHLTPathsByName.push_back(trigNames.triggerName(k));
				//cout << "\t matching paths = " << trigNames.triggerName(k) << endl;
			}
		}
	}

}


float TrigPrescaleWeightProducer::GetHTthreshold(const string trigName)
{
	int ndigits = 0; //max should be 4
	
	vector<size_t> pos;
	const size_t p1 = trigName.find_first_of('_');
	size_t p2 = string::npos;
	if (p1 != string::npos) p2 = trigName.find('_',p1+1);
	if (p2 == string::npos) cout << " p2 is not found! " << endl;
	
	string part = trigName.substr(p1+3, p2-6);

	float num = atof(part.c_str());
	//cout << trigName <<  " :: part = " << part  << ":: cut = " << num << endl;

	return num;

}

vector<string> TrigPrescaleWeightProducer::TokenizeByAsterix(const string str)
{
	vector<size_t> asterix_positions;
	for (unsigned p=0; p < str.length(); ++p)
	{
//		cout << p << " = " << str[p] << endl;
		if (str[p] == '*') 
		{
//			cout << "found * at = " << p << endl;
			asterix_positions.push_back(p);
		}
	}

	vector<string> parts;

	if (asterix_positions.size() > 0) //for
	{
		string part1 = str.substr(0, asterix_positions.at(0));
		parts.push_back(part1);
//		cout << "part 1 = " << part1 << endl;
	}

	if (asterix_positions.size() > 1) //for
	{
		const string str_part2 = str.substr(asterix_positions.at(0)+1, str.length());
		//cout << str_part2 << endl;
		string part2 = str_part2.substr(0,str_part2.length()-1); 
		parts.push_back(part2);
//		cout << "part 2 = " << part2 << endl;
	}

	if (asterix_positions.size() > 2)
	{
		cout << __FUNCTION__ << ":: WARNING! More then 2 wild cards found in the trigger names!" << endl;
		assert(false);
	}


	return parts;
}

#include "FWCore/Framework/interface/MakerMacros.h"
//define this as a plug-in
DEFINE_FWK_MODULE(TrigPrescaleWeightProducer);

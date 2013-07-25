
import FWCore.ParameterSet.Config as cms

triggerFilter = cms.EDFilter(
		"TriggerFilter",
		trigWeight = cms.InputTag('trigWgtProd:prescaleWeight'),
		highestPrescaledTrigName = cms.InputTag('trigWgtProd:highestPrescaledTriggerName'),
		triggerResults = cms.InputTag("TriggerResults::HLT"),
		HLTProcess = cms.string("HLT"),
		hltPaths = cms.vstring('HLT_HT100_v*','HLT_HT150_v*','HLT_HT160_v*','HLT_HT200_v*','HLT_HT240_v*','HLT_HT250_v*','HLT_HT260_v*','HLT_HT300_v*','HLT_HT350_v*','HLT_HT400_v*','HLT_HT450_v*','HLT_HT500_v*','HLT_HT550_v*','HLT_HT600_v*','HLT_HT650_v*','HLT_HT700_v*','HLT_HT750_v*','HLT_HT800_v*','HLT_HT850_v*','HLT_HT2000_v*','HLT_HT*_L1FastJet_v*'),
		debug = cms.int32(0)
 )


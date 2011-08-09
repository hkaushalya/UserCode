#from __future__ import print_function
import sys, string
import FWCore.ParameterSet.Config as cms
from JetMETAnalysis.ecalDeadCellTools.EcalDeadCellEventFilter_cfi import *
import sys, string

runningDATA = False
evts2Process = 5000
outputFile = 'METFilter_TTbarMC.root'
skimFile = 'TrkFiltOutEventsFromTTbarMC.root'
print 'runningDATA=', runningDATA, ' Event to be processed = ', evts2Process 
print 'Output files will be : ', outputFile , ', ' , skimFile 
process = cms.Process("anomtrk")

## Options and Output Report
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.PyReleaseValidation.autoCond import autoCond
process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )
process.load("Configuration.StandardSequences.MagneticField_cff")

#Global tag picks up necessay calibration etc for some of the filters
#like the TPEcal filter. The tag is different for MC and DATA.
#https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions
process.GlobalTag.globaltag = "MC_42_V13::All"
print 'Setting Global Tag to ', process.GlobalTag.globaltag 


process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(evts2Process) )

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
		'/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/FEEC7BB2-E497-E011-8F47-003048678B30.root',
		'/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/FEB5DA32-E297-E011-9901-0030486792B6.root',
		'/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/FE7166D2-7798-E011-9A7A-0026189438C0.root',
		'/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/FE525894-E897-E011-A74F-0018F3D0961E.root',
		'/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/FCD5AD8C-E797-E011-8AE2-00261894393E.root',
		'/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/FCBCF980-7998-E011-9676-00261894393C.root',
		'/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/FCB97880-A398-E011-9D3F-0018F3D096AE.root',
		'/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/FCB78C13-F597-E011-91AC-0018F3D0963C.root',
		'/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/FC87DB1C-7C98-E011-944A-00261894387E.root',
		'/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/FC39ABAB-7A98-E011-9D78-00261894385A.root',
		'/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/FC1C32F7-9C98-E011-B8A8-0018F34D0D62.root',
		'/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/FC0955A5-EE97-E011-952C-003048679236.root',
		'/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/FAE74503-E797-E011-87FB-001A92971BB8.root',
		'/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/FA55A8A8-F397-E011-8E23-003048678FF4.root',
		'/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/FA535A2A-7698-E011-BA91-0018F34D0D62.root',
		'/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/FA19149D-7498-E011-87E8-00248C0BE01E.root',
		'/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/F8DA3948-E597-E011-AC89-002618943832.root',
		'/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/F8C08759-E197-E011-A02E-0018F3D0967A.root',
		'/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/F8B610F1-7998-E011-9CB1-0018F3D09624.root',
		'/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/F8B252AF-A098-E011-AB17-0018F3D09608.root',
		'/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/F853C7FA-1298-E011-997C-003048678FFE.root',
		'/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/F6DC3E54-7498-E011-915C-0018F3D09678.root',
		'/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/F656870C-0998-E011-B394-001A92971B82.root',
		'/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/F64376FE-8F98-E011-9A98-002618943856.root',
		'/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/F6380E65-E897-E011-A7FC-002618943981.root',
		'/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/F624F7D7-8198-E011-9D59-00261894393A.root',
		'/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/F610EDD1-EC97-E011-B481-002618943836.root',
		'/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/F4EC37F5-7698-E011-849F-003048679180.root'
		)
)

#secFiles.extend( [
#] )


process.TFileService = cms.Service("TFileService",
		fileName = cms.string(outputFile)
		)

process.out = cms.OutputModule("PoolOutputModule",
		fileName = cms.untracked.string(skimFile),
		# save only events passing the full path
		SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
		# save PAT Layer 1 output; you need a '*' to
		# unpack the list of commands 'patEventContent'
		# outputCommands = cms.untracked.vstring('drop *', *patEventContent ) #To further strip down the eventsize
)

#### remove beam scraping events
process.noScraping= cms.EDFilter("FilterOutScraping",
						applyfilter = cms.untracked.bool(True),
						debugOn = cms.untracked.bool(False), ## Or 'True' to get some per-event info
						numtrack = cms.untracked.uint32(10),
						thresh = cms.untracked.double(0.25)
)
process.goodVertices = cms.EDFilter(
  "VertexSelector",
  filter = cms.bool(False),
  src = cms.InputTag("offlinePrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
)
process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')
process.HBHENoiseFilter.minIsolatedNoiseSumE = cms.double(999999.) 
process.HBHENoiseFilter.minNumIsolatedNoiseChannels = cms.int32(999999) 
process.HBHENoiseFilter.minIsolatedNoiseSumEt = cms.double(999999.)
process.load('RecoMET.METAnalyzers.CSCHaloFilter_cfi') 
process.load('JetMETAnalysis.ecalDeadCellTools.EcalDeadCellEventFilter_cfi')
process.ecalDeadCellTPfilter = EcalDeadCellEventFilter.clone()
process.ecalDeadCellTPfilter.tpDigiCollection = cms.InputTag("ecalTPSkim")
process.ecalDeadCellTPfilter.etValToBeFlagged = cms.double(63.75)
process.ecalDeadCellTPfilter.ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB")
process.ecalDeadCellTPfilter.eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE")

process.anomtrk = cms.EDFilter('AnomTrkMETFilter',
							trackCollection  = cms.InputTag('generalTracks'),
							caloJetInputTag_ = cms.InputTag("ak5CaloJets"),
							pfJetInputTag_   = cms.InputTag("ak5PFJets"),
							HltTriggerResults = cms.InputTag("TriggerResults::HLT"),
							req_trigger = cms.untracked.bool(False),
							TriggerPathsToStore = cms.vstring("HLT_MET65_v3"),
							nMinVtx = cms.untracked.int32(1),	
							ndofVtx = cms.untracked.int32(4),	
							maxDelzTrkVtx = cms.untracked.double(24.0),	
							minMet = cms.untracked.double(0.0),	
							print_hists = cms.untracked.bool(True),
							verbose = cms.untracked.int32(0)	
)

process.myfilter_path = cms.Path(process.noScraping 
										* process.goodVertices
										* process.HBHENoiseFilter 
										* process.CSCTightHaloFilter 
										* process.ecalDeadCellTPfilter
										* process.anomtrk
)

process.out.SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('myfilter_path') )

process.outpath = cms.EndPath(process.out)


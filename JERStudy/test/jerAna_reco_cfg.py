import FWCore.ParameterSet.Config as cms
from JetMETAnalysis.ecalDeadCellTools.EcalDeadCellEventFilter_cfi import *
import sys, string

evts2Process = 10000
outputFile = 'JER.root'
print 'Output files will be : ', outputFile 
process = cms.Process("jerana")

## Options and Output Report
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

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

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
	'/store/mc/Summer11/QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6/GEN-SIM-RECO/PU_S3_START42_V11-v2/0004/FEE60FB5-457E-E011-98A5-0026189438B8.root'
	)
)

#secFiles.extend( [
#] )


process.TFileService = cms.Service("TFileService",
		fileName = cms.string(outputFile)
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
#std. RA2 halo filter
process.load('SandBox.Skims.beamHaloFilter_cfi')
#process.load('RecoMET.METAnalyzers.CSCHaloFilter_cfi') 

process.load('JetMETAnalysis.ecalDeadCellTools.EcalDeadCellEventFilter_cfi')
process.ecalDeadCellTPfilter = EcalDeadCellEventFilter.clone()
process.ecalDeadCellTPfilter.tpDigiCollection = cms.InputTag("ecalTPSkim")
process.ecalDeadCellTPfilter.etValToBeFlagged = cms.double(63.75)
process.ecalDeadCellTPfilter.ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB")
process.ecalDeadCellTPfilter.eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE")

process.jerana = cms.EDFilter('recoJERStudy',
							detJetsInputTag   = cms.InputTag("ak5PFJets"),
							genJetsInputTag   = cms.InputTag("ak5GenJets"),
							#HltTriggerResults = cms.InputTag("TriggerResults::HLT"),
							#req_trigger = cms.untracked.bool(False),
							#TriggerPathsToStore = cms.vstring("HLT_MET65_v3"),
							#nMinVtx = cms.untracked.int32(1),	
							#ndofVtx = cms.untracked.int32(4),	
							#maxDelzTrkVtx = cms.untracked.double(24.0),	
							#minMet = cms.untracked.double(0),	
							verbose = cms.untracked.int32(0)
)

process.myfilter_path = cms.Path(
										#process.noScraping 
										#* process.goodVertices
										#* process.HBHENoiseFilter 
										#* process.beamHaloFilter 
										#* process.ecalDeadCellTPfilter
										#* 
										process.jerana
)

#process.out.SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('myfilter_path') )

#process.outpath = cms.EndPath(process.out)

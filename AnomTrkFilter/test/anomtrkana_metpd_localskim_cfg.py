import FWCore.ParameterSet.Config as cms
from JetMETAnalysis.ecalDeadCellTools.EcalDeadCellEventFilter_cfi import *
import sys, string

runningDATA = True
evts2Process = -1

badEvents = True

if badEvents: 
	outputFile = 'METFilter_JETPDSkim_Bad.root'
	skimFile = 'TrkFiltOutEventsFromJETPDSkim_Bad.root'
else:
	outputFile = 'METFilter_METPDSkim.root'
	skimFile = 'TrkFiltOutEventsFromMETPDSkim.root'

print 'Output files will be : ', outputFile , ', ' , skimFile 

process = cms.Process("anomtrk")

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
#if runningDATA:  process.GlobalTag.globaltag = "GR_R_42_V19:All"
#if runningDATA:  process.GlobalTag.globaltag = "GR_R_42_V6::All"
if runningDATA:  process.GlobalTag.globaltag = "GR_R_42_V12::All"
else: process.GlobalTag.globaltag = "MC_42_V13::All"
print 'Setting Global Tag to ', process.GlobalTag.globaltag 


process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(evts2Process) )

if badEvents:
	process.source = cms.Source("PoolSource",
   	fileNames = cms.untracked.vstring(
					#'file:./output/TrkFiltOutEventsFromMETPD.root',
					#'file:/uscms_data/d3/samantha/TrkFiltOutEventsFromJETPD.root'
#this has one of Tai's bad events: 167151:20:22174344
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/151/94223088-409B-E011-8FB5-003048F11C58.root',
#163270:409:252238770 from METB-tag sample
			'/store/data/Run2011A/METBTag/RECO/May10ReReco-v1/0001/40321B37-807E-E011-806C-001A92811720.root',
#165415:103:91591777 , /Jet/Run2011A-PromptReco-v4/RECO
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/165/415/A44BE14A-BD85-E011-964F-001D09F24024.root'
   	)
	)
else:
	process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring('/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/151/94223088-409B-E011-8FB5-003048F11C58.root'
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


process.anomtrk = cms.EDFilter('AnomTrkMETFilter',
							trackCollection  = cms.InputTag('generalTracks'),
							caloJetInputTag_ = cms.InputTag("ak5CaloJets"),
							pfJetInputTag_   = cms.InputTag("ak5PFJets"),
							HltTriggerResults = cms.InputTag("TriggerResults::HLT"),
							print_hists = cms.untracked.bool(True),
							req_trigger = cms.untracked.bool(False),
							TriggerPathsToStore = cms.vstring("HLT_MET65_v3"),
							nMinVtx = cms.untracked.int32(1),	
							ndofVtx = cms.untracked.int32(4),	
							maxDelzTrkVtx = cms.untracked.double(24.0),	
							minMet = cms.untracked.double(0.0),	
							verbose = cms.untracked.int32(0),
							printHists = cms.untracked.bool(False),
							processBadOnly = cms.untracked.bool(badEvents)
							#processBadOnly = cms.untracked.bool(False)

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


process.myfilter_path = cms.Path(process.noScraping 
										* process.goodVertices
										* process.HBHENoiseFilter 
										* process.beamHaloFilter
										* process.ecalDeadCellTPfilter
										* process.anomtrk
)

process.out.SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('myfilter_path') )

process.outpath = cms.EndPath(process.out)


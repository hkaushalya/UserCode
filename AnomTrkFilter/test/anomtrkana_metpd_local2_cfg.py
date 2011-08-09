import FWCore.ParameterSet.Config as cms
from JetMETAnalysis.ecalDeadCellTools.EcalDeadCellEventFilter_cfi import *
import sys, string

runningDATA = True
evts2Process = 5000
outputFile = 'METFilter_METPD_local2.root'
skimFile = 'TrkFiltOutEventsFromMETPDSkim_local2.root'
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
process.GlobalTag.globaltag = "GR_R_42_V6::All"
print 'Setting Global Tag to ', process.GlobalTag.globaltag 


process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(evts2Process) )

process.source = cms.Source("PoolSource",
		fileNames = cms.untracked.vstring(
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/830/10173618-19A3-E011-A8B1-BCAEC5329702.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/807/FAA1BE64-EDA2-E011-8D5D-003048F11C28.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/807/F40207C3-FFA2-E011-8B59-001D09F251B8.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/807/F27FA7E1-F7A2-E011-92DE-003048F024E0.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/807/E8E4FEB8-FFA2-E011-B5AB-001D09F2A690.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/807/D83B0C84-EAA2-E011-AC15-003048D2BC42.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/807/CE500F85-EAA2-E011-BAE6-0030487CAF5E.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/807/CC6F714B-B9A2-E011-ACD5-003048D375AA.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/807/C6A76336-BEA2-E011-B7FE-003048D2C174.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/807/B47CACB7-AEA2-E011-913E-0030487CD716.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/807/A8F9C766-E1A2-E011-9F4F-003048F118DE.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/807/A0FCA8C3-BCA2-E011-9BB4-003048F118C2.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/807/6AFED713-7CA2-E011-9E6E-BCAEC5364C93.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/807/680E4528-D8A2-E011-A57D-BCAEC518FF56.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/807/54C2EFBB-AEA2-E011-A917-003048D2C0F0.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/807/4A400F20-8AA2-E011-BEEE-003048F1BF68.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/807/40CE1EEF-CCA2-E011-BAC4-003048F0258C.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/807/3A5F8D38-B2A2-E011-8CD3-003048F1C832.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/807/1AACFA32-BEA2-E011-810A-003048D2C16E.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/807/10037ECB-C8A2-E011-8468-003048D375AA.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/807/0CC1C2F0-A4A2-E011-A0F4-0030487CD6D8.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/786/F89808F6-F5A1-E011-BDE0-003048F118C2.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/786/CA4075D3-FAA1-E011-8F1F-003048F1C420.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/786/4AA476D5-FAA1-E011-8396-003048D3750A.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/785/4EC9660C-C4A1-E011-8C2F-0030487CD77E.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/784/BC00F24F-F0A1-E011-BA6B-0030486780B8.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/784/BA40B54E-F0A1-E011-B876-003048F118C2.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/754/50F0BDC0-A3A1-E011-897F-003048F118C6.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/746/FEA137D3-A3A1-E011-A845-0030486730C6.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/746/A2675B87-85A1-E011-AC8F-003048F1110E.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/746/728076D8-A3A1-E011-9884-003048D2C108.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/746/641FEB68-7CA1-E011-93EE-001D09F24F65.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/746/184F73BA-82A1-E011-81B4-003048F024C2.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/740/E40FD6A7-5AA1-E011-9D4C-001D09F2516D.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/740/B8F95FCA-63A1-E011-B504-BCAEC532971F.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/715/5CA8CFEF-26A1-E011-9232-0030487CD6E8.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/696/802096D3-26A1-E011-98E1-003048D2BB58.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/676/D81C5175-1BA1-E011-A708-003048D3756A.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/676/CE9A8CC9-21A1-E011-94E1-0019B9F704D6.root',
			'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/676/CC91642D-50A1-E011-BECF-003048D37514.root'
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

process.outpath = cms.EndPath(process.out
)


#
#  SUSY-PAT configuration file adapted for RA2 workflow
#
#  PAT configuration for the SUSY group - 41X series
#  More information here:
#  https://twiki.cern.ch/twiki/bin/view/CMS/SusyPatLayer1DefV10
#

# Starting with a skeleton process which gets imported with the following line
#from PhysicsTools.PatAlgos.patTemplate_cfg import *

import FWCore.ParameterSet.Config as cms

process = cms.Process("MySimplePATAnalysis")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True)
)


evts2Process = 100000
#runningOnMC = True 
runningOnMC = False

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "START42_V12::All"
if runningOnMC == False:
    process.GlobalTag.globaltag = "GR_R_42_V12::All"

#-- Message Logger ------------------------------------------------------------
process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
    limit = cms.untracked.int32(10),
    reportEvery = cms.untracked.int32(1)
    )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


#-- Input Source --------------------------------------------------------------
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(evts2Process))

#process.load('SusyAnalysis.MyAnalysis.HTRun2011AMay10ReReco_160404to163869_cfi')
#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
# "dcap://pnfs/cms/WAX/resilient/seema/SusyRA2Analysis/12June2011_HTRun2011A_May10ReReco_160404to163869_AOD_lpcA/susypat_1000_1_jHl.root"
#,"dcap://pnfs/cms/WAX/resilient/seema/SusyRA2Analysis/12June2011_HTRun2011A_May10ReReco_160404to163869_AOD_lpcA/susypat_1001_1_jZd.root"
#,"dcap://pnfs/cms/WAX/resilient/seema/SusyRA2Analysis/12June2011_HTRun2011A_May10ReReco_160404to163869_AOD_lpcA/susypat_1002_1_yUP.root"
#,"dcap://pnfs/cms/WAX/resilient/seema/SusyRA2Analysis/12June2011_HTRun2011A_May10ReReco_160404to163869_AOD_lpcA/susypat_1003_1_GQj.root"
#,"dcap://pnfs/cms/WAX/resilient/seema/SusyRA2Analysis/12June2011_HTRun2011A_May10ReReco_160404to163869_AOD_lpcA/susypat_1004_1_94m.root"
#,"dcap://pnfs/cms/WAX/resilient/seema/SusyRA2Analysis/12June2011_HTRun2011A_May10ReReco_160404to163869_AOD_lpcA/susypat_1005_1_QpZ.root"
#    )
#)
if evts2Process != -1 : 
	print "Using short list of files for testing."
	process.load("UserCode.RA2QCDvetoAna.HTRun2011APromptRecoV4_160404to167151_shortlist_cfi")
	#process.load("UserCode.RA2QCDvetoAna.2011RA2_Jun09_forrestp_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_cfi")
else : process.load("UserCode.RA2QCDvetoAna.HTRun2011APromptRecoV4_160404to167151_cfi")
#-- check RA2 recipe here ------------------------------------------------------------
process.prefilterCounter        = cms.EDProducer("EventCountProducer")
process.postStdCleaningCounter  = cms.EDProducer("EventCountProducer")

#-- Output module configuration -----------------------------------------------
process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')

#process.load('SandBox.Skims.puWeightProducer_cfi')

#process.susyRA2Analysis = cms.EDAnalyzer("SusyRA2Analysis",
#                                         VertexSource = cms.InputTag("goodVerticesRA2"),
#                                         MHTSource = cms.InputTag("mhtPF"),
#                                         HTSource  = cms.InputTag("htPF"),
#                                         JetSource = cms.InputTag("patJetsAK5PFPt50Eta25"),
#                                         METSource = cms.InputTag("patMETsPF"),
#                                         DoPUReweight = cms.bool(False),
#                                         PUWeigthSource = cms.InputTag("puWeight"),
#                                         DoOptimizePlots = cms.bool(True),
#                                         Debug       = cms.bool(False)
#)

process.QCDvetoAna = cms.EDFilter('RA2QCDvetoAna',
							#caloJetInputTag_ = cms.InputTag("ak5CaloJets"),
							#pfJetInputTag_   = cms.InputTag("ak5PFJets"),
							#HltTriggerResults = cms.InputTag("TriggerResults::HLT"),
							#req_trigger = cms.untracked.bool(False),
							#TriggerPathsToStore = cms.vstring("HLT_MET65_v3"),
							patJetsPFPt30InputTag = cms.InputTag("patJetsAK5PFPt30"),
							patJetsPFPt50Eta25InputTag = cms.InputTag("patJetsAK5PFPt50Eta25"),
							mhtInputTag = cms.InputTag("mhtPF"),
							htInputTag = cms.InputTag("htPF"),
							nMinVtx = cms.untracked.int32(1),	
							ndofVtx = cms.untracked.int32(4),	
							maxDelzTrkVtx = cms.untracked.double(24.0),	
							minMet = cms.untracked.double(0.0),	
							verbose = cms.untracked.int32(0)
)


process.fullseq = cms.Sequence(
                      process.ra2StdCleaning *
                      process.ra2PostCleaning *
                      process.ra2EcalTPFilter *
                      process.vetoBEFilterHTRun2011AMay10ReRecov1 *
                      process.ra2FullPFSelection *
                      process.QCDvetoAna
)

process.fullseq.remove(process.jetMHTPFDPhiFilter)

process.ppf = cms.Path( process.fullseq )

##-- Output module configuration ---------------------------------------
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('ra2Histograms.root')
)

#process.out.SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('ppf') )
#process.out.fileName = cms.untracked.string('susypat.root')
##process.out.outputCommands = cms.untracked.vstring('keep *')
#process.outpath = cms.EndPath( process.out )
#
#del process.outpath

###-- Dump config ------------------------------------------------------------
##file = open('SusyPAT_RA2414_cfg.py','w')
##file.write(str(process.dumpPython()))
##file.close()


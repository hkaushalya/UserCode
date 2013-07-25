#
#  SUSY-PAT configuration file adapted for RA2 workflow
#
#  PAT configuration for the SUSY group - 41X series
#  More information here:
#  https://twiki.cern.ch/twiki/bin/view/CMS/SusyPatLayer1DefV10
#

# Starting with a skeleton process which gets imported with the following line
#from PhysicsTools.PatAlgos.patTemplate_cfg import *

outputFile="Default.root"
#evts = 5000
evts = -1
runningOnMC = True 
#runningOnMC = False

import FWCore.ParameterSet.Config as cms

process = cms.Process("ra2b")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True)
)

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
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(evts))
#load test input files
#process.load("UserCode.RA2QCDvetoAna.QCD_Pt_120_1800_cfi")
process.load("UserCode.RA2QCDvetoAna.LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_cfi")

#quit()
#-- check RA2 recipe here ------------------------------------------------------------
process.prefilterCounter        = cms.EDProducer("EventCountProducer")
process.postStdCleaningCounter  = cms.EDProducer("EventCountProducer")

#-- Output module configuration -----------------------------------------------
process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')
process.mhtPFFilter.MinMHT = 0

process.load('SandBox.Skims.electronSelector_cfi')
process.electronSelector.DoElectronVeto = True

process.QCDvetoAna = cms.EDFilter('RA2bQCDvetoAna',
							patJetsPFPt30InputTag = cms.InputTag("patJetsAK5PFPt30"),
							patJetsPFPt50Eta25InputTag = cms.InputTag("patJetsAK5PFPt50Eta25"),
							mhtInputTag = cms.InputTag("mhtPF"),
							htInputTag = cms.InputTag("htPF"),
							nMinVtx = cms.untracked.int32(1),	
							ndofVtx = cms.untracked.int32(4),	
							maxDelzVtx = cms.untracked.double(24.0),	
							maxDelRho = cms.untracked.double(15.0),
							minJetEt4MHt = cms.untracked.double(30.0),	
							maxJetEta4MHt = cms.untracked.double(5.0),	
							dMinHT = cms.untracked.double(350.0),	
							dMinMHT = cms.untracked.double(0.0),	
							verbose = cms.untracked.int32(0)
)

process.fullseq = cms.Sequence(
                      process.ra2StdCleaning *
                      process.ra2PostCleaning *
                      process.ra2EcalTPFilter *
                      process.vetoBEFilterHTRun2011AMay10ReRecov1 *
                      process.ra2FullPFSelection *
							 process.electronSelector *
                      process.QCDvetoAna
)

process.fullseq.remove(process.jetMHTPFDPhiFilter)

process.ppf = cms.Path( process.fullseq )

##-- Output module configuration ---------------------------------------
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outputFile)
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


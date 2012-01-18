runningOnMC = True 
evts = -1
outputFile="Data_1.root"

import FWCore.ParameterSet.Config as cms

process = cms.Process("factorization")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True)
)

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('standard')

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

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
#'/store/user/tarasc/2011RA2/Sept28/tarasc/QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6/QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11_ra2StdCleaned_Inclusive2Jets_PROD/49cc45d3507c9435f30a6d7aecca3838/ra2SUSYPAT_468_1_nZx.root'
#'file:dcap://dcache-cms-dcap.desy.de//pnfs/desy.de/cms/tier2/store/user/mschrode/HT/RA2PreSelectionOnData_Run2011A_05Aug2011-v1_V5/231eedad1f0346b0c473749caab737bd/RA2SkimsOnData_8_1_7vv.root'

#'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_HTRun2011A_PromptRecoV4_165088to167913/susypat_477_1_4Kb.root',
#'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_HTRun2011A_PromptRecoV4_165088to167913/susypat_9_1_1dH.root'
#'/store/user/lhx/2011RA2/Jul01/lhx/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/WJetsToLNu_TuneZ2-madgraph_Summer11_ra2StdCleaned_Inclusive2Jets_PROD/a4709de23b27e18cf9a5b5ad2ea85731/ra2SUSYPAT_96_1_5ny.root'
#,'/store/user/lhx/2011RA2/Jul01/lhx/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/WJetsToLNu_TuneZ2-madgraph_Summer11_ra2StdCleaned_Inclusive2Jets_PROD/a4709de23b27e18cf9a5b5ad2ea85731/ra2SUSYPAT_97_1_BDX.root'
)
)
#load test input files
#process.load("UserCode.RA2QCDvetoAna.QCD_Pt_120_1800_cfi")

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

process.dummyCounter = cms.EDProducer("EventCountProducer")

from UserCode.RA2QCDvetoAna.factorization_mcdebug_cfi import *
process.factorization_debug = factorization_mcdebug.clone()
	

process.preseq = cms.Sequence(
                      process.ra2StdCleaning *
                      process.ra2PostCleaning *
                      process.ra2EcalTPFilter *
                      process.ra2FullPFSelection *
							 process.electronSelector *
							 process.factorization_debug
)

process.preseq.remove(process.jetMHTPFDPhiFilter)

process.ppf = cms.Path( process.preseq )

##-- Output module configuration ---------------------------------------
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outputFile)
)

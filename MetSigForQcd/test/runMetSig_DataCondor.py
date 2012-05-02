evts2Process = -1
outputFile="Data_1.root"

import FWCore.ParameterSet.Config as cms
import sys

fileNameContainer = cms.untracked.vstring()
file = sys.argv[3]

infile = open (file, "r")
for line in infile.readlines():
       line = line.rstrip('\n') # equivalent to Perl's chomp
       fileNameContainer.append (line)

print fileNameContainer

dataset = sys.argv[6]
print 'runMetSig_MCcondor dataset = ' , dataset

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True)
)
#-- Message Logger ------------------------------------------------------------
process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
    limit = cms.untracked.int32(10),
    reportEvery = cms.untracked.int32(1)
    )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "GR_R_42_V12::All"


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(evts2Process) )

process.source = cms.Source("PoolSource", 
                            fileNames = fileNameContainer           
)


process.trigWgtProd = cms.EDProducer('TrigPrescaleWeightProducer',
	 debug = cms.untracked.int32(0),
	 triggerResults = cms.InputTag("TriggerResults::HLT"),
    HLTProcess = cms.string("HLT"),
	 hltPaths = cms.vstring('HLT_HT100_v*','HLT_HT150_v*','HLT_HT160_v*','HLT_HT200_v*','HLT_HT240_v*','HLT_HT250_v*','HLT_HT260_v*','HLT_HT300_v*','HLT_HT350_v*','HLT_HT400_v*','HLT_HT450_v*','HLT_HT500_v*','HLT_HT550_v*','HLT_HT600_v*','HLT_HT650_v*','HLT_HT700_v*','HLT_HT750_v*','HLT_HT800_v*','HLT_HT850_v*','HLT_HT2000_v*','HLT_HT*_L1FastJet_v*')
)

process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')
process.mhtPFFilter.MinMHT = 0

process.load('SandBox.Skims.electronSelector_cfi')
process.electronSelector.DoElectronVeto = True

#========== relevant snippet to run over Seema's data skim ======
process.load("UserCode.RA2QCDvetoAna.filterBoolean_cfi")
process.RA2_HBHENoiseFilterRA2                 = process.booleanFilter.clone()
process.RA2_HBHENoiseFilterRA2.ResultSource    = cms.InputTag("HBHENoiseFilterRA2","Result","PAT")
process.RA2_beamHaloFilter                     = process.booleanFilter.clone()
process.RA2_beamHaloFilter.ResultSource        = cms.InputTag("beamHaloFilter","Result","PAT")
process.RA2_eeNoiseFilter                      = process.booleanFilter.clone()
process.RA2_eeNoiseFilter.ResultSource         = cms.InputTag("eeNoiseFilter","Result","PAT")
process.RA2_trackingFailureFilter              = process.booleanFilter.clone()
process.RA2_trackingFailureFilter.ResultSource = cms.InputTag("trackingFailureFilter","Result","PAT")
process.RA2_inconsistentMuons                  = process.booleanFilter.clone()
process.RA2_inconsistentMuons.ResultSource     = cms.InputTag("inconsistentMuons","Result","PAT")
process.RA2_greedyMuons                        = process.booleanFilter.clone()
process.RA2_greedyMuons.ResultSource           = cms.InputTag("greedyMuons","Result","PAT")

process.postProcessingSeq = cms.Sequence(process.RA2_HBHENoiseFilterRA2 *
			process.RA2_beamHaloFilter *
			process.RA2_eeNoiseFilter *
			process.RA2_trackingFailureFilter *
			process.RA2_inconsistentMuons *
			process.RA2_greedyMuons )

print 'Loading PBNE filter, Processing DATA.'
process.load('SandBox.Skims.jetIDSelector_cfi')

from UserCode.MetSigForQcd.metsig_cfi import *

usePrescaleWeight = 0
applyLumiweighing = 0
applyEventweighing = 0

print 'usePrescaleWeight is set to  = ', usePrescaleWeight
print 'applyLumiweighing is set to  = ', applyLumiweighing
print 'applyEventweighing is set to = ', applyEventweighing

njet = 3
verbose = False
njetCut = 1
dPhiCut = 1
dMinMetSig = 0



process.metsignomht = metsig.clone(
                      dMinMHT = cms.untracked.double(0.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig),   
							 dMinNjet = cms.untracked.double(njet),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

process.metsigmht200 = metsig.clone(
                      dMinMHT = cms.untracked.double(200.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig),   
							 dMinNjet = cms.untracked.double(njet),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

process.metsigmht350 = metsig.clone(
                      dMinMHT = cms.untracked.double(350.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig),   
							 dMinNjet = cms.untracked.double(njet),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

process.metsigmht500 = metsig.clone(
                      dMinMHT = cms.untracked.double(500.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig),   
							 dMinNjet = cms.untracked.double(njet),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

process.preseq = cms.Sequence(
							 process.postProcessingSeq *	 #only for Seema's datasets
                      process.ra2EcalTPFilter *
                      process.ra2FullPFSelection *
							 process.electronSelector *
							 process.trigWgtProd * 
							 process.metsignomht +
							 process.metsigmht200 +
							 process.metsigmht350 +
							 process.metsigmht500
)


process.p = cms.Path(process.preseq)

##-- Output module configuration ---------------------------------------
process.TFileService = cms.Service("TFileService",
  			fileName = cms.string(sys.argv[4]+'_'+sys.argv[5]+'.root')
)

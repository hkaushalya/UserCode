evts2Process = -1
outputFile="Data_1.root"

import FWCore.ParameterSet.Config as cms

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
process.GlobalTag.globaltag = "START42_V12::All"


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(evts2Process) )

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(

	)
)

process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')
process.mhtPFFilter.MinMHT = 0
from SandBox.Skims.RA2Leptons_cff import *

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
dMinMetSig = 300

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
                      process.ra2StdCleaning *
                      process.ra2PostCleaning *
                      process.ra2EcalTPFilter *
#                      process.ra2FullPFSelection *
							 ra2PFMuonVeto *
							 ra2PFElectronVeto *
							 process.metsignomht +
							 process.metsigmht200 +
							 process.metsigmht350 +
							 process.metsigmht500
)


process.p = cms.Path(process.preseq)

##-- Output module configuration ---------------------------------------
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outputFile)
)

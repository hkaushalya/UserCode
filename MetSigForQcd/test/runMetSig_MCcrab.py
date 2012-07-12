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


process.load("UserCode.MetSigForQcd.metsigall_cfi")

usePrescaleWeight = 0
applyLumiweighing = 0
applyEventweighing = 0

print 'runMetSig_MCcrab: usePrescaleWeight is set to  = ', usePrescaleWeight
print 'runMetSig_MCcrab: applyLumiweighing is set to  = ', applyLumiweighing
print 'runMetSig_MCcrab: applyEventweighing is set to = ', applyEventweighing

njetmin = 6
njetmax = 1000
verbose = False
njetCut = 1
dPhiCut = 1

print 'runMetSig_MCcrab:  njetmin = ', njetmin
print 'runMetSig_MCcrab:  njetmax = ', njetmax
print 'runMetSig_MCcrab:  njetcut = ', njetCut
print 'runMetSig_MCcrab:  dPhiCut = ', dPhiCut


process.metsignomht.dMinNjet = njetmin
process.metsigmht25.dMinNjet = njetmin
process.metsigmht50.dMinNjet = njetmin
process.metsigmht75.dMinNjet = njetmin
process.metsigmht100.dMinNjet = njetmin

process.metsignomht.applyDphiCut = dPhiCut
process.metsigmht25.applyDphiCut = dPhiCut
process.metsigmht50.applyDphiCut = dPhiCut
process.metsigmht75.applyDphiCut = dPhiCut
process.metsigmht100.applyDphiCut = dPhiCut

process.metsig50nomht.dMinNjet = njetmin
process.metsig50mht25.dMinNjet = njetmin
process.metsig50mht50.dMinNjet = njetmin
process.metsig50mht75.dMinNjet = njetmin
process.metsig50mht100.dMinNjet = njetmin

process.metsig50nomht.applyDphiCut = dPhiCut
process.metsig50mht25.applyDphiCut = dPhiCut
process.metsig50mht50.applyDphiCut = dPhiCut
process.metsig50mht75.applyDphiCut = dPhiCut
process.metsig50mht100.applyDphiCut = dPhiCut

process.metsig75nomht.dMinNjet = njetmin
process.metsig75mht25.dMinNjet = njetmin
process.metsig75mht50.dMinNjet = njetmin
process.metsig75mht75.dMinNjet = njetmin
process.metsig75mht100.dMinNjet = njetmin

process.metsig75nomht.applyDphiCut = dPhiCut
process.metsig75mht25.applyDphiCut = dPhiCut
process.metsig75mht50.applyDphiCut = dPhiCut
process.metsig75mht75.applyDphiCut = dPhiCut
process.metsig75mht100.applyDphiCut = dPhiCut

process.metsig100nomht.dMinNjet = njetmin
process.metsig100mht25.dMinNjet = njetmin
process.metsig100mht50.dMinNjet = njetmin
process.metsig100mht75.dMinNjet = njetmin
process.metsig100mht100.dMinNjet = njetmin

process.metsig100nomht.applyDphiCut = dPhiCut
process.metsig100mht25.applyDphiCut = dPhiCut
process.metsig100mht50.applyDphiCut = dPhiCut
process.metsig100mht75.applyDphiCut = dPhiCut
process.metsig100mht100.applyDphiCut = dPhiCut



#process.metsignomht.dMinNjet = njetmin
#process.metsigmht100.dMinNjet = njetmin
#process.metsigmht150.dMinNjet = njetmin
#process.metsigmht175.dMinNjet = njetmin
#process.metsigmht200.dMinNjet = njetmin
#process.metsigmht350.dMinNjet = njetmin
#process.metsigmht500.dMinNjet = njetmin
#process.metsig100nomht.dMinNjet = njetmin
#process.metsig100mht100.dMinNjet = njetmin
#process.metsig100mht150.dMinNjet = njetmin
#process.metsig100mht175.dMinNjet = njetmin
#process.metsig100mht200.dMinNjet = njetmin
#process.metsig100mht350.dMinNjet = njetmin
#process.metsig100mht500.dMinNjet = njetmin
#process.metsig200nomht.dMinNjet = njetmin
#process.metsig200mht100.dMinNjet = njetmin
#process.metsig200mht150.dMinNjet = njetmin
#process.metsig200mht175.dMinNjet = njetmin
#process.metsig200mht200.dMinNjet = njetmin
#process.metsig200mht350.dMinNjet = njetmin
#process.metsig200mht500.dMinNjet = njetmin
#process.metsig300nomht.dMinNjet = njetmin
#process.metsig300mht100.dMinNjet = njetmin
#process.metsig300mht150.dMinNjet = njetmin
#process.metsig300mht175.dMinNjet = njetmin
#process.metsig300mht200.dMinNjet = njetmin
#process.metsig300mht350.dMinNjet = njetmin
#process.metsig300mht500.dMinNjet = njetmin
#
#process.metsignomht.dMaxNjet = njetmax
#process.metsigmht100.dMaxNjet = njetmax
#process.metsigmht150.dMaxNjet = njetmax
#process.metsigmht175.dMaxNjet = njetmax
#process.metsigmht200.dMaxNjet = njetmax
#process.metsigmht350.dMaxNjet = njetmax
#process.metsigmht500.dMaxNjet = njetmax
#process.metsig100nomht.dMaxNjet = njetmax
#process.metsig100mht100.dMaxNjet = njetmax
#process.metsig100mht150.dMaxNjet = njetmax
#process.metsig100mht175.dMaxNjet = njetmax
#process.metsig100mht200.dMaxNjet = njetmax
#process.metsig100mht350.dMaxNjet = njetmax
#process.metsig100mht500.dMaxNjet = njetmax
#process.metsig200nomht.dMaxNjet = njetmax
#process.metsig200mht100.dMaxNjet = njetmax
#process.metsig200mht150.dMaxNjet = njetmax
#process.metsig200mht175.dMaxNjet = njetmax
#process.metsig200mht200.dMaxNjet = njetmax
#process.metsig200mht350.dMaxNjet = njetmax
#process.metsig200mht500.dMaxNjet = njetmax
#process.metsig300nomht.dMaxNjet = njetmax
#process.metsig300mht100.dMaxNjet = njetmax
#process.metsig300mht150.dMaxNjet = njetmax
#process.metsig300mht175.dMaxNjet = njetmax
#process.metsig300mht200.dMaxNjet = njetmax
#process.metsig300mht350.dMaxNjet = njetmax
#process.metsig300mht500.dMaxNjet = njetmax
#
#
#process.metsignomht.applyDphiCut = dPhiCut
#process.metsigmht100.applyDphiCut = dPhiCut
#process.metsigmht150.applyDphiCut = dPhiCut
#process.metsigmht175.applyDphiCut = dPhiCut
#process.metsigmht200.applyDphiCut = dPhiCut
#process.metsigmht350.applyDphiCut = dPhiCut
#process.metsigmht500.applyDphiCut = dPhiCut
#process.metsig100nomht.applyDphiCut = dPhiCut
#process.metsig100mht100.applyDphiCut = dPhiCut
#process.metsig100mht150.applyDphiCut = dPhiCut
#process.metsig100mht175.applyDphiCut = dPhiCut
#process.metsig100mht200.applyDphiCut = dPhiCut
#process.metsig100mht350.applyDphiCut = dPhiCut
#process.metsig100mht500.applyDphiCut = dPhiCut
#process.metsig200nomht.applyDphiCut = dPhiCut
#process.metsig200mht100.applyDphiCut = dPhiCut
#process.metsig200mht150.applyDphiCut = dPhiCut
#process.metsig200mht175.applyDphiCut = dPhiCut
#process.metsig200mht200.applyDphiCut = dPhiCut
#process.metsig200mht350.applyDphiCut = dPhiCut
#process.metsig200mht500.applyDphiCut = dPhiCut
#process.metsig300nomht.applyDphiCut = dPhiCut
#process.metsig300mht100.applyDphiCut = dPhiCut
#process.metsig300mht150.applyDphiCut = dPhiCut
#process.metsig300mht175.applyDphiCut = dPhiCut
#process.metsig300mht200.applyDphiCut = dPhiCut
#process.metsig300mht350.applyDphiCut = dPhiCut
#process.metsig300mht500.applyDphiCut = dPhiCut

process.preseq = cms.Sequence(
                      process.ra2StdCleaning *
                      process.ra2PostCleaning *
                      process.ra2EcalTPFilter *
#                      process.ra2FullPFSelection *
							 ra2PFMuonVeto *
							 ra2PFElectronVeto *
#							 process.metsigseq 
							 process.metsignomht *
							 process.metsigmht25 *
							 process.metsigmht50 *
							 process.metsigmht75 *
							 process.metsigmht100 *
							 process.metsig25nomht *
							 process.metsig25mht25 *
							 process.metsig25mht50 *
							 process.metsig25mht75 *
							 process.metsig25mht100 *
							 process.metsig50nomht *
							 process.metsig50mht25 *
							 process.metsig50mht50 *
							 process.metsig50mht75 *
							 process.metsig50mht100 *
							 process.metsig75nomht *
							 process.metsig75mht25 *
							 process.metsig75mht50 *
							 process.metsig75mht75 *
							 process.metsig75mht100 
)


process.p = cms.Path(process.preseq)

##-- Output module configuration ---------------------------------------
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outputFile)
)

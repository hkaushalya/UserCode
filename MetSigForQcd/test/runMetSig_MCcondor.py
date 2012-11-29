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
process.GlobalTag.globaltag = "START42_V12::All"


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(evts2Process) )

process.source = cms.Source("PoolSource",
                            fileNames = fileNameContainer           
)


process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')
process.mhtPFFilter.MinMHT = 0


from SandBox.Skims.RA2Leptons_cff import *

process.load("UserCode.MetSigForQcd.metsigall_cfi")
process.load("UserCode.MetSigForQcd.signifMHTProducer_cfi")

process.mymhtPFforSgnf.JetCollection = cms.InputTag("patJetsAK5PF")
process.mymhtPFforSgnf.MinJetPt      = cms.double(15)



usePrescaleWeight = 0
applyLumiweighing = 0
applyEventweighing = 0

print 'runMetSig_MCcrab: usePrescaleWeight is set to  = ', usePrescaleWeight
print 'runMetSig_MCcrab: applyLumiweighing is set to  = ', applyLumiweighing
print 'runMetSig_MCcrab: applyEventweighing is set to = ', applyEventweighing


njetmin = 3
njetmax = 1000
verbose = False
njetCut = 1
dPhiCut = 0

process.metsignomht.dMinNjet = njetmin
process.metsigmht25.dMinNjet = njetmin
process.metsigmht50.dMinNjet = njetmin
process.metsignomht.applyDphiCut = dPhiCut
process.metsigmht25.applyDphiCut = dPhiCut
process.metsigmht50.applyDphiCut = dPhiCut

process.preseq = cms.Sequence(
                      #process.ra2StdCleaning *
                      #process.ra2PostCleaning *
                      #process.ra2EcalTPFilter *
                      #process.ra2FullPFSelection *
							 ra2PFMuonVeto *
							 ra2PFElectronVeto *
							 process.mymhtPFforSgnf *
							 process.metsignomht *
							 process.metsigmht25 *
							 process.metsigmht50 
)


process.p = cms.Path(process.preseq)

##-- Output module configuration ---------------------------------------
process.TFileService = cms.Service("TFileService",
  			fileName = cms.string(sys.argv[4]+'_'+sys.argv[5]+'.root')
)

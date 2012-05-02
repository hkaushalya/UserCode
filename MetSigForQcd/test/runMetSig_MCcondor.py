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


process.dummyCounter = cms.EDProducer("EventCountProducer")

if dataset == 0:
	print 'Loading BEFilter 0: Using FLAT QCD Sample'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pt_15to3000_TuneZ2_Flat_Summer11_cfi')
	process.beFilterSeq = cms.Sequence(
			process.vetoBEFilterQCDPt15to3000TuneZ2FlatSummer11
			)
elif dataset == 1:
	print 'Loading BEFilter 1'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_120to170_cfi')
	process.beFilterSeq = cms.Sequence(
			process.vetoBEFilterQCDPythia120To170
			)
elif dataset == 2:
	print 'Loading BEFilter 2'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_170to300_cfi')
	process.beFilterSeq = cms.Sequence(
			process.vetoBEFilterQCDPythia170To300
			)
elif dataset == 3:
	print 'Loading BEFilter 3'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_300to470_cfi')
	process.beFilterSeq = cms.Sequence(
			process.vetoBEFilterQCDPythia300To470
			)
elif dataset == 4:
	print 'Loading BEFilter 4'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_470to600_cfi')
	process.beFilterSeq = cms.Sequence(
			process.vetoBEFilterQCDPythia470To600
			)
elif dataset == 5:
	print 'Loading BEFilter 5'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_600to800_cfi')
	process.beFilterSeq = cms.Sequence(
			process.vetoBEFilterQCDPythia600To800
			)
elif dataset == 6:
	print 'Loading BEFilter 6'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_800to1000_cfi')
	process.beFilterSeq = cms.Sequence(
			process.vetoBEFilterQCDPythia800To1000
			)
elif dataset == 7:
	print 'Loading BEFilter 7'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_1000to1400_cfi')
	process.beFilterSeq = cms.Sequence(
			process.vetoBEFilterQCDPythia1000To1400
			)
elif dataset == 8:
	print 'Loading BEFilter 8'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_1400to1800_cfi')
	process.beFilterSeq = cms.Sequence(
			process.vetoBEFilterQCDPythia1400To1800
			)
elif dataset == 9:
	print 'Loading BEFilter 9'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_1800_cfi')
	process.beFilterSeq = cms.Sequence(
			process.vetoBEFilterQCDPythia1800
			)
elif dataset == 10:
	print 'Loading BEFilter 10 for MAD Graph 100to250'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_TuneD6T_HT_100To250_cfi')
	process.beFilterSeq = cms.Sequence(
			process.vetoBEFilterQCDTuneD6THT100To250
			)
elif dataset == 11:
	print 'Loading BEFilter 11 for MAD Graph 250to500'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_TuneD6T_HT_250To500_cfi')
	process.beFilterSeq = cms.Sequence(
			process.vetoBEFilterQCDTuneD6THT250To500
			)
elif dataset == 12:
	print 'Loading BEFilter 12 for MAD Graph 500to1000'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_TuneD6T_HT_500To1000_cfi')
	process.beFilterSeq = cms.Sequence(
			process.vetoBEFilterQCDTuneD6THT500To1000
			)
elif dataset == 13:
	print 'Loading BEFilter 13 for MAD Graph 100'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_TuneD6T_HT_1000_cfi')
	process.beFilterSeq = cms.Sequence(
			process.vetoBEFilterQCDTuneD6THT1000
			)
else:
	print 'No BE Filter for the given dataset', dataset
	process.beFilterSeq = cms.Sequence(
			process.dummyCounter
			)


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
                      #process.ra2StdCleaning *
                      #process.ra2PostCleaning *
                      process.ra2EcalTPFilter *
                      #process.ra2FullPFSelection *
							 ra2PFMuonVeto *
							 ra2PFElectronVeto *
							 process.beFilterSeq *
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

outputFile="Default.root"
#evts = 1000
evts = -1
runningOnMC = True 
#runningOnMC = False

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
print 'myRa2bAna_condor: dataset = ' , dataset


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
process.source = cms.Source("PoolSource",
                            fileNames = fileNameContainer           
                            #skipEvents = cms.untracked.uint32(10)
)

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
							dMinHT = cms.untracked.double(500.0),	
							dMinMHT = cms.untracked.double(0.0),	
							verbose = cms.untracked.int32(0)
)

process.dummyCounter = cms.EDProducer("EventCountProducer")

process.preseq = cms.Sequence( process.dummyCounter *
                      process.ra2StdCleaning *
                      process.ra2PostCleaning *
                      process.ra2EcalTPFilter *
                      process.vetoBEFilterHTRun2011AMay10ReRecov1 *
                      process.ra2FullPFSelection *
							 process.electronSelector *
                      process.QCDvetoAna
)

process.preseq.remove(process.jetMHTPFDPhiFilter)

if dataset == 0:
	print 'Loading BEFilter 0: Using FLAT QCD Sample'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pt_15to3000_TuneZ2_Flat_Summer11_cfi')
	process.preseq.replace(process.dummyCounter,process.vetoBEFilterQCDPt15to3000TuneZ2FlatSummer11)
elif dataset == 1:
	print 'Loading BEFilter 1'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_120to170_cfi')
	process.preseq.replace(process.dummyCounter,process.vetoBEFilterQCDPythia120To170)
elif dataset == 2:
	print 'Loading BEFilter 2'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_170to300_cfi')
	process.preseq.replace(process.dummyCounter,process.vetoBEFilterQCDPythia170To300)
elif dataset == 3:
	print 'Loading BEFilter 3'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_300to470_cfi')
	process.preseq.replace(process.dummyCounter,process.vetoBEFilterQCDPythia300To470)
elif dataset == 4:
	print 'Loading BEFilter 4'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_470to600_cfi')
	process.preseq.replace(process.dummyCounter,process.vetoBEFilterQCDPythia470To600)
elif dataset == 5:
	print 'Loading BEFilter 5'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_600to800_cfi')
	process.preseq.replace(process.dummyCounter,process.vetoBEFilterQCDPythia600To800)
elif dataset == 6:
	print 'Loading BEFilter 6'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_800to1000_cfi')
	process.preseq.replace(process.dummyCounter,process.vetoBEFilterQCDPythia800To1000)
elif dataset == 7:
	print 'Loading BEFilter 7'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_1000to1400_cfi')
	process.preseq.replace(process.dummyCounter,process.vetoBEFilterQCDPythia1000To1400)
elif dataset == 8:
	print 'Loading BEFilter 8'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_1400to1800_cfi')
	process.preseq.replace(process.dummyCounter,process.vetoBEFilterQCDPythia1400To1800)
elif dataset == 9:
	print 'Loading BEFilter 9'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_1800_cfi')
	process.preseq.replace(process.dummyCounter,process.vetoBEFilterQCDPythia1800)
else:
	print 'No BE Filter for the given dataset', dataset

##-- Output module configuration ---------------------------------------
process.TFileService = cms.Service("TFileService",
  #fileName = cms.string('file_j'+sys.argv[2]+'.root')
  fileName = cms.string(sys.argv[4]+'_'+sys.argv[5]+'.root')
)

process.ppf = cms.Path( process.preseq )

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


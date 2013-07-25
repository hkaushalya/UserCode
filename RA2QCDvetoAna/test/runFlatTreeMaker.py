runningOnMC = True 
evts = 1000
outputFile="FlatTree.root"

import FWCore.ParameterSet.Config as cms

process = cms.Process("treeMaker")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True)
)

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('standard')

options.register('dataset', True, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "qcd flat/binned dataset")
options.parseArguments()
options._tagOrder =[]
print options

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
'/store/user/tarasc/2011RA2/Sept28/tarasc/QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6/QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11_ra2StdCleaned_Inclusive2Jets_PROD/49cc45d3507c9435f30a6d7aecca3838/ra2SUSYPAT_468_1_nZx.root'
#'file:dcap://dcache-cms-dcap.desy.de//pnfs/desy.de/cms/tier2/store/user/mschrode/HT/RA2PreSelectionOnData_Run2011A_May10ReReco-v1_V5/231eedad1f0346b0c473749caab737bd/RA2SkimsOnData_358_1_pxG.root'
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

process.treeMaker = cms.EDFilter('FlatTreeMaker',
							jetSrc = cms.InputTag("patJetsAK5PFPt50Eta25"),
							mhtSrc = cms.InputTag("mhtPF"),
							metSrc = cms.InputTag("pfMet"),
							htSrc = cms.InputTag("htPF"),
							muonSrc = cms.InputTag("cleanPatMuons"),
							eleSrc = cms.InputTag("cleanPatElectrons"),
							forVetoMuonSrc     = cms.InputTag("patMuonsPFIDIso"),
							forVetoElectronSrc = cms.InputTag("patElectronsPFIDIso"),
							doFillTree         = cms.bool(True),
							outRootName        = cms.string("StdRA2.root"),
							genJetSrc = cms.InputTag("ak5GenJets"),
							genMETSrc = cms.InputTag("genMetTrue"),
							genParticleSrc = cms.InputTag("genParticles"),
							recoGenJetsDR = cms.double(0.5),
							vtxSrc = cms.InputTag("goodVerticesRA2"),
							nVtxPUcut = cms.uint32(1),
							debug = cms.bool(True),
							doSgnf = cms.bool(False),
							verbose = cms.untracked.int32(0)
)

process.preseq = cms.Sequence(
                      process.ra2StdCleaning *
                      process.ra2PostCleaning *
                      process.ra2EcalTPFilter *
                      #process.vetoBEFilterHTRun2011AMay10ReRecov1 *
							 #process.vetoBEFilterQCDPt15to3000TuneZ2FlatSummer11 *
                      process.ra2FullPFSelection *
							 process.electronSelector *
                      process.treeMaker
)

process.preseq.remove(process.jetMHTPFDPhiFilter)

if options.dataset == -1:
	print 'Loading PBNE filter, Processing DATA.'
	process.load('SandBox.Skims.jetIDSelector_cfi')
	process.fullseq = cms.Sequence(process.preseq
			)
elif options.dataset == 0:
	print 'Loading BEFilter 0: Using FLAT QCD Sample'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pt_15to3000_TuneZ2_Flat_Summer11_cfi')
	process.fullseq = cms.Sequence(
			process.vetoBEFilterQCDPt15to3000TuneZ2FlatSummer11
			* process.preseq
			)
elif options.dataset == 1:
	print 'Loading BEFilter 1'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_120to170_cfi')
	process.fullseq = cms.Sequence(
			process.vetoBEFilterQCDPythia120To170
			* process.preseq
			)
elif options.dataset == 2:
	print 'Loading BEFilter 2'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_170to300_cfi')
	process.fullseq = cms.Sequence(
			process.vetoBEFilterQCDPythia170To300
			* process.preseq
			)
elif options.dataset == 3:
	print 'Loading BEFilter 3'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_300to470_cfi')
	process.fullseq = cms.Sequence(
			process.vetoBEFilterQCDPythia300To470
			* process.preseq
			)
elif options.dataset == 4:
	print 'Loading BEFilter 4'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_470to600_cfi')
	process.fullseq = cms.Sequence(
			process.vetoBEFilterQCDPythia470To600
			* process.preseq
			)
elif options.dataset == 5:
	print 'Loading BEFilter 5'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_600to800_cfi')
	process.fullseq = cms.Sequence(
			process.vetoBEFilterQCDPythia600To800
			* process.preseq
			)
elif options.dataset == 6:
	print 'Loading BEFilter 6'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_800to1000_cfi')
	process.fullseq = cms.Sequence(
			process.vetoBEFilterQCDPythia800To1000
			* process.preseq
			)
elif options.dataset == 7:
	print 'Loading BEFilter 7'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_1000to1400_cfi')
	process.fullseq = cms.Sequence(
			process.vetoBEFilterQCDPythia1000To1400
			* process.preseq
			)
elif options.dataset == 8:
	print 'Loading BEFilter 8'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_1400to1800_cfi')
	process.fullseq = cms.Sequence(
			process.vetoBEFilterQCDPythia1400To1800
			* process.preseq
			)
elif options.dataset == 9:
	print 'Loading BEFilter 9'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_1800_cfi')
	process.fullseq = cms.Sequence(
			process.preseq *
			process.vetoBEFilterQCDPythia1800
			)
else:
	print 'No BE Filter for the given dataset', dataset



process.ppf = cms.Path( process.fullseq )

##-- Output module configuration ---------------------------------------
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outputFile)
)


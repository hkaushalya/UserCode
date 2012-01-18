runningOnMC = False
evts = -1
outputFile="SmearedJets.root"

import FWCore.ParameterSet.Config as cms


process = cms.Process("factorizationSmearedJets")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True)
)

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('standard')

options.register('dataset', -1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "Specify correct dataset to load appropriate filters.")
options.parseArguments()
options._tagOrder =[]
print options

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "START42_V12::All"

if runningOnMC == False:
    process.GlobalTag.globaltag = "GR_R_42_V12::All"
#    process.GlobalTag.globaltag = "GR_R_42_V20::All"

print '          |**********************************************|'
if runningOnMC == True:
	print '          | Running on MC with GT',  process.GlobalTag.globaltag
else:
	print '          | Running on DATA with GT',  process.GlobalTag.globaltag
print '          |**********************************************|'

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
#'file:dcap://dcache-cms-dcap.desy.de//pnfs/desy.de/cms/tier2/store/user/mschrode/HT/RA2PreSelectionOnData_Run2011A_05Aug2011-v1_V5/231eedad1f0346b0c473749caab737bd/RA2SkimsOnData_8_1_7vv.root'

#'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_HTRun2011A_PromptRecoV4_165088to167913/susypat_477_1_4Kb.root',
#'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_HTRun2011A_PromptRecoV4_165088to167913/susypat_9_1_1dH.root'
#'file:/uscms/home/samantha/SandBox/susypat.root'
#'/store/data/Run2011A/HT/AOD/May10ReReco-v1/0003/A0E0B56A-0A7C-E011-BB2B-003048678FD6.root'
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

import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
process.hltFilter1         = hlt.hltHighLevel.clone(
		#HLTPaths = cms.vstring('HLT_HT350_v*'), # provide list of HLT paths (or patterns) you want		
		HLTPaths = cms.vstring('HLT_HT100_v*','HLT_HT150_v*','HLT_HT160_v*','HLT_HT200_v*','HLT_HT240_v*','HLT_HT250_v*','HLT_HT260_v*','HLT_HT300_v*','HLT_HT350_v*','HLT_HT400_v*','HLT_HT450_v*','HLT_HT500_v*','HLT_HT550_v*','HLT_HT600_v*','HLT_HT650_v*','HLT_HT700_v*','HLT_HT750_v*','HLT_HT800_v*','HLT_HT850_v*','HLT_HT2000_v*','HLT_HT*_L1FastJet_v*'), # provide list of HLT paths (or patterns) you want		
		TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"), 
		andOr = cms.bool(True), # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
		throw = False # throw exception on unknown path names
		)

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
#===========================================================================================

###############################################################################
# Weight producer for prescaling effects
###############################################################################
process.load("RA2.WeightProducer.prescaleweightproducer_cfi")
#process.prescaleweightProducer.PrescaleCut = -1
###############################################################################

###############################################################################
process.load("RA2.QCDSmearProd.qcdsmearprod_cfi")
###############################################################################

###############################################################################
# Rebalancing and Smearing configuration
###############################################################################
process.QCDfromSmearing.SmearingFile    = '/uscms_data/d3/samantha/SUSYRA2work/CMSSW_4_2_5/src/RA2/QCDSmearProd/test/MCJetResolution_Summer11_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_fine_v2.root'
#process.QCDfromSmearing.jetCollection  = process.PatTagNames.jetTag
process.QCDfromSmearing.jetCollection   = 'patJetsAK5PFPt30'
process.QCDfromSmearing.genjetCollection= 'ak5GenJets'
process.QCDfromSmearing.uncertaintyName = ''
process.QCDfromSmearing.InputHisto1     = 'h_tot_JetAll_ResponsePt'
process.QCDfromSmearing.InputHisto2     = 'h_tot_JetAll_ResponsePt'
process.QCDfromSmearing.InputHisto3p    = 'h_tot_JetAll_ResponsePt'
process.QCDfromSmearing.NRebin          = 1
#process.QCDfromSmearing.SmearCollection = 'Reco'
process.QCDfromSmearing.SmearCollection = 'Gen'

### Here you can do systematic variations of the response templates
process.QCDfromSmearing.PtBinEdges_scaling  = cms.vdouble(0.,5000.)
process.QCDfromSmearing.EtaBinEdges_scaling = cms.vdouble(0.,5.)
process.QCDfromSmearing.LowerTailScaling    = cms.vdouble(1.0)
process.QCDfromSmearing.UpperTailScaling    = cms.vdouble(1.0)
process.QCDfromSmearing.AdditionalSmearing  = cms.vdouble(1.0)

process.QCDfromSmearing.absoluteTailScaling = False
process.QCDfromSmearing.SmearedJetPt        = 0.
process.QCDfromSmearing.RebalanceJetPt      = 10.
process.QCDfromSmearing.RebalanceMode       = 'MHThigh'
process.QCDfromSmearing.weightName          = 'prescaleweightProducer:weight'
process.QCDfromSmearing.ControlPlots        = False
process.QCDfromSmearing.Ntries              = 1
process.QCDfromSmearing.MHTcut_low          = cms.double(200.)
process.QCDfromSmearing.MHTcut_medium       = cms.double(350.)
process.QCDfromSmearing.MHTcut_high         = cms.double(500.)
process.QCDfromSmearing.HTcut_low           = cms.double(500.)
process.QCDfromSmearing.HTcut_medium        = cms.double(800.)
process.QCDfromSmearing.HTcut_high          = cms.double(1000.)
process.QCDfromSmearing.HTcut_veryhigh      = cms.double(1200.)
process.QCDfromSmearing.HTcut_extremehigh   = cms.double(1400.)
#process.QCDfromSmearing.probExtreme        = cms.double(3.e-5)
###############################################################################


doEventWeighing = False
if options.dataset == 0:
	doEventWeighing = True

print 'Event Weighing is set to doEventWeighing = ', doEventWeighing

from UserCode.RA2QCDvetoAna.factorizationSmearedJets_cfi import *

process.factorization_ht350 = factorizationSmearedJets.clone(
							dMinHT = cms.untracked.double(350.0),
							usePrescaleWeight = cms.untracked.int32(1) 
)
process.factorization_ht500 = factorizationSmearedJets.clone(
							dMinHT = cms.untracked.double(500.0)	
)
process.factorization_ht600 = factorizationSmearedJets.clone(
							dMinHT = cms.untracked.double(600.0)	
)
process.factorization_ht700 = factorizationSmearedJets.clone(
							dMinHT = cms.untracked.double(700.0)	
)
process.factorization_ht800 = factorizationSmearedJets.clone(
							dMinHT = cms.untracked.double(800.0)	
)
process.factorization_ht1000 = factorizationSmearedJets.clone(
							dMinHT = cms.untracked.double(1000.0)	
)
process.factorization_ht1200 = factorizationSmearedJets.clone(
							dMinHT = cms.untracked.double(1200.0)	
)
process.factorization_ht1400 = factorizationSmearedJets.clone(
							dMinHT = cms.untracked.double(1400.0)	
)


from UserCode.RA2QCDvetoAna.flatTree_cfi import *
process.treeMaker = treeMaker.clone(
					jetSrc = cms.InputTag("patJetsAK5PFPt30"),
					outRootName = cms.string("StdRA2FlatTree.root")
)

from UserCode.RA2QCDvetoAna.genPtHatFilter_cfi import *
process.ptHatFilter = genPtHatFilter.clone(
									dMinPtHat = cms.untracked.double(150)
)

print 'dMinPtHat = ', process.ptHatFilter.dMinPtHat

process.dumpcont = cms.EDAnalyzer("EventContentAnalyzer")
#process.pppp = cms.Path(process.dumpcont)
#process.Tracer = cms.Service("Tracer")



process.dummyCounter = cms.EDProducer("EventCountProducer")

process.preseq = cms.Sequence(
							 #process.ptHatFilter *
				          #process.hltFilter1 *
							 #process.postProcessingSeq *	 #only for Seema's datasets
                      #process.ra2StdCleaning *    #remove this when using Christian's Dataset. Already applied.
                      #process.ra2PostCleaning *    #remove this when using Christian's Dataset. Already applied.
                      process.ra2EcalTPFilter *
                      process.ra2FullPFSelection *
							 process.electronSelector *
			 	          #process.treeMaker *
							 process.prescaleweightProducer * #only for data
							 process.QCDfromSmearing *
							 #process.dumpcont *
                      process.factorization_ht350 *
                      process.factorization_ht500 *
                      process.factorization_ht600 *
                      process.factorization_ht700 *
                      process.factorization_ht800 *
                      process.factorization_ht1000 *
                      process.factorization_ht1200 *
							 process.factorization_ht1400
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
			process.vetoBEFilterQCDPythia1800
			* process.preseq
			)
elif options.dataset == 10:
	print 'Loading BEFilter 10 for MAD Graph 100to250'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_TuneD6T_HT_100To250_cfi')
	process.fullseq = cms.Sequence(
			process.vetoBEFilterQCDTuneD6THT100To250
			* process.preseq
			)
elif options.dataset == 11:
	print 'Loading BEFilter 11 for MAD Graph 250to500'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_TuneD6T_HT_250To500_cfi')
	process.fullseq = cms.Sequence(
			process.vetoBEFilterQCDTuneD6THT250To500
			* process.preseq
			)
elif options.dataset == 12:
	print 'Loading BEFilter 12 for MAD Graph 500to1000'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_TuneD6T_HT_500To1000_cfi')
	process.fullseq = cms.Sequence(
			process.vetoBEFilterQCDTuneD6THT500To1000
			* process.preseq
			)
elif options.dataset == 13:
	print 'Loading BEFilter 13 for MAD Graph 100'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_TuneD6T_HT_1000_cfi')
	process.fullseq = cms.Sequence(
			process.vetoBEFilterQCDTuneD6THT1000
			* process.preseq
			)
else:
	print 'No BE Filter for the given dataset', options.dataset
	process.fullseq = cms.Sequence(
			process.preseq
			)


process.ppf = cms.Path( process.fullseq )

##-- Output module configuration ---------------------------------------
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outputFile)
)

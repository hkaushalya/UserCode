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
print 'runFactorization: dataset = ' , dataset
process = cms.Process("factorization")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True)
)

#import FWCore.ParameterSet.VarParsing as VarParsing
#options = VarParsing.VarParsing ('standard')

#options.register('dataset', True, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "qcd flat/binned dataset")
#options.parseArguments()
#options._tagOrder =[]
#print options

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

import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
process.hltFilter1         = hlt.hltHighLevel.clone(
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

doEventWeighing = False
if dataset == 0:
	doEventWeighing = True

print 'Event Weighing is set to doEventWeighing = ', doEventWeighing

from UserCode.RA2QCDvetoAna.factorization_cfi import *

process.factorization_ht350 = factorization.clone(
							dMinHT = cms.untracked.double(350.0)	
)
process.factorization_ht500 = factorization.clone(
							dMinHT = cms.untracked.double(500.0)	
)
process.factorization_ht600 = factorization.clone(
							dMinHT = cms.untracked.double(600.0)	
)
process.factorization_ht700 = factorization.clone(
							dMinHT = cms.untracked.double(700.0)	
)
process.factorization_ht800 = factorization.clone(
							dMinHT = cms.untracked.double(800.0)	
)
process.factorization_ht1000 = factorization.clone(
							dMinHT = cms.untracked.double(1000.0)	
)
process.factorization_ht1200 = factorization.clone(
							dMinHT = cms.untracked.double(1200.0)	
)
process.factorization_ht1400 = factorization.clone(
							dMinHT = cms.untracked.double(1400.0)	
)


from UserCode.RA2QCDvetoAna.flatTree_cfi import *
process.treeMaker = treeMaker.clone(
					jetSrc = cms.InputTag("patJetsAK5PFPt30"),
					outRootName  = cms.string('StdRA2FlatTree'+sys.argv[4]+'_'+sys.argv[5]+'.root')
)

from UserCode.RA2QCDvetoAna.genPtHatFilter_cfi import *
process.ptHatFilter = genPtHatFilter.clone(
									dMinPtHat = cms.untracked.double(150)
)

print 'dMinPtHat = ', process.ptHatFilter.dMinPtHat


process.dummyCounter = cms.EDProducer("EventCountProducer")

process.preseq = cms.Sequence(
							 #process.ptHatFilter *
							 #process.hltFilter1 *
							 #process.postProcessingSeq * #only when running over Seema's data skim's & Anwar's MG pats
                      #process.ra2StdCleaning *    #remove this when using Christian's Dataset. Already applied.
                      #process.ra2PostCleaning *    #remove this when using Christian's Dataset. Already applied.
                      process.ra2EcalTPFilter *
                      process.ra2FullPFSelection *
							 process.electronSelector *
							 #process.treeMaker *
                      process.factorization_ht350 *
                      process.factorization_ht500 *
                      process.factorization_ht600 *
                      process.factorization_ht700 *
                      process.factorization_ht800 *
                      process.factorization_ht1000 *
                      process.factorization_ht1200 *
                      process.factorization_ht1400
                      #process.Factorization
)
process.preseq.remove(process.jetMHTPFDPhiFilter)
#process.preseq.remove(process.countJetsAK5PFPt50Eta25)

if dataset == -1:
	print 'Loading PBNE filter, Processing DATA.'
	process.load('SandBox.Skims.jetIDSelector_cfi')
	process.fullseq = cms.Sequence(process.preseq
			)
elif dataset == 0:
	print 'Loading BEFilter 0: Using FLAT QCD Sample'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pt_15to3000_TuneZ2_Flat_Summer11_cfi')
	process.fullseq = cms.Sequence(
			process.vetoBEFilterQCDPt15to3000TuneZ2FlatSummer11
			* process.preseq
			)
elif dataset == 1:
	print 'Loading BEFilter 1'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_120to170_cfi')
	process.fullseq = cms.Sequence(
			process.vetoBEFilterQCDPythia120To170
			* process.preseq
			)
elif dataset == 2:
	print 'Loading BEFilter 2'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_170to300_cfi')
	process.fullseq = cms.Sequence(
			process.vetoBEFilterQCDPythia170To300
			* process.preseq
			)
elif dataset == 3:
	print 'Loading BEFilter 3'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_300to470_cfi')
	process.fullseq = cms.Sequence(
			process.vetoBEFilterQCDPythia300To470
			* process.preseq
			)
elif dataset == 4:
	print 'Loading BEFilter 4'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_470to600_cfi')
	process.fullseq = cms.Sequence(
			process.vetoBEFilterQCDPythia470To600
			* process.preseq
			)
elif dataset == 5:
	print 'Loading BEFilter 5'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_600to800_cfi')
	process.fullseq = cms.Sequence(
			process.vetoBEFilterQCDPythia600To800
			* process.preseq
			)
elif dataset == 6:
	print 'Loading BEFilter 6'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_800to1000_cfi')
	process.fullseq = cms.Sequence(
			process.vetoBEFilterQCDPythia800To1000
			* process.preseq
			)
elif dataset == 7:
	print 'Loading BEFilter 7'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_1000to1400_cfi')
	process.fullseq = cms.Sequence(
			process.vetoBEFilterQCDPythia1000To1400
			* process.preseq
			)
elif dataset == 8:
	print 'Loading BEFilter 8'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_1400to1800_cfi')
	process.fullseq = cms.Sequence(
			process.vetoBEFilterQCDPythia1400To1800
			* process.preseq
			)
elif dataset == 9:
	print 'Loading BEFilter 9'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_1800_cfi')
	process.fullseq = cms.Sequence(
			process.vetoBEFilterQCDPythia1800
			* process.preseq
			)
elif dataset == 10:
	print 'Loading BEFilter 10 for MAD Graph 100to250'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_TuneD6T_HT_100To250_cfi')
	process.fullseq = cms.Sequence(
			process.vetoBEFilterQCDTuneD6THT100To250
			* process.preseq
			)
elif dataset == 11:
	print 'Loading BEFilter 11 for MAD Graph 250to500'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_TuneD6T_HT_250To500_cfi')
	process.fullseq = cms.Sequence(
			process.vetoBEFilterQCDTuneD6THT250To500
			* process.preseq
			)
elif dataset == 12:
	print 'Loading BEFilter 12 for MAD Graph 500to1000'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_TuneD6T_HT_500To1000_cfi')
	process.fullseq = cms.Sequence(
			process.vetoBEFilterQCDTuneD6THT500To1000
			* process.preseq
			)
elif dataset == 13:
	print 'Loading BEFilter 13 for MAD Graph 100'
	process.load('PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_TuneD6T_HT_1000_cfi')
	process.fullseq = cms.Sequence(
			process.vetoBEFilterQCDTuneD6THT1000
			* process.preseq
			)
else:
	print 'No BE Filter for the given dataset', dataset
	process.fullseq = cms.Sequence(
			process.preseq
			)


##-- Output module configuration ---------------------------------------
process.TFileService = cms.Service("TFileService",
  fileName = cms.string(sys.argv[4]+'_'+sys.argv[5]+'.root')
)

process.ppf = cms.Path( process.preseq )


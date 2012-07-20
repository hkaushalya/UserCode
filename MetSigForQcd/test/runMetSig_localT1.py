runningOnMC = True
outputFile="SUSYT1_5Jets.root"

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import sys
options = VarParsing.VarParsing ('standard')

options.register('dataset', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "dataset selection, default = 0")
options.register('njetmin', 3, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "Minimum Number of jets (default = 3)")
options.register('dphicut', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "Apply dphi cut? yes=1, no=0 (default)")
options.register('maxevts', -1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "events to process (default = -1)")
options.parseArguments()
options._tagOrder =[]
print options



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
if runningOnMC == False:
    process.GlobalTag.globaltag = "GR_R_42_V12::All"

print "Using GT ",  process.GlobalTag.globaltag

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxevts) )

if options.dataset==0:
	process.load("UserCode.MetSigForQcd.SMS_T1tttt_Mgluino_450to1200_mLSP_50to800_7TeV_Pythia6Z_PU_START42_V11_FSIM_v2_cfi");
	outputFile = cms.string("SUSY_T1_njetmin"+ `options.njetmin` + "_dphicut" + `options.dphicut` + ".root")
else: 
	sys.exit("Invalid dataset!")

print "Output file will be:", outputFile


process.trigWgtProd = cms.EDProducer('TrigPrescaleWeightProducer',
	 debug = cms.untracked.int32(0),
	 triggerResults = cms.InputTag("TriggerResults::HLT"),
    HLTProcess = cms.string("HLT"),
	 hltPaths = cms.vstring('HLT_HT100_v*','HLT_HT150_v*','HLT_HT160_v*','HLT_HT200_v*','HLT_HT240_v*','HLT_HT250_v*','HLT_HT260_v*','HLT_HT300_v*','HLT_HT350_v*','HLT_HT400_v*','HLT_HT450_v*','HLT_HT500_v*','HLT_HT550_v*','HLT_HT600_v*','HLT_HT650_v*','HLT_HT700_v*','HLT_HT750_v*','HLT_HT800_v*','HLT_HT850_v*','HLT_HT2000_v*','HLT_HT*_L1FastJet_v*')
)

process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')
process.mhtPFFilter.MinMHT = 0

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


from SandBox.Skims.RA2Leptons_cff import *
from SandBox.Skims.mhtProducer_cfi import *
# MHT using PF Jets
process.mhtPF = mht.clone()
process.mhtPF.JetCollection = cms.InputTag('patJetsAK5PFPt30')


process.load("UserCode.MetSigForQcd.metsigall_cfi")
process.load("UserCode.MetSigForQcd.signifMHTProducer_cfi")

usePrescaleWeight = 0
applyLumiweighing = 0
applyEventweighing = 0

print '[runMetSig_localT1] usePrescaleWeight is set to  = ', usePrescaleWeight
print '[runMetSig_localT1] applyLumiweighing is set to  = ', applyLumiweighing
print '[runMetSig_localT1] applyEventweighing is set to = ', applyEventweighing

#njetmax = 1000 #need to include in the settings 
verbose = False
njetCut = 1

#print 'runMetSig_MCcrab:  njetmax = ', njetmax
print '[runMetSig_localT1] runMetSig_MCcrab:  njetcut = ', njetCut

process.metsignomht.dMinNjet = options.njetmin
process.metsigmht25.dMinNjet = options.njetmin
process.metsigmht50.dMinNjet = options.njetmin
process.metsignomht.applyDphiCut = options.dphicut
process.metsigmht25.applyDphiCut = options.dphicut
process.metsigmht50.applyDphiCut = options.dphicut

process.preseq = cms.Sequence(
							 #process.postProcessingSeq *	 #only for Seema's datasets
                      #process.ra2StdCleaning * #exclude this when running on Seema's skims
                      #process.ra2PostCleaning * #exclude this when running on Seema's skims
                     #process.ra2EcalTPFilter *
#                     process.ra2FullPFSelection *
							 ra2PFMuonVeto *
							 ra2PFElectronVeto *
							 #process.mhtPF *
							 #process.trigWgtProd * 
							 #process.metsigseq 
							 process.mymhtPFforSgnf *
							 process.metsignomht *
							 process.metsigmht25 *
							 process.metsigmht50 
)


process.p = cms.Path(process.preseq)

##-- Output module configuration ---------------------------------------
process.TFileService = cms.Service("TFileService",
                                   #fileName = cms.string(outputFile)
                                   fileName = outputFile
)

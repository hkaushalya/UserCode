from PhysicsTools.PatAlgos.patTemplate_cfg import *

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

import os
import sys
import re

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('standard')

options.register('GlobalTag', "START53_V7G::All", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "GlobaTTag to use (otherwise default Pat GT is used)")
options.register('mcInfo', True, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "process MonteCarlo data, default is data")
options.register('jetCorrections', 'L2Relative', VarParsing.VarParsing.multiplicity.list, VarParsing.VarParsing.varType.string, "Level of jet corrections to use: Note the factors are read from DB via GlobalTag")
options.jetCorrections.append('L3Absolute')

options.register('hltName', 'HLT', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "HLT menu to use for trigger matching") 
options.register('mcVersion', '', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "'36X' for example. Used for specific MC fix")
options.register('jetTypes', 'AK5PF', VarParsing.VarParsing.multiplicity.list, VarParsing.VarParsing.varType.string, "Additional jet types that will be produced (AK5Calo and AK5PF, cross cleaned in PF2PAT, are included anyway)")
options.register('hltSelection', '*', VarParsing.VarParsing.multiplicity.list, VarParsing.VarParsing.varType.string, "hlTriggers (OR) used to filter events. for data: ''HLT_Mu9', 'HLT_IsoMu9', 'HLT_IsoMu13_v*''; for MC, HLT_Mu9")
options.register('addKeep', '', VarParsing.VarParsing.multiplicity.list, VarParsing.VarParsing.varType.string, "Additional keep and drop statements to trim the event content")

options.register('dataVersion', '', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "'36X' for example. Used for specific DATA fix")

options.register('debug', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "switch on/off debug mode")

options.register('test', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "switch on/off debug mode")

options.register('doPtHatWeighting', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "PtHat weighting for QCD flat samples, default is False")

options.register('fileslist', '', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "name of a file with input source list")

options.register('fastsim', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "fastsim sample or not, default is False")

options.register('usePhiCorrMET', True, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "use phi corrected MET or not, default is False")

options.parseArguments()
options._tagOrder =[]

print options
 
print 'USING GT : ', process.GlobalTag.globaltag


#-- Message Logger ------------------------------------------------------------
process.MessageLogger.cerr.FwkReport.reportEvery = 100
if options.debug:
   process.MessageLogger.cerr.FwkReport.reportEvery = 1

#-- Input Source --------------------------------------------------------------
inputfiles = cms.untracked.vstring()
if options.fileslist:
   if os.path.exists(options.fileslist) == False or os.path.isfile(options.fileslist) == False:
      print 'fileslist ',options.fileslist,' does not exist\n'
      sys.exit(5)
   else:
      ifile = open(options.fileslist, 'r')
      for line in ifile.readlines():
         inputfiles.append(line)

print "inputfiles : \n", inputfiles, "\n"

if options.files:
   process.source.fileNames = options.files
elif options.fileslist:
   process.source.fileNames = inputfiles
else:
   process.source.fileNames = [
      'file:susypat.root'
   ]

process.source.inputCommands = cms.untracked.vstring( "keep *", "drop *_MEtoEDMConverter_*_*" )
process.maxEvents.input = options.maxEvents

#-- Calibration tag -----------------------------------------------------------
if options.GlobalTag:
   process.GlobalTag.globaltag = options.GlobalTag

if options.mcInfo == False: options.jetCorrections.append('L2L3Residual')
options.jetCorrections.insert(0, 'L1FastJet')

print "jetCorrections: "
print options.jetCorrections

# Default hltFilter with path "*"
import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
if options.hltSelection:
   process.hltFilter = hlt.hltHighLevel.clone(
      TriggerResultsTag = cms.InputTag("TriggerResults","",options.hltName),
      HLTPaths = cms.vstring(options.hltSelection),
      throw = True, # Don't throw?!
      andOr = True
   )

# Load RA2 related configuration
process.load('SandBox.Skims.RA2Cleaning_cff')
process.load('SandBox.Skims.RA2Jets_cff')
process.load('SandBox.Skims.RA2HT_cff')
process.load('SandBox.Skims.RA2MHT_cff')
process.load('SandBox.Skims.RA2DeltaR_cff')

process.load('RA2CleaningFilterResults_cff')

##____________________________________________________________________________||
process.dummyCounter = cms.EDProducer("EventCountProducer")


# make default RA2Jet collections
#from SandBox.Skims.basicJetSelector_cfi import *
#from PhysicsTools.PatAlgos.selectionLayer1.jetCountFilter_cfi import *

process.load('RecoMET.METFilters.jetIDFailureFilter_cfi')
process.jetIDFailure.JetSource = cms.InputTag("stopJetsPFchsPt30")
process.jetIDFailure.MinJetPt  = cms.double(30.0)
process.jetIDFailure.MaxJetEta = cms.double(999.0)

process.ra2AllCleaning = cms.Sequence(
      process.cleaningOnFilterResults * process.jetIDFailure * process.dummyCounter
   )

if options.fastsim:
   process.ra2AllCleaning.remove(process.RA2_HBHENoiseFilterRA2)
   process.ra2AllCleaning.remove(process.RA2_beamHaloFilter)

process.load('SandBox.Stop.StopJets_cff')

# Define CHS correctors
from JetMETCorrections.Configuration.DefaultJEC_cff import *
process.ak5PFchsL1Fastjet  = ak5PFL1Fastjet.clone ( algorithm = cms.string('AK5PFchs') )
process.ak5PFchsL2Relative = ak5PFL2Relative.clone( algorithm = cms.string('AK5PFchs') )
process.ak5PFchsL3Absolute = ak5PFL3Absolute.clone( algorithm = cms.string('AK5PFchs') )
process.ak5PFchsResidual   = ak5PFResidual.clone  ( algorithm = cms.string('AK5PFchs') )

process.ak5PFchsL1FastL2L3         = ak5PFL1FastL2L3Residual.clone( correctors = cms.vstring('ak5PFchsL1Fastjet', 'ak5PFchsL2Relative', 'ak5PFchsL3Absolute') )
process.ak5PFchsL1FastL2L3Residual = ak5PFL1FastL2L3Residual.clone( correctors = cms.vstring('ak5PFchsL1Fastjet', 'ak5PFchsL2Relative', 'ak5PFchsL3Absolute', 'ak5PFchsResidual') )
process.ak5PFchsL2L3               = ak5PFL1FastL2L3Residual.clone( correctors = cms.vstring('ak5PFchsL2Relative', 'ak5PFchsL3Absolute') )


#for new JECs
process.load("SandBox.Skims.jesChange_cfi")

process.load("SandBox.Stop.trackIsolationMaker_cfi")

process.load("SandBox.Skims.RA2Leptons_cff")

process.load("SandBox.Stop.StopLeptons_cff")

process.load("SandBox.Stop.StopMETPhiCorr_cff")

if options.mcInfo == False:
   process.newMETwPhiCorr = process.newMETwPhiCorrData.clone()
else:
   process.newMETwPhiCorr = process.newMETwPhiCorrMC.clone()

# Other sequence
process.comb_seq = cms.Sequence(
# All cleaning && all basic variables, e.g., mht, ht...     
#   process.ppf *
# hlt requirement
   process.hltFilter *
   process.ra2PFMuons * process.ra2Electrons * process.stopPFMuons * process.stopElectrons * #process.trackIsolation * 
   process.metphi_seq * process.newMETwPhiCorr *
   process.stopPFJets
)

# Fix a bug if load PhysicsTools.PatAlgos.patTemplate_cfg, but no output path used
process.out.SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('*') )

process.load("SandBox.Stop.StopDPhiSelection_cff")
process.jetsMETDPhiFilter.jetSrc = cms.InputTag("stopJetsPFchsPt30")
if options.usePhiCorrMET == True: 
   process.jetsMETDPhiFilter.metSrc = cms.InputTag("newMETwPhiCorr")
else:
   process.jetsMETDPhiFilter.metSrc = cms.InputTag("patMETsPF")
process.jetsMETDPhiFilter.dPhiCuts = cms.untracked.vdouble(0.5, 0.5, 0.3)

process.load("SandBox.Stop.StopBTagJets_cff")
process.stopBJets.src = cms.InputTag("stopJetsPFchsPt30")

process.load("SandBox.Stop.StopType3TopTagger_cff")
if options.usePhiCorrMET == True: 
   process.type3topTagger.metSrc = cms.InputTag("newMETwPhiCorr")
else:
   process.type3topTagger.metSrc = cms.InputTag("patMETsPF")
process.type3topTagger.taggingMode = cms.untracked.bool(False)
process.type3topTagger.jetSrc = cms.InputTag("stopJetsPFchsPt30")

process.metPFchsFilter = process.mhtPFchsFilter.clone()
if options.usePhiCorrMET == True: 
   process.metPFchsFilter.MHTSource = cms.InputTag("newMETwPhiCorr")
else:
   process.metPFchsFilter.MHTSource = cms.InputTag("patMETsPF")

process.metPFchsFilter.MinMHT = cms.double(175)

process.met350PFchsFilter = process.mhtPFchsFilter.clone()
if options.usePhiCorrMET == True:
   process.met350PFchsFilter.MHTSource = cms.InputTag("newMETwPhiCorr")
else:
   process.met350PFchsFilter.MHTSource = cms.InputTag("patMETsPF")
process.met350PFchsFilter.MinMHT = cms.double(350)


from UserCode.LostLeptonTree.lostleptontree_cfi import *
process.treeMaker = demo.clone(
	#FileNamePUDataDistribution = cms.string("UserCode/LostLeptonTree/DataPileupHistogram_RA2Summer12_190456-208686_ABCD.root"),   
	FileNamePUDataDistribution = cms.string("/uscms_data/d3/samantha/SUSYRA2work/CMSSW_5_3_5/src/UserCode/LostLeptonTree/test/DataPileupHistogram_RA2Summer12_190456-208686_ABCD.root"),   
)

runningOnMC = options.mcInfo
outputFile="Ntuple.root"
process.treeMaker.MCflag = runningOnMC
process.treeMaker.DoPUReweight = runningOnMC
process.treeMaker.ApplyEventWeighing = False  #for flat sample 
#process.treeMaker.JetAllSource = cms.InputTag("patJetsPF")
#process.treeMaker.JetAllSource = cms.InputTag("stopJetsPFchsPt30")
process.treeMaker.JetAllSource = cms.InputTag("stopJetsPFchsPt0")
process.treeMaker.PFMetSource = cms.InputTag("newMETwPhiCorr")


#process.StdStop_RA2LeptonVeto_path = cms.Path(
#                                   process.comb_seq * process.ra2AllCleaning *
#                                   process.ra2MuonVeto * process.ra2ElectronVeto *
#                                   process.count2JetsPFchsPt70Eta24Std *
#                                   process.count4JetsPFchsPt50Eta24Std *
#                                   process.count5JetsPFchsPt30Eta24Std *
#                                   process.count5JetsPFchsPt30Std *
#                                   process.jetsMETDPhiFilter *
#                                   process.stopBJets * process.stopCountBJets *
#                                   process.type3topTagger #*
#                                   process.metPFchsFilter #*
#                                   process.met350PFchsFilter
#)

process.StdStop_DirectionalLeptonVeto_path = cms.Path(
                                   process.comb_seq * process.ra2AllCleaning * 
                                   #process.count2JetsPFchsPt70Eta24Std *
                                   #process.count4JetsPFchsPt50Eta24Std *
                                   #process.count5JetsPFchsPt30Eta24Std *
                                   #process.count5JetsPFchsPt30Std *
                                   process.stopPFMuonVeto * process.stopElectronVeto *
                                   #process.jetsMETDPhiFilter *
                                   #process.stopBJets * process.stopCountBJets *
                                   #process.type3topTagger ##*
                                   #process.metPFchsFilter #*
#                                   process.met350PFchsFilter
												process.newpatJetsPFchsPt0
											*	process.treeMaker 
)

process.out.fileName = cms.untracked.string('filterTags.root')

##____________________________________________________________________________||
process.hltTriggerSummaryAOD = cms.EDProducer(
    "TriggerSummaryProducerAOD",
    processName = cms.string( "@" )
    )
   
##____________________________________________________________________________||
from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger(process)

##____________________________________________________________________________||
process.patTriggerFilter = process.patTrigger.clone() 
process.patTriggerFilter.processName = cms.string('SUSY')
process.MessageLogger.suppressWarning += ['patTriggerFilter']

process.outpath = cms.EndPath(
    process.patTrigger *
    process.hltTriggerSummaryAOD *
    process.patTriggerFilter *
    process.out
)

##____________________________________________________________________________||
from PhysicsTools.PatAlgos.patEventContent_cff import *
process.out.outputCommands = cms.untracked.vstring(
    'drop *',
    'keep patTriggerPaths_patTriggerFilter_*_*',
)   
process.out.SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('*', '!*'))

##-- Output module configuration ---------------------------------------
process.TFileService = cms.Service("TFileService",
					 fileName = cms.string(outputFile)
)

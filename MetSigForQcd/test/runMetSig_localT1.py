runningOnMC = True
evts2Process = -1
outputFile="SUSYT1_6JetsInclNoDphi.root"

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
if runningOnMC == False:
    process.GlobalTag.globaltag = "GR_R_42_V12::All"

print "Using GT ",  process.GlobalTag.globaltag


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(evts2Process) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
		fileNames = cms.untracked.vstring(

	##########T1->ttttt sample
	'file:/uscms_data/d3/samantha/SUSYRA2work/CMSSW_4_2_5/src/SandBox/Skims/T1pats/susy_pat1.root'
	,'file:/uscms_data/d3/samantha/SUSYRA2work/CMSSW_4_2_5/src/SandBox/Skims/T1pats/susy_pat2.root'
	,'file:/uscms_data/d3/samantha/SUSYRA2work/CMSSW_4_2_5/src/SandBox/Skims/T1pats/susy_pat3.root'
	,'file:/uscms_data/d3/samantha/SUSYRA2work/CMSSW_4_2_5/src/SandBox/Skims/T1pats/susy_pat4.root'
	,'file:/uscms_data/d3/samantha/SUSYRA2work/CMSSW_4_2_5/src/SandBox/Skims/T1pats/susy_pat5.root'

		)
)

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

usePrescaleWeight = 0
applyLumiweighing = 0
applyEventweighing = 0

print 'runMetSig_localLM5: usePrescaleWeight is set to  = ', usePrescaleWeight
print 'runMetSig_localLM5: applyLumiweighing is set to  = ', applyLumiweighing
print 'runMetSig_localLM5: applyEventweighing is set to = ', applyEventweighing

njetmin = 6
#njetmax = 1000 #need to include in the settings 
verbose = False
njetCut = 1
dPhiCut = 0

print 'runMetSig_MCcrab:  njetmin = ', njetmin
#print 'runMetSig_MCcrab:  njetmax = ', njetmax
print 'runMetSig_MCcrab:  njetcut = ', njetCut
print 'runMetSig_MCcrab:  dPhiCut = ', dPhiCut

process.metsignomht.dMinNjet = njetmin
process.metsigmht25.dMinNjet = njetmin
process.metsigmht50.dMinNjet = njetmin
process.metsigmht75.dMinNjet = njetmin
process.metsigmht100.dMinNjet = njetmin
process.metsigmht125.dMinNjet = njetmin
process.metsigmht150.dMinNjet = njetmin
process.metsigmht175.dMinNjet = njetmin
process.metsigmht200.dMinNjet = njetmin
process.metsigmht350.dMinNjet = njetmin
process.metsigmht500.dMinNjet = njetmin

process.metsignomht.applyDphiCut = dPhiCut
process.metsigmht25.applyDphiCut = dPhiCut
process.metsigmht50.applyDphiCut = dPhiCut
process.metsigmht75.applyDphiCut = dPhiCut
process.metsigmht100.applyDphiCut = dPhiCut
process.metsigmht125.applyDphiCut = dPhiCut
process.metsigmht150.applyDphiCut = dPhiCut
process.metsigmht175.applyDphiCut = dPhiCut
process.metsigmht200.applyDphiCut = dPhiCut
process.metsigmht350.applyDphiCut = dPhiCut
process.metsigmht500.applyDphiCut = dPhiCut


process.metsig25nomht.dMinNjet = njetmin
process.metsig25mht25.dMinNjet = njetmin
process.metsig25mht50.dMinNjet = njetmin
process.metsig25mht75.dMinNjet = njetmin
process.metsig25mht100.dMinNjet = njetmin
process.metsig25mht125.dMinNjet = njetmin
process.metsig25mht150.dMinNjet = njetmin
process.metsig25mht175.dMinNjet = njetmin
process.metsig25mht200.dMinNjet = njetmin
process.metsig25mht350.dMinNjet = njetmin
process.metsig25mht500.dMinNjet = njetmin

process.metsig25nomht.applyDphiCut = dPhiCut
process.metsig25mht25.applyDphiCut = dPhiCut
process.metsig25mht50.applyDphiCut = dPhiCut
process.metsig25mht75.applyDphiCut = dPhiCut
process.metsig25mht100.applyDphiCut = dPhiCut
process.metsig25mht125.applyDphiCut = dPhiCut
process.metsig25mht150.applyDphiCut = dPhiCut
process.metsig25mht175.applyDphiCut = dPhiCut
process.metsig25mht200.applyDphiCut = dPhiCut
process.metsig25mht350.applyDphiCut = dPhiCut
process.metsig25mht500.applyDphiCut = dPhiCut


process.metsig50nomht.dMinNjet = njetmin
process.metsig50mht25.dMinNjet = njetmin
process.metsig50mht50.dMinNjet = njetmin
process.metsig50mht75.dMinNjet = njetmin
process.metsig50mht100.dMinNjet = njetmin
process.metsig50mht125.dMinNjet = njetmin
process.metsig50mht150.dMinNjet = njetmin
process.metsig50mht175.dMinNjet = njetmin
process.metsig50mht200.dMinNjet = njetmin
process.metsig50mht350.dMinNjet = njetmin
process.metsig50mht500.dMinNjet = njetmin

process.metsig50nomht.applyDphiCut = dPhiCut
process.metsig50mht25.applyDphiCut = dPhiCut
process.metsig50mht50.applyDphiCut = dPhiCut
process.metsig50mht75.applyDphiCut = dPhiCut
process.metsig50mht100.applyDphiCut = dPhiCut
process.metsig50mht125.applyDphiCut = dPhiCut
process.metsig50mht150.applyDphiCut = dPhiCut
process.metsig50mht175.applyDphiCut = dPhiCut
process.metsig50mht200.applyDphiCut = dPhiCut
process.metsig50mht350.applyDphiCut = dPhiCut
process.metsig50mht500.applyDphiCut = dPhiCut


process.metsig75nomht.dMinNjet = njetmin
process.metsig75mht25.dMinNjet = njetmin
process.metsig75mht50.dMinNjet = njetmin
process.metsig75mht75.dMinNjet = njetmin
process.metsig75mht100.dMinNjet = njetmin
process.metsig75mht125.dMinNjet = njetmin
process.metsig75mht150.dMinNjet = njetmin
process.metsig75mht175.dMinNjet = njetmin
process.metsig75mht200.dMinNjet = njetmin
process.metsig75mht350.dMinNjet = njetmin
process.metsig75mht500.dMinNjet = njetmin

process.metsig75nomht.applyDphiCut = dPhiCut
process.metsig75mht25.applyDphiCut = dPhiCut
process.metsig75mht50.applyDphiCut = dPhiCut
process.metsig75mht75.applyDphiCut = dPhiCut
process.metsig75mht100.applyDphiCut = dPhiCut
process.metsig75mht125.applyDphiCut = dPhiCut
process.metsig75mht150.applyDphiCut = dPhiCut
process.metsig75mht175.applyDphiCut = dPhiCut
process.metsig75mht200.applyDphiCut = dPhiCut
process.metsig75mht350.applyDphiCut = dPhiCut
process.metsig75mht500.applyDphiCut = dPhiCut


process.metsig100nomht.dMinNjet = njetmin
process.metsig100mht25.dMinNjet = njetmin
process.metsig100mht50.dMinNjet = njetmin
process.metsig100mht75.dMinNjet = njetmin
process.metsig100mht100.dMinNjet = njetmin
process.metsig100mht125.dMinNjet = njetmin
process.metsig100mht150.dMinNjet = njetmin
process.metsig100mht175.dMinNjet = njetmin
process.metsig100mht200.dMinNjet = njetmin
process.metsig100mht350.dMinNjet = njetmin
process.metsig100mht500.dMinNjet = njetmin

process.metsig100nomht.applyDphiCut = dPhiCut
process.metsig100mht25.applyDphiCut = dPhiCut
process.metsig100mht50.applyDphiCut = dPhiCut
process.metsig100mht75.applyDphiCut = dPhiCut
process.metsig100mht100.applyDphiCut = dPhiCut
process.metsig100mht125.applyDphiCut = dPhiCut
process.metsig100mht150.applyDphiCut = dPhiCut
process.metsig100mht175.applyDphiCut = dPhiCut
process.metsig100mht200.applyDphiCut = dPhiCut
process.metsig100mht350.applyDphiCut = dPhiCut
process.metsig100mht500.applyDphiCut = dPhiCut


process.metsig200nomht.dMinNjet = njetmin
process.metsig200mht25.dMinNjet = njetmin
process.metsig200mht50.dMinNjet = njetmin
process.metsig200mht75.dMinNjet = njetmin
process.metsig200mht100.dMinNjet = njetmin
process.metsig200mht125.dMinNjet = njetmin
process.metsig200mht150.dMinNjet = njetmin
process.metsig200mht175.dMinNjet = njetmin
process.metsig200mht200.dMinNjet = njetmin
process.metsig200mht350.dMinNjet = njetmin
process.metsig200mht500.dMinNjet = njetmin

process.metsig200nomht.applyDphiCut = dPhiCut
process.metsig200mht25.applyDphiCut = dPhiCut
process.metsig200mht50.applyDphiCut = dPhiCut
process.metsig200mht75.applyDphiCut = dPhiCut
process.metsig200mht100.applyDphiCut = dPhiCut
process.metsig200mht125.applyDphiCut = dPhiCut
process.metsig200mht150.applyDphiCut = dPhiCut
process.metsig200mht175.applyDphiCut = dPhiCut
process.metsig200mht200.applyDphiCut = dPhiCut
process.metsig200mht350.applyDphiCut = dPhiCut
process.metsig200mht500.applyDphiCut = dPhiCut


process.metsig300nomht.dMinNjet = njetmin
process.metsig300mht25.dMinNjet = njetmin
process.metsig300mht50.dMinNjet = njetmin
process.metsig300mht75.dMinNjet = njetmin
process.metsig300mht100.dMinNjet = njetmin
process.metsig300mht125.dMinNjet = njetmin
process.metsig300mht150.dMinNjet = njetmin
process.metsig300mht175.dMinNjet = njetmin
process.metsig300mht200.dMinNjet = njetmin
process.metsig300mht350.dMinNjet = njetmin
process.metsig300mht500.dMinNjet = njetmin

process.metsig300nomht.applyDphiCut = dPhiCut
process.metsig300mht25.applyDphiCut = dPhiCut
process.metsig300mht50.applyDphiCut = dPhiCut
process.metsig300mht75.applyDphiCut = dPhiCut
process.metsig300mht100.applyDphiCut = dPhiCut
process.metsig300mht125.applyDphiCut = dPhiCut
process.metsig300mht150.applyDphiCut = dPhiCut
process.metsig300mht175.applyDphiCut = dPhiCut
process.metsig300mht200.applyDphiCut = dPhiCut
process.metsig300mht350.applyDphiCut = dPhiCut
process.metsig300mht500.applyDphiCut = dPhiCut




process.preseq = cms.Sequence(
							 #process.postProcessingSeq *	 #only for Seema's datasets
                      process.ra2StdCleaning * #exclude this when running on Seema's skims
                      #process.ra2PostCleaning * #exclude this when running on Seema's skims
                     #process.ra2EcalTPFilter *
#                     process.ra2FullPFSelection *
							 ra2PFMuonVeto *
							 ra2PFElectronVeto *
							 #process.mhtPF *
							 #process.trigWgtProd * 
							 process.metsigseq 
)


process.p = cms.Path(process.preseq)

##-- Output module configuration ---------------------------------------
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outputFile)
)

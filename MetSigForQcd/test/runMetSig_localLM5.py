runningOnMC = True
evts2Process = -1
#outputFile="Data_1.root"
outputFile="SUSYLM5_2JetsExclWithDphi.root"

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


##########Seema's LM5 sample
'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM5SUSYsftshtpythia6_Summer11/susypat_10_1_SBb.root',
'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM5SUSYsftshtpythia6_Summer11/susypat_11_1_Hkw.root',
'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM5SUSYsftshtpythia6_Summer11/susypat_1_1_gzT.root',
'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM5SUSYsftshtpythia6_Summer11/susypat_12_1_Stc.root',
'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM5SUSYsftshtpythia6_Summer11/susypat_13_1_se3.root',
'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM5SUSYsftshtpythia6_Summer11/susypat_14_1_Xzg.root',
'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM5SUSYsftshtpythia6_Summer11/susypat_15_1_ApF.root',
'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM5SUSYsftshtpythia6_Summer11/susypat_16_1_N1c.root',
'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM5SUSYsftshtpythia6_Summer11/susypat_17_1_pgW.root',
'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM5SUSYsftshtpythia6_Summer11/susypat_2_1_rGh.root',
'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM5SUSYsftshtpythia6_Summer11/susypat_3_1_SZr.root',
'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM5SUSYsftshtpythia6_Summer11/susypat_4_1_PZW.root',
'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM5SUSYsftshtpythia6_Summer11/susypat_5_1_B8G.root',
'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM5SUSYsftshtpythia6_Summer11/susypat_6_1_xKn.root',
'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM5SUSYsftshtpythia6_Summer11/susypat_7_1_zpK.root',
'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM5SUSYsftshtpythia6_Summer11/susypat_8_1_Aly.root',
'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM5SUSYsftshtpythia6_Summer11/susypat_9_1_rPo.root'

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

njetmin = 2
njetmax = 2
verbose = False
njetCut = 1
dPhiCut = 1

process.metsignomht.dMinNjet = njetmin
process.metsigmht100.dMinNjet = njetmin
process.metsigmht150.dMinNjet = njetmin
process.metsigmht175.dMinNjet = njetmin
process.metsigmht200.dMinNjet = njetmin
process.metsigmht350.dMinNjet = njetmin
process.metsigmht500.dMinNjet = njetmin
process.metsig100nomht.dMinNjet = njetmin
process.metsig100mht100.dMinNjet = njetmin
process.metsig100mht150.dMinNjet = njetmin
process.metsig100mht175.dMinNjet = njetmin
process.metsig100mht200.dMinNjet = njetmin
process.metsig100mht350.dMinNjet = njetmin
process.metsig100mht500.dMinNjet = njetmin
process.metsig200nomht.dMinNjet = njetmin
process.metsig200mht100.dMinNjet = njetmin
process.metsig200mht150.dMinNjet = njetmin
process.metsig200mht175.dMinNjet = njetmin
process.metsig200mht200.dMinNjet = njetmin
process.metsig200mht350.dMinNjet = njetmin
process.metsig200mht500.dMinNjet = njetmin
process.metsig300nomht.dMinNjet = njetmin
process.metsig300mht100.dMinNjet = njetmin
process.metsig300mht150.dMinNjet = njetmin
process.metsig300mht175.dMinNjet = njetmin
process.metsig300mht200.dMinNjet = njetmin
process.metsig300mht350.dMinNjet = njetmin
process.metsig300mht500.dMinNjet = njetmin

process.metsignomht.dMaxNjet = njetmax
process.metsigmht100.dMaxNjet = njetmax
process.metsigmht150.dMaxNjet = njetmax
process.metsigmht175.dMaxNjet = njetmax
process.metsigmht200.dMaxNjet = njetmax
process.metsigmht350.dMaxNjet = njetmax
process.metsigmht500.dMaxNjet = njetmax
process.metsig100nomht.dMaxNjet = njetmax
process.metsig100mht100.dMaxNjet = njetmax
process.metsig100mht150.dMaxNjet = njetmax
process.metsig100mht175.dMaxNjet = njetmax
process.metsig100mht200.dMaxNjet = njetmax
process.metsig100mht350.dMaxNjet = njetmax
process.metsig100mht500.dMaxNjet = njetmax
process.metsig200nomht.dMaxNjet = njetmax
process.metsig200mht100.dMaxNjet = njetmax
process.metsig200mht150.dMaxNjet = njetmax
process.metsig200mht175.dMaxNjet = njetmax
process.metsig200mht200.dMaxNjet = njetmax
process.metsig200mht350.dMaxNjet = njetmax
process.metsig200mht500.dMaxNjet = njetmax
process.metsig300nomht.dMaxNjet = njetmax
process.metsig300mht100.dMaxNjet = njetmax
process.metsig300mht150.dMaxNjet = njetmax
process.metsig300mht175.dMaxNjet = njetmax
process.metsig300mht200.dMaxNjet = njetmax
process.metsig300mht350.dMaxNjet = njetmax
process.metsig300mht500.dMaxNjet = njetmax


process.metsignomht.applyDphiCut = dPhiCut
process.metsigmht100.applyDphiCut = dPhiCut
process.metsigmht150.applyDphiCut = dPhiCut
process.metsigmht175.applyDphiCut = dPhiCut
process.metsigmht200.applyDphiCut = dPhiCut
process.metsigmht350.applyDphiCut = dPhiCut
process.metsigmht500.applyDphiCut = dPhiCut
process.metsig100nomht.applyDphiCut = dPhiCut
process.metsig100mht100.applyDphiCut = dPhiCut
process.metsig100mht150.applyDphiCut = dPhiCut
process.metsig100mht175.applyDphiCut = dPhiCut
process.metsig100mht200.applyDphiCut = dPhiCut
process.metsig100mht350.applyDphiCut = dPhiCut
process.metsig100mht500.applyDphiCut = dPhiCut
process.metsig200nomht.applyDphiCut = dPhiCut
process.metsig200mht100.applyDphiCut = dPhiCut
process.metsig200mht150.applyDphiCut = dPhiCut
process.metsig200mht175.applyDphiCut = dPhiCut
process.metsig200mht200.applyDphiCut = dPhiCut
process.metsig200mht350.applyDphiCut = dPhiCut
process.metsig200mht500.applyDphiCut = dPhiCut
process.metsig300nomht.applyDphiCut = dPhiCut
process.metsig300mht100.applyDphiCut = dPhiCut
process.metsig300mht150.applyDphiCut = dPhiCut
process.metsig300mht175.applyDphiCut = dPhiCut
process.metsig300mht200.applyDphiCut = dPhiCut
process.metsig300mht350.applyDphiCut = dPhiCut
process.metsig300mht500.applyDphiCut = dPhiCut
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
							 process.metsigseq 
)


process.p = cms.Path(process.preseq)

##-- Output module configuration ---------------------------------------
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outputFile)
)

runningOnMC = True
evts2Process = -1
outputFile="SUSYLMX_5JetsExclWithDphi.root"

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

#Anwar's new SUSY MC sample
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_10_1_WdZ.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_11_1_NAm.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_12_1_GLv.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_13_1_iPo.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_14_1_VJA.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_15_1_vbQ.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_16_1_NK0.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_17_1_IhE.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_18_1_Qdr.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_19_1_EFk.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_1_1_7Gk.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_20_1_a07.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_21_1_Mfl.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_22_1_e0A.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_23_1_oBS.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_24_1_bC4.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_25_1_xGz.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_26_1_JBb.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_27_1_0Et.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_28_1_MDe.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_29_1_6oy.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_2_1_ahD.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_30_1_Y7F.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_31_1_4Je.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_32_1_si9.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_33_1_Hzo.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_34_1_ZQZ.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_35_1_7c2.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_36_1_Vb9.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_37_1_wFv.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_38_1_7W6.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_39_1_K0c.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_3_1_UQR.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_40_2_sKu.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_4_1_FL3.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_5_1_9Jw.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_6_1_dre.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_7_1_txf.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_8_1_kTd.root',
'/store/user/lpcjm/samantha/SUSY_1500_260_bhatti/samantha/SUSY_1500_260-AODSIM/samantha_SUSY_1500_260_bhatti_2Jets_noLeptonVeto/ccb8981a4f7138e417de052ad96dfa62/susy_pat_9_1_3Le.root'

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

njetmin = 5
njetmax = 5
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

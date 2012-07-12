runningOnMC = True
evts2Process = -1
outputFile="Data_1.root"
#outputFile="SUSYLM9_2jetsDphi.root"
#outputFile="SUSYLM9_2jetsOnly_lumiWgtd.root"
#outputFile="SUSYLM9_2jetsOnly.root"
#outputFile="SUSY_anwar_3jetsDphiMetSig300.root"

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
#		 '/store/user/lhx/2011RA2/Jul01/lhx/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/WJetsToLNu_TuneZ2-madgraph_Summer11_ra2StdCleaned_Inclusive2Jets_PROD/a4709de23b27e18cf9a5b5ad2ea85731/ra2SUSYPAT_96_1_5ny.root'
#		 ,'/store/user/lhx/2011RA2/Jul01/lhx/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/WJetsToLNu_TuneZ2-madgraph_Summer11_ra2StdCleaned_Inclusive2Jets_PROD/a4709de23b27e18cf9a5b5ad2ea85731/ra2SUSYPAT_97_1_BDX.root'
#		'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_HTRun2011B_PromptRecoV1_177516to178677/susypat_42_1_4Bd.root',
#		'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_HTRun2011A_PromptRecoV4_165088to167913/susypat_100_1_lvB.root',
#		'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_HTRun2011A_PromptRecoV4_165088to167913/susypat_101_1_kON.root'

#'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_10_1_TXG.root',
#'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_11_1_x1a.root',
#'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_1_1_QO3.root',
#'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_12_1_pfP.root',
#'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_13_1_7OQ.root',
#'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_14_1_wDF.root',
#'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_15_1_Wmj.root',
#'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_16_1_sE4.root',
#'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_17_1_Pqr.root',
#'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_18_1_5Za.root',
#'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_2_1_Kz9.root',
#'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_3_1_ij9.root',
#'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_4_1_1Cj.root',
#'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_5_1_MNO.root',
#'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_6_1_Eq3.root',
#'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_7_1_uE4.root',
#'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_8_1_uWQ.root',
#'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_9_1_qXc.root'

##########Seema's LM9 sample
#'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM9SUSYsftshtpythia6_Summer11/susypat_10_1_tFc.root',
#'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM9SUSYsftshtpythia6_Summer11/susypat_11_1_979.root',
#'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM9SUSYsftshtpythia6_Summer11/susypat_1_1_ptU.root',
#'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM9SUSYsftshtpythia6_Summer11/susypat_21_1_uuC.root',
#'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM9SUSYsftshtpythia6_Summer11/susypat_22_1_Hgh.root',
#'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM9SUSYsftshtpythia6_Summer11/susypat_23_1_AdH.root',
#'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM9SUSYsftshtpythia6_Summer11/susypat_24_1_myX.root',
#'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM9SUSYsftshtpythia6_Summer11/susypat_8_1_43T.root',
#'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM9SUSYsftshtpythia6_Summer11/susypat_9_1_z6d.root'
#

##########Seema's LM5 sample
#'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM5SUSYsftshtpythia6_Summer11/susypat_10_1_SBb.root',
#'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM5SUSYsftshtpythia6_Summer11/susypat_11_1_Hkw.root',
#'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM5SUSYsftshtpythia6_Summer11/susypat_1_1_gzT.root',
#'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM5SUSYsftshtpythia6_Summer11/susypat_12_1_Stc.root',
#'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM5SUSYsftshtpythia6_Summer11/susypat_13_1_se3.root',
#'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM5SUSYsftshtpythia6_Summer11/susypat_14_1_Xzg.root',
#'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM5SUSYsftshtpythia6_Summer11/susypat_15_1_ApF.root',
#'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM5SUSYsftshtpythia6_Summer11/susypat_16_1_N1c.root',
#'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM5SUSYsftshtpythia6_Summer11/susypat_17_1_pgW.root',
#'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM5SUSYsftshtpythia6_Summer11/susypat_2_1_rGh.root',
#'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM5SUSYsftshtpythia6_Summer11/susypat_3_1_SZr.root',
#'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM5SUSYsftshtpythia6_Summer11/susypat_4_1_PZW.root',
#'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM5SUSYsftshtpythia6_Summer11/susypat_5_1_B8G.root',
#'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM5SUSYsftshtpythia6_Summer11/susypat_6_1_xKn.root',
#'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM5SUSYsftshtpythia6_Summer11/susypat_7_1_zpK.root',
#'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM5SUSYsftshtpythia6_Summer11/susypat_8_1_Aly.root',
#'file:/eos/uscms/store/user/seema/SusyRA2Analysis/19September2011_2JetSkim_LM5SUSYsftshtpythia6_Summer11/susypat_9_1_rPo.root'

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


from UserCode.MetSigForQcd.metsig_cfi import *

usePrescaleWeight = 0
applyLumiweighing = 0
applyEventweighing = 0

print 'usePrescaleWeight is set to  = ', usePrescaleWeight
print 'applyLumiweighing is set to  = ', applyLumiweighing
print 'applyEventweighing is set to = ', applyEventweighing

from SandBox.Skims.RA2Leptons_cff import *
from SandBox.Skims.mhtProducer_cfi import *
# MHT using PF Jets
process.mhtPF = mht.clone()
process.mhtPF.JetCollection = cms.InputTag('patJetsAK5PFPt30')

njet = 2
njetmax = 1000
verbose = False
njetCut = 1
dPhiCut = 1
dMinMetSig = 0

process.metsignomht = metsig.clone(
                      dMinMHT = cms.untracked.double(0.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig),   
							 dMinNjet = cms.untracked.double(njet),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing),
                      verbose = cms.untracked.int32(verbose)
)

process.metsigmht200 = metsig.clone(
                      dMinMHT = cms.untracked.double(200.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig),   
							 dMinNjet = cms.untracked.double(njet),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing),
                      verbose = cms.untracked.int32(verbose)

)

process.metsigmht350 = metsig.clone(
                      dMinMHT = cms.untracked.double(350.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig),   
							 dMinNjet = cms.untracked.double(njet),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing),
                      verbose = cms.untracked.int32(verbose)

)

process.metsigmht500 = metsig.clone(
                      dMinMHT = cms.untracked.double(500.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig),   
							 dMinNjet = cms.untracked.double(njet),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing),
                      verbose = cms.untracked.int32(verbose)

)


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
							 process.metsignomht +
							 process.metsigmht200 +
							 process.metsigmht350 +
							 process.metsigmht500
)


process.p = cms.Path(process.preseq)

##-- Output module configuration ---------------------------------------
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outputFile)
)

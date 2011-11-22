import FWCore.ParameterSet.Config as cms

#maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_10_1_TXG.root'
,'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_11_1_x1a.root'
,'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_12_1_pfP.root'
,'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_13_1_7OQ.root'
,'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_14_1_wDF.root'
,'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_15_1_Wmj.root'
,'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_16_1_sE4.root'
,'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_17_1_Pqr.root'
,'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_18_1_5Za.root'
,'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_1_1_QO3.root'
,'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_2_1_Kz9.root'
,'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_3_1_ij9.root'
,'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_4_1_1Cj.root'
,'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_5_1_MNO.root'
,'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_6_1_Eq3.root'
,'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_7_1_uE4.root'
,'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_8_1_uWQ.root'
,'/store/user/samantha/2011RA2/LM9_SUSY_7TeV_pythia6/samantha/LM9_SUSY_sftsht_7TeV-pythia6/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_ra2StdCleaned_NoMHTDelPhiCuts/9973ee9e540802ff8e80fdb7484f97b5/LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_9_1_qXc.root'

] );


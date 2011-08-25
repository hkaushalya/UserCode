import FWCore.ParameterSet.Config as cms

#maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
 '/store/user/samantha/2011RA2/QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6/2011Aug_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_01.root'
 ,'/store/user/samantha/2011RA2/QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6/2011Aug_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_02.root'
 ,'/store/user/samantha/2011RA2/QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6/2011Aug_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_03.root'
 ,'/store/user/samantha/2011RA2/QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6/2011Aug_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_04.root'
 ,'/store/user/samantha/2011RA2/QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6/2011Aug_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_05.root'
 ,'/store/user/samantha/2011RA2/QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6/2011Aug_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_06.root'
 ,'/store/user/samantha/2011RA2/QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6/2011Aug_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_07.root'
 ,'/store/user/samantha/2011RA2/QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6/2011Aug_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_08.root'
 ,'/store/user/samantha/2011RA2/QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6/2011Aug_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_09.root'
 ,'/store/user/samantha/2011RA2/QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6/2011Aug_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_10.root'
 ,'/store/user/samantha/2011RA2/QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6/2011Aug_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_11.root'
] );


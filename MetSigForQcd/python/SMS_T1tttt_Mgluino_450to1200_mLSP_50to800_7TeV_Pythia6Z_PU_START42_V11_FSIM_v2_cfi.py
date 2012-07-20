import FWCore.ParameterSet.Config as cms

print '[MetSigForQCD] loading SMS_T1tttt_Mgluino_450to1200_mLSP_50to800_7TeV_Pythia6Z_PU_START42_V11_FSIM_v2_cfi'

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()

source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
	#Total no. of events = 49997
	#Cross-section = 0.01 pb
	'file:/uscms_data/d3/samantha/SUSYRA2work/CMSSW_4_2_5/src/SandBox/Skims/T1pats/susy_pat1.root'
	,'file:/uscms_data/d3/samantha/SUSYRA2work/CMSSW_4_2_5/src/SandBox/Skims/T1pats/susy_pat2.root'
	,'file:/uscms_data/d3/samantha/SUSYRA2work/CMSSW_4_2_5/src/SandBox/Skims/T1pats/susy_pat3.root'
	,'file:/uscms_data/d3/samantha/SUSYRA2work/CMSSW_4_2_5/src/SandBox/Skims/T1pats/susy_pat4.root'
	,'file:/uscms_data/d3/samantha/SUSYRA2work/CMSSW_4_2_5/src/SandBox/Skims/T1pats/susy_pat5.root'

] );

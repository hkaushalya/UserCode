evts = -1

import FWCore.ParameterSet.Config as cms
import sys

fileNameContainer = cms.untracked.vstring()
file = sys.argv[3]

print 'arg 4=' , sys.argv[4]
print 'arg 5=' , sys.argv[5]


infile = open (file, "r")
for line in infile.readlines():
       line = line.rstrip('\n') # equivalent to Perl's chomp
       fileNameContainer.append (line)

print fileNameContainer

runningOnMC = True

process = cms.Process("treemaker")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True)
)

process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# ############ Database ##############
process.EcalLaserCorrectionService = cms.ESProducer("EcalLaserCorrectionService")

process.GlobalTag.globaltag = "START53_V7G::All"

print 'USING GT : ', process.GlobalTag.globaltag 

#-- Message Logger ------------------------------------------------------------
process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
    limit = cms.untracked.int32(10),
    reportEvery = cms.untracked.int32(1000)
    )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


#-- Input Source --------------------------------------------------------------
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(evts))

process.source = cms.Source("PoolSource", 
			fileNames = fileNameContainer
			#fileNames = cms.untracked.vstring(
       #'/store/user/lpcsusyhad/QCD_Pt_300to470_TuneZ2star_8TeV_pythia6_Summer12/samantha/QCD_Pt-300to470_TuneZ2star_8TeV_pythia6/Summer12_PU_S7_START52_V9_RA2_NoCuts_09Aug2012V1/f42b6bbcc7ea08f6809cc21efa861dac/susypat_261_1_RZ9.root'
   #)
)

process.load('SandBox.Skims.RA2Leptons_cff')
#process.LeptonVeto = cms.Sequence(
#		process.ra2PFMuonVeto *
#		process.ra2ElectronVeto
#		)


#for new JECs
process.load("SandBox.Skims.jesChange_cfi")

#Please use correct GlobalTag. 
#Please be careful to correct JECLevels and met-phi-corrections for data and MC
if runningOnMC == True :
	process.newJetsMET.JECLevel = cms.string('ak5PFchsL1FastL2L3')
	process.patMETPhiCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runABCvsNvtx_mc
else :
	process.newJetsMET.JECLevel = cms.string('ak5PFchsL1FastL2L3Residual')
	process.patMETPhiCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runABCvsNvtx_data

# update missing noise filters
process.load("SandBox.Skims.RA2Cleaning_cff")
process.ra2PBNR.JetSource   = "newpatJetsPFchsPt30"
process.ra2PBNR.taggingMode = False

#process.load("SandBox.Skims.ecalLaserCorrFilter_cfi")
process.load("RecoMET.METFilters.ecalLaserCorrFilter_cfi")
process.ecalLaserCorrFilter.Debug = False
process.ecalLaserCorrFilter.taggingMode = False

from UserCode.MuonEffForHadTau.muoneff_cfi import *
process.muoneff = mueff.clone(
	#MuonSource        = cms.InputTag("patMuonsPF")
	MuonSource        = cms.InputTag("patMuonsPFIDIso")
	#FileNamePUDataDistribution = cms.string("UserCode/LostLeptonTree/DataPileupHistogram_RA2Summer12_190456-208686_ABCD.root"),   
#	FileNamePUDataDistribution = cms.string("DataPileupHistogram_RA2Summer12_190456-208686_ABCD.root"),   
	
)

process.preseq = cms.Sequence(
					#process.HCALLaserEvtFilterList2012 *
					process.ra2Electrons *   #for new JECs
#					process.LeptonVeto *
#					process.puWeight *
					process.newra2PFchsJets *  #new JECs. must use the new correct GTs
					process.ra2PBNR *
              	process.ecalLaserCorrFilter *
					process.muoneff
)

##-- Output module configuration ---------------------------------------
process.TFileService = cms.Service("TFileService",
  fileName = cms.string(sys.argv[4]+'_'+sys.argv[5]+'.root')
)

process.ppf = cms.Path( process.preseq )


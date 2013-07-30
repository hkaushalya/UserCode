evts = 100

import FWCore.ParameterSet.Config as cms
import sys

print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

if (len(sys.argv) != 5): 
	print ' ********************************************** '
	print ' Need 5 arguments'
	str = ' ' + sys.argv[0] + ' ' + sys.argv[1] +' filelist.txt rootfile_namepart1 rootfile_namepart2'
	print str
	print ' ********************************************** '
	sys.exit(0)


fileNameContainer = cms.untracked.vstring()
file = sys.argv[2]

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


#This global tag is good for 53X MC samples in lpcsusyhad area
#for collision data (or other MC, need to use appropriate GT.
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

from UserCode.RA2Example.ra2example_cfi import *
process.ra2example = demo.clone(
)

process.ra2example.MCflag = runningOnMC

process.load('SandBox.Skims.RA2Leptons_cff')
process.LeptonVeto = cms.Sequence(
		process.ra2PFMuonVeto *
		process.ra2ElectronVeto
		)


from RA2Classic.WeightProducer.weightProducer_cfi import *

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


if runningOnMC == True :
	process.preseq = cms.Sequence(
					#process.HCALLaserEvtFilterList2012 *
					process.ra2Electrons *   #for new JECs
					process.LeptonVeto *
#					process.puWeight *
					#process.dumpcont *
					process.newra2PFchsJets *  #new JECs. must use the new correct GTs
					process.ra2PBNR *
              	process.ecalLaserCorrFilter *
					process.ra2example 
	)
else :
	process.preseq = cms.Sequence(
					#process.HCALLaserEvtFilterList2012 *
					process.ra2Electrons *   #for new JECs
					process.ra2NoiseCleaning *
					process.LeptonVeto *
					process.trigProd * 
					process.newra2PFchsJets *  #new JECs. must use the new correct GTs
					process.ra2PBNR *
              	process.ecalLaserCorrFilter *
					process.ra2example
	)

##-- Output module configuration ---------------------------------------
process.TFileService = cms.Service("TFileService",
  fileName = cms.string(sys.argv[3]+'_'+sys.argv[4]+'.root')
)

process.ppf = cms.Path( process.preseq )


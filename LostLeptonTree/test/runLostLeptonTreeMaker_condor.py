evts = -1
outputFile="Ntuple.root"

import FWCore.ParameterSet.Config as cms
import sys

fileNameContainer = cms.untracked.vstring()
file = sys.argv[3]

infile = open (file, "r")
for line in infile.readlines():
       line = line.rstrip('\n') # equivalent to Perl's chomp
       fileNameContainer.append (line)

print fileNameContainer


runningOnMC = False

process = cms.Process("treemaker")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True)
)

process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#GT to use from RA2 2012 page
#2012A+B (Jul13 rereco, 53X) 	FT_53_V6_AN2::All
#2012A (Aug06 rereco, 53X) 	FT_53_V6C_AN2::All
#2012Cv1 (Aug24 rereco, 53X) 	FT_53_V10_AN2::All
#2012Cv2 (prompt reco, 53X) 	GR_P_V41_AN2::All
#Summer12_DR53X 53X 	START53_V7F::All

if runningOnMC:
	process.GlobalTag.globaltag = "START52_V11C::All"
else:
	process.GlobalTag.globaltag = "FT_53_V6_AN2::All"
print 'USING GT : ', process.GlobalTag.globaltag 

#-- Message Logger ------------------------------------------------------------
process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
    limit = cms.untracked.int32(10),
    reportEvery = cms.untracked.int32(1)
    )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


#-- Input Source --------------------------------------------------------------
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(evts))

process.source = cms.Source("PoolSource", 
			fileNames = fileNameContainer
			#fileNames = cms.untracked.vstring(
       #'/store/user/lpcsusyhad/QCD_Pt_300to470_TuneZ2star_8TeV_pythia6_Summer12/samantha/QCD_Pt-300to470_TuneZ2star_8TeV_pythia6/Summer12_PU_S7_START52_V9_RA2_NoCuts_09Aug2012V1/f42b6bbcc7ea08f6809cc21efa861dac/susypat_261_1_RZ9.root'
   )
)

from UserCode.LostLeptonTree.lostleptontree_cfi import *
process.treeMaker = demo.clone(
	no_beamHaloFilter = cms.bool(True),
	no_HBHENoiseFilter = cms.bool(True)

)
process.treeMaker.MCflag = runningOnMC
process.treeMaker.DoPUReweight = False
process.treeMaker.ApplyEventWeighing = False  #for flat sample 


process.preseq = cms.Sequence(
		      process.treeMaker 
)

##-- Output module configuration ---------------------------------------
process.TFileService = cms.Service("TFileService",
  fileName = cms.string(sys.argv[4]+'_'+sys.argv[5]+'.root')
)

process.ppf = cms.Path( process.preseq )


import FWCore.ParameterSet.Config as cms

outputFile="Ntuple.root"
evts = -1

process = cms.Process("treemaker")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True)
)

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "START52_V11C::All"
print 'USING GT for 523 : ', process.GlobalTag.globaltag 

#-- Message Logger ------------------------------------------------------------
process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
    limit = cms.untracked.int32(10),
    reportEvery = cms.untracked.int32(1)
    )
process.MessageLogger.cerr.FwkReport.reportEvery = -1


#-- Input Source --------------------------------------------------------------
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(evts))

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
#       '/store/user/lpcsusyhad/QCD_Pt_300to470_TuneZ2star_8TeV_pythia6_Summer12/samantha/QCD_Pt-300to470_TuneZ2star_8TeV_pythia6/Summer12_PU_S7_START52_V9_RA2_NoCuts_09Aug2012V1/f42b6bbcc7ea08f6809cc21efa861dac/susypat_261_1_RZ9.root'
   )
)

from UserCode.RA2QCDvetoAna.lostleptontree_cfi import *
process.treeMaker = demo.clone(
)

from SandBox.Skims.mhtProducer_cfi import *

#process.mhtprod = mht.clone()
#	process.mhtprod.JetCollection = cms.InputTag("patJetsAK5PF")

process.preseq = cms.Sequence(
			#	process.mhtprod *
		      process.treeMaker 
)

##-- Output module configuration ---------------------------------------
process.TFileService = cms.Service("TFileService",
					 fileName = cms.string(outputFile)
)

process.ppf = cms.Path( process.preseq )


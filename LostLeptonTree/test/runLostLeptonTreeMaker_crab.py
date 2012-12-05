import FWCore.ParameterSet.Config as cms

runningOnMC = False
outputFile="Ntuple.root"
evts = -1

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
#'/store/user/lpcsusyhad/QCD_Pt_15to3000_TuneZ2_Flat_8TeV_pythia6_Summer12/samantha/QCD_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6/Summer12_PU_S7_START52_V9_RA2_NoCuts_09Aug2012V1/f42b6bbcc7ea08f6809cc21efa861dac/susypat_1_1_Zts.root'
'/store/user/lpcsusyhad/53X_ntuples/HT_Run2012A_13Jul2012_v1_lpc1/seema/HT/HT_Run2012A-13Jul2012-v1_NOCUTS_HLTPFHTInc_12Oct2012V3/32ed558ccb6d6c1922414263922f006c/susypat_99_1_N7w.root'
   )
)

#process.load("UserCode.RA2QCDvetoAna.T1tttt_mG800_mLSP200_NOCUTS_09Aug2012V1_lacroix_cfi")

from UserCode.LostLeptonTree.lostleptontree_cfi import *
process.treeMaker = demo.clone(
	no_beamHaloFilter = cms.bool(True),
	no_HBHENoiseFilter = cms.bool(True)

)
process.treeMaker.MCflag = runningOnMC
process.treeMaker.DoPUReweight = False
process.treeMaker.ApplyEventWeighing = False  #for flat sample 



from SandBox.Skims.mhtProducer_cfi import *

#process.mhtprod = mht.clone()
#	process.mhtprod.JetCollection = cms.InputTag("patJetsAK5PF")


#from UserCode.RA2QCDvetoAna.trigprescaleWeight_cfi import *
#process.trgProd = trigPrescaleWgtProd.clone()

process.preseq = cms.Sequence(
			#	process.mhtprod *
				#process.trgProd *
		      process.treeMaker 
)

##-- Output module configuration ---------------------------------------
process.TFileService = cms.Service("TFileService",
					 fileName = cms.string(outputFile)
)

process.ppf = cms.Path( process.preseq )


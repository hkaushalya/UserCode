import FWCore.ParameterSet.Config as cms
import os

runningOnMC = False
outputFile="Ntuple.root"
evts = 10000

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
#2012D (prompt reco, 53X) 	GR_P_V42_AN2::All
#Summer12_DR53X 53X 	START53_V7F::All

if runningOnMC:
	process.GlobalTag.globaltag = "START53_V7F::All"
else:
	#process.GlobalTag.globaltag = "FT_53_V6_AN2::All"
	process.GlobalTag.globaltag = "GR_P_V42_AN2::All"
	#process.GlobalTag.globaltag = "FT_P_V42_AN3::All"


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

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
#       '/store/user/lpcsusyhad/QCD_Pt_300to470_TuneZ2star_8TeV_pythia6_Summer12/samantha/QCD_Pt-300to470_TuneZ2star_8TeV_pythia6/Summer12_PU_S7_START52_V9_RA2_NoCuts_09Aug2012V1/f42b6bbcc7ea08f6809cc21efa861dac/susypat_261_1_RZ9.root'
#'/store/user/lpcsusyhad/QCD_Pt_15to3000_TuneZ2_Flat_8TeV_pythia6_Summer12/samantha/QCD_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6/Summer12_PU_S7_START52_V9_RA2_NoCuts_09Aug2012V1/f42b6bbcc7ea08f6809cc21efa861dac/susypat_1_1_Zts.root'
#'/store/user/lpcsusyhad/53X_ntuples/HT_Run2012A_13Jul2012_v1_lpc1/seema/HT/HT_Run2012A-13Jul2012-v1_NOCUTS_HLTPFHTInc_12Oct2012V3/32ed558ccb6d6c1922414263922f006c/susypat_99_1_N7w.root'
#'/store/user/lpcsusyhad/53X_ntuples/HT_Run2012A_13Jul2012_v1_lpc1/seema/HT/HT_Run2012A-13Jul2012-v1_NOCUTS_HLTPFHTInc_12Oct2012V3/32ed558ccb6d6c1922414263922f006c/susypat_98_1_C1E.root'
#'/store/user/lpcsusyhad/53X_ntuples/QCD_HT_250To500_MGPythia_v1_lpc1/seema/QCD_HT-250To500_TuneZ2star_8TeV-madgraph-pythia6/QCD_HT-250To500_TuneZ2star_8TeV-madgraph-pythia_v1_NOCUTS_12Oct2012V3/c31a0db70b73f0b9a355af58227b92dd/susypat_1000_1_nDc.root'

#2012A_July13
#'/store/user/lpcsusyhad/53X_ntuples/HT_Run2012A_13Jul2012_v1_lpc1/seema/HT/HT_Run2012A-13Jul2012-v1_NOCUTS_HLTPFHTInc_12Oct2012V3/32ed558ccb6d6c1922414263922f006c/susypat_99_1_N7w.root'
#'/store/user/lpcsusyhad/53X_ntuples/HT_Run2012A_13Jul2012_v1_lpc1/seema/HT/HT_Run2012A-13Jul2012-v1_NOCUTS_HLTPFHTInc_12Oct2012V3/32ed558ccb6d6c1922414263922f006c/susypat_1_1_RjB.root'

#2012A_Aug06
#'/store/user/lpcsusyhad/53X_ntuples/HT_Run2012Arecover_06Aug2012_v1_lpc1/seema/HT/HT_Run2012A-recover-06Aug2012-v1_NOCUTS_HLTPFHTInc_12Oct2012V3/8e90e164b25bf8451f454d8c6a3799cd/susypat_1_1_mPZ.root'
#'/store/user/lpcsusyhad/53X_ntuples/HT_Run2012Arecover_06Aug2012_v1_lpc1/seema/HT/HT_Run2012A-recover-06Aug2012-v1_NOCUTS_HLTPFHTInc_12Oct2012V3/8e90e164b25bf8451f454d8c6a3799cd/susypat_390_1_0GD.root'

#2012B_July13
#	'/store/user/lpcsusyhad/53X_ntuples/HTMHT_Run2012B_13Jul2012_v1_lpc1/seema/HTMHT/HTMHT_Run2012B-13Jul2012-v1_NOCUTS_HLTPFHTInc_12Oct2012V3/32ed558ccb6d6c1922414263922f006c/susypat_1_1_yIg.root'
#	'/store/user/lpcsusyhad/53X_ntuples/HTMHT_Run2012B_13Jul2012_v1_lpc1/seema/HTMHT/HTMHT_Run2012B-13Jul2012-v1_NOCUTS_HLTPFHTInc_12Oct2012V3/32ed558ccb6d6c1922414263922f006c/susypat_999_1_jfv.root'

#2012C promptReco
	#'/store/user/lpcsusyhad/53X_ntuples/HTMHT_Run2012C_PromptReco_v2_lpc1/seema/HTMHT/HTMHT_Run2012C-PromptReco-v2_NoHepTopTagger_NOCUTS_HLTPFHTInc_12Oct2012V3/993e6ea01263feda15a347a7daafedb6/susypat_1_1_3FP.root'
	#'/store/user/lpcsusyhad/53X_ntuples/HTMHT_Run2012C_PromptReco_v2_lpc1/seema/HTMHT/HTMHT_Run2012C-PromptReco-v2_NoHepTopTagger_NOCUTS_HLTPFHTInc_12Oct2012V3/993e6ea01263feda15a347a7daafedb6/susypat_999_1_thS.root'

#HTMHT_Run2012C_24Aug2012_v1_lpc1
#	'/store/user/lpcsusyhad/53X_ntuples/HTMHT_Run2012C_24Aug2012_v1_lpc1/seema/HTMHT/HTMHT_Run2012C-24Aug2012-v1_NoHepTopTagger_NOCUTS_HLTPFHTInc_12Oct2012V3/5752ae8db176959272f91d699d678612/susypat_1_1_Op8.root',
#	'/store/user/lpcsusyhad/53X_ntuples/HTMHT_Run2012C_24Aug2012_v1_lpc1/seema/HTMHT/HTMHT_Run2012C-24Aug2012-v1_NoHepTopTagger_NOCUTS_HLTPFHTInc_12Oct2012V3/5752ae8db176959272f91d699d678612/susypat_111_1_NF3.root'

#HTMHT_Run2012D_PromptReco_v1_lpc1
#'/store/user/lpcsusyhad/53X_ntuples/HTMHT_Run2012D_PromptReco_v1_lpc1/seema/HTMHT/HTMHT_Run2012D-PromptReco-v1_NoHepTopTagger_NOCUTS_HLTPFHTInc_12Oct2012V3/b187b7cb1a0ab43ecf500fbb58f6aeb9/susypat_1088_1_LtJ.root'
#,'/store/user/lpcsusyhad/53X_ntuples/HTMHT_Run2012D_PromptReco_v1_lpc1/seema/HTMHT/HTMHT_Run2012D-PromptReco-v1_NoHepTopTagger_NOCUTS_HLTPFHTInc_12Oct2012V3/b187b7cb1a0ab43ecf500fbb58f6aeb9/susypat_1089_1_Scm.root'

'dcap://pnfs/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/HTMHT_Run2012D_PromptReco_v1_lpc1/seema/HTMHT/HTMHT_Run2012D-PromptReco-v1_NoHepTopTagger_NOCUTS_HLTPFHTInc_12Oct2012V3/b187b7cb1a0ab43ecf500fbb58f6aeb9/susypat_821_1_SP3.root'
,'dcap://pnfs/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/HTMHT_Run2012D_PromptReco_v1_lpc1/seema/HTMHT/HTMHT_Run2012D-PromptReco-v1_NoHepTopTagger_NOCUTS_HLTPFHTInc_12Oct2012V3/b187b7cb1a0ab43ecf500fbb58f6aeb9/susypat_82_1_loz.root'

   )
)

#process.load("UserCode.RA2QCDvetoAna.T1tttt_mG800_mLSP200_NOCUTS_09Aug2012V1_lacroix_cfi")

from UserCode.LostLeptonTree.lostleptontree_cfi import *
process.treeMaker = demo.clone(
	FileNamePUDataDistribution = cms.string("UserCode/LostLeptonTree/DataPileupHistogram_RA2Summer12_190456-208686_ABCD.root"),   
	no_beamHaloFilter = cms.bool(False),
	no_HBHENoiseFilter = cms.bool(False)
)

process.treeMaker.MCflag = runningOnMC
process.treeMaker.DoPUReweight = runningOnMC
process.treeMaker.ApplyEventWeighing = False  #for flat sample 


#from SandBox.Skims.mhtProducer_cfi import *
#process.mhtprod = mht.clone()
#	process.mhtprod.JetCollection = cms.InputTag("patJetsAK5PF")


from UserCode.TrigFilter.trigprescaleWeight_cfi import *
process.trigProd = trigPrescaleWgtProd.clone(
		#hltPaths = cms.vstring('HLT_HT*', 'HLT_PFHT*'),
		hltPaths = cms.vstring('HLT_HT*', 'HLT_PF*'),
		debug =cms.int32(1)
)

Tracer = cms.Service(
     "Tracer" 
     )


 ## --- Selection sequences ---------------------------------------------

     # Filter-related selection
#from RecoMET.METFilters.jetIDFailureFilter_cfi import jetIDFailure
#process.PBNRFilter = jetIDFailure.clone(
#		JetSource = cms.InputTag('MHTJets'),
#		MinJetPt      = cms.double(30.0),
#		taggingMode   = cms.bool(True)
#		)
#from RecoMET.METFilters.multiEventFilter_cfi import multiEventFilter
#os.environ['CMSSW_SEARCH_PATH'] = '/pnfs/cms/WAX/11/store/user/samantha/2011RA2/'
#process.HCALLaserEvtFilterList2012 = multiEventFilter.clone(
#		file        = cms.FileInPath('AllBadHCALLaser_Nov192012.txt'),
#		taggingMode = cms.bool(True)
#		)

process.load('SandBox.Skims.RA2Leptons_cff')
process.LeptonVeto = cms.Sequence(
		process.ra2PFMuonVeto *
		process.ra2ElectronVeto
		)


from RA2Classic.WeightProducer.weightProducer_cfi import *

#for this to work the CMMSW_SEARCH_PATH should be updated using 
# setenv CMSSW_SEARCH_PATH mylocation:$CMSSW_SEARCH_PATH 
process.puWeight = weightProducer.clone(
#	FileNamePUDataDistribution = cms.string("DataPileupHistogram_RA2Summer12_190456-208686_ABCD.root"),
	FileNamePUDataDistribution = cms.string("RA2Classic/AdditionalInputFiles/DataPileupHistogram_RA2Summer12_190456-208686_ABCD.root"),   
 # for summer 12 MC samples
	PU = cms.int32(3)
)

process.dumpcont = cms.EDAnalyzer("EventContentAnalyzer")
#process.pppp = cms.Path(process.dumpcont)
#process.Tracer = cms.Service("Tracer")
process.preseq = cms.Sequence(
			#	process.mhtprod *
				#process.PBNRFilter *
				#process.HCALLaserEvtFilterList2012 *
				process.LeptonVeto *
				process.trigProd * 
				#process.puWeight *
				#process.dumpcont *
		      process.treeMaker 
)

##-- Output module configuration ---------------------------------------
process.TFileService = cms.Service("TFileService",
					 fileName = cms.string(outputFile)
)

process.ppf = cms.Path( process.preseq )


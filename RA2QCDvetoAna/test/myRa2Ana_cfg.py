#
#  SUSY-PAT configuration file adapted for RA2 workflow
#
#  PAT configuration for the SUSY group - 41X series
#  More information here:
#  https://twiki.cern.ch/twiki/bin/view/CMS/SusyPatLayer1DefV10
#

# Starting with a skeleton process which gets imported with the following line
#from PhysicsTools.PatAlgos.patTemplate_cfg import *

outputFile="Default.root"
evts = 6000
#evts = -1
runningOnMC = True 
#runningOnMC = False
myDataset = 1

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('standard')

# setup any defaults you want
#options.myEvents = 200
#options.outputFile="default.root"
#options.myDataset=2
#options.register('GlobalTag', "MC_42_V10::All", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "GlobaTTag to use (otherwise default Pat GT is used)")
#options.register('mcInfo', True, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "process MonteCarlo data, default is data")

options.register('myDataset', 10, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "dataset QCD or SUSY.")
options.register('outputFile', "default.root", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Output file name. Default will be choosen according to dataset.")
options.register('evts', 5000, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "Events to process. Default 1000 events")
#"Dataset (QCD ==1, SUSY LM9==2)")
#options.myDataset=1
#options.maxEvents = 200

# get and parse the command line arguments
options.parseArguments()


import FWCore.ParameterSet.Config as cms

process = cms.Process("MySimplePATAnalysis")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True)
)

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "START42_V12::All"
if runningOnMC == False:
    process.GlobalTag.globaltag = "GR_R_42_V12::All"

#-- Message Logger ------------------------------------------------------------
process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
    limit = cms.untracked.int32(10),
    reportEvery = cms.untracked.int32(1)
    )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


#-- Input Source --------------------------------------------------------------
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(evts))

#process.load("UserCode.JERStudy.JER_2011Aug_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_shortlist_cfi")
#print "myDataset=", myDataset
if (options.myDataset == 1):
	outputFile = "DelPhiMin_QCD.root"
	#process.load("UserCode.RA2QCDvetoAna.JER_2011Aug_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_MERGED_cfi")
	print "Running over QCD sample, output file ", outputFile
elif (options.myDataset == 11):
	outputFile = "DelPhiMin_QCD_1.root"
	#process.load("UserCode.RA2QCDvetoAna.JER_2011Aug_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_MERGED_1of11_cfi")
	process.load("UserCode.RA2QCDvetoAna.QCD_Pt_120to170_TuneZ2_7TeV_pythia6_Summer11_PU_S3_START42_V11_v2_AODSIM_cfi")
	print "Running over 1st QCD sample, output file ", outputFile
elif (options.myDataset == 12):
	outputFile = "DelPhiMin_QCD_2.root"
	#process.load("UserCode.RA2QCDvetoAna.JER_2011Aug_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_MERGED_2of11_cfi")
	process.load("UserCode.RA2QCDvetoAna.QCD_Pt_170to300_TuneZ2_7TeV_pythia6_Summer11_PU_S3_START42_V11_v2_AODSIM_cfi")
	print "Running over 2nd QCD sample, output file ", outputFile
elif (options.myDataset == 13):
	outputFile = "DelPhiMin_QCD_3.root"
	#process.load("UserCode.RA2QCDvetoAna.JER_2011Aug_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_MERGED_3of11_cfi")
	process.load("UserCode.RA2QCDvetoAna.QCD_Pt_300to470_TuneZ2_7TeV_pythia6_Summer11_PU_S3_START42_V11_v2_AODSIM_cfi")
	print "Running over 3rd QCD sample, output file ", outputFile
elif (options.myDataset == 14):
	outputFile = "DelPhiMin_QCD_4.root"
	#process.load("UserCode.RA2QCDvetoAna.JER_2011Aug_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_MERGED_4of11_cfi")
	process.load("UserCode.RA2QCDvetoAna.QCD_Pt_470to600_TuneZ2_7TeV_pythia6_Summer11_PU_S3_START42_V11_v2_AODSIM_cfi")
	print "Running over 4th QCD sample, output file ", outputFile
elif (options.myDataset == 15):
	outputFile = "DelPhiMin_QCD_5.root"
	#process.load("UserCode.RA2QCDvetoAna.JER_2011Aug_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_MERGED_5of11_cfi")
	process.load("UserCode.RA2QCDvetoAna.QCD_Pt_600to800_TuneZ2_7TeV_pythia6_Summer11_PU_S3_START42_V11_v2_AODSIM_cfi")
	print "Running over 4th QCD sample, output file ", outputFile
elif (options.myDataset == 16):
	outputFile = "DelPhiMin_QCD_6.root"
	#process.load("UserCode.RA2QCDvetoAna.JER_2011Aug_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_MERGED_6of11_cfi")
	process.load("UserCode.RA2QCDvetoAna.QCD_Pt_800to1000_TuneZ2_7TeV_pythia6_Summer11_PU_S3_START42_V11_v2_AODSIM_cfi")
	print "Running over 4th QCD sample, output file ", outputFile
elif (options.myDataset == 17):
	outputFile = "DelPhiMin_QCD_7.root"
	#process.load("UserCode.RA2QCDvetoAna.JER_2011Aug_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_MERGED_7of11_cfi")
	process.load("UserCode.RA2QCDvetoAna.QCD_Pt_1000to1400_TuneZ2_7TeV_pythia6_Summer11_PU_S3_START42_V11_v2_AODSIM_cfi")
	print "Running over 4th QCD sample, output file ", outputFile

elif (options.myDataset == 18):
	outputFile = "DelPhiMin_QCD_8.root"
	#process.load("UserCode.RA2QCDvetoAna.JER_2011Aug_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_MERGED_8of11_cfi")
	process.load("UserCode.RA2QCDvetoAna.QCD_Pt_1400to1800_TuneZ2_7TeV_pythia6_Summer11_PU_S3_START42_V11_v2_AODSIM_cfi")
	print "Running over 4th QCD sample, output file ", outputFile

elif (options.myDataset == 19):
	outputFile = "DelPhiMin_QCD_9.root"
	#process.load("UserCode.RA2QCDvetoAna.JER_2011Aug_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_MERGED_9of11_cfi")
	process.load("UserCode.RA2QCDvetoAna.QCD_Pt_1800_TuneZ2_7TeV_pythia6_Summer11_PU_S3_START42_V11_v2_AODSIM_cfi")
	print "Running over 4th QCD sample, output file ", outputFile

elif (options.myDataset == 20):
#	outputFile = "DelPhiMin_QCD_10.root"
#	process.load("UserCode.RA2QCDvetoAna.JER_2011Aug_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_MERGED_10of11_cfi")
	print "Running over 4th QCD sample, output file ", outputFile

elif (options.myDataset == 21):
#	outputFile = "DelPhiMin_QCD_11.root"
#	process.load("UserCode.RA2QCDvetoAna.JER_2011Aug_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_MERGED_11of11_cfi")
	print "Running over 4th QCD sample, output file ", outputFile

elif (options.myDataset == 2):
	outputFile = "DelPhiMin_SUSYLM9.root"
	process.load("UserCode.RA2QCDvetoAna.LM9_SUSY_sftsht_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1_cfi");
	print "Running over SUSY LM9 sample, output file ", outputFile

#if evts != -1 : 
#	print "Using short list of files for testing."
#	process.load("UserCode.RA2QCDvetoAna.HTRun2011APromptRecoV4_160404to167151_shortlist_cfi")
	#process.load("UserCode.RA2QCDvetoAna.2011RA2_Jun09_forrestp_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_cfi")
#else : process.load("UserCode.RA2QCDvetoAna.HTRun2011APromptRecoV4_160404to167151_cfi")

print options
#quit()
#-- check RA2 recipe here ------------------------------------------------------------
process.prefilterCounter        = cms.EDProducer("EventCountProducer")
process.postStdCleaningCounter  = cms.EDProducer("EventCountProducer")

#-- Output module configuration -----------------------------------------------
process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')
process.mhtPFFilter.MinMHT = 0

#process.load('SandBox.Skims.puWeightProducer_cfi')

#process.susyRA2Analysis = cms.EDAnalyzer("SusyRA2Analysis",
#                                         VertexSource = cms.InputTag("goodVerticesRA2"),
#                                         MHTSource = cms.InputTag("mhtPF"),
#                                         HTSource  = cms.InputTag("htPF"),
#                                         JetSource = cms.InputTag("patJetsAK5PFPt50Eta25"),
#                                         METSource = cms.InputTag("patMETsPF"),
#                                         DoPUReweight = cms.bool(False),
#                                         PUWeigthSource = cms.InputTag("puWeight"),
#                                         DoOptimizePlots = cms.bool(True),
#                                         Debug       = cms.bool(False)
#)


process.load('SandBox.Skims.electronSelector_cfi')
process.electronSelector.DoElectronVeto = True

process.QCDvetoAna = cms.EDFilter('RA2QCDvetoAna',
							patJetsPFPt30InputTag = cms.InputTag("patJetsAK5PFPt30"),
							patJetsPFPt50Eta25InputTag = cms.InputTag("patJetsAK5PFPt50Eta25"),
							mhtInputTag = cms.InputTag("mhtPF"),
							htInputTag = cms.InputTag("htPF"),
							nMinVtx = cms.untracked.int32(1),	
							ndofVtx = cms.untracked.int32(4),	
							maxDelzVtx = cms.untracked.double(24.0),	
							maxDelRho = cms.untracked.double(15.0),
							minJetEt4MHt = cms.untracked.double(30.0),	
							maxJetEta4MHt = cms.untracked.double(5.0),	
							dMinHT = cms.untracked.double(350.0),	
							dMinMHt = cms.untracked.double(0.0),	
							verbose = cms.untracked.int32(0)
)


process.fullseq = cms.Sequence(
                      process.ra2StdCleaning *
                      process.ra2PostCleaning *
                      process.ra2EcalTPFilter *
                      process.vetoBEFilterHTRun2011AMay10ReRecov1 *
                      process.ra2FullPFSelection *
							 process.electronSelector *
                      process.QCDvetoAna
)

process.fullseq.remove(process.jetMHTPFDPhiFilter)

process.ppf = cms.Path( process.fullseq )

##-- Output module configuration ---------------------------------------
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outputFile)
)

#process.out.SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('ppf') )
#process.out.fileName = cms.untracked.string('susypat.root')
##process.out.outputCommands = cms.untracked.vstring('keep *')
#process.outpath = cms.EndPath( process.out )
#
#del process.outpath

###-- Dump config ------------------------------------------------------------
##file = open('SusyPAT_RA2414_cfg.py','w')
##file.write(str(process.dumpPython()))
##file.close()


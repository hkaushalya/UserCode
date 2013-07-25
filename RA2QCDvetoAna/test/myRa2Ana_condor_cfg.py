import FWCore.ParameterSet.Config as cms
import sys

process = cms.Process("MySimplePATAnalysis")

fileNameContainer = cms.untracked.vstring()
file = sys.argv[3]

infile = open (file, "r")
for line in infile.readlines():
       line = line.rstrip('\n') # equivalent to Perl's chomp
       fileNameContainer.append (line)

print fileNameContainer
process.source = cms.Source("PoolSource",
                            fileNames = fileNameContainer           
                            #skipEvents = cms.untracked.uint32(10)
)

process.source.duplicateCheckMode = cms.untracked.string('checkAllFilesOpened')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")



process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True)
)
#-- Message Logger ------------------------------------------------------------
process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
    limit = cms.untracked.int32(10),
    reportEvery = cms.untracked.int32(1)
    )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


runningOnMC = True 
process.GlobalTag.globaltag = "START42_V12::All"
if runningOnMC == False:
    process.GlobalTag.globaltag = "GR_R_42_V12::All"



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


##-- Output module configuration ---------------------------------------
process.TFileService = cms.Service("TFileService",
  #fileName = cms.string('file_j'+sys.argv[2]+'.root')
  fileName = cms.string(sys.argv[4]+'_'+sys.argv[5]+'.root')
)

process.ppf = cms.Path(process.fullseq)


###-- Dump config ------------------------------------------------------------
print "---------------- Config file info --------------------------"
#file = open('myRa2Ana_condor_cfg','w')
#file.write(str(process.dumpPython()))
#file.close()


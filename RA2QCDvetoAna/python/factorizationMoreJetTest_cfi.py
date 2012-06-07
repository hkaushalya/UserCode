
import FWCore.ParameterSet.Config as cms

factorizationMoreJetTest = cms.EDFilter(
							"FactorizationMoreJetTest",
                      patJetsPFPt50Eta25InputTag = cms.InputTag("patJetsAK5PFPt50Eta25"),
                      patJetsPFPt30Eta50InputTag = cms.InputTag("patJetsAK5PFPt30"),
                      mhtInputTag = cms.InputTag("mhtPF"),
                      htInputTag = cms.InputTag("htPF"),
                      dMinHT = cms.untracked.double(0.0),   
                      dMinMHT = cms.untracked.double(0.0),   
                      dMaxMHT = cms.untracked.double(99999.0),   
                      dMinNjet50 = cms.untracked.double(3.0),   
                      dMinNjetUsedForDphi = cms.untracked.double(3.0),   
							 usePrescaleWeight = cms.untracked.int32(0),
							 #prescaleWeight = cms.InputTag("prescaleweightProducer:weight"),
							 prescaleWeight = cms.InputTag('trigWgtProd:prescaleWeight'),
                      verbose = cms.untracked.int32(0),
                      ApplyLumiWeighing = cms.untracked.int32(0),
                      ApplyEventWeighing = cms.untracked.int32(0),
							 htBins = cms.vdouble(0,500,800,1000,1200,1400,7000)
 )


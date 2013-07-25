
import FWCore.ParameterSet.Config as cms

factorization = cms.EDFilter(
							"Factorization",
                      patJetsPFPt50Eta25InputTag = cms.InputTag("patJetsAK5PFPt50Eta25"),
                      patJetsPFPt30Eta50InputTag = cms.InputTag("patJetsAK5PFPt30"),
                      mhtInputTag = cms.InputTag("mhtPF"),
                      htInputTag = cms.InputTag("htPF"),
                    #  outputFile = cms.string("Factorization"),
                      dMinHT = cms.untracked.double(0.0),   
                      dMinMHT = cms.untracked.double(0.0),   
                      dMaxMHT = cms.untracked.double(99999.0),   
							 usePrescaleWeight = cms.untracked.int32(0),
							 #prescaleWeight = cms.InputTag("prescaleweightProducer:weight"),
							 prescaleWeight = cms.InputTag('trigWgtProd:prescaleWeight'),
                      verbose = cms.untracked.int32(0),
                      ApplyLumiWeighing = cms.untracked.int32(0),
                      ApplyEventWeighing = cms.untracked.int32(0),
							 HltTriggerResults = cms.InputTag("TriggerResults::HLT"),
							 TriggerPathsToStore = cms.vstring('HLT_HT*'),
							 htBins = cms.vdouble(0,500,800,1000,1200,1400,7000)
 )


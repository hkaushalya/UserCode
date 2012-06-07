
import FWCore.ParameterSet.Config as cms

factorization_mcdebug = cms.EDFilter(
							"Factorization_MCdebug",
                      patJetsPFPt50Eta25InputTag = cms.InputTag("patJetsAK5PFPt50Eta25"),
                      mhtInputTag = cms.InputTag("mhtPF"),
                      htInputTag = cms.InputTag("htPF"),
                      dMinHT = cms.untracked.double(0.0),   
                      dMinMHT = cms.untracked.double(0.0),   
                      verbose = cms.untracked.int32(0)
)

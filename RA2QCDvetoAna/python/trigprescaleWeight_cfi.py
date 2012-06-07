
import FWCore.ParameterSet.Config as cms

factorization = cms.EDFilter(
							"Factorization",
                      patJetsPFPt50Eta25InputTag = cms.InputTag("patJetsAK5PFPt50Eta25"),
                      mhtInputTag = cms.InputTag("mhtPF"),
                      htInputTag = cms.InputTag("htPF"),
                    #  outputFile = cms.string("Factorization"),
                      dMinHT = cms.untracked.double(0.0),   
                      dMinMHT = cms.untracked.double(0.0),   
                      verbose = cms.untracked.int32(0),
                      ApplyLumiWeighing = cms.untracked.int32(0),
                      ApplyEventWeighing = cms.untracked.int32(0)
 )


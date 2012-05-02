import FWCore.ParameterSet.Config as cms

metsig = cms.EDAnalyzer('MetSigForQcd',
                      patJetsPFPt50Eta25InputTag = cms.InputTag("patJetsAK5PFPt50Eta25"),
                      patJetsPFPt30Eta50InputTag = cms.InputTag("patJetsAK5PFPt30"),
                      mhtInputTag = cms.InputTag("mhtPF"),
                      htInputTag = cms.InputTag("htPF"),
							 dMinNjet = cms.untracked.double(3.0),
                      dMinHT = cms.untracked.double(0.0),   
                      dMinMHT = cms.untracked.double(0.0),   
                      dMinMetSig = cms.untracked.double(0.0),   
							 usePrescaleWeight = cms.untracked.int32(0),
							 prescaleWeight = cms.InputTag('trigWgtProd:prescaleWeight'),
                      verbose = cms.untracked.int32(0),
                      ApplyLumiWeighing = cms.untracked.int32(0),
                      ApplyEventWeighing = cms.untracked.int32(0),
							 htBins = cms.vdouble(0,350, 500,800,1000,1200,1400,7000),
							 applyNjetCut = cms.untracked.int32(1),
							 applyHtCut   = cms.untracked.int32(1),
							 applyMhtCut  = cms.untracked.int32(1),
							 applyDphiCut = cms.untracked.int32(1)

)

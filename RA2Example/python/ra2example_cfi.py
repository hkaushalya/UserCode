import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('RA2Example',
                                        Debug           = cms.bool    (False),
                                        VertexSource    = cms.InputTag('goodVertices'),
  													 ApplyEventWeighing = cms.untracked.int32(0),
                                        DoPUReweight    = cms.bool    (True),
                                        PUWeigthSource  = cms.InputTag("puWeight"),
                                        PUWeigthABSource  = cms.InputTag("puWeightAB:weight"),
                                        PUWeigthABCSource = cms.InputTag("puWeightABC:weight"),
                                        PFMetSource     = cms.InputTag("patMETsPF"),
                                        JetAllSource    = cms.InputTag("patJetsPF"),
                                        genJetAllSource = cms.InputTag("ak5GenJets"),
                                        MHTSource       = cms.InputTag("mhtPFchs"),
                                        HTSource        = cms.InputTag("htPFchs"),
                                        MinPFJetPt      = cms.untracked.double(0.0),
                                        MinGENJetPt     = cms.untracked.double(0.0),
                                        HBHENoiseFiltSrc = cms.InputTag('HBHENoiseFilterRA2'),
													 MCflag          = cms.bool (False)
)

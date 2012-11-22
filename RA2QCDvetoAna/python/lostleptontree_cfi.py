import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('LostLeptonTree',
                                        Debug           = cms.bool    (False),
                                        VertexSource    = cms.InputTag('goodVertices'),
                                        DoPUReweight    = cms.bool    (False),
                                        PUWeigthSource  = cms.InputTag("puWeight"),
                                        PFMetSource     = cms.InputTag("patMETsPF"),
                                        JetAllSource    = cms.InputTag("patJetsAK5PF"),
                                        genJetAllSource = cms.InputTag("ak5GenJets"),
                                        MHTSource       = cms.InputTag("mhtPFchs"),
                                        HTSource        = cms.InputTag("htPFchs"),
                                        MinPFJetPt      = cms.untracked.double(0.0),
                                        MinGENJetPt     = cms.untracked.double(0.0)
)

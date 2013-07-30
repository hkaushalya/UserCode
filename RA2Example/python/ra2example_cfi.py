import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('RA2Example',
                                        Debug           = cms.bool    (False),
                                        VertexSource    = cms.InputTag('goodVertices'),
                                        PFMetSource     = cms.InputTag("patMETsPF"),
                                        JetAllSource    = cms.InputTag("patJetsPF"),
                                        chsJetAllSource = cms.InputTag("newpatJetsPFchsPt30"),
                                        genJetAllSource = cms.InputTag("ak5GenJets"),
                                        MHTSource       = cms.InputTag("mhtPFchs"),
                                        HTSource        = cms.InputTag("htPFchs"),
                                        MinPFJetPt      = cms.untracked.double(0.0),
                                        MinGENJetPt     = cms.untracked.double(0.0),
													 genParticleSource = cms.InputTag("genParticles"),
													 #std. model quarks + muon+neutrinos +SUSY, quarks + gluiono,chargino,neutralino,mmuin
													 genParticleList   = cms.vdouble(1,2,3,4,5,6,7,8,12,13,14,16,18,1000001,1000002,1000003,1000004, 1000005, 1000006, 1000012,1000013,1000014,1000016, 2000001,2000002,2000003,2000004,2000005, 2000006,2000013,1000021,1000022,1000023,1000024, 1000025, 1000035,1000037,1000039),
													 MCflag          = cms.bool (False)
)

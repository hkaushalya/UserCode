import FWCore.ParameterSet.Config as cms

mueff= cms.EDAnalyzer('MuonEff',
								#Debug              = cms.bool(False),
								VertexSource       = cms.InputTag('goodVertices'),
								DoPUReweight       = cms.bool(True),
								#FileNamePUDataDistribution = cms.string("RA2Classic/AdditionalInputFiles/DataPileupHistogram_RA2Summer12_190456-208686_ABCD.root"),   
								MuonSource        = cms.InputTag("patMuons"),
								recoJetAllSource      = cms.InputTag("newpatJetsPFchsPt30"),
								genJetAllSource   = cms.InputTag("ak5GenJets"),
								genParticleSource = cms.InputTag("genParticles"),
								MHTSource         = cms.InputTag("mhtPFchs"),
								HTSource          = cms.InputTag("htPFchs"),
                      	MinHT             = cms.untracked.double(0.0)   
)

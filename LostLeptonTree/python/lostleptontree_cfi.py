import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('LostLeptonTree',
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
									MCflag          = cms.bool (False),
									no_beamHaloFilter       = cms.bool(False),
									no_eeBadScFilter        = cms.bool(False),
									no_eeNoiseFilter        = cms.bool(False),
									no_greedyMuonsFilter    = cms.bool(False),
									no_hcalLaserEventFilter = cms.bool(False),
									no_inconsistentFilter   = cms.bool(False),
									no_ra2EcalBEFilter      = cms.bool(False),
									no_ra2EcalTPFilter      = cms.bool(False),
									no_trackingFailureFilter= cms.bool(False),
									no_HBHENoiseFilter      = cms.bool(False)
)
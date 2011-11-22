
import FWCore.ParameterSet.Config as cms

treeMaker = cms.EDFilter('FlatTreeMaker',
							jetSrc = cms.InputTag("patJetsAK5PF"),
							mhtSrc = cms.InputTag("mhtPF"),
							metSrc = cms.InputTag("pfMet"),
							htSrc = cms.InputTag("htPF"),
							muonSrc = cms.InputTag("cleanPatMuons"),
							eleSrc = cms.InputTag("cleanPatElectrons"),
							forVetoMuonSrc     = cms.InputTag("patMuonsPFIDIso"),
							forVetoElectronSrc = cms.InputTag("patElectronsPFIDIso"),
							doFillTree         = cms.bool(True),
							outRootName        = cms.string("StdRA2FlatTree.root"),
							genJetSrc = cms.InputTag("ak5GenJets"),
							genMETSrc = cms.InputTag("genMetTrue"),
							genParticleSrc = cms.InputTag("genParticles"),
							recoGenJetsDR = cms.double(0.5),
							vtxSrc = cms.InputTag("goodVerticesRA2"),
							nVtxPUcut = cms.uint32(1),
							debug = cms.bool(True),
							doSgnf = cms.bool(False)
)


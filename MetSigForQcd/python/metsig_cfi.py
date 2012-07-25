import FWCore.ParameterSet.Config as cms

metsig = cms.EDAnalyzer('MetSigForQcd',
                      patJetsPFPt50Eta25InputTag = cms.InputTag("patJetsAK5PFPt50Eta25"),
                      patJetsPFPt30Eta50InputTag = cms.InputTag("patJetsAK5PFPt30"),
                      mhtInputTag = cms.InputTag("mhtPF"),
                      htInputTag = cms.InputTag("htPF"),
							 dMinNjet = cms.untracked.double(3.0),
							 dMaxNjet = cms.untracked.double(100.0),
                      dMinHT = cms.untracked.double(0.0),   
                      dMinMHT = cms.untracked.double(0.0),   
                      dMinMetSig = cms.untracked.double(0.0),   
							 usePrescaleWeight = cms.untracked.int32(0),
							 prescaleWeight = cms.InputTag('trigWgtProd:prescaleWeight'),
							 mhtSigProducer = cms.InputTag('mymhtPFforSgnf'),
                      verbose = cms.untracked.int32(0),
                      ApplyLumiWeighing = cms.untracked.int32(0),
                      ApplyEventWeighing = cms.untracked.int32(0),
							 htBins = cms.vdouble(0,350, 500,800,1000,1200,1400,7000),
							 mhtBins = cms.vdouble(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,100,110,120,130,140,150,160,170,180,190,200,220,240,260,280,300,350,400,1000),
							 msigBins = cms.vdouble(0,10, 20, 30,40,50,60,70,80,90,100,110,120,130,140,160,180,200,220,240,260,280, 300, 320, 340, 360, 380,400, 420, 440, 460, 480, 500,1000),
							 msigAltBins = cms.vdouble(0,1, 2, 3,4,5,6,7,8,9,10,11,12,13,14,16,18,20,21,22,23,24,25,26,27,28,29,30,32, 34, 36, 38,40,45,50,60),
							 applyNjetCut = cms.untracked.int32(1),
							 applyHtCut   = cms.untracked.int32(1),
							 applyMhtCut  = cms.untracked.int32(1),
							 applyDphiCut = cms.untracked.int32(1)

)

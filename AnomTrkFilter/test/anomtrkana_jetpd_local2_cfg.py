import FWCore.ParameterSet.Config as cms

process = cms.Process("AnomTrk")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
		fileNames = cms.untracked.vstring(
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/830/AE8FE4D3-FCA2-E011-963F-0030487CD14E.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/830/A6200FF2-0CA3-E011-AE1C-003048F117EA.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/830/80D5E6C8-FCA2-E011-B3C6-003048F024C2.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/830/6050A3CC-FAA2-E011-94BC-0019B9F70468.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/830/56F35DCE-FAA2-E011-A2B7-0019B9F72BFF.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/830/540A2238-30A3-E011-B09B-001D09F2305C.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/830/508C2067-07A3-E011-A1C8-003048F118D4.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/830/4EAFAF79-2BA3-E011-8857-0030487C5CFA.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/830/4A0E584A-2AA3-E011-B110-003048F118D2.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/830/40714AC1-3FA3-E011-9334-BCAEC518FF6B.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/830/3087C39B-27A3-E011-A934-0019B9F72F97.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/830/2ECCF68F-1CA3-E011-912E-001D09F297EF.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/830/1CFEDB63-20A3-E011-B964-003048F11942.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/830/12BA2C61-20A3-E011-88D2-BCAEC5329701.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/830/08947ADD-03A3-E011-B49A-003048F1BF68.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/807/F20BAD10-B5A2-E011-AFC4-003048D2BBF0.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/807/D2E2FFAA-A7A2-E011-9B18-003048D2C0F2.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/807/C47E5B45-D3A2-E011-8EB1-BCAEC53296F7.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/807/B6D0B667-07A3-E011-B983-003048F024DE.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/807/AE35BEA0-BFA2-E011-AC85-003048D373AE.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/807/A8E72E76-AAA2-E011-A771-003048D375AA.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/807/A4CAFDF8-D3A2-E011-9F01-003048F11CF0.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/807/A044234C-9AA2-E011-9752-003048D2C01E.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/807/9A007639-D1A2-E011-9281-003048F1110E.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/807/92B99615-B5A2-E011-B517-003048D375AA.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/807/8E42FFDE-E4A2-E011-B0B5-003048F1183E.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/807/761EF3B7-DBA2-E011-8815-003048F1C420.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/807/6EF7D968-E1A2-E011-8B43-0030487C8E00.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/807/6EB6FC49-05A3-E011-9167-BCAEC5329713.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/807/6A4B3066-EDA2-E011-8B09-0030487CD76A.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/807/66AF4EB4-C6A2-E011-BFAA-0030487CD77E.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/807/66802C31-96A2-E011-9CA1-003048F11114.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/807/6053A2BD-F3A2-E011-B00D-0030487C5CFA.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/807/5A72CDE8-98A2-E011-ADC7-0030486780E6.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/807/44A5A38C-DEA2-E011-9181-003048F118C4.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/807/28CC5B30-96A2-E011-96D8-003048F11C5C.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/807/22D8A233-BEA2-E011-AFCF-003048D374F2.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/807/1C710F67-FEA2-E011-9480-003048F1BF66.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/807/0C8F66A0-BFA2-E011-B411-003048D2C108.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/786/64E79786-11A2-E011-BB43-003048F117EA.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/786/5AF28097-02A2-E011-BE56-003048F24A04.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/786/3EE3C285-11A2-E011-9257-003048F0258C.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/785/C0BBA30B-C4A1-E011-8BBE-001D09F29114.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/784/225722F5-03A2-E011-B400-BCAEC5329719.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/784/18EF24F6-03A2-E011-9016-BCAEC5329730.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/754/52E8FBBF-A3A1-E011-8FE4-003048F1BF66.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/746/F6450E5D-A2A1-E011-8908-BCAEC518FF50.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/746/7CDE79CC-A3A1-E011-B98D-003048D2C1C4.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/746/5A1014CD-A3A1-E011-A954-003048D3756A.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/746/3AB8EFBF-DEA1-E011-8111-003048F118D4.root'
			)
)

#secFiles.extend( [
#] )


process.TFileService = cms.Service("TFileService",
		fileName = cms.string('METFilter_JETPD2.root')
		)

process.out = cms.OutputModule("PoolOutputModule",
		fileName = cms.untracked.string('TrkFiltOutEventsFromJETPD2.root'),
		# save only events passing the full path
		SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
		# save PAT Layer 1 output; you need a '*' to
		# unpack the list of commands 'patEventContent'
		# outputCommands = cms.untracked.vstring('drop *', *patEventContent ) #To further strip down the eventsize
)

#### remove beam scraping events
process.noScraping= cms.EDFilter("FilterOutScraping",
						applyfilter = cms.untracked.bool(True),
						debugOn = cms.untracked.bool(False), ## Or 'True' to get some per-event info
						numtrack = cms.untracked.uint32(10),
						thresh = cms.untracked.double(0.20)
)



process.anomtrk = cms.EDFilter('AnomTrkMETFilter',
							trackCollection = cms.InputTag('generalTracks'),
							caloJetInputTag_ = cms.InputTag("ak5CaloJets"),
							pfJetInputTag_   = cms.InputTag("ak5PFJets"),
							HltTriggerResults         = cms.InputTag("TriggerResults::HLT"),
							req_trigger = cms.untracked.bool(True),
							TriggerPathsToStore = cms.vstring("HLT_Jet300_v5"),
							nMinVtx = cms.untracked.int32(1),	
							ndofVtx = cms.untracked.int32(4),	
							maxDelzTrkVtx = cms.untracked.double(12.0),	
							minMet = cms.untracked.double(0.0)	
)

process.myfilter_path = cms.Path(process.noScraping +
											process.anomtrk
)

#process.out.SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('myfilter_path') )

#process.outpath = cms.EndPath(process.out)


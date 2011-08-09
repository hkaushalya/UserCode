import FWCore.ParameterSet.Config as cms

process = cms.Process("AnomTrk")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
		fileNames = cms.untracked.vstring(
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/168/437/7A9E0623-12A6-E011-A77E-BCAEC518FF89.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/168/425/C650DEF3-1AA6-E011-98F7-BCAEC53296F9.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/168/423/D2F46030-1AA6-E011-A2B6-003048D373AE.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/168/229/0E1D9C74-C9A4-E011-9AC5-003048D37580.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/168/207/1A1BEF33-BCA4-E011-A0A2-003048F1BF66.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/168/162/7EAD54E5-AEA4-E011-A56E-003048F1183E.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/168/155/B8D474B9-90A4-E011-B764-BCAEC53296F3.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/168/152/4ED883CA-90A4-E011-A2F0-BCAEC518FF80.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/969/BA183411-99A3-E011-A6E8-003048F1182E.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/913/DAA6D645-96A3-E011-9C36-003048F118C6.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/913/C62B0C0A-8BA3-E011-B0BD-003048F024FE.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/913/761EC12D-8DA3-E011-86D2-003048F11C58.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/913/7241F93A-88A3-E011-9F02-BCAEC5329729.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/913/66B3F987-8EA3-E011-8D81-003048F24A04.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/913/28A392F5-8FA3-E011-95CB-003048F11114.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/913/18058AD8-86A3-E011-AE2B-0030486780E6.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/898/F80508D9-86A3-E011-92B2-BCAEC53296FB.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/898/EC9857C3-57A3-E011-AD17-485B3989721B.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/898/E4D7A499-8EA3-E011-BEB4-003048F1C832.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/898/E488D61A-6CA3-E011-9FA0-0030487C90EE.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/898/E01EABC5-78A3-E011-9EBF-001D09F2423B.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/898/DC5E7753-56A3-E011-9323-BCAEC5329727.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/898/DA2D1BD0-59A3-E011-8DF5-BCAEC532971B.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/898/BC85AC9D-61A3-E011-BD71-BCAEC53296F6.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/898/A6E6A2E6-81A3-E011-96EF-003048F1183E.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/898/98A86F94-53A3-E011-9347-E0CB4E4408D2.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/898/92826E65-69A3-E011-BA26-BCAEC53296F8.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/898/84CB343B-88A3-E011-B8FC-003048F1BF68.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/898/6C55155F-67A3-E011-845E-001D09F24FEC.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/898/504D9C0A-8BA3-E011-A8EC-003048F11DE2.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/898/448F2F60-67A3-E011-929B-001D09F25217.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/898/34CB4B9F-53A3-E011-B480-BCAEC5329728.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/898/34C38B50-5BA3-E011-9DA7-BCAEC532971D.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/898/2CED2B6E-85A3-E011-B2AE-003048D3750A.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/898/2662251B-4DA3-E011-9C69-BCAEC5364C93.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/898/1EBDE412-A7A3-E011-BF5E-003048F11942.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/898/1C37BF9C-61A3-E011-99CD-E0CB4E5536AE.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/898/0A5CF485-4EA3-E011-8488-BCAEC5329729.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/830/FA3915CD-FAA2-E011-B9C3-001D09F25438.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/830/E88E6330-2EA3-E011-9827-001D09F254CE.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/830/E442F451-F9A2-E011-BE87-003048D2C01A.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/830/D2EFA62E-2EA3-E011-94CF-001D09F2462D.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/830/CA9E5204-01A3-E011-B85B-003048F118C6.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/830/CA38297C-2CA3-E011-81E6-001D09F2AF1E.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/830/C4921E22-2EA3-E011-8B89-003048D2BB58.root',
			'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/830/BEE38DDD-03A3-E011-80FF-003048F118AC.root'
			)
)

#secFiles.extend( [
#] )


process.TFileService = cms.Service("TFileService",
		fileName = cms.string('METFilter_JETPD1.root')
		)

process.out = cms.OutputModule("PoolOutputModule",
		fileName = cms.untracked.string('TrkFiltOutEventsFromJETPD1.root'),
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


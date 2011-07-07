import FWCore.ParameterSet.Config as cms

process = cms.Process("AnomTrk")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
   fileNames = cms.untracked.vstring(
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/168/437/6C9AB0EA-11A6-E011-A2B4-003048D374F2.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/168/425/1E6CC011-1BA6-E011-88D6-003048F1110E.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/168/423/9AD30D23-1AA6-E011-A34B-BCAEC53296F6.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/168/229/8295E773-C9A4-E011-87D4-BCAEC518FF68.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/168/207/4C83C333-BCA4-E011-BABE-003048F1C832.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/168/162/E640B7DC-AEA4-E011-82AC-BCAEC518FF91.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/168/155/58C7EEC0-90A4-E011-A987-BCAEC5329702.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/168/152/6ED60EA9-90A4-E011-933B-003048F024C2.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/969/B0B59111-99A3-E011-87D2-BCAEC518FF5F.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/913/EC5F2845-8EA3-E011-AD99-BCAEC518FF69.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/913/CAED0406-9EA3-E011-85A3-BCAEC5329730.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/913/AEE86FE7-81A3-E011-905B-003048F024FA.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/913/86754BD9-86A3-E011-89D0-003048678098.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/898/FEAF6255-69A3-E011-88E8-BCAEC53296F4.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/898/FA21D047-6EA3-E011-8E64-001D09F24EE3.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/898/E83CE1D2-46A3-E011-96AE-E0CB4E4408E3.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/898/D487DB47-6EA3-E011-9D0A-001D09F2545B.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/898/C8F0D009-50A3-E011-B840-BCAEC5329716.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/898/C86AD181-80A3-E011-97F9-BCAEC518FF63.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/898/BE6EC607-84A3-E011-8BC3-003048F1C420.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/898/94932348-6EA3-E011-8E28-001D09F295A1.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/898/804E1078-64A3-E011-AB50-BCAEC518FF8F.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/898/7ED776BC-7DA3-E011-AEF6-0030487C8E00.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/898/7A0C70E7-41A3-E011-8782-BCAEC5329702.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/898/7203A11B-4DA3-E011-BF88-BCAEC5329707.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/898/6857B4CF-59A3-E011-A72B-BCAEC518FF8E.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/898/5E5E5509-50A3-E011-8AFA-BCAEC518FF8D.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/898/5E5CEC14-4BA3-E011-876B-BCAEC5364CFB.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/898/24F354B7-6AA3-E011-A80F-003048F117B6.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/898/1407F677-64A3-E011-8206-BCAEC5329721.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/898/048D43C6-57A3-E011-8D7C-BCAEC518FF80.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/830/E42A8D8F-1CA3-E011-8138-E0CB4E55367F.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/830/C6503150-05A3-E011-9573-003048F01E88.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/830/969DDE49-05A3-E011-AC10-003048F1C836.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/830/8A361B86-13A3-E011-9BE3-BCAEC518FF69.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/830/7C59C7D8-F0A2-E011-A974-003048D37560.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/830/68DC7E10-15A3-E011-B11A-485B39897227.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/830/582FD079-12A3-E011-BBCA-BCAEC5329700.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/830/4C30692A-29A3-E011-A9FD-001D09F2AD84.root',
		'/store/data/Run2011A/MET/RECO/PromptReco-v4/000/167/830/1C60BB41-F2A2-E011-835D-0030487CD6D2.root'
    )
)

#secFiles.extend( [
#] )


process.TFileService = cms.Service("TFileService",
		fileName = cms.string('METFilter_METPD.root')
		)

process.out = cms.OutputModule("PoolOutputModule",
		fileName = cms.untracked.string('TrkFiltOutEventsFromMETPD.root'),
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
							HltTriggerResults         = cms.InputTag("TriggerResults::HLT"),
							req_trigger = cms.untracked.bool(False),
							TriggerPathsToStore = cms.vstring("HLT_MET65_v3"),
							nMinVtx = cms.untracked.int32(1),	
							minMet = cms.untracked.double(1000.0)	
)

process.myfilter_path = cms.Path(process.noScraping +
											process.anomtrk
)

process.out.SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('myfilter_path') )

process.outpath = cms.EndPath(process.out
)


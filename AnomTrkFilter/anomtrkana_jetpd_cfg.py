import FWCore.ParameterSet.Config as cms

process = cms.Process("AnomTrk")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
   fileNames = cms.untracked.vstring(
		'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/785/C0BBA30B-C4A1-E011-8BBE-001D09F29114.root',
		'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/754/52E8FBBF-A3A1-E011-8FE4-003048F1BF66.root',
		'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/740/60347DAF-74A1-E011-B510-485B3989721B.root',
		'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/740/385F3474-69A1-E011-8CE8-BCAEC5364C42.root',
		'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/715/6C7500E2-26A1-E011-A6D8-003048D3750A.root',
		'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/696/56FBBAD0-26A1-E011-AEA9-003048F1C58C.root',
		'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/676/DC36AF3E-3FA1-E011-9C35-001D09F25401.root',
		'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/676/6474B07E-47A1-E011-B30D-001D09F2516D.root',
		'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/676/481490D9-3DA1-E011-BEAD-001D09F24F65.root',
		'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/676/3E65DBDA-31A1-E011-BA66-003048D374F2.root',
		'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/676/1829ECD9-3DA1-E011-BC0F-001D09F2426D.root',
		'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/675/F811DEC6-E1A0-E011-B9D6-0030487CD716.root',
		'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/675/EE8A93E3-DCA0-E011-897A-003048F024E0.root',
		'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/675/C2E752CA-FBA0-E011-B703-00304879BAB2.root',
		'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/675/BCA05C06-F9A0-E011-8B22-0030487A3DE0.root',
		'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/675/B6C7A28B-16A1-E011-82E2-003048D2C108.root',
		'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/675/B47FF50A-DAA0-E011-8B81-0030487C60AE.root',
		'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/675/A8101BB2-D8A0-E011-A826-0030487CD6B4.root',
		'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/675/92495DEA-EAA0-E011-9432-001D09F2441B.root',
		'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/675/7E5C5E87-16A1-E011-85C0-003048D2C16E.root',
		'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/675/70B6E6A6-18A1-E011-A212-003048D2BC42.root',
		'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/675/4C1E9808-F9A0-E011-8B19-0030487CD14E.root',
		'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/675/1EC443E3-DCA0-E011-90F1-003048F024F6.root',
		'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/675/0CE9C2F9-DEA0-E011-BD86-001D09F2527B.root',
		'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/674/F83EBEF5-A5A0-E011-95DA-001D09F28F11.root',
		'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/674/D217F9B3-CCA0-E011-8FCD-003048D375AA.root',
		'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/674/68BF54C6-A8A0-E011-9203-00304879FBB2.root',
		'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/674/541C79B2-9FA0-E011-BC89-003048F118E0.root',
		'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/674/1401708D-ABA0-E011-BA46-001D09F24399.root',
		'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/673/D80EE6AE-9CA0-E011-9051-001D09F231B0.root',
		'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/673/C6A958AE-9CA0-E011-80A0-003048F118C6.root',
		'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/673/90121D17-97A0-E011-90D9-001D09F244BB.root',
		'/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/167/673/46E1CC02-95A0-E011-BE30-001D09F241B9.root'
    )
)

#secFiles.extend( [
#] )


process.TFileService = cms.Service("TFileService",
		fileName = cms.string('METFilter_JETPD.root')
		)

process.out = cms.OutputModule("PoolOutputModule",
		fileName = cms.untracked.string('TrkFiltOutEventsFromJETPD.root'),
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
							req_trigger = cms.untracked.bool(True),
							TriggerPathsToStore = cms.vstring("HLT_Jet300_v5"),
							nMinVtx = cms.untracked.int32(1),	
							minMet = cms.untracked.double(0.0)	
)

process.myfilter_path = cms.Path(process.noScraping +
											process.anomtrk
)

process.out.SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('myfilter_path') )

process.outpath = cms.EndPath(process.out
)


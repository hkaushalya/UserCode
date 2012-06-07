
import FWCore.ParameterSet.Config as cms

genPtHatFilter = cms.EDFilter('GenPtHatFilter',
							dMinPtHat = cms.untracked.double(0.0)
)



import FWCore.ParameterSet.Config as cms

runfilter = cms.EDFilter(
							"RunFilter",
                      lowRunNumber = cms.double(-1.0),   
                      hiRunNumber  = cms.double(99999999)
 )


import FWCore.ParameterSet.Config as cms

process = cms.EDAnalyzer('RPCInclusionAnalyzer')


process.RPCInclusionAnalyzer = cms.EDAnalyzer('RPCInclusionAnalyzer',
	dtag = cms.string('csctf'),
	isMC = cms.untracked.bool(False)
)

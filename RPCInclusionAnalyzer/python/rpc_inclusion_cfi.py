import FWCore.ParameterSet.Config as cms

process = cms.Process('RPCInclusionAnalyzer')


process.RPCInclusionAnalyzer = cms.EDAnalyzer('RPCInclusionAnalyzer',
	dtag = cms.string('csctf'),
	isMC = cms.bool(False)
)

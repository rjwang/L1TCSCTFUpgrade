import FWCore.ParameterSet.Config as cms

process = cms.Process('RPCInclusionAnalyzer')


process.RPCInclusionAnalyzer = cms.EDAnalyzer('RPCInclusionAnalyzer',
	dtag = cms.string('csctf'),
	RPCTPTag = cms.InputTag("L1TMuonTriggerPrimitives","RPC"),
	CSCTFTag = cms.InputTag("L1CSCTFTrackConverter"),
	isMC = cms.bool(False)
)

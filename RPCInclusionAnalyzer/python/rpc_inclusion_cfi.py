import FWCore.ParameterSet.Config as cms

process = cms.Process('RPCInclusionAnalyzer')


process.RPCInclusionAnalyzer = cms.EDAnalyzer('RPCInclusionAnalyzer',
	dtag = cms.string('csctf'),
	GenParticles = cms.untracked.InputTag("genParticles"),
	RPCTPTag = cms.InputTag("L1TMuonTriggerPrimitives","RPC"),
	CSCTFTag = cms.InputTag("L1CSCTFTrackConverter"),
	lutParam     = cms.PSet(
		isBeamStartConf = cms.untracked.bool(True),
		ReadPtLUT = cms.bool(False),
		PtMethod = cms.untracked.uint32(32)
	),
	isMC = cms.bool(False)
)

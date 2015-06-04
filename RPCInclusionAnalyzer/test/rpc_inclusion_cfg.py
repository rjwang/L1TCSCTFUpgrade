import FWCore.ParameterSet.Config as cms


from L1TCSCTFUpgrade.RPCInclusionAnalyzer.rpc_inclusion_cfi import *



process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True),
                                        SkipEvent = cms.untracked.vstring('ProductNotFound')
                                        )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	'/store/data/Run2015A/Commissioning/RAW/v1/000/246/926/00000/6A496D1A-F709-E511-B2B9-02163E0145D2.root',
    )
)



process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("csctf_Run246926.root")
                                  )



process.p = cms.Path(process.RPCInclusionAnalyzer)

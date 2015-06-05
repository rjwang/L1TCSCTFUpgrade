import FWCore.ParameterSet.Config as cms


from L1TCSCTFUpgrade.RPCInclusionAnalyzer.rpc_inclusion_cfi import *



process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True),
                                        SkipEvent = cms.untracked.vstring('ProductNotFound')
                                        )

#load run conditions
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')


# Global Tag to specify the "good" events to run over
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	'/store/data/Run2012D/SingleMu/RAW-RECO/ZMu-22Jan2013-v1/10000/1010C4E3-A1A7-E211-B288-00259074AE9A.root',
    )
)



process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("csctf_test.root")
                                  )



process.p = cms.Path(process.RPCInclusionAnalyzer)

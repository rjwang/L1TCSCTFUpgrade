import FWCore.ParameterSet.Config as cms

from L1TCSCTFUpgrade.RPCInclusionAnalyzer.rpc_inclusion_cfi import *

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('L1TriggerDPGUpgrade.L1TMuon.L1TMuonTriggerPrimitiveProducer_cfi')
process.load('L1TriggerDPGUpgrade.L1TMuon.L1CSCTFTrackConverter_cfi')
process.load('L1TriggerDPGUpgrade.L1TMuon.L1DTTFTrackConverter_cfi')
process.load('L1TriggerDPGUpgrade.L1TMuon.L1RPCTFTrackConverter_cfi')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')

# load reco process
process.load("RecoMuon.TrackingTools.MuonServiceProxy_cff")
process.load("RecoMuon.TrackingTools.MuonTrackLoader_cff")
process.load('Configuration/StandardSequences/SimL1Emulator_cff')
process.load("Configuration.StandardSequences.RawToDigi_Data_cff")
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load("Configuration.StandardSequences.Generator_cff")

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


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	'/store/user/rewang/csctf/SingleMuLowPt_5GeVto200GeV_GEN_SIM_DIGI_L1_4_9_2015.root',
    )
)



process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("csctf_mc.root")
                                  )




#-----------Final Process -------------------------------------------------------------------------

process.L1TMuonSeq = cms.Sequence(
                                   process.L1TMuonTriggerPrimitives +
                                   process.L1CSCTFTrackConverter    +
				   process.L1DTTFTrackConverter     +
				   process.L1RPCTFTrackConverters   +
                                   process.RPCInclusionAnalyzer
                                   )
#---------------------------------------------------------------------------------------------------

process.L1TMuonPath = cms.Path(process.L1TMuonSeq)
process.Schedule = cms.Schedule(process.L1TMuonPath)

#process.p = cms.Path(process.RPCInclusionAnalyzer)

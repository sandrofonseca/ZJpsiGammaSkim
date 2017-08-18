#input_filename = '/store/data/Run2017A/MuonEG/MINIAOD/PromptReco-v3/000/296/888/00000/446CB6BD-6A55-E711-89DF-02163E011A9C.root'
input_filename = '/store/data/Run2017B/MuonEG/MINIAOD/PromptReco-v1/000/297/101/00000/F6C53CDD-2557-E711-8F55-02163E01A6D8.root' 
ouput_filename = 'rootuple.root'

import FWCore.ParameterSet.Config as cms
process = cms.Process("TagProbeZJpsiGamma")

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '92X_dataRun2_Prompt_v8', '')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 2

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(20))
process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(input_filename))
process.TFileService = cms.Service("TFileService",fileName = cms.string(ouput_filename))
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False))
######## HLT Trigger ########################################################################### 
process.triggerSelection = cms.EDFilter("TriggerResultsFilter",
                                        triggerConditions = cms.vstring('HLT_DoubleMu20_7_Mass0to30_Photon23_v*'
                                                                       ),
                                        hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                        l1tResults = cms.InputTag( "" ),
                                        throw = cms.bool(False)
)
############################################################################################





process.TPSequence = cms.Sequence(
                                   process.triggerSelection #*
)

process.p = cms.Path(process.TPSequence)




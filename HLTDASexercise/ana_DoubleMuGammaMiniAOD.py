import FWCore.ParameterSet.Config as cms

process = cms.Process('HLTANALYZER')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
#   fileNames = cms.untracked.vstring('/store/data/Run2016G/SingleMuon/MINIAOD/23Sep2016-v1/1110000/A2C0F697-B19C-E611-A4D8-F04DA275BFF2.root'),
   fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/s/sfonseca/ZJpsi_Gamma/TagProbeTrigger2017/CMSSW_9_2_7/src/ZJpsiGammaSkim/PhotonSelectionProducer/dimugamma_skim10k.root'), # skimmed file on EOS at LPC
)


process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('RelVal nevts:100'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('histos_DoubleMuGammaTrigAnalyzer.root')
                                   )

### analyzer configuration

process.doubleMuGammaTrigAnalyzerMiniAOD = cms.EDAnalyzer("DoubleMuGammaTrigAnalyzerMiniAOD",
        JPsiFiltered = cms.InputTag("ZmmgDimuons"), 
        photonInputTag = cms.InputTag("ZmmgDimuonGamma"),
        JPsiGammaFiltered = cms.InputTag("ZmmgDimuonGamma")

)

process.doubleMuGammaTrigAnalyzerMiniAOD.triggerName = cms.untracked.string("HLT_DoubleMu20_7_Mass0to30_Photon23_v1")
process.doubleMuGammaTrigAnalyzerMiniAOD.verbose = cms.untracked.bool(True)

process.GlobalTag.globaltag = "92X_dataRun2_Prompt_v8"

# Path and EndPath definitions
process.HLTanalyzers = cms.Path(process.doubleMuGammaTrigAnalyzerMiniAOD)

import FWCore.ParameterSet.Config as cms

process = cms.Process('SKIM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)
"""
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                  calibratedPatPhotons = cms.PSet(
    initialSeed = cms.untracked.uint32(12345),
    engineName = cms.untracked.string('TRandom3')
    )
"""
# Input source
process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
      # 'file:446CB6BD-6A55-E711-89DF-02163E011A9C.root'
        '/store/data/Run2016G/SingleMuon/MINIAOD/23Sep2016-v1/1110000/A2C0F697-B19C-E611-A4D8-F04DA275BFF2.root'
   ),
)


process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('RelVal nevts:100'),
    name = cms.untracked.string('Applications')
)

# JSON

import FWCore.PythonUtilities.LumiList as LumiList
#process.source.lumisToProcess = LumiList.LumiList(filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-297723_13TeV_PromptReco_Collisions17_JSON.txt').getVLuminosityBlockRange()

process.source.lumisToProcess = LumiList.LumiList(filename = 'Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt').getVLuminosityBlockRange()



# skim definitions
###____________________________________________________________________________
###
###  HLT Filter
###____________________________________________________________________________

import copy
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
process.ZmmgHLTFilter = copy.deepcopy(hltHighLevel)
process.ZmmgHLTFilter.throw = cms.bool(False)
process.ZmmgHLTFilter.HLTPaths = ['HLT_*']


### Get muons of needed quality for Z -> mumugamma 
process.ZmmgTrailingMuons = cms.EDFilter('MuonSelector',
    src = cms.InputTag('slimmedMuons'),
    cut = cms.string('''pt > 7.0 && 
                        abs(eta) < 2.4 && 
                        isGlobalMuon = 1 && 
                        isTrackerMuon = 1 '''),
    filter = cms.bool(True)                                
    )

### Require a harder pt cut on the leading leg
process.ZmmgLeadingMuons = cms.EDFilter('MuonSelector',
    src = cms.InputTag('ZmmgTrailingMuons'),
    cut = cms.string('pt > 20'),
    filter = cms.bool(True)                                
    )


### Build dimuon candidates
process.ZmmgDimuons = cms.EDProducer('CandViewShallowCloneCombiner',
    decay = cms.string('ZmmgLeadingMuons@+ ZmmgTrailingMuons@-'),
    checkCharge = cms.bool(True),
    cut = cms.string('mass > 3.0 && mass < 3.2')
    )

### Require at least one dimuon candidate
process.ZmmgDimuonFilter = cms.EDFilter('CandViewCountFilter',
    src = cms.InputTag('ZmmgDimuons'),
    minNumber = cms.uint32(1)
)



# PAT Photon #######################################################
#https://github.com/cmkuo/ggAnalysis/blob/a28309dbabc16b9d27bb5c015d15c369426e1bab/ggNtuplizer/plugins/ggNtuplizer.cc#L79
process.calibratedPatPhotons = cms.EDProducer("CalibratedPatPhotonProducerRun2",
    correctionFile = cms.string('ScalesSmearings/Moriond17_23Jan_ele'),# Moriond17
    photons = cms.InputTag("slimmedPhotons"),
    isMC = cms.bool(True),
    isSynchronization = cms.bool(False)
)


### Select photon candidates with Et > 23.0 GeV
process.ZmmgPhotons = cms.EDFilter('CandViewSelector',
    src = cms.InputTag('calibratedPatPhotons'),
    cut = cms.string('et > 23.0'),
    filter = cms.bool(True)
)

process.dimugamma_filter_step = cms.Path(process.ZmmgHLTFilter + process.ZmmgTrailingMuons + process.ZmmgLeadingMuons + process.ZmmgDimuons + process.ZmmgDimuonFilter)#  + process.calibratedPatPhotons + process.ZmmgPhotons )

#process.dimugamma_filter_step = cms.Path(process.ZmmgTrailingMuons + process.ZmmgLeadingMuons + process.ZmmgDimuons + process.ZmmgDimuonFilter)#  + process.calibratedPatPhotons + process.ZmmgPhotons )


# Output definition

process.Out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string ("dimugamma_skim.root"),
                               outputCommands = cms.untracked.vstring('keep *'),
                               SelectEvents = cms.untracked.PSet(
                                   SelectEvents = cms.vstring('dimugamma_filter_step')
                               ),
)
  
process.GlobalTag.globaltag = "80X_dataRun2_2016SeptRepro_v3"


# Path and EndPath definitions
process.output_step = cms.EndPath(process.Out)

#process.load('FWCore.MessageService.MessageLogger_cfi')
#process.MessageLogger.cerr.FwkReport.reportEvery = 1000000

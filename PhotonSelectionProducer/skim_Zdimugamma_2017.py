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
    input = cms.untracked.int32(1000)
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
     #  'file:446CB6BD-6A55-E711-89DF-02163E011A9C.root'
        '/store/data/Run2017B/MuonEG/MINIAOD/PromptReco-v1/000/297/101/00000/F6C53CDD-2557-E711-8F55-02163E01A6D8.root'
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
process.source.lumisToProcess = LumiList.LumiList(filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-299420_13TeV_PromptReco_Collisions17_JSON.txt').getVLuminosityBlockRange()

#process.source.lumisToProcess = LumiList.LumiList(filename = 'Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt').getVLuminosityBlockRange()



# skim definitions
###____________________________________________________________________________
###
###  HLT Filter
###____________________________________________________________________________

import copy
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
process.ZmmgHLTFilter = copy.deepcopy(hltHighLevel)
process.ZmmgHLTFilter.throw = cms.bool(False)
process.ZmmgHLTFilter.HLTPaths = ['HLT_DoubleMu20_7_Mass0to30_Photon23_v*']


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

#process.calibratedPatPhotons = cms.EDProducer("CalibratedPatPhotonProducerRun2",
#    correctionFile = cms.string('ScalesSmearings/Moriond17_23Jan_ele'),# Moriond17
#    photons = cms.InputTag("slimmedPhotons"),
#    isMC = cms.bool(True),
#    isSynchronization = cms.bool(False)
#)


process.ZmmgDimuonGamma = cms.EDProducer("PhotonSelectionProducer",
        JPsiFiltered = cms.InputTag("ZmmgDimuons"),        
        photonSrc      = cms.InputTag("slimmedPhotons"),
        #calibphotonSrc = cms.InputTag("calibratedPatPhotons"),
        #ebReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedEBRecHits"),
        #eeReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedEERecHits"),
        #esReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedESRecHits"),
        #recoPhotonSrc = cms.InputTag("reducedEgamma", "reducedGedPhotonCores"),
        triggerMatch = cms.bool(False) # trigger match is performed in ZmmgHLTFilter 

)
'''
process.ZmmgDimuonGammaFilter = cms.EDFilter('CandViewCountFilter',
    src = cms.InputTag('PatPhotons'),
    minNumber = cms.uint32(1)
)
'''


### Select photon candidates with Et > 23.0 GeV
#process.ZmmgPhotons = cms.EDFilter('CandViewSelector',
#    src = cms.InputTag('process.PatPhotons'),
#    cut = cms.string('et > 23.0'),
#    filter = cms.bool(True)
#)

process.dimugamma_filter_step = cms.Path(process.ZmmgHLTFilter + process.ZmmgTrailingMuons + process.ZmmgLeadingMuons + process.ZmmgDimuons + process.ZmmgDimuonFilter + process.ZmmgDimuonGamma)# + process.ZmmgDimuonGammaFilter)

#process.dimugamma_filter_step = cms.Path(process.ZmmgTrailingMuons + process.ZmmgLeadingMuons + process.ZmmgDimuons + process.ZmmgDimuonFilter)


# Output definition

process.Out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string ("dimugamma_skim.root"),
                               outputCommands = cms.untracked.vstring('keep *'),
                               SelectEvents = cms.untracked.PSet(
                                   SelectEvents = cms.vstring('dimugamma_filter_step')
                               ),
)
  
process.GlobalTag.globaltag = "92X_dataRun2_Prompt_v8" # ??????


# Path and EndPath definitions
process.output_step = cms.EndPath(process.Out)

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

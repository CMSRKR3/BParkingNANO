import FWCore.ParameterSet.Config as cms

paths=['HLT_DoubleEle10_eta1p22_mMax6',
       'HLT_DoubleEle9p5_eta1p22_mMax6',
       'HLT_DoubleEle9_eta1p22_mMax6',
       'HLT_DoubleEle8p5_eta1p22_mMax6',
       'HLT_DoubleEle8_eta1p22_mMax6',
       'HLT_DoubleEle7p5_eta1p22_mMax6',
       'HLT_DoubleEle7_eta1p22_mMax6',
       'HLT_DoubleEle6p5_eta1p22_mMax6',
       'HLT_DoubleEle6_eta1p22_mMax6',
       'HLT_DoubleEle5p5_eta1p22_mMax6',
       'HLT_DoubleEle5_eta1p22_mMax6',
       'HLT_DoubleEle4p5_eta1p22_mMax6',
       'HLT_DoubleEle4_eta1p22_mMax6'
]
paths_OR = " || ".join([ 'path( "{:s}_v*" )'.format(path) for path in paths])

seeds = ['L1_DoubleEG11_er1p2_dR_Max0p6',
         'L1_DoubleEG10p5_er1p2_dR_Max0p6',
         'L1_DoubleEG10_er1p2_dR_Max0p6',
         'L1_DoubleEG9p5_er1p2_dR_Max0p6',
         'L1_DoubleEG9_er1p2_dR_Max0p7',
         'L1_DoubleEG8p5_er1p2_dR_Max0p7',
         'L1_DoubleEG8_er1p2_dR_Max0p7',
         'L1_DoubleEG7p5_er1p2_dR_Max0p7',
         'L1_DoubleEG7_er1p2_dR_Max0p8',
         'L1_DoubleEG6p5_er1p2_dR_Max0p8',
         'L1_DoubleEG6_er1p2_dR_Max0p8',
         'L1_DoubleEG5p5_er1p2_dR_Max0p8',
         'L1_DoubleEG5_er1p2_dR_Max0p9',
         'L1_DoubleEG4p5_er1p2_dR_Max0p9',
         'L1_DoubleEG4_er1p2_dR_Max0p9',
]

# https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatAlgos/plugins/PATTriggerObjectStandAloneUnpacker.cc
myUnpackedPatTrigger = cms.EDProducer(
    "PATTriggerObjectStandAloneUnpacker",
    patTriggerObjectsStandAlone = cms.InputTag("slimmedPatTrigger"),
    triggerResults = cms.InputTag("TriggerResults::HLT"),
    unpackFilterLabels = cms.bool(True),
)

# https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatAlgos/python/triggerLayer1/triggerMatcherExamples_cfi.py
# https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatAlgos/plugins/PATTriggerMatcher.cc
myTriggerMatches = cms.EDProducer(
    "PATTriggerMatcherDEtaLessByDR", # match by DeltaEta only, best match by DeltaR
    #"PATTriggerMatcherDEtaLessByDEta", # match by DeltaEta only, best match by DeltaEta
    #"PATTriggerMatcherDRDPtLessByR", # match by DeltaR only, best match by DeltaR
    src = cms.InputTag("slimmedElectrons"),
    matched = cms.InputTag("myUnpackedPatTrigger"),
    matchedCuts = cms.string(paths_OR), # e.g. 'path("HLT_DoubleEle6_eta1p22_mMax6_v*")'
    maxDeltaR = cms.double(2.0),
    maxDeltaEta = cms.double(0.5),
    #maxDPtRel = cms.double(0.5),
    resolveAmbiguities    = cms.bool( True ), # only one match per trigger object
    resolveByMatchQuality = cms.bool( True ), # take best match found per reco object (e.g. by DeltaR)
)

# https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatAlgos/plugins/PATTriggerMatchEmbedder.cc
mySlimmedElectronsWithEmbeddedTrigger = cms.EDProducer(
    "PATTriggerMatchElectronEmbedder",
    src = cms.InputTag("slimmedElectrons"),
    matches = cms.VInputTag('myTriggerMatches'),
)

electronTrgSelector = cms.EDProducer(
    "ElectronTriggerSelector",
    electronCollection = cms.InputTag("mySlimmedElectronsWithEmbeddedTrigger"),
    bits = cms.InputTag("TriggerResults","","HLT"),
    prescales = cms.InputTag("patTrigger"),
    objects = cms.InputTag("slimmedPatTrigger"),
    vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
    maxdR_matching = cms.double(10.), # not used
    dzForCleaning_wrtTrgElectron = cms.double(1.),
    filterElectron = cms.bool(True),
    ptMin = cms.double(2.),
    absEtaMax = cms.double(1.25),
    HLTPaths=cms.vstring(paths),
    L1seeds=cms.vstring(seeds),
)

countTrgElectrons = cms.EDFilter(
    "PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("electronTrgSelector", "trgElectrons"),
)

from FWCore.ParameterSet.VarParsing import VarParsing
import FWCore.ParameterSet.Config as cms

options = VarParsing('python')

options.register('isMC', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('globalTag', 'NOTSET',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Set global tag"
)
options.register('wantSummary', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('wantFullRECO', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('reportEvery', 10,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "report every N events"
)
options.register('skip', 0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "skip first N events"
)
options.register('lhcRun', 3,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "LHC Run 2 or 3 (default)"
)

options.setDefault('maxEvents', 1000)
options.setDefault('tag', '124X')
options.parseArguments()
print(options)

globaltag = None
if   options.lhcRun == 3: globaltag = '124X_mcRun3_2022_realistic_v11' if options.isMC else '124X_dataRun3_Prompt_v4'
elif options.lhcRun == 2: globaltag = '102X_upgrade2018_realistic_v15' if options.isMC else '102X_dataRun2_v11'
if options._beenSet['globalTag']: globaltag = options.globalTag

extension = {False : 'data', True : 'mc'}
outputFileNANO = cms.untracked.string('_'.join(['BParkNANO', extension[options.isMC], options.tag])+'.root')
outputFileFEVT = cms.untracked.string('_'.join(['BParkFullEvt', extension[options.isMC], options.tag])+'.root')
if not options.inputFiles:
    if options.lhcRun == 2:
        options.inputFiles = [
            'file:///eos/cms/store/cmst3/group/bpark/BToKmumu_1000Events_MINIAOD.root'
        ] if options.isMC else [
            'root://cms-xrd-global.cern.ch//store/data/Run2018B/ParkingBPH4/MINIAOD/05May2019-v2/230000/6B5A24B1-0E6E-504B-8331-BD899EB60110.root'
        ]
    elif options.lhcRun == 3:
        options.inputFiles = [
            'root://cms-xrd-global.cern.ch//store/user/jodedra/BuTOjpsiKEE20221103FIFTYMminiaod/BuTOjpsiKEE20221103FIFTYM/SUMMER22_MINIAOD/221106_001759/0000/step1_inMINIAODSIM_1.root'
        ] if options.isMC else [
            'root://cms-xrd-global.cern.ch//store/data/Run2022C/ParkingDoubleElectronLowMass0/MINIAOD/PromptReco-v1/000/356/309/00000/2e62494a-6ef9-4262-9afc-95c2b8951970.root'
        ]
annotation = '%s nevts:%d' % (outputFileNANO, options.maxEvents)

# Process
from Configuration.StandardSequences.Eras import eras
from PhysicsTools.BParkingNano.modifiers_cff import *
process = None
if   options.lhcRun == 3: process = cms.Process('BParkNANO',eras.Run3,BToKEE_DiEle)
elif options.lhcRun == 2: process = cms.Process('BParkNANO',eras.Run2_2018)

# import of standard configurations
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('PhysicsTools.BParkingNano.nanoBPark_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Input source
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    secondaryFileNames = cms.untracked.vstring(),
    skipEvents=cms.untracked.uint32(options.skip),
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(options.wantSummary),
)

process.nanoMetadata.strings.tag = annotation
# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string(annotation),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition
process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = outputFileFEVT,
    outputCommands = (cms.untracked.vstring('keep *',
                                            'drop *_*_SelectedTransient*_*',
                     )),
    splitLevel = cms.untracked.int32(0)
)

process.NANOAODoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAOD'),
        filterName = cms.untracked.string('')
    ),
    fileName = outputFileNANO,
    outputCommands = cms.untracked.vstring(
      'drop *',
      "keep nanoaodFlatTable_*Table_*_*",     # event data
      "keep nanoaodUniqueString_nanoMetadata_*_*",   # basic metadata
    )

)


# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globaltag, '')

from PhysicsTools.BParkingNano.nanoBPark_cff import *
if options.lhcRun == 2:
    process = nanoAOD_customizeMuonTriggerBPark(process)
    process = nanoAOD_customizeElectronFilteredBPark(process)
    process = nanoAOD_customizeTrackFilteredBPark(process)
    process = nanoAOD_customizeBToKLL(process)
    process = nanoAOD_customizeBToKstarEE(process)
    process = nanoAOD_customizeBToKstarMuMu(process)
elif options.lhcRun == 3:
    from PhysicsTools.BParkingNano.electronsTrigger_cff import *
    process = nanoAOD_customizeDiEle(process)
    process = nanoAOD_customizeElectronFilteredBPark(process)
    process = nanoAOD_customizeTriggerBitsBPark(process)
    process = nanoAOD_customizeTrackFilteredBPark(process)
    process = nanoAOD_customizeBToKLL(process)

# Path and EndPath definitions
if options.lhcRun == 2:
    process.nanoAOD_KMuMu_step = cms.Path(process.nanoSequence + process.nanoTracksSequence + process.nanoBKMuMuSequence + CountBToKmumu )
    process.nanoAOD_Kee_step   = cms.Path(process.nanoSequence + process.nanoTracksSequence + process.nanoBKeeSequence   + CountBToKee   )
    process.nanoAOD_KstarMuMu_step = cms.Path(process.nanoSequence + process.nanoTracksSequence + process.KstarToKPiSequence + process.nanoBKstarMuMuSequence + CountBToKstarMuMu )
    process.nanoAOD_KstarEE_step  = cms.Path(process.nanoSequence + process.nanoTracksSequence + process.KstarToKPiSequence + process.nanoBKstarEESequence + CountBToKstarEE  )
elif options.lhcRun == 3:
    process.nanoAOD_DiEle_step = cms.Path(process.nanoSequence
                                          +process.nanoDiEleSequence
                                          +process.nanoTracksSequence
                                          +process.nanoBKeeSequence
                                          +CountBToKee)

# customisation of the process.
if options.isMC:
    from PhysicsTools.BParkingNano.nanoBPark_cff import nanoAOD_customizeMC
    nanoAOD_customizeMC(process)

process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)
process.NANOAODoutput_step = cms.EndPath(process.NANOAODoutput)

# Schedule definition
if options.lhcRun == 3:

    process.schedule = cms.Schedule(process.nanoAOD_DiEle_step,
                                    process.endjob_step,
                                    process.NANOAODoutput_step)
    if options.wantFullRECO:
        process.schedule = cms.Schedule(process.nanoAOD_DiEle_step,
                                        process.endjob_step,
                                        process.FEVTDEBUGHLToutput_step,
                                        process.NANOAODoutput_step)
    from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
    associatePatAlgosToolsTask(process)
    process.NANOAODoutput.SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('nanoAOD_DiEle_step')
    )

elif options.lhcRun == 2:

    process.schedule = cms.Schedule(
        process.nanoAOD_KMuMu_step,
        process.nanoAOD_Kee_step,
        #process.nanoAOD_KstarMuMu_step,
        #process.nanoAOD_KstarEE_step,
        process.endjob_step,
        process.NANOAODoutput_step
    )
    if options.wantFullRECO:
        process.schedule = cms.Schedule(
            process.nanoAOD_KMuMu_step,
            process.nanoAOD_Kee_step,
            #process.nanoAOD_KstarMuMu_step,
            #process.nanoAOD_KstarEE_step,
            process.endjob_step,
            process.FEVTDEBUGHLToutput_step,
            process.NANOAODoutput_step
        )
    from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
    associatePatAlgosToolsTask(process)
    process.NANOAODoutput.SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring(
            'nanoAOD_KMuMu_step',
            'nanoAOD_Kee_step',
            #'nanoAOD_KstarMuMu_step',
            #'nanoAOD_KstarEE_step',
        )
    )

### from https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/3287/1/1/1/1/1.html
process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))
process.NANOAODoutput.fakeNameForCrab=cms.untracked.bool(True)    

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)

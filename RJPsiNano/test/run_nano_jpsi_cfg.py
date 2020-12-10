from FWCore.ParameterSet.VarParsing import VarParsing
import FWCore.ParameterSet.Config as cms

options = VarParsing('python')

options.register('isMC', True,
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

options.register('reportEvery', 1000,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "report every N events"
)
options.register('skip',1861,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "skip first N events"
)

options.setDefault('maxEvents',1)
options.setDefault('tag', '10614')
options.parseArguments()

#globaltag = '102X_dataRun2_v11' if not options.isMC else '102X_upgrade2018_realistic_v15'
globaltag = '106X_dataRun2_v28' if not options.isMC else '106X_upgrade2018_realistic_v11_L1v1'


if options._beenSet['globalTag']:
    globaltag = options.globalTag

extension = {False : 'data', True : 'mc'}
outputFileNANO = cms.untracked.string('_'.join(['RJPsi', extension[options.isMC], options.tag])+'.root')

#input files (it can be a list of files)
if not options.inputFiles:
    options.inputFiles = ['root://cms-xrd-global.cern.ch//store/data/Run2018C/Charmonium/MINIAOD/17Sep2018-v1/60000/5865682B-91E0-F047-969B-C52A2FCB241F.root'] if not options.isMC else \
        ['root://cms-xrd-global.cern.ch///store/user/manzoni/RJPsi_Bc_PMX_HLT_RECO_MINI_28oct20_v5/RJpsi-BcToXToJpsiMuMuSelected-RunIISummer19UL18MiniAOD_0.root']
                         #['file:/afs/cern.ch/user/g/garamire/work/private/CMSPisa/RJPsiAnalysis/ulBcReconstruction/CMSSW_10_6_14/src/PhysicsTools/RJPsiNano/test/02F13381-1D94-CC43-948A-2EFFB8572949.root']
                         #['root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18MiniAOD/OniaAndX_ToMuMu_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/260000/AACF29A5-3793-7843-80B9-3396383C32EC.root']
                         #['root://cms-xrd-global.cern.ch//store/user/manzoni/RJPsi_Bc_PMX_HLT_RECO_MINI_28oct20_v5/RJpsi-BcToXToJpsiMuMuSelected-RunIISummer19UL18MiniAOD_1000.root']
                         #['root://cms-xrd-global.cern.ch//store/mc/RunIISummer19UL18MiniAOD/BcToJPsiTauNu_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/MINIAODSIM/106X_upgrade2018_realistic_v11_L1v1_ext1-v2/100000/02F13381-1D94-CC43-948A-2EFFB8572949.root']


#what's the purpose of this 
annotation = '%s nevts:%d' % (outputFileNANO, options.maxEvents)

from Configuration.StandardSequences.Eras import eras
process = cms.Process('RJPsiNANO',eras.Run2_2018)

# import of standard configurations (do we need all of them?)
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('PhysicsTools.RJPsiNano.nanoRJPsi_3Mu_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


#prints the time report
process.Timing = cms.Service("Timing",
                             summaryOnly = cms.untracked.bool(True),
                             useJobReport = cms.untracked.bool(True)
)

#load all the chosen options
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

#purpose?
process.nanoMetadata.strings.tag = annotation
# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string(annotation),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition
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


from PhysicsTools.RJPsiNano.nanoRJPsi_3Mu_cff import *
process = nanoAOD_customizeMuonTriggerBPark(process)  #delete?
process = nanoAOD_customizeTrackFilteredBPark(process)
process = nanoAOD_customizeBTo3Mu(process)
process = nanoAOD_customizeBTo2MuK(process)
process = nanoAOD_customizeBTo2MuP(process)
process = nanoAOD_customizeTriggerBitsBPark(process)  #? delete?




# Path and EndPath definitions
process.nanoAOD_3Mu_step = cms.Path(process.nanoSequence + process.nanoBTo3MuSequence + CountBTo3Mu )
process.nanoAOD_2MuK_step = cms.Path(process.nanoSequence + process.nanoBTo2MuKSequence + CountBTo2MuK )
process.nanoAOD_2MuP_step = cms.Path(process.nanoSequence + process.nanoBTo2MuPSequence + CountBTo2MuP )

# customisation of the process.
if options.isMC:
    from PhysicsTools.RJPsiNano.nanoRJPsi_3Mu_cff import nanoAOD_customizeMC
    nanoAOD_customizeMC(process)

process.endjob_step = cms.EndPath(process.endOfProcess)

process.NANOAODoutput_step = cms.EndPath(process.NANOAODoutput)

# Schedule definition
process.schedule = cms.Schedule(
                                process.nanoAOD_3Mu_step,
                                process.nanoAOD_2MuK_step,
                                process.nanoAOD_2MuP_step,
                                process.endjob_step, 
                                process.NANOAODoutput_step
                               )

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

process.NANOAODoutput.SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring(
                                   'nanoAOD_3Mu_step', 
                                   'nanoAOD_2MuK_step', 
                                   'nanoAOD_2MuP_step', 
                                   )
)

# ?? 
### from https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/3287/1/1/1/1/1.html
process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))
process.NANOAODoutput.fakeNameForCrab=cms.untracked.bool(True)    

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)

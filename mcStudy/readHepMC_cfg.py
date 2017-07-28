# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: readHepMC_cff.py -n 10 --conditions 75X_mcRun2_HeavyIon_v13 --eventcontent AODSIM -s GEN,SIM,DIGI,L1,DIGI2RAW,HLT:HIon,RAW2DIGI,L1Reco,RECO --beamspot Realistic50ns13TeVCollision --datatier GEN-SIM-DIGI-RECO --pileup NoPileUp --era Run2_25ns --no_exec --mc --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1_HI
import FWCore.ParameterSet.Config as cms
import sys

from Configuration.StandardSequences.Eras import eras

process = cms.Process('HLT',eras.Run2_25ns)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.Geometry.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
#process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic50ns13TeVCollision_cfi')
process.load('IOMC.EventVertexGenerators.VtxSmearedNominalHICollision2015_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('HLTrigger.Configuration.HLT_HIon_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.VtxSmeared.src = 'source'
process.genParticles.src = 'source'
process.g4SimHits.Generator.HepMCProductLabel = cms.string('source')
# Input source
process.source = cms.Source("MCFileSource",
    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
    fileNames = cms.untracked.vstring(sys.argv[2]),
    #fileNames = cms.untracked.vstring("hepmc/PbPb2015_jet_hepmc_b12/rapgap_98.dat"),
    inputCommands = cms.untracked.vstring('keep *')
)


process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('readHepMC_cff.py nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.AODSIMoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RECO'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    fileName = cms.untracked.string(sys.argv[3]),
#    fileName = cms.untracked.string("rapgapReco.root"),
    outputCommands = process.AODSIMEventContent.outputCommands
    #outputCommands = process.FEVTSIMEventContent.outputCommands,
)
process.AODSIMoutput.outputCommands.extend(process.FEVTSIMEventContent.outputCommands)

#process.AODSIMoutput.outputCommands.extend(process.FEVTContent.outputCommands)
# Additional output definition
# Other statements
#########################################################################################################
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
process.mix.digitizers = cms.PSet(process.theDigitizersValid)
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '75X_mcRun2_HeavyIon_v13', '')
########################################################################################################
process.generator = cms.EDFilter("inputHepMC")
process.ProductionFilterSequence = cms.Sequence(process.generator)
#########################################################################################################
process.l1GtTriggerMenuXml = cms.ESProducer("L1GtTriggerMenuXmlProducer",
  TriggerMenuLuminosity = cms.string('startup'),
  DefXmlFile = cms.string('L1Menu_CollisionsHeavyIons2015_v5_L1T_Scales_20141121.xml'),
  VmeXmlFile = cms.string('')
)
process.L1GtTriggerMenuRcdSource = cms.ESSource("EmptyESSource",
  recordName = cms.string('L1GtTriggerMenuRcd'),
  iovIsRunNotTime = cms.bool(True),
  firstValid = cms.vuint32(1)
)
process.es_prefer_l1GtParameters = cms.ESPrefer('L1GtTriggerMenuXmlProducer','l1GtTriggerMenuXml')
#########################################################################################################
# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.AODSIMoutput_step = cms.EndPath(process.AODSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.L1simulation_step,process.digi2raw_step)

process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.endjob_step,process.AODSIMoutput_step])

# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 

# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.postLS1Customs
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1_HI 

#call to customisation function customisePostLS1_HI imported from SLHCUpgradeSimulations.Configuration.postLS1Customs
process = customisePostLS1_HI(process)

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforFullSim 

#call to customisation function customizeHLTforFullSim imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforFullSim(process)

# End of customisation functions


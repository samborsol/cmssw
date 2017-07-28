import FWCore.ParameterSet.Config as cms
import sys


process = cms.Process("Forward")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load("CondCore.DBCommon.CondDBSetup_cfi")
process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
#process.load('Configuration.StandardSequences.SkimsHeavyIons_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.MessageLogger.cerr.FwkReport.reportEvery = 1000


process.source = cms.Source("PoolSource",
#                            fileNames = cms.untracked.vstring('file:/home/jgomez2/EAFE4330-EB64-E211-896F-BCAEC5329719.root')
                            fileNames = cms.untracked.vstring('file:/home/jgomez2/CMSSW_5_3_8_HI_patch2/src/Analyzers/ForwardAnalyzer/2011run182398_10_1_6Lr.root')
                            )

process.TFileService = cms.Service("TFileService",
#                                   fileName = cms.string('blah.root')
                                   fileName = cms.string('yay3.root')
                                   )


process.GlobalTag.globaltag = 'GR_R_44_V12::All'


process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFtowers"),
    nonDefaultGlauberModel = cms.string(""),
    centralitySrc = cms.InputTag("hiCentrality")
    )



#########################FILTER##################################################
import FWCore.ParameterSet.Config as cms

# Coincidence of HF towers above threshold
from HeavyIonsAnalysis.Configuration.hfCoincFilter_cff import *

# Selection of at least a two-track fitted vertex
primaryVertexFilter = cms.EDFilter("VertexSelector",
                                   src = cms.InputTag("hiSelectedVertex"),
                                   cut = cms.string("!isFake && abs(z) <= 10 &&  tracksSize >= 2"),
                                   filter = cms.bool(True),   # otherwise it won't filter the events
                                   )

# Cluster-shape filter re-run offline
from RecoLocalTracker.SiPixelRecHits.SiPixelRecHits_cfi import *
from HLTrigger.special.hltPixelClusterShapeFilter_cfi import *
hltPixelClusterShapeFilter.inputTag = "siPixelRecHits"

# Reject BSC beam halo L1 technical bits
from L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff import *
from HLTrigger.HLTfilters.hltLevel1GTSeed_cfi import hltLevel1GTSeed
noBSChalo = hltLevel1GTSeed.clone(
    L1TechTriggerSeeding = cms.bool(True),
    L1SeedsLogicalExpression = cms.string('NOT (36 OR 37 OR 38 OR 39)')
    )

collisionEventSelection = cms.Sequence(noBSChalo *
                                       hfCoincFilter3 *
                                       primaryVertexFilter *
                                       siPixelRecHits *
                                       hltPixelClusterShapeFilter)

#collisionEventSelection = cms.Sequence(primaryVertexFilter)

############################################################################### 

###############################################################################

process.triggerSelection = cms.EDFilter( "TriggerResultsFilter",
                                         triggerConditions = cms.vstring("HLT_HIMinBiasHfOrBSC_*"),
                                         hltResults = cms.InputTag("TriggerResults","","HLT"),
                                         l1tResults = cms.InputTag("gtDigis","","RECO"),
                                         daqPartitions = cms.uint32( 0x01 ),
                                         l1tIgnoreMask = cms.bool( False ),
                                         l1techIgnorePrescales = cms.bool( False ),
                                         throw = cms.bool( True )
                                         )


process.hltbitanalysis = cms.EDAnalyzer("HLTBitAnalyzer",
                                        ### Trigger objects
                                        l1GctHFBitCounts                = cms.InputTag("gctDigis"),
                                        l1GctHFRingSums                 = cms.InputTag("gctDigis"),
                                        l1GtObjectMapRecord             = cms.InputTag("hltL1GtObjectMap::HLT"),
                                        l1GtReadoutRecord               = cms.InputTag("gtDigis::RECO"),

                                        l1extramc                       = cms.string('l1extraParticles'),
                                        l1extramu                       = cms.string('l1extraParticles'),
                                        hltresults                      = cms.InputTag("TriggerResults::HLT"),
                                        HLTProcessName                  = cms.string("HLT"),
                                        UseTFileService                 = cms.untracked.bool(True),

                                        ### Run parameters
                                        RunParameters = cms.PSet(
    HistogramFile = cms.untracked.string('hltbitanalysis.root')
    )
                                        )##end of HLT


process.fwdana = cms.EDAnalyzer('ForwardAnalyzer')

process.upcvertexana = cms.EDAnalyzer('UPCVertexAnalyzer',
                                      vertexCollection=cms.string("hiSelectedVertex")
                                      )

process.upcselectedtrackana = cms.EDAnalyzer('UPCTrackAnalyzer',
                                          trackCollection=cms.string("hiSelectedTracks")
                                             )

process.upccentralityana = cms.EDAnalyzer('UPCCentralityAnalyzer',
                                          centralityVariable=process.HeavyIonGlobalParameters.centralityVariable
                                          )

process.calotowerana = cms.EDAnalyzer('CaloTowerAnalyzer',
                                      towerCollection=cms.string("CaloTower")
                                      )


process.trackSequence = cms.Sequence(process.upcvertexana*process.upcselectedtrackana)
process.forwardSequence = cms.Sequence(process.fwdana)
process.triggerSequence = cms.Sequence(process.hltbitanalysis)
process.centralitySequence = cms.Sequence(process.upccentralityana)
process.caloSequence = cms.Sequence(process.calotowerana)
process.filterSequence = cms.Sequence(process.triggerSelection*process.collisionEventSelection)

#process.path = cms.Path(process.triggerSelection+
#                      process.trackSequence+
#                     process.forwardSequence+
#                    process.triggerSequence+
#                   process.centralitySequence+
#                  process.caloSequence
# )

process.path = cms.Path(process.filterSequence+
                        process.trackSequence+
                        process.forwardSequence+
                        process.triggerSequence+
                        process.centralitySequence+
                        process.caloSequence
                        )

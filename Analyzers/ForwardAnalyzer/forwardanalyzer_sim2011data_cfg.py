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
#process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
#process.load('Configuration.StandardSequences.SkimsHeavyIons_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load('Analyzers/ForwardAnalyzer/FakeTrackAnalyzer_cfi')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.MessageLogger.cerr.FwkReport.reportEvery = 1000


process.source = cms.Source("PoolSource",
#                            fileNames = cms.untracked.vstring('file:/home/jgomez2/EAFE4330-EB64-E211-896F-BCAEC5329719.root')
                            fileNames = cms.untracked.vstring('file:/home/jgomez2/Simulation/CMSSW_4_4_2_patch5/src/amptDefault_cfi_py_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO.root')
                            )

process.TFileService = cms.Service("TFileService",
#                                   fileName = cms.string('castorfixed.root')
                                   fileName = cms.string('ForwardTrees_skim_alltracks_2010.root')
                                   )


process.GlobalTag.globaltag = 'GR_R_44_V15::All'



#from CmsHi.Analysis2010.CommonFunctions_cff import *
#overrideCentrality(process)
#process.HeavyIonGlobalParameters = cms.PSet(
#    centralityVariable = cms.string("HFhits"),
#    nonDefaultGlauberModel = cms.string(""),
#    centralitySrc = cms.InputTag("hiCentrality")
#    )
process.HeavyIonGlobalParameters = cms.PSet(
        centralityVariable = cms.string("HFtowers"),
            nonDefaultGlauberModel = cms.string(""),
            centralitySrc = cms.InputTag("hiCentrality")
            )


#########################FILTER##################################################
import FWCore.ParameterSet.Config as cms


process.fwdana = cms.EDAnalyzer('ForwardAnalyzer')

process.upcselectedtrackana = cms.EDAnalyzer('UPCTrackAnalyzer',
                                             trackCollection=cms.string("hiGoodTightMergedTracks")
                                             )

process.calotowerana = cms.EDAnalyzer('CaloTowerAnalyzer',
                                      towerCollection=cms.string("CaloTower")
                                      )

process.upccentralityana = cms.EDAnalyzer('UPCCentralityAnalyzer',
                                          centralityVariable=process.HeavyIonGlobalParameters.centralityVariable
                                          )
#process.fakeana = cms.EDAnalyzer('FakeTrackAnalyzer',
 #                                trackCollection=cms.string("hiGoodTightMergedTracks")
  #                               )



process.trackSequence = cms.Sequence(process.upcselectedtrackana)
process.centralitySequence = cms.Sequence(process.upccentralityana)
#process.trackSequence = cms.Sequence(process.upcselectedtrackana*process.FakeTrkCorr)
process.caloSequence = cms.Sequence(process.calotowerana)
process.path = cms.Path(process.trackSequence+
                        process.centralitySequence+
                        process.caloSequence
                        )

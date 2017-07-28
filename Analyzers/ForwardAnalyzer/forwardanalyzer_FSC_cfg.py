import FWCore.ParameterSet.Config as cms
import sys

if len(sys.argv) > 2:
    file=open(sys.argv[2])
    outfile=sys.argv[3]

process = cms.Process("Forward")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load("RecoLocalCalo.HcalRecProducers.HcalHitReconstructor_hf_cfi")
process.load("RecoLocalCalo.HcalRecProducers.HcalHitReconstructor_zdc_cfi")
process.load("EventFilter.CastorRawToDigi.CastorRawToDigi_cfi")
process.load("RecoLocalCalo.CastorReco.CastorSimpleReconstructor_cfi")

process.load("CondCore.DBCommon.CondDBSetup_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.castorDigis = cms.EDProducer("CastorRawToDigi",
                                     CastorFirstFED = cms.untracked.int32(690),
                                     FilterDataQuality = cms.bool(True),
                                     ExceptionEmptyData = cms.untracked.bool(True),
                                     InputLabel = cms.InputTag("rawDataCollector"),
                                     UnpackCalib = cms.untracked.bool(False),
                                     FEDs = cms.untracked.vint32(690,691,692),
                                     lastSample = cms.int32(9),
                                     firstSample = cms.int32(0),
                                     )

process.hcalDigis = cms.EDProducer("HcalRawToDigi",
                                   UnpackZDC = cms.untracked.bool(True),
                                   UnpackTTP = cms.untracked.bool(True),
                                   FilterDataQuality = cms.bool(True),
                                   HcalFirstFED = cms.untracked.int32(700),
                                   InputLabel = cms.InputTag("rawDataCollector"),
                                   ComplainEmptyData = cms.untracked.bool(True),
                                   UnpackCalib = cms.untracked.bool(True),
                                   FEDs = cms.untracked.vint32(722),
                                   streams = cms.untracked.vstring(
                  'HCAL_Trigger','HCAL_SlowData','HCAL_QADCTDC'),
                                   lastSample = cms.int32(9),
                                   firstSample = cms.int32(0),
                                   )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )


process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(file.readlines()) 
                            #    'file:/home/ferraioc/CMSSW_4_4_2_patch5/src/crabjob/res/2011run182398_10_1_6Lr.root'
                            )


process.hcal_db_producer = cms.ESProducer("HcalDbProducer",
                                          dump = cms.untracked.vstring(''),
                                          file = cms.untracked.string('')
                                          )

#process.GlobalTag.globaltag = 'GR_R_44_V10::All'
#process.GlobalTag.globaltag = 'GR_R_53_V19::All'
process.GlobalTag.globaltag = 'GR_E_V33::All'

#########################HLT Stuff...Seems to Only Work for Reco##############

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
                                        )###############END OF HLT STUFF


###L1 Stuff
process.GtDigis = cms.EDProducer( "L1GlobalTriggerRawToDigi",
                                  #DaqGtInputTag = cms.InputTag( "rawDataRepacker" ),
                                  DaqGtInputTag = cms.InputTag("rawDataCollector"),
                                  DaqGtFedId = cms.untracked.int32( 813 ),
                                  ActiveBoardsMask = cms.uint32( 0xffff ),
                                  UnpackBxInEvent = cms.int32( 5 ),
                                  Verbosity = cms.untracked.int32( 0 )
                                  )
#####END OF L1 Stuff 


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outfile)
                                   )

process.hfreco.firstSample = 1
process.hfreco.samplesToAdd = 6

process.fwdana = cms.EDAnalyzer('ForwardAnalyzer_FSC',
                                l1GtRR=cms.InputTag("GtDigis")
                                )

#process.p = cms.Path(process.zdcreco*process.fwdana)
#process.p = cms.Path(process.GtDigis*process.hltbitanalysis*process.hcalDigis*process.castorDigis*process.castorreco*process.hfreco*process.fwdana)
process.p = cms.Path(process.hcalDigis*process.fwdana)

import FWCore.ParameterSet.Config as cms
import sys

process = cms.Process('HiForest')
process.options = cms.untracked.PSet()

process.load("HeavyIonsAnalysis.JetAnalysis.HiForest_cff")
process.HiForest.inputLines = cms.vstring("HiForest V3",)
import subprocess
version = subprocess.Popen(["(cd $CMSSW_BASE/src && git describe --tags)"], stdout=subprocess.PIPE, shell=True).stdout.read()
if version == '':
    version = 'no git info'
process.HiForest.HiForestVersion = cms.string(version)

process.source = cms.Source("PoolSource",
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
                            fileNames = cms.untracked.vstring(
#                                "file:aod/PbPb2015_jet_aod/aodMC1.root"
				sys.argv[2]
                            )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1))


process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '75X_mcRun2_HeavyIon_v13', '')
process.HiForest.GlobalTagLabel = process.GlobalTag.globaltag

from HeavyIonsAnalysis.Configuration.CommonFunctions_cff import overrideJEC_pp5020
process = overrideJEC_pp5020(process)

process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string(
#					"tree/PbPb2015_jet_tree/hiforest1.root"
					sys.argv[3]
					))

process.load('RecoJets.JetProducers.ak5PFJets_cfi')
process.ak5PFJets.doAreaFastjet = True
process.ak3PFJets = process.ak5PFJets.clone(rParam = 0.3)
process.ak4PFJets = process.ak5PFJets.clone(rParam = 0.4)

process.load('RecoJets.JetProducers.ak5GenJets_cfi')
process.ak3GenJets = process.ak5GenJets.clone(rParam = 0.3)
process.ak4GenJets = process.ak5GenJets.clone(rParam = 0.4)

process.load('RecoJets.Configuration.GenJetParticles_cff')
process.load('RecoHI.HiJetAlgos.HiGenJets_cff')
process.load('HeavyIonsAnalysis.JetAnalysis.makePartons_cff')
process.highPurityTracks = cms.EDFilter("TrackSelector",
                                src = cms.InputTag("generalTracks"),
                                cut = cms.string('quality("highPurity")')
)
process.load('HeavyIonsAnalysis.JetAnalysis.jets.ak3PFJetSequence_pp_mc_cff')
process.load('HeavyIonsAnalysis.JetAnalysis.jets.ak4PFJetSequence_pp_mc_cff')
process.load('HeavyIonsAnalysis.JetAnalysis.jets.ak5PFJetSequence_pp_mc_cff')

process.ak3PFJetAnalyzer.useHepMC = cms.untracked.bool(True)
process.ak3PFJetAnalyzer.eventInfoTag = cms.InputTag("source")
process.ak4PFJetAnalyzer.useHepMC = cms.untracked.bool(True)
process.ak4PFJetAnalyzer.eventInfoTag = cms.InputTag("source")
process.ak5PFJetAnalyzer.useHepMC = cms.untracked.bool(True)
process.ak5PFJetAnalyzer.eventInfoTag = cms.InputTag("source")

process.jetSequences = cms.Sequence(
    process.myPartons +
    process.genParticlesForJets +
    process.ak3GenJets +
    process.ak4GenJets +
    process.ak5GenJets +
    process.ak3PFJets +
    process.ak4PFJets +
    process.ak5PFJets +
    process.highPurityTracks +
    process.ak3PFJetSequence +
    process.ak4PFJetSequence +
    process.ak5PFJetSequence 
)

process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cff')

process.hiEvtAnalyzer = cms.EDAnalyzer("HiEvtAnalyzer",
    CentralityBinSrc = cms.InputTag("centralityBin","HFtowers"),
    CentralitySrc = cms.InputTag("hiCentrality"),
    EvtPlane = cms.InputTag("hiEvtPlane"),
    EvtPlaneFlat = cms.InputTag("hiEvtPlaneFlat"),
    HiMC = cms.InputTag("heavyIon"),
    Vertex = cms.InputTag("offlinePrimaryVertices"),
    doCentrality = cms.bool(False),
    doEvtPlane = cms.bool(False),
    doEvtPlaneFlat = cms.bool(False),
    doHiMC = cms.bool(False),
    doMC = cms.bool(False),
    doVertex = cms.bool(True),
    evtPlaneLevel = cms.int32(0),
    useHepMC = cms.bool(False)
)

process.load('HeavyIonsAnalysis.JetAnalysis.HiGenAnalyzer_cfi')
process.HiGenParticleAna.ptMin = -9999
process.HiGenParticleAna.genParticleSrc = cms.untracked.InputTag("genParticles")
process.load("GeneratorInterface.HiGenCommon.HeavyIon_cff")
process.load('HeavyIonsAnalysis.EventAnalysis.HiMixAnalyzerRECO_cff')
process.load('GeneratorInterface.HiGenCommon.HeavyIon_cff')
process.load('HeavyIonsAnalysis.JetAnalysis.HiGenAnalyzer_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.runanalyzer_cff')
process.HiGenParticleAna.doHI = False


process.load('HeavyIonsAnalysis.EventAnalysis.runanalyzer_cff')
process.load("HeavyIonsAnalysis.JetAnalysis.pfcandAnalyzer_cfi")
process.pfcandAnalyzer.skipCharged = False
process.pfcandAnalyzer.pfPtMin = 0
process.pfcandAnalyzer.pfCandidateLabel = cms.InputTag("particleFlow")
process.pfcandAnalyzer.doVS = cms.untracked.bool(False)
process.pfcandAnalyzer.doUEraw_ = cms.untracked.bool(False)
process.pfcandAnalyzer.genLabel = cms.InputTag("genParticles")
process.load("HeavyIonsAnalysis.JetAnalysis.hcalNoise_cff")

process.load('HeavyIonsAnalysis.JetAnalysis.ExtraTrackReco_cff')
process.load('HeavyIonsAnalysis.JetAnalysis.TrkAnalyzers_cff')

process.load('HeavyIonsAnalysis.PhotonAnalysis.ggHiNtuplizer_cfi')
process.ggHiNtuplizer.gsfElectronLabel   = cms.InputTag("gedGsfElectrons")
process.ggHiNtuplizer.recoPhotonHiIsolationMap = cms.InputTag('photonIsolationHIProducerpp')
process.ggHiNtuplizer.VtxLabel           = cms.InputTag("offlinePrimaryVertices")
process.ggHiNtuplizer.particleFlowCollection = cms.InputTag("particleFlow")
process.ggHiNtuplizer.doVsIso            = cms.bool(False)
process.ggHiNtuplizer.doElectronVID      = cms.bool(True)
process.ggHiNtuplizerGED = process.ggHiNtuplizer.clone(recoPhotonSrc = cms.InputTag('gedPhotons'),
                                                       recoPhotonHiIsolationMap = cms.InputTag('photonIsolationHIProducerppGED'))

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.AOD
switchOnVIDElectronIdProducer(process, dataFormat)

my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff']

for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.load("HeavyIonsAnalysis.VectorBosonAnalysis.tupelSequence_pp_mc_cff")
process.load('HeavyIonsAnalysis.JetAnalysis.rechitanalyzer_aod_cfi')

process.rechitAna = cms.Sequence(process.rechitanalyzer)

process.load("HeavyIonsAnalysis.MuonAnalysis.hltMuTree_cfi")
process.hltMuTree.vertices = cms.InputTag("offlinePrimaryVertices")

process.UPCTrigger = cms.EDFilter("HLTHighLevel",
     HLTPaths = cms.vstring('HLT_HIUPCSingleEG5NotHF2Pixel_SingleTrack_v1'),
     TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
     andOr = cms.bool(True),
     eventSetupPathsKey = cms.string(''),
     throw = cms.bool(True)
)

process.fwdana = cms.EDAnalyzer("ForwardAnalyzer",
     trackSrc = cms.InputTag("generalTracks"),
     vtxCollection_ = cms.InputTag("offlinePrimaryVertices")
#     CentralitySrc = cms.InputTag("hiCentrality"),
#     CentralityBinSrc = cms.InputTag("centralityBin","HFtowers")
)

process.patPhotonSequence = cms.Sequence(process.patPhotons)

process.hltobject = cms.EDAnalyzer("TriggerObjectAnalyzer",
    processName = cms.string('HLT'),
    treeName = cms.string('JetTriggers'),
    triggerEvent = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    triggerNames = cms.vstring('HLT_AK4CaloJet100_Eta5p1_v',
        'HLT_AK4CaloJet100_Jet35_Eta0p7_v',
        'HLT_AK4CaloJet100_Jet35_Eta1p1_v',
        'HLT_AK4CaloJet110_Eta5p1_v',
        'HLT_AK4CaloJet120_Eta5p1_v',
        'HLT_AK4CaloJet150_v',
        'HLT_AK4CaloJet40_Eta5p1_v',
        'HLT_AK4CaloJet60_Eta5p1_v',
        'HLT_AK4CaloJet80_45_45_Eta2p1_v',
        'HLT_AK4CaloJet80_Eta5p1_v',
        'HLT_AK4CaloJet80_Jet35_Eta0p7_v',
        'HLT_AK4CaloJet80_Jet35_Eta1p1_v',
        'HLT_AK4PFBJetBCSV60_Eta2p1_v',
        'HLT_AK4PFBJetBCSV80_Eta2p1_v',
        'HLT_AK4PFBJetBSSV60_Eta2p1_v',
        'HLT_AK4PFBJetBSSV80_Eta2p1_v',
        'HLT_AK4PFDJet60_Eta2p1_v',
        'HLT_AK4PFDJet80_Eta2p1_v',
        'HLT_AK4PFJet100_Eta5p1_v',
        'HLT_AK4PFJet110_Eta5p1_v',
        'HLT_AK4PFJet120_Eta5p1_v',
        'HLT_AK4PFJet40_Eta5p1_v',
        'HLT_AK4PFJet60_Eta5p1_v',
        'HLT_AK4PFJet80_Eta5p1_v',
        'HLT_HIL1DoubleMu0_v',
        'HLT_HIL1DoubleMu10_v',
        'HLT_HIL2DoubleMu0_NHitQ_v',
        'HLT_HIL2Mu15_v',
        'HLT_HIL2Mu20_v',
        'HLT_HIL2Mu3Eta2p5_AK4CaloJet100Eta2p1_v',
        'HLT_HIL2Mu3Eta2p5_AK4CaloJet40Eta2p1_v',
        'HLT_HIL2Mu3Eta2p5_AK4CaloJet60Eta2p1_v',
        'HLT_HIL2Mu3Eta2p5_AK4CaloJet80Eta2p1_v',
        'HLT_HIL2Mu3Eta2p5_HIPhoton10Eta1p5_v',
        'HLT_HIL2Mu3Eta2p5_HIPhoton15Eta1p5_v',
        'HLT_HIL2Mu3Eta2p5_HIPhoton20Eta1p5_v',
        'HLT_HIL2Mu3Eta2p5_HIPhoton30Eta1p5_v',
        'HLT_HIL2Mu3Eta2p5_HIPhoton40Eta1p5_v',
        'HLT_HIL2Mu3_NHitQ10_v',
        'HLT_HIL2Mu5_NHitQ10_v',
        'HLT_HIL2Mu7_NHitQ10_v',
        'HLT_HIL3DoubleMu0_OS_m2p5to4p5_v',
        'HLT_HIL3DoubleMu0_OS_m7to14_v',
        'HLT_HIL3Mu15_v',
        'HLT_HIL3Mu20_v',
        'HLT_HIL3Mu3_NHitQ15_v',
        'HLT_HIL3Mu5_NHitQ15_v',
        'HLT_HIL3Mu7_NHitQ15_v',
        'HLT_HISinglePhoton10_Eta1p5_v',
        'HLT_HISinglePhoton10_Eta3p1_v',
        'HLT_HISinglePhoton15_Eta1p5_v',
        'HLT_HISinglePhoton15_Eta3p1_v',
        'HLT_HISinglePhoton20_Eta1p5_v',
        'HLT_HISinglePhoton20_Eta3p1_v',
        'HLT_HISinglePhoton30_Eta1p5_v',
        'HLT_HISinglePhoton30_Eta3p1_v',
        'HLT_HISinglePhoton40_Eta1p5_v',
        'HLT_HISinglePhoton40_Eta3p1_v',
        'HLT_HISinglePhoton50_Eta1p5_v',
        'HLT_HISinglePhoton50_Eta3p1_v',
        'HLT_HISinglePhoton60_Eta1p5_v',
        'HLT_HISinglePhoton60_Eta3p1_v',
        'HLT_FullTrack18ForPPRef_v',
        'HLT_FullTrack24ForPPRef_v',
        'HLT_FullTrack34ForPPRef_v',
        'HLT_FullTrack45ForPPRef_v',
        'HLT_DmesonPPTrackingGlobal_Dpt8_v',
        'HLT_DmesonPPTrackingGlobal_Dpt15_v',
        'HLT_DmesonPPTrackingGlobal_Dpt20_v',
        'HLT_DmesonPPTrackingGlobal_Dpt30_v',
        'HLT_DmesonPPTrackingGlobal_Dpt40_v',
        'HLT_DmesonPPTrackingGlobal_Dpt50_v',
        'HLT_DmesonPPTrackingGlobal_Dpt60_v'),
    triggerResults = cms.InputTag("TriggerResults","","HLT")
)

process.castorDigis = cms.EDProducer("CastorRawToDigi",
    CastorCtdc = cms.bool(False),
    CastorFirstFED = cms.int32(690),
    ComplainEmptyData = cms.untracked.bool(False),
    ExceptionEmptyData = cms.untracked.bool(False),
    ExpectedOrbitMessageTime = cms.int32(-1),
    FEDs = cms.untracked.vint32(690, 691, 692, 693, 722),
    FilterDataQuality = cms.bool(True),
    InputLabel = cms.InputTag("rawDataCollector"),
    UnpackTTP = cms.bool(True),
    UnpackZDC = cms.bool(True),
    UseNominalOrbitMessageTime = cms.bool(True),
    ZDCFirstFED = cms.int32(693),
    firstSample = cms.int32(0),
    lastSample = cms.int32(9),
    silent = cms.untracked.bool(False)
)

process.castorreco = cms.EDProducer("CastorSimpleReconstructor",
    Subdetector = cms.string('CASTOR'),
    correctForPhaseContainment = cms.bool(False),
    correctForTimeslew = cms.bool(False),
    correctionPhaseNS = cms.double(0.0),
    digiLabel = cms.InputTag("castorDigis"),
    doSaturationCorr = cms.bool(True),
    firstSample = cms.int32(4),
    maxADCvalue = cms.int32(127),
    samplesToAdd = cms.int32(2),
    setSaturationFlag = cms.bool(True),
    tsFromDB = cms.bool(True)
)

process.hfreco = cms.EDProducer("HcalHitReconstructor",
    HFInWindowStat = cms.PSet(
        hflongEthresh = cms.double(40.0),
        hflongMaxWindowTime = cms.vdouble(10),
        hflongMinWindowTime = cms.vdouble(-10),
        hfshortEthresh = cms.double(40.0),
        hfshortMaxWindowTime = cms.vdouble(10),
        hfshortMinWindowTime = cms.vdouble(-12)
    ),
    PETstat = cms.PSet(
        HcalAcceptSeverityLevel = cms.int32(9),
        longETParams = cms.vdouble(0, 0, 0, 0, 0,
            0, 0, 0, 0, 0,
            0, 0, 0),
        longEnergyParams = cms.vdouble(43.5, 45.7, 48.32, 51.36, 54.82,
            58.7, 63.0, 67.72, 72.86, 78.42,
            84.4, 90.8, 97.62),
        long_R = cms.vdouble(0.98),
        long_R_29 = cms.vdouble(0.8),
        shortETParams = cms.vdouble(0, 0, 0, 0, 0,
            0, 0, 0, 0, 0,
            0, 0, 0),
        shortEnergyParams = cms.vdouble(35.1773, 35.37, 35.7933, 36.4472, 37.3317,
            38.4468, 39.7925, 41.3688, 43.1757, 45.2132,
            47.4813, 49.98, 52.7093),
        short_R = cms.vdouble(0.8),
        short_R_29 = cms.vdouble(0.8)
    ),
    S8S1stat = cms.PSet(
        HcalAcceptSeverityLevel = cms.int32(9),
        isS8S1 = cms.bool(True),
        longETParams = cms.vdouble(0, 0, 0, 0, 0,
            0, 0, 0, 0, 0,
            0, 0, 0),
        longEnergyParams = cms.vdouble(40, 100, 100, 100, 100,
            100, 100, 100, 100, 100,
            100, 100, 100),
        long_optimumSlope = cms.vdouble(0.3, 0.1, 0.1, 0.1, 0.1,
            0.1, 0.1, 0.1, 0.1, 0.1,
            0.1, 0.1, 0.1),
        shortETParams = cms.vdouble(0, 0, 0, 0, 0,
            0, 0, 0, 0, 0,
            0, 0, 0),
        shortEnergyParams = cms.vdouble(40, 100, 100, 100, 100,
            100, 100, 100, 100, 100,
            100, 100, 100),
        short_optimumSlope = cms.vdouble(0.3, 0.1, 0.1, 0.1, 0.1,
            0.1, 0.1, 0.1, 0.1, 0.1,
            0.1, 0.1, 0.1)
    ),
    S9S1stat = cms.PSet(
        HcalAcceptSeverityLevel = cms.int32(9),
        isS8S1 = cms.bool(False),
        longETParams = cms.vdouble(0, 0, 0, 0, 0,
            0, 0, 0, 0, 0,
            0, 0, 0),
        longEnergyParams = cms.vdouble(43.5, 45.7, 48.32, 51.36, 54.82,
            58.7, 63.0, 67.72, 72.86, 78.42,
            84.4, 90.8, 97.62),
        long_optimumSlope = cms.vdouble(-99999, 0.0164905, 0.0238698, 0.0321383, 0.041296,
            0.0513428, 0.0622789, 0.0741041, 0.0868186, 0.100422,
            0.135313, 0.136289, 0.0589927),
        shortETParams = cms.vdouble(0, 0, 0, 0, 0,
            0, 0, 0, 0, 0,
            0, 0, 0),
        shortEnergyParams = cms.vdouble(35.1773, 35.37, 35.7933, 36.4472, 37.3317,
            38.4468, 39.7925, 41.3688, 43.1757, 45.2132,
            47.4813, 49.98, 52.7093),
        short_optimumSlope = cms.vdouble(-99999, 0.0164905, 0.0238698, 0.0321383, 0.041296,
            0.0513428, 0.0622789, 0.0741041, 0.0868186, 0.100422,
            0.135313, 0.136289, 0.0589927)
    ),
    Subdetector = cms.string('HF'),
    correctForPhaseContainment = cms.bool(False),
    correctForTimeslew = cms.bool(False),
    correctTiming = cms.bool(True),
    correctionPhaseNS = cms.double(13.0),
    dataOOTCorrectionCategory = cms.string('Data'),
    dataOOTCorrectionName = cms.string(''),
    digiLabel = cms.InputTag("hcalDigis"),
    digiTimeFromDB = cms.bool(True),
    digistat = cms.PSet(
        HFdigiflagCoef = cms.vdouble(0.93, -0.38275, -0.012667),
        HFdigiflagExpectedPeak = cms.int32(2),
        HFdigiflagFirstSample = cms.int32(1),
        HFdigiflagMinEthreshold = cms.double(40),
        HFdigiflagSamplesToAdd = cms.int32(3)
    ),
    dropZSmarkedPassed = cms.bool(True),
    firstAuxTS = cms.int32(1),
    firstSample = cms.int32(2),
    hfTimingTrustParameters = cms.PSet(
        hfTimingTrustLevel1 = cms.int32(1),
        hfTimingTrustLevel2 = cms.int32(4)
    ),
    mcOOTCorrectionCategory = cms.string('MC'),
    mcOOTCorrectionName = cms.string(''),
    puCorrMethod = cms.int32(0),
    recoParamsFromDB = cms.bool(True),
    samplesToAdd = cms.int32(1),
    saturationParameters = cms.PSet(
        maxADCvalue = cms.int32(127)
    ),
    setHSCPFlags = cms.bool(False),
    setNegativeFlags = cms.bool(False),
    setNoiseFlags = cms.bool(True),
    setPulseShapeFlags = cms.bool(False),
    setSaturationFlags = cms.bool(True),
    setTimingTrustFlags = cms.bool(True),
    tsFromDB = cms.bool(True),
    useLeakCorrection = cms.bool(False)
)

process.PixelCPEGenericESProducer = cms.ESProducer("PixelCPEGenericESProducer",
    Alpha2Order = cms.bool(True),
    ClusterProbComputationFlag = cms.int32(0),
    ComponentName = cms.string('PixelCPEGeneric'),
    DoCosmics = cms.bool(False),
    EdgeClusterErrorX = cms.double(50.0),
    EdgeClusterErrorY = cms.double(85.0),
    IrradiationBiasCorrection = cms.bool(False),
    LoadTemplatesFromDB = cms.bool(True),
    MagneticFieldRecord = cms.ESInputTag(""),
    PixelErrorParametrization = cms.string('NOTcmsim'),
    TruncatePixelCharge = cms.bool(True),
    UseErrorsFromTemplates = cms.bool(True),
    eff_charge_cut_highX = cms.double(1.0),
    eff_charge_cut_highY = cms.double(1.0),
    eff_charge_cut_lowX = cms.double(0.0),
    eff_charge_cut_lowY = cms.double(0.0),
    inflate_all_errors_no_trk_angle = cms.bool(False),
    inflate_errors = cms.bool(False),
    size_cutX = cms.double(3.0),
    size_cutY = cms.double(3.0),
    useLAAlignmentOffsets = cms.bool(False),
    useLAWidthFromDB = cms.bool(True)
)

process.hltanalysis = cms.EDAnalyzer("HLTBitAnalyzer",
    HLTProcessName = cms.string('HLT'),
    RunParameters = cms.PSet(
        HistogramFile = cms.untracked.string('hltbitanalysis.root')
    ),
    UseTFileService = cms.untracked.bool(True),
    dummyBranches = cms.untracked.vstring('HLT_HIUPCL1SingleMuOpenNotHF2_v1',
        'HLT_HIUPCSingleMuNotHF2Pixel_SingleTrack_v1',
        'HLT_HIUPCL1DoubleMuOpenNotHF2_v1',
        'HLT_HIUPCDoubleMuNotHF2Pixel_SingleTrack_v1',
        'HLT_HIUPCL1SingleEG2NotHF2_v1',
        'HLT_HIUPCSingleEG2NotHF2Pixel_SingleTrack_v1',
        'HLT_HIUPCL1DoubleEG2NotHF2_v1',
        'HLT_HIUPCDoubleEG2NotHF2Pixel_SingleTrack_v1',
        'HLT_HIUPCL1SingleEG5NotHF2_v1',
        'HLT_HIUPCSingleEG5NotHF2Pixel_SingleTrack_v1',
        'HLT_HIUPCL1DoubleMuOpenNotHF1_v1',
        'HLT_HIUPCDoubleMuNotHF1Pixel_SingleTrack_v1',
        'HLT_HIUPCL1DoubleEG2NotZDCAND_v1',
        'HLT_HIUPCL1DoubleEG2NotZDCANDPixel_SingleTrack_v1',
        'HLT_HIUPCL1DoubleMuOpenNotZDCAND_v1',
        'HLT_HIUPCL1DoubleMuOpenNotZDCANDPixel_SingleTrack_v1',
        'HLT_HIUPCL1EG2NotZDCAND_v1',
        'HLT_HIUPCEG2NotZDCANDPixel_SingleTrack_v1',
        'HLT_HIUPCL1MuOpenNotZDCAND_v1',
        'HLT_HIUPCL1MuOpenNotZDCANDPixel_SingleTrack_v1',
        'HLT_HIUPCL1NotHFplusANDminusTH0BptxAND_v1',
        'HLT_HIUPCL1NotHFplusANDminusTH0BptxANDPixel_SingleTrack_v1',
        'HLT_HIUPCL1NotHFMinimumbiasHFplusANDminustTH0_v1',
        'HLT_HIUPCL1NotHFMinimumbiasHFplusANDminustTH0Pixel_SingleTrack_v1',
        'HLT_HIUPCL1DoubleMuOpenNotHFMinimumbiasHFplusANDminustTH0_v1',
        'HLT_HIUPCL1DoubleMuOpenNotHFMinimumbiasHFplusANDminustTH0Pixel_SingleTrack_v1',
        'HLT_HIUPCL1NotMinimumBiasHF2_AND_v1',
        'HLT_HIUPCL1NotMinimumBiasHF2_ANDPixel_SingleTrack_v1',
        'HLT_HIUPCL1ZdcOR_BptxAND_v1',
        'HLT_HIUPCL1ZdcOR_BptxANDPixel_SingleTrack_v1',
        'HLT_HIUPCL1ZdcXOR_BptxAND_v1',
        'HLT_HIUPCL1ZdcXOR_BptxANDPixel_SingleTrack_v1',
        'HLT_HIUPCL1NotZdcOR_BptxAND_v1',
        'HLT_HIUPCL1NotZdcOR_BptxANDPixel_SingleTrack_v1',
        'HLT_HIPuAK4CaloJet40_Eta5p1_v1',
        'HLT_HIPuAK4CaloJet60_Eta5p1_v1',
        'HLT_HIPuAK4CaloJet80_Eta5p1_v1',
        'HLT_HIPuAK4CaloJet100_Eta5p1_v1',
        'HLT_HIZeroBias_v1',
        'HLT_HIZeroBiasPixel_SingleTrack_v1',
        'HLT_HIRandom_v1'),
    hltresults = cms.InputTag("TriggerResults","","HLT"),
    l1GctHFBitCounts = cms.InputTag("gctDigis"),
    l1GctHFRingSums = cms.InputTag("gctDigis"),
    l1GtObjectMapRecord = cms.InputTag("hltL1GtObjectMap","","HLT"),
    l1GtReadoutRecord = cms.InputTag("gtDigis"),
    l1extramc = cms.string('l1extraParticles'),
    l1extramu = cms.string('l1extraParticles')
)

process.siPixelTemplateDBObjectESProducer = cms.ESProducer("SiPixelTemplateDBObjectESProducer")

process.stripCPEESProducer = cms.ESProducer("StripCPEESProducer",
    ComponentName = cms.string('stripCPE'),
    ComponentType = cms.string('SimpleStripCPE'),
    parameters = cms.PSet(

    )
)

process.templates = cms.ESProducer("PixelCPETemplateRecoESProducer",
    Alpha2Order = cms.bool(True),
    ClusterProbComputationFlag = cms.int32(0),
    ComponentName = cms.string('PixelCPETemplateReco'),
    DoCosmics = cms.bool(False),
    DoLorentz = cms.bool(True),
    LoadTemplatesFromDB = cms.bool(True),
    UseClusterSplitter = cms.bool(False),
    speed = cms.int32(-2)
)

process.ak4PFcorr = cms.EDProducer("JetCorrFactorsProducer",
    emf = cms.bool(False),
    extraJPTOffset = cms.string('L1FastJet'),
    flavorType = cms.string('J'),
    levels = cms.vstring('L2Relative',
        'L3Absolute',
        'L2L3Residual'),
    payload = cms.string('AK4PF_offline'),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    src = cms.InputTag("ak4PFJets"),
    useNPV = cms.bool(False),
    useRho = cms.bool(False)
)

process.zdcreco = cms.EDProducer("ZdcHitReconstructor",
    AuxTSvec = cms.vint32(4, 5, 6, 7),
    Subdetector = cms.string('ZDC'),
    correctForPhaseContainment = cms.bool(False),
    correctForTimeslew = cms.bool(False),
    correctTiming = cms.bool(True),
    correctionPhaseNS = cms.double(0.0),
    digiLabel = cms.InputTag("hcalDigis"),
    dropZSmarkedPassed = cms.bool(True),
    lowGainFrac = cms.double(8.15),
    lowGainOffset = cms.int32(1),
    recoMethod = cms.int32(2),
    saturationParameters = cms.PSet(
        maxADCvalue = cms.int32(127)
    ),
    setHSCPFlags = cms.bool(True),
    setNoiseFlags = cms.bool(True),
    setSaturationFlags = cms.bool(True),
    setTimingTrustFlags = cms.bool(False)
)

process.pfTowers = cms.EDAnalyzer("RecHitTreeProducer",
    BasicClusterSrc1 = cms.untracked.InputTag("islandBasicClusters","islandBarrelBasicClusters"),
    EBRecHitSrc = cms.untracked.InputTag("ecalRecHit","EcalRecHitsEB"),
    EERecHitSrc = cms.untracked.InputTag("ecalRecHit","EcalRecHitsEE"),
    FastJetTag = cms.untracked.InputTag("kt4CaloJets"),
    HFlongMin = cms.untracked.double(0.5),
    HFshortMin = cms.untracked.double(0.85),
    HFtowerMin = cms.untracked.double(3.0),
    JetSrc = cms.untracked.InputTag("ak4CaloJets"),
    TowerTreePtMin = cms.untracked.double(-99),
    bkg = cms.InputTag("voronoiBackgroundCalo"),
    doBasicClusters = cms.untracked.bool(False),
    doEbyEonly = cms.untracked.bool(False),
    doEcal = cms.untracked.bool(False),
    doFastJet = cms.untracked.bool(False),
    doHcal = cms.untracked.bool(False),
    doTowers = cms.untracked.bool(True),
    doUEraw_ = cms.untracked.bool(False),
    doVS = cms.untracked.bool(False),
    etaBins = cms.int32(15),
    fourierOrder = cms.int32(5),
    hasVtx = cms.untracked.bool(False),
    hcalHBHERecHitSrc = cms.untracked.InputTag("hbhereco","reducedHcalRecHits"),
    hcalHFRecHitSrc = cms.untracked.InputTag("hfreco","reducedHcalRecHits"),
    towersSrc = cms.untracked.InputTag("PFTowers"),
    useJets = cms.untracked.bool(False),
    vtxSrc = cms.untracked.InputTag("offlinePrimaryVerticesWithBS")
)
process.rechitanalyzer = cms.EDAnalyzer("RecHitTreeProducer",
    BasicClusterSrc1 = cms.untracked.InputTag("islandBasicClusters","islandBarrelBasicClusters"),
    EBRecHitSrc = cms.untracked.InputTag("reducedEcalRecHitsEB"),
    EBTreePtMin = cms.untracked.double(0),
    EERecHitSrc = cms.untracked.InputTag("reducedEcalRecHitsEE"),
    EETreePtMin = cms.untracked.double(0),
    FastJetTag = cms.untracked.InputTag("kt4CaloJets"),
    HBHETreePtMin = cms.untracked.double(0),
    HFTreePtMin = cms.untracked.double(0),
    HFlongMin = cms.untracked.double(0.5),
    HFshortMin = cms.untracked.double(0.85),
    HFtowerMin = cms.untracked.double(3.0),
    JetSrc = cms.untracked.InputTag("ak4CaloJets"),
    TowerTreePtMin = cms.untracked.double(-9999),
    bkg = cms.InputTag("voronoiBackgroundCalo"),
    calZDCDigi = cms.untracked.bool(False),
    doBasicClusters = cms.untracked.bool(False),
    doCASTOR = cms.untracked.bool(False),
    doEbyEonly = cms.untracked.bool(False),
    doEcal = cms.untracked.bool(True),
    doFastJet = cms.untracked.bool(False),
    doHF = cms.untracked.bool(True),
    doHcal = cms.untracked.bool(True),
    doTowers = cms.untracked.bool(True),
    doUEraw = cms.untracked.bool(False),
    doUEraw_ = cms.untracked.bool(False),
    doVS = cms.untracked.bool(False),
    doZDCDigi = cms.untracked.bool(False),
    doZDCRecHit = cms.untracked.bool(False),
    etaBins = cms.int32(15),
    fourierOrder = cms.int32(5),
    hasVtx = cms.untracked.bool(False),
    hcalHBHERecHitSrc = cms.untracked.InputTag("reducedHcalRecHits","hbhereco"),
    hcalHFRecHitSrc = cms.untracked.InputTag("reducedHcalRecHits","hfreco"),
    towersSrc = cms.untracked.InputTag("towerMaker"),
    useJets = cms.untracked.bool(False),
    vtxSrc = cms.untracked.InputTag("offlinePrimaryVerticesWithBS")
)

process.patElectronSequence = cms.Sequence(process.patElectrons)

process.ana_step = cms.Path(    process.UPCTrigger*
				process.hltanalysis+
				process.hltobject+
				process.hiEvtAnalyzer+
				process.fwdana+
				process.HiGenParticleAna+
				process.jetSequences+
				process.egmGsfElectronIDSequence+
				process.ggHiNtuplizer+
				process.ggHiNtuplizerGED+
				process.pfcandAnalyzer+
				process.rechitAna+
				process.HiForest+
				process.trackSequencesPP
)

process.load('HeavyIonsAnalysis.JetAnalysis.EventSelection_cff')
process.pHBHENoiseFilterResultProducer = cms.Path(process.HBHENoiseFilterResultProducer)
#process.UPCTrigger + process.HBHENoiseFilterResultProducer )
process.HBHENoiseFilterResult = cms.Path(process.fHBHENoiseFilterResult)
#process.UPCTrigger + process.fHBHENoiseFilterResult)
process.HBHENoiseFilterResultRun1 = cms.Path(process.fHBHENoiseFilterResultRun1) 
#process.UPCTrigger + process.fHBHENoiseFilterResultRun1)
process.HBHENoiseFilterResultRun2Loose = cms.Path(process.fHBHENoiseFilterResultRun2Loose) 
#process.UPCTrigger + process.fHBHENoiseFilterResultRun2Loose)
process.HBHENoiseFilterResultRun2Tight = cms.Path(process.fHBHENoiseFilterResultRun2Tight) 
#process.UPCTrigger + process.fHBHENoiseFilterResultRun2Tight)
process.HBHEIsoNoiseFilterResult = cms.Path(process.fHBHEIsoNoiseFilterResult) 
#process.UPCTrigger + process.fHBHEIsoNoiseFilterResult)

process.PAprimaryVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && abs(z) <= 25 && position.Rho <= 2 && tracksSize >= 2"),
    filter = cms.bool(True), # otherwise it won't filter the events
)

process.NoScraping = cms.EDFilter("FilterOutScraping",
 applyfilter = cms.untracked.bool(True),
 debugOn = cms.untracked.bool(False),
 numtrack = cms.untracked.uint32(10),
 thresh = cms.untracked.double(0.25)
)

process.pPAprimaryVertexFilter = cms.Path(process.PAprimaryVertexFilter) 
#process.UPCTrigger + process.PAprimaryVertexFilter)
process.pBeamScrapingFilter=cms.Path(process.NoScraping)
#process.UPCTrigger + process.NoScraping)
process.load("HeavyIonsAnalysis.VertexAnalysis.PAPileUpVertexFilter_cff")
process.pVertexFilterCutG = cms.Path(process.pileupVertexFilterCutG) 
#process.UPCTrigger + process.pileupVertexFilterCutG)
process.pVertexFilterCutGloose = cms.Path(process.pileupVertexFilterCutGloose) 
#process.UPCTrigger + process.pileupVertexFilterCutGloose)
process.pVertexFilterCutGtight = cms.Path(process.pileupVertexFilterCutGtight) 
#process.UPCTrigger + process.pileupVertexFilterCutGtight)
process.pVertexFilterCutGplus = cms.Path(process.pileupVertexFilterCutGplus) 
#process.UPCTrigger + process.pileupVertexFilterCutGplus)
process.pVertexFilterCutE = cms.Path(process.pileupVertexFilterCutE) 
#process.UPCTrigger + process.pileupVertexFilterCutE)
process.pVertexFilterCutEandG = cms.Path(process.pileupVertexFilterCutEandG) 
#process.UPCTrigger + process.pileupVertexFilterCutEandG)

process.pAna = cms.EndPath(process.skimanalysis) 
#process.UPCTrigger + process.skimanalysis)


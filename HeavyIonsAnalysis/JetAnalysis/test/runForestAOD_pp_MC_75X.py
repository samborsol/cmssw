### HiForest Configuration
# Collisions: pp
# Type: MC
# Input: AOD

import FWCore.ParameterSet.Config as cms
process = cms.Process('HiForest')
process.options = cms.untracked.PSet()

#####################################################################################
# HiForest labelling info
#####################################################################################

process.load("HeavyIonsAnalysis.JetAnalysis.HiForest_cff")
process.HiForest.inputLines = cms.vstring("HiForest V3",)
import subprocess
version = subprocess.Popen(["(cd $CMSSW_BASE/src && git describe --tags)"], stdout=subprocess.PIPE, shell=True).stdout.read()
if version == '':
    version = 'no git info'
process.HiForest.HiForestVersion = cms.string(version)

#####################################################################################
# Input source
#####################################################################################

process.source = cms.Source("PoolSource",
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
                            fileNames = cms.untracked.vstring(
                                "file:SingleGammaFlatPt0To50_pythia8_GENSIMRAWRECO.root"
                            )
)

# Number of events we want to process, -1 = all events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10))

import HLTrigger.HLTfilters.hltHighLevel_cfi

process.UPCTrigger = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.UPCTrigger.HLTPaths = [
                                "HLT_HIUPCL1DoubleEG2NotHF2_v1",
                                "HLT_HIUPCDoubleEG2NotHF2Pixel_SingleTrack_v1",
                                "HLT_HIUPCL1SingleEG5NotHF2_v1",
                                "HLT_HIUPCSingleEG5NotHF2Pixel_SingleTrack_v1"
]

#####################################################################################
# Load Global Tag, Geometry, etc.
#####################################################################################

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

##### library for forward analyzer :

process.load("RecoLocalCalo.HcalRecProducers.HcalHitReconstructor_hf_cfi")
process.load("RecoLocalCalo.HcalRecProducers.HcalHitReconstructor_zdc_cfi")
process.load("EventFilter.CastorRawToDigi.CastorRawToDigi_cfi")
process.load("RecoLocalCalo.CastorReco.CastorSimpleReconstructor_cfi")
process.load('RecoLocalTracker.SiPixelRecHits.PixelCPEESProducers_cff')
#process.fwdana = cms.EDAnalyzer('ForwardAnalyzer'
#                                , CentralitySrc = cms.InputTag("hiCentrality")
#                                , CentralityBinSrc = cms.InputTag("centralityBin","HFtowers")
#                                , trackSrc = cms.InputTag("hiGeneralTracks")
#                                , vtxCollection_ = cms.InputTag("hiSelectedVertex")
#)

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '75X_mcRun2_asymptotic_ppAt5TeV_v3', '')
process.HiForest.GlobalTagLabel = process.GlobalTag.globaltag


# Customization
from HeavyIonsAnalysis.Configuration.CommonFunctions_cff import overrideJEC_pp5020
process = overrideJEC_pp5020(process)

#####################################################################################
# Define tree output
#####################################################################################

process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string("HiForestAOD.root"))

#####################################################################################
# Additional Reconstruction and Analysis: Main Body
#####################################################################################

####################################################################################

#############################
# Jets
#############################

process.load("HeavyIonsAnalysis.JetAnalysis.FullJetSequence_nominalPP")

#process.load('HeavyIonsAnalysis.JetAnalysis.jets.HiReRecoJets_pp_cff')

process.load('HeavyIonsAnalysis.JetAnalysis.jets.akPu2CaloJetSequence_pp_mc_cff')
process.load('HeavyIonsAnalysis.JetAnalysis.jets.akVs2CaloJetSequence_pp_mc_cff')
process.load('HeavyIonsAnalysis.JetAnalysis.jets.akVs2PFJetSequence_pp_mc_cff')
process.load('HeavyIonsAnalysis.JetAnalysis.jets.akPu2PFJetSequence_pp_mc_cff')

process.load('HeavyIonsAnalysis.JetAnalysis.jets.akPu3CaloJetSequence_pp_mc_cff')
process.load('HeavyIonsAnalysis.JetAnalysis.jets.akVs3CaloJetSequence_pp_mc_cff')
process.load('HeavyIonsAnalysis.JetAnalysis.jets.akVs3PFJetSequence_pp_mc_cff')
process.load('HeavyIonsAnalysis.JetAnalysis.jets.akPu3PFJetSequence_pp_mc_cff')

process.load('HeavyIonsAnalysis.JetAnalysis.jets.akPu4CaloJetSequence_pp_mc_cff')
process.load('HeavyIonsAnalysis.JetAnalysis.jets.akVs4CaloJetSequence_pp_mc_cff')
process.load('HeavyIonsAnalysis.JetAnalysis.jets.akVs4PFJetSequence_pp_mc_cff')
process.load('HeavyIonsAnalysis.JetAnalysis.jets.akPu4PFJetSequence_pp_mc_cff')

process.load('HeavyIonsAnalysis.JetAnalysis.jets.akPu5CaloJetSequence_pp_mc_cff')
process.load('HeavyIonsAnalysis.JetAnalysis.jets.akVs5CaloJetSequence_pp_mc_cff')
process.load('HeavyIonsAnalysis.JetAnalysis.jets.akVs5PFJetSequence_pp_mc_cff')
process.load('HeavyIonsAnalysis.JetAnalysis.jets.akPu5PFJetSequence_pp_mc_cff')

process.highPurityTracks = cms.EDFilter("TrackSelector",
                                        src = cms.InputTag("hiGeneralTracks"),
                                        cut = cms.string('quality("highPurity")'))

highPurityTracks = cms.EDFilter("TrackSelector",
                                src = cms.InputTag("generalTracks"),
                                cut = cms.string('quality("highPurity")')
)

process.jetSequences = cms.Sequence(
     process.myPartons +
     process.genParticlesForJets +
     process.ak3GenJets +
     process.ak4GenJets +
     process.ak5GenJets +
     process.ak3PFJets +
     process.ak5PFJets +
     process.akSoftDrop4PFJets +
     process.akSoftDrop5PFJets +
     process.akFilter4PFJets +
     process.akFilter5PFJets +
     process.akSoftDrop4GenJets +
     process.akSoftDrop5GenJets +
     process.highPurityTracks +
     process.ak3PFJetSequence +
     process.ak4PFJetSequence +
     process.ak5PFJetSequence +
     process.ak4CaloJetSequence +
     process.akSoftDrop4PFJetSequence +
     process.akSoftDrop5PFJetSequence
)

# Use this version for JEC
#process.load("HeavyIonsAnalysis.JetAnalysis.FullJetSequence_JECPP")

#####################################################################################

############################
# Event Analysis
############################
process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cff')
from HeavyIonsAnalysis.EventAnalysis.dummybranches_cff import addHLTdummybranchesForPP
addHLTdummybranchesForPP(process)
#process.load('HeavyIonsAnalysis.EventAnalysis.hltobject_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_data_cfi') #use data version to avoid PbPb MC
process.hiEvtAnalyzer.Vertex = cms.InputTag("offlinePrimaryVertices")
process.hiEvtAnalyzer.doCentrality = cms.bool(False)
process.hiEvtAnalyzer.doEvtPlane = cms.bool(False)
process.hiEvtAnalyzer.doMC = cms.bool(True) #general MC info
process.hiEvtAnalyzer.doHiMC = cms.bool(False) #HI specific MC info

process.load('HeavyIonsAnalysis.JetAnalysis.HiGenAnalyzer_cfi')
process.HiGenParticleAna.genParticleSrc = cms.untracked.InputTag("genParticles")
process.HiGenParticleAna.doHI = False
process.load('HeavyIonsAnalysis.EventAnalysis.runanalyzer_cff')
process.load("HeavyIonsAnalysis.JetAnalysis.pfcandAnalyzer_pp_cfi")
process.pfcandAnalyzer.skipCharged = False
process.pfcandAnalyzer.pfPtMin = 0
process.pfcandAnalyzer.pfCandidateLabel = cms.InputTag("particleFlow")
process.pfcandAnalyzer.doVS = cms.untracked.bool(False)
process.pfcandAnalyzer.doUEraw_ = cms.untracked.bool(False)
process.pfcandAnalyzer.genLabel = cms.InputTag("genParticles")

process.load('HeavyIonsAnalysis.JetAnalysis.rechitanalyzer_cfi')
process.rechitanalyzer.EBTreePtMin = -9999
process.rechitanalyzer.EETreePtMin = -9999
process.rechitanalyzer.HBHETreePtMin = -9999
process.rechitanalyzer.HFTreePtMin = -9999
process.rechitanalyzer.HFlongMin = -9999
process.rechitanalyzer.HFshortMin = -9999
process.rechitanalyzer.HFtowerMin = -9999

process.rechitAna = cms.Sequence(process.rechitanalyzer)

#####################################################################################

#########################
# Track Analyzer
#########################
process.load('HeavyIonsAnalysis.JetAnalysis.ExtraTrackReco_cff')
process.load('HeavyIonsAnalysis.JetAnalysis.TrkAnalyzers_cff')

# Use this instead for track corrections
## process.load('HeavyIonsAnalysis.JetAnalysis.TrkAnalyzers_Corr_cff')

#####################################################################################

#####################
# photons
######################
process.load('HeavyIonsAnalysis.PhotonAnalysis.ggHiNtuplizer_cfi')
process.ggHiNtuplizer.gsfElectronLabel   = cms.InputTag("gedGsfElectrons")
process.ggHiNtuplizer.recoPhotonHiIsolationMap = cms.InputTag('photonIsolationHIProducerpp')
process.ggHiNtuplizer.VtxLabel           = cms.InputTag("offlinePrimaryVertices")
process.ggHiNtuplizer.particleFlowCollection = cms.InputTag("particleFlow")
process.ggHiNtuplizer.doVsIso            = cms.bool(False)
process.ggHiNtuplizer.doElectronVID      = cms.bool(True)
process.ggHiNtuplizerGED = process.ggHiNtuplizer.clone(recoPhotonSrc = cms.InputTag('gedPhotons'),
                                                       recoPhotonHiIsolationMap = cms.InputTag('photonIsolationHIProducerppGED'))

####################################################################################
#####################
# Electron ID
#####################

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format to be processed
# DataFormat.AOD or DataFormat.MiniAOD
dataFormat = DataFormat.AOD
switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce. Check here https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2#Recipe_for_regular_users_for_7_4
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
####################################################################################


#####################
# tupel and necessary PAT sequences
#####################

process.load("HeavyIonsAnalysis.VectorBosonAnalysis.tupelSequence_pp_mc_cff")

#####################################################################################

#########################
# Main analysis list
#########################
process.ana_step = cms.Path(process.hltanalysis *
                            #process.hltobject *
                            process.hiEvtAnalyzer *
                            #process.fwdana *
                            process.HiGenParticleAna *
                            process.UPCTrigger *
                            process.jetSequences +
                            process.egmGsfElectronIDSequence + #Should be added in the path for VID module
                            process.ggHiNtuplizer +
                            process.ggHiNtuplizerGED +
                            process.pfcandAnalyzer +
                            process.rechitAna +
                            process.HiForest +
			    process.trackSequencesPP +
                            process.runAnalyzer +
                            process.tupelPatSequence
)

#####################################################################################

#########################
# Event Selection
#########################

process.load('HeavyIonsAnalysis.JetAnalysis.EventSelection_cff')
process.pHBHENoiseFilterResultProducer = cms.Path( process.UPCTrigger * process.HBHENoiseFilterResultProducer )
process.HBHENoiseFilterResult = cms.Path(process.UPCTrigger * process.fHBHENoiseFilterResult)
process.HBHENoiseFilterResultRun1 = cms.Path(process.UPCTrigger * process.fHBHENoiseFilterResultRun1)
process.HBHENoiseFilterResultRun2Loose = cms.Path(process.UPCTrigger * process.fHBHENoiseFilterResultRun2Loose)
process.HBHENoiseFilterResultRun2Tight = cms.Path(process.UPCTrigger * process.fHBHENoiseFilterResultRun2Tight)
process.HBHEIsoNoiseFilterResult = cms.Path(process.UPCTrigger * process.fHBHEIsoNoiseFilterResult)

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

process.pPAprimaryVertexFilter = cms.Path(process.UPCTrigger * process.PAprimaryVertexFilter)
process.pBeamScrapingFilter=cms.Path(process.UPCTrigger * process.NoScraping)

process.load("HeavyIonsAnalysis.VertexAnalysis.PAPileUpVertexFilter_cff")

process.pVertexFilterCutG = cms.Path(process.UPCTrigger * process.pileupVertexFilterCutG)
process.pVertexFilterCutGloose = cms.Path(process.UPCTrigger * process.pileupVertexFilterCutGloose)
process.pVertexFilterCutGtight = cms.Path(process.UPCTrigger * process.pileupVertexFilterCutGtight)
process.pVertexFilterCutGplus = cms.Path(process.UPCTrigger * process.pileupVertexFilterCutGplus)
process.pVertexFilterCutE = cms.Path(process.UPCTrigger * process.pileupVertexFilterCutE)
process.pVertexFilterCutEandG = cms.Path(process.UPCTrigger * process.pileupVertexFilterCutEandG)

process.pAna = cms.EndPath(process.UPCTrigger * process.skimanalysis)

# Customization

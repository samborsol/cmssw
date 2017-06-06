### HiForest Configuration
# Collisions: pp
# Type: MC
# Input: AOD

import FWCore.ParameterSet.Config as cms
import sys

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
                                sys.argv[2]
                            )
)

# Number of events we want to process, -1 = all events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1))


#####################################################################################
# Load Global Tag, Geometry, etc.
#####################################################################################

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '75X_mcRun2_HeavyIon_v13', '')
process.HiForest.GlobalTagLabel = process.GlobalTag.globaltag


# Customization
from HeavyIonsAnalysis.Configuration.CommonFunctions_cff import overrideJEC_pp5020
process = overrideJEC_pp5020(process)

#####################################################################################
# Define tree output
#####################################################################################

process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string(sys.argv[3]))

#####################################################################################
# Additional Reconstruction and Analysis: Main Body
#####################################################################################

####################################################################################

#############################
# Jets
#############################
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



#process.schedule.extend([process.raw2digi_step,process.reconstruction_step,process.endjob_step,process.AODSIMoutput_step
# Include this to turn on storing the jet constituents and new jet variables for q/g separation
#process.ak4PFJetAnalyzer.doJetConstituents = cms.untracked.bool(True)
#process.ak4PFJetAnalyzer.doNewJetVars = cms.untracked.bool(True)
# Use this version for JEC
#process.load("HeavyIonsAnalysis.JetAnalysis.FullJetSequence_JECPP")

#####################################################################################

############################
# Event Analysis
############################
process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cff')
process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_data_cfi') #use data version to avoid PbPb MC
process.hiEvtAnalyzer.Vertex = cms.InputTag("offlinePrimaryVertices")
process.hiEvtAnalyzer.doCentrality = cms.bool(False)
process.hiEvtAnalyzer.doEvtPlane = cms.bool(False)
process.hiEvtAnalyzer.doMC = cms.bool(True) #general MC info
process.hiEvtAnalyzer.doHiMC = cms.bool(False) #HI specific MC info
'''
process.load('HeavyIonsAnalysis.JetAnalysis.HiGenAnalyzer_cfi')
process.HiGenParticleAna.genParticleSrc = cms.untracked.InputTag("genParticles")
process.HiGenParticleAna.doHI = False
'''
#Gen particles
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

#process.load('HeavyIonsAnalysis.JetAnalysis.rechitanalyzer_cfi')
process.load('HeavyIonsAnalysis.JetAnalysis.rechitanalyzer_aod_cfi')

##### modifications to rechitanalyzer inputs #####

process.load('HeavyIonsAnalysis.JetAnalysis.rechitanalyzer_cfi')
process.rechitanalyzer.EBRecHitSrc = cms.untracked.InputTag("reducedEcalRecHitsEB")
process.rechitanalyzer.EERecHitSrc = cms.untracked.InputTag("reducedEcalRecHitsEE")
process.rechitanalyzer.hcalHFRecHitSrc = cms.untracked.InputTag("reducedHcalRecHits", "hfreco")
process.rechitanalyzer.hcalHBHERecHitSrc = cms.untracked.InputTag("reducedHcalRecHits", "hbhereco")
process.rechitanalyzer.hasVtx = cms.untracked.bool(True)
process.rechitanalyzer.useJets = cms.untracked.bool(True)
process.rechitanalyzer.doBasicClusters = cms.untracked.bool(True)
process.rechitanalyzer.doTowers = cms.untracked.bool(True)
process.rechitanalyzer.doCASTOR = cms.untracked.bool(False)
process.rechitanalyzer.doZDCDigi = cms.untracked.bool(False)
process.rechitanalyzer.doZDCRecHit = cms.untracked.bool(False)
process.rechitanalyzer.doVS = cms.untracked.bool(False)
process.rechitanalyzer.doUEraw = cms.untracked.bool(True)
process.rechitanalyzer.doFastJet = cms.untracked.bool(False)
process.rechitanalyzer.calZDCDigi = cms.untracked.bool(False)
process.rechitanalyzer.doHF = cms.untracked.bool(True)
process.rechitanalyzer.vtxSrc = cms.untracked.InputTag("offlinePrimaryVerticesWithBS")
process.rechitanalyzer.JetSrc = cms.untracked.InputTag("ak4CaloJets")
process.rechitAna = cms.Sequence(process.rechitanalyzer)

##### modifications to rechitanalyzer inputs - END #####

process.load("HeavyIonsAnalysis.MuonAnalysis.hltMuTree_cfi")
process.hltMuTree.vertices = cms.InputTag("offlinePrimaryVertices")


process.demo = cms.EDAnalyzer('TestSolution')

#########################
# Main analysis list
#########################
process.ana_step = cms.Path(process.hltanalysis *
			    process.demo *
                            process.hiEvtAnalyzer *
                            process.HiGenParticleAna*
                            process.jetSequences +
                            #process.egmGsfElectronIDSequence + #Should be added in the path for VID module
                            process.ggHiNtuplizer +
                            process.ggHiNtuplizerGED +
                            process.pfcandAnalyzer +
                            process.HiForest +
			    process.trackSequencesPP +
                            process.runAnalyzer +
                            process.hltMuTree + 
			    process.rechitanalyzer
                            #process.tupelPatSequence
)

#####################################################################################

#########################
# Event Selection
#########################
'''
process.load('HeavyIonsAnalysis.JetAnalysis.EventSelection_cff')
process.pHBHENoiseFilterResultProducer = cms.Path( process.HBHENoiseFilterResultProducer )
process.HBHENoiseFilterResult = cms.Path(process.fHBHENoiseFilterResult)
process.HBHENoiseFilterResultRun1 = cms.Path(process.fHBHENoiseFilterResultRun1)
process.HBHENoiseFilterResultRun2Loose = cms.Path(process.fHBHENoiseFilterResultRun2Loose)
process.HBHENoiseFilterResultRun2Tight = cms.Path(process.fHBHENoiseFilterResultRun2Tight)
process.HBHEIsoNoiseFilterResult = cms.Path(process.fHBHEIsoNoiseFilterResult)

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
process.pBeamScrapingFilter=cms.Path(process.NoScraping)

process.load("HeavyIonsAnalysis.VertexAnalysis.PAPileUpVertexFilter_cff")

process.pVertexFilterCutG = cms.Path(process.pileupVertexFilterCutG)
process.pVertexFilterCutGloose = cms.Path(process.pileupVertexFilterCutGloose)
process.pVertexFilterCutGtight = cms.Path(process.pileupVertexFilterCutGtight)
process.pVertexFilterCutGplus = cms.Path(process.pileupVertexFilterCutGplus)
process.pVertexFilterCutE = cms.Path(process.pileupVertexFilterCutE)
process.pVertexFilterCutEandG = cms.Path(process.pileupVertexFilterCutEandG)

process.pAna = cms.EndPath(process.skimanalysis)
'''
# Customization

import FWCore.ParameterSet.Config as cms

### PP RECO does not include R=3 or R=5 jets.
### re-RECO is only possible for PF, RECO is missing calotowers
from RecoJets.JetProducers.ak5PFJets_cfi import ak5PFJets
ak5PFJets.doAreaFastjet = True
ak3PFJets = ak5PFJets.clone(rParam = 0.3)
from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets
ak3GenJets = ak5GenJets.clone(rParam = 0.3)
ak4GenJets = ak5GenJets.clone(rParam = 0.4)

#SoftDrop PF jets
from RecoJets.JetProducers.PFJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
akSoftDrop4PFJets = cms.EDProducer(
    "FastjetJetProducer",
    PFJetParameters,
    AnomalousCellParameters,
    jetAlgorithm = cms.string("AntiKt"),
    rParam       = cms.double(0.4),
    useSoftDrop = cms.bool(True),
    zcut = cms.double(0.1),
    beta = cms.double(0.0),
    R0   = cms.double(0.4),
    useExplicitGhosts = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
akSoftDrop5PFJets = akSoftDrop4PFJets.clone(rParam = cms.double(0.5), R0 = cms.double(0.5))

from HeavyIonsAnalysis.JetAnalysis.akSoftDrop4GenJets_cfi import akSoftDrop4GenJets
akSoftDrop5GenJets = akSoftDrop4GenJets.clone(rParam = 0.5)

#Filter PF jets
akFilter4PFJets = cms.EDProducer(
    "FastjetJetProducer",
    PFJetParameters,
    AnomalousCellParameters,
    jetAlgorithm = cms.string("AntiKt"),
    rParam       = cms.double(0.4),
    useFiltering = cms.bool(True),
    nFilt = cms.int32(4),
    rFilt = cms.double(0.15),
    useExplicitGhosts = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
akFilter5PFJets = akFilter4PFJets.clone(rParam = cms.double(0.5))

from RecoJets.Configuration.GenJetParticles_cff import *
from RecoHI.HiJetAlgos.HiGenJets_cff import *
from HeavyIonsAnalysis.JetAnalysis.makePartons_cff import myPartons

from HeavyIonsAnalysis.JetAnalysis.jets.ak3PFJetSequence_pp_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.ak4PFJetSequence_pp_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.ak5PFJetSequence_pp_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.ak4CaloJetSequence_pp_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akSoftDrop4PFJetSequence_pp_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akSoftDrop5PFJetSequence_pp_mc_cff import *

highPurityTracks = cms.EDFilter("TrackSelector",
                                src = cms.InputTag("generalTracks"),
                                cut = cms.string('quality("highPurity")')
)

# Other radii jets and calo jets need to be reconstructed
jetSequences = cms.Sequence(
    myPartons +
    genParticlesForJets +
    ak3GenJets +
    ak4GenJets +
    ak5GenJets +
    ak3PFJets +
    ak5PFJets +
    akSoftDrop4PFJets +
    akSoftDrop5PFJets +
    akFilter4PFJets +
    akFilter5PFJets +
    akSoftDrop4GenJets +
    akSoftDrop5GenJets +
    highPurityTracks +
    ak3PFJetSequence +
    ak4PFJetSequence +
    ak5PFJetSequence +
    ak4CaloJetSequence +
    akSoftDrop4PFJetSequence +
    akSoftDrop5PFJetSequence 
)

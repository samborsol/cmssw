import FWCore.ParameterSet.Config as cms

source = cms.Source("MCFileSource",
    fileNames = cms.untracked.vstring(
       'file:FPMC_WW_Inclusive_13TeV.hepmc'
    )
)

generator = cms.EDFilter("inputHepMC")

ProductionFilterSequence = cms.Sequence(generator)

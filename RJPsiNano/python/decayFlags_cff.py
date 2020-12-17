import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

decayFlagsTable = cms.EDProducer( "DecayFlagsTableProducer",
                                  srcGen = cms.InputTag("prunedGenParticles"),
                                  objName = cms.string("DecayFlag"),
                                  objType = cms.string("DecayFlag"),
                                  branchName = cms.string("GenDecay"),
                                  src = cms.string("GenDecay"),
                                )

decayFlagsTables = cms.Sequence(decayFlagsTable)


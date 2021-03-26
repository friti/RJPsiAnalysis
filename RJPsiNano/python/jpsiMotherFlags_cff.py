import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

jpsiMotherFlagsTable = cms.EDProducer( "jpsiMotherFlagsTableProducer",
                                  prunedGen = cms.InputTag("prunedGenParticles"),
                                  packedGen = cms.InputTag("packedGenParticles"),
                                  objName = cms.string("jpsiMotherFlag"),
                                  objType = cms.string("jpsiMotherFlag"),
                                  branchName = cms.string("GenDecay"),
                                  src = cms.string("GenDecay"),
                                )

jpsiMotherFlagsTables = cms.Sequence(jpsiMotherFlagsTable)


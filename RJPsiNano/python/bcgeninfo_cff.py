import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

BcGenInfoTable = cms.EDProducer( "BcGenInfoTableProducer",
                                  srcGen = cms.InputTag("prunedGenParticles"),
                                  objName = cms.string("BcGenInfo"),
                                  objType = cms.string("BcGenInfo"),
                                  branchName = cms.string("BcGenInfo"),
                                  src = cms.string("BcGenInfo"),
                                )

BcGenInfoTables = cms.Sequence(BcGenInfoTable)


import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.RJPsiNano.common_cff import RJpsiCandVars

HighMassLowMassFlags = cms.EDProducer("HighMassLowMassbkgFlagsProducer",
                                  prunedGen = cms.InputTag("prunedGenParticles"),
                                  packedGen = cms.InputTag("packedGenParticles"),
                                  #objName = cms.string("HighMassLowMassFlag"),
                                  #objType = cms.string("HighMassLowMassFlag"),
                                  #branchName = cms.string("GenDecay"),
                                  #src = cms.string("GenDecay"),
                                )

HighMassLowMassFlagsTable = cms.EDProducer("SimpleCompositeCandidateFlatTableProducer",
                                           src = cms.InputTag("HighMassLowMassFlags","AncestorFlags"),
                                           cut = cms.string(""),
                                           name = cms.string("HighMassLowMassFlags"),
                                           doc = cms.string("doc"),
                                           singleton = cms.bool(False), 
                                           extension = cms.bool(False), 
                                           variables = cms.PSet(RJpsiCandVars,
                                               mu1_idx = Var("userInt('mu1_idx')", int, doc = 'gen index of the mu1 from the jpsi'),
                                               mu2_idx = Var("userInt('mu2_idx')", int, doc = 'gen index of the mu2 from the jpsi'),
                                               mu3_idx = Var("userInt('mu3_idx')", int, doc = 'gen index of the mu3 from the jpsi'),
                                               jpsi_idx = Var("userInt('jpsi_idx')", int, doc = 'gen index of  the jpsi'),
                                               jpsi_ancestor_idx = Var("userInt('jpsi_ancestor_idx')", int, doc = 'gen index of the ancestor of the jpsi'),
                                               mu3_ancestor_idx = Var("userInt('mu3_ancestor_idx')", int, doc = 'gen index of the ancestor of the mu3'),),)


HighMassLowMassFlagsTables = cms.Sequence(HighMassLowMassFlags * HighMassLowMassFlagsTable)
#HighMassLowMassFlagsTables = cms.Sequence(HighMassLowMassFlags)

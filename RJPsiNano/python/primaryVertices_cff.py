import FWCore.ParameterSet.Config as cms
from PhysicsTools.RJPsiNano.common_cff import *

pvSelector = cms.EDProducer("PrimaryVertexSelector",
                            vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
                            dimuonCollection = cms.InputTag("JpsiMuonPairs", "muonPairsForBTo3Mu"),
                            dimuonTTCollection = cms.InputTag("JpsiMuonPairs", "dimuonTransientTracks"),
                            dzForCleaning_wrtSV = cms.double(0.0)
                            )

pvSelectorSequence = cms.Sequence(pvSelector)


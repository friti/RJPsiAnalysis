import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var, CandVars
from PhysicsTools.RJPsiNano.common_cff import RJpsiCandVars, ufloat, uint, ubool
from PhysicsTools.RJPsiNano.rjpsi_common_cff import JpsiMuonPairs, BuilderDefaultCfg, TableDefaultVariables, TableDefault, Final3PartTableVariables
from PhysicsTools.RJPsiNano.primaryVertices_cff import *

BTo2MuKCfg = BuilderDefaultCfg.clone()
BTo2MuKCfg.dileptons = cms.InputTag('JpsiMuonPairs')
BTo2MuKCfg.leptonTransientTracks = JpsiMuonPairs.transientTracksSrc
BTo2MuKCfg.postVtxSelection = cms.string(' && '.join([
        BuilderDefaultCfg.postVtxSelection.value(),
        'mass > 4.5',
        ])
)

BTo2MuK = cms.EDProducer(
    'BTo2MuTkBuilder',
    BTo2MuKCfg,
    srcGen = cms.InputTag("prunedGenParticles"),
    particle_mass            = cms.double(0.493677),
    particles = cms.InputTag('tracksBPark', 'SelectedTracks'),
    particlesTransientTracks  = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    particleSelection         = cms.string(''),
)

BTo2MuKTableVariables = Final3PartTableVariables.clone()

BTo2MuKTable = TableDefault.clone()
BTo2MuKTable.src       = cms.InputTag("BTo2MuK")
BTo2MuKTable.name      = cms.string("BTo2MuK")
BTo2MuKTable.doc       = cms.string("BTo2MuK Variable")
BTo2MuKTable.variables = BTo2MuKTableVariables

BTo2MuKSequence = cms.Sequence(
    (JpsiMuonPairs*pvSelector * BTo2MuK)
)


CountBTo2MuK = cms.EDFilter(
    "PATCandViewCountFilter",
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("BTo2MuK")
)

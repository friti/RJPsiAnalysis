import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var, CandVars
from PhysicsTools.RJPsiNano.common_cff import RJpsiCandVars, ufloat, uint, ubool
from PhysicsTools.RJPsiNano.rjpsi_common_cff import JpsiMuonPairs, BuilderDefaultCfg, TableDefaultVariables, TableDefault, Final3PartTableVariables
from PhysicsTools.RJPsiNano.primaryVertices_cff import *

BTo2MuPCfg = BuilderDefaultCfg.clone()
#BTo2MuPCfg.dimuons = cms.InputTag('JpsiMuonPairs')
BTo2MuPCfg.muonTransientTracks = JpsiMuonPairs.transientTracksSrc
BTo2MuPCfg.postVtxSelection = cms.string(' && '.join([
        BuilderDefaultCfg.postVtxSelection.value(),
        'mass > 4.5',
        ])
)

BTo2MuP = cms.EDProducer(
    'BTo2MuTkBuilder',
    BTo2MuPCfg,
    srcGen = cms.InputTag("prunedGenParticles"),
    particle_mass            = cms.double(0.139571),
    particles = cms.InputTag('tracksBPark', 'SelectedTracks'),
    particlesTransientTracks  = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    particleSelection         = cms.string(''),
)

BTo2MuPTableVariables = Final3PartTableVariables.clone()

BTo2MuPTable = TableDefault.clone()
BTo2MuPTable.src       = cms.InputTag("BTo2MuP")
BTo2MuPTable.name      = cms.string("BTo2MuP")
BTo2MuPTable.doc       = cms.string("BTo2MuP Variable")
BTo2MuPTable.variables = BTo2MuPTableVariables

BTo2MuPSequence = cms.Sequence(
    (JpsiMuonPairs * pvSelector * BTo2MuP)
)


CountBTo2MuP = cms.EDFilter(
    "PATCandViewCountFilter",
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("BTo2MuP")
)  

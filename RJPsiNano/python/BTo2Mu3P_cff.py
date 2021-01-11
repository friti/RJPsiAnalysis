import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var, CandVars
from PhysicsTools.RJPsiNano.common_cff import RJpsiCandVars, ufloat, uint, ubool
from PhysicsTools.RJPsiNano.rjpsi_common_cff import JpsiMuonPairs, BuilderDefaultCfg, TableDefaultVariables, TableDefault
from PhysicsTools.RJPsiNano.primaryVertices_cff import *

BTo2Mu3PCfg = BuilderDefaultCfg.clone()
#BTo2Mu3PCfg.dimuons = cms.InputTag('JpsiMuonPairs')
BTo2Mu3PCfg.muonTransientTracks = JpsiMuonPairs.transientTracksSrc
BTo2Mu3PCfg.postVtxSelection = cms.string(' && '.join([
        BuilderDefaultCfg.postVtxSelection.value(),
        'mass > 4.5',
        ])
)

BTo2Mu3PCfg.kaonSelection = cms.string(' && '.join([
    #    BuilderDefaultCfg.kaonSelection.value(),
    'pt > 1',
    'eta < 2.5',
        ])
)

BTo2Mu3P = cms.EDProducer(
    'BTo2Mu3PiBuilder',
    BTo2Mu3PCfg,
    srcGen = cms.InputTag("prunedGenParticles"),
    particles = cms.InputTag('tracksBPark', 'SelectedTracks'),
    particlesTransientTracks  = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    particleSelection         = cms.string(''),
)

BTo2Mu3PTableVariables = TableDefaultVariables.clone(
    pi1_dxy     = ufloat('pi1_dxy'),
    pi1_dz     = ufloat('pi1_dz'),
    pi2_dxy     = ufloat('pi2_dxy'),
    pi2_dz     = ufloat('pi2_dz'),
    pi3_dxy     = ufloat('pi3_dxy'),
    pi3_dz     = ufloat('pi3_dz'),
    pi1_dxyErr     = ufloat('pi1_dxyErr'),
    pi1_dzErr     = ufloat('pi1_dzErr'),
    pi2_dxyErr     = ufloat('pi2_dxyErr'),
    pi2_dzErr     = ufloat('pi2_dzErr'),
    pi3_dxyErr     = ufloat('pi3_dxyErr'),
    pi3_dzErr     = ufloat('pi3_dzErr'),

    pi1Idx     = uint('pi1_idx'),
    pi2Idx     = uint('pi2_idx'),
    pi3Idx     = uint('pi3_idx'),
    bodies3_fit_pi1_pt    = ufloat('fitted_pi1_pt'),
    bodies3_fit_pi1_eta   = ufloat('fitted_pi1_eta'),
    bodies3_fit_pi1_phi   = ufloat('fitted_pi1_phi'),
    bodies3_fit_pi2_pt    = ufloat('fitted_pi2_pt'),
    bodies3_fit_pi2_eta   = ufloat('fitted_pi2_eta'),
    bodies3_fit_pi2_phi   = ufloat('fitted_pi2_phi'),
    bodies3_fit_pi3_pt    = ufloat('fitted_pi3_pt'),
    bodies3_fit_pi3_eta   = ufloat('fitted_pi3_eta'),
    bodies3_fit_pi3_phi   = ufloat('fitted_pi3_phi'),
    pi1_iso03     = ufloat('pi1_iso03'),
    pi1_iso04     = ufloat('pi1_iso04'),
    pi2_iso03     = ufloat('pi2_iso03'),
    pi2_iso04     = ufloat('pi2_iso04'),
    pi3_iso03     = ufloat('pi3_iso03'),
    pi3_iso04     = ufloat('pi3_iso04'),

)

BTo2Mu3PTable = TableDefault.clone()
BTo2Mu3PTable.src       = cms.InputTag("BTo2Mu3P")
BTo2Mu3PTable.name      = cms.string("BTo2Mu3P")
BTo2Mu3PTable.doc       = cms.string("BTo2Mu3P Variable")
BTo2Mu3PTable.variables = BTo2Mu3PTableVariables

BTo2Mu3PSequence = cms.Sequence(
    (JpsiMuonPairs * pvSelector * BTo2Mu3P)
)


CountBTo2Mu3P = cms.EDFilter(
    "PATCandViewCountFilter",
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("BTo2Mu3P")
)  

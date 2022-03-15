import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var, CandVars
from PhysicsTools.RJPsiNano.common_cff import RJpsiCandVars, ufloat, uint, ubool
from PhysicsTools.RJPsiNano.rjpsi_common_cff import JpsiMuonPairs, BuilderDefaultCfg, TableDefaultVariables, TableDefault,Final3PartTableVariables, Final3MuonsTableVariables
from PhysicsTools.RJPsiNano.primaryVertices_cff import *

BTo3MuCfg = BuilderDefaultCfg.clone()
#BTo3MuCfg.dileptons             = cms.InputTag('JpsiMuonPairs')
#BTo3MuCfg.leptonTransientTracks = JpsiMuonPairs.transientTracksSrc

BTo3Mu = cms.EDProducer(
    'BTo3MuBuilder',
    BTo3MuCfg,
    muons = cms.InputTag('muonTrgSelector', 'SelectedMuons'),
    muonSelection = cms.string(''),
    srcGen = cms.InputTag("prunedGenParticles"),
    dimuons_fitter = cms.InputTag('JpsiMuonPairs','dimuonFitter')
)

BTo3MuTableVariables = Final3MuonsTableVariables.clone()

BTo3MuTable = TableDefault.clone()
BTo3MuTable.src       = cms.InputTag("BTo3Mu")
BTo3MuTable.cut       = cms.string("")
BTo3MuTable.name      = cms.string("BTo3Mu")
BTo3MuTable.doc       = cms.string("BTo3Mu Variable")
BTo3MuTable.singleton = cms.bool(False)
BTo3MuTable.extension = cms.bool(False)
BTo3MuTable.variables = BTo3MuTableVariables

BTo3MuSequence = cms.Sequence(
    (JpsiMuonPairs * pvSelector * BTo3Mu)
)

CountBTo3Mu = cms.EDFilter(
    "PATCandViewCountFilter",
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("BTo3Mu")
)

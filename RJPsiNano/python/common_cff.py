import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var, CandVars

def ufloat(expr, precision = -1, doc = ''):
  return Var('userFloat("%s")' % expr, 
             float, precision = precision, doc = doc)

def uint(expr, doc = ''):
  return Var('userInt("%s")' % expr, int, doc = doc)

def ubool(expr, precision = -1, doc = ''):
  return Var('userInt("%s") == 1' % expr, bool, doc = doc)

RJpsiCandVars = CandVars.clone()
RJpsiCandVars.charge.precision = cms.int32(-1)
RJpsiCandVars.eta   .precision = cms.int32(-1)
RJpsiCandVars.mass  .precision = cms.int32(-1)
RJpsiCandVars.pdgId .precision = cms.int32(-1)
RJpsiCandVars.phi   .precision = cms.int32(-1)
RJpsiCandVars.pt    .precision = cms.int32(-1)

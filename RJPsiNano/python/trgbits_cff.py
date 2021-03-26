import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *



trgTable = cms.EDProducer( "TrgBitTableProducer",
                          hltresults = cms.InputTag("TriggerResults::HLT"),
                          l1results  = cms.InputTag("gtStage2Digis::RECO"),
                          #add interesting paths
                          paths      = cms.vstring(
                                             "HLT_DoubleMu4_JpsiTrk_Displaced",
                                             "HLT_Dimuon0_Jpsi3p5_Muon2",
                                             "HLT_DoubleMu4_LowMassNonResonantTrk_Displaced",
                                             "HLT_DoubleMu4_PsiPrimeTrk_Displaced"
                                              ),
                           #add interesting seeds
                           seeds     = cms.vstring(
                                              ),
                            
)

trgTables = cms.Sequence(trgTable)




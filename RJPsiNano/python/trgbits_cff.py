import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *



trgTable = cms.EDProducer( "TrgBitTableProducer",
                          hltresults = cms.InputTag("TriggerResults::HLT"),
                          l1results  = cms.InputTag("gtStage2Digis::RECO"),
                          #add interesting paths
                          paths      = cms.vstring(
                                             "HLT_DoubleMu4_JpsiTrk_Displaced"
                                              ),
                           #add interesting seeds
                           seeds     = cms.vstring(
                                              ),
                            
)

trgTables = cms.Sequence(trgTable)




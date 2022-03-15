from __future__ import print_function
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *   #why do we need all these nanoaod functions?
from PhysicsTools.NanoAOD.globals_cff import *
from PhysicsTools.NanoAOD.nano_cff import *
from PhysicsTools.NanoAOD.vertices_cff import *
from PhysicsTools.NanoAOD.NanoAODEDMEventContent_cff import *
from PhysicsTools.RJPsiNano.trgbits_cff import *



##for gen and trigger muon
from PhysicsTools.RJPsiNano.genparticlesBPark_cff import *
from PhysicsTools.RJPsiNano.particlelevelBPark_cff import *
from PhysicsTools.RJPsiNano.decayFlags_cff import *
from PhysicsTools.RJPsiNano.bcgeninfo_cff import *
from PhysicsTools.RJPsiNano.jpsiMotherFlags_cff import *
from PhysicsTools.RJPsiNano.HighMassLowMassbkg_cff import *
#G: from PhysicsTools.RJPsiNano.triggerObjectsBPark_cff import *
from PhysicsTools.RJPsiNano.muonsBPark_cff import * 

## filtered input collections
#from PhysicsTools.RJPsiNano.electronsBPark_cff import * 
from PhysicsTools.RJPsiNano.tracksBPark_cff import *

## B collections
from PhysicsTools.RJPsiNano.BTo3Mu_cff import *
from PhysicsTools.RJPsiNano.BTo2MuK_cff import *
from PhysicsTools.RJPsiNano.BTo2MuP_cff import *
from PhysicsTools.RJPsiNano.BTo2Mu3P_cff import *


#G: nanoSequenceOnlyFullSim = cms.Sequence(triggerObjectBParkTables + l1bits)  #purpose?

# from PhysiscsTools.NanoAOD
nanoSequence = cms.Sequence(nanoMetadata + globalTables)# + 
                            #vertexSequence)# +           
                            #vertexTables)# + 
                            #globalTables + vertexTables)# + 
                            #l1bits)
                            #triggerObjectBParkTables + l1bits)

nanoSequenceMC = cms.Sequence(particleLevelBParkSequence + genParticleBParkSequence + 
                              genParticleBParkTables + decayFlagsTables + jpsiMotherFlagsTables + BcGenInfoTables + HighMassLowMassFlagsTables +# + lheInfoTable) 
                              globalTablesMC) # + genWeightsTable + genParticleBParkTables + lheInfoTable) 


def nanoAOD_customizeMuonTriggerBPark(process):  #delete trigger inside
    #process.nanoSequence = cms.Sequence( process.nanoSequence + muonBParkSequence + muonTriggerMatchedTables + muonBParkTables)
    process.nanoSequence = cms.Sequence( process.nanoSequence + muonBParkSequence + muonBParkTables)
    return process

def nanoAOD_customizeTrackFilteredBPark(process): #not needed now, but useful for Bc->jpsi pi
    process.nanoSequence = cms.Sequence( process.nanoSequence + tracksBParkSequence + tracksBParkTables)
    return process

#def nanoAOD_customizeElectronFilteredBPark(process): #can be useful
#    process.nanoBKeeSequence = cms.Sequence( electronsBParkSequence + electronBParkTables)
#    return process

def nanoAOD_customizeTriggerBitsBPark(process): #needed??
    process.nanoSequence = cms.Sequence( process.nanoSequence + trgTables)
    return process

def nanoAOD_customizeBTo3Mu(process):
    process.nanoBTo3MuSequence = cms.Sequence( BTo3MuSequence + BTo3MuTable )#+HighMassLowMassFlagsTables )
    return process

def nanoAOD_customizeBTo2MuK(process):
    process.nanoBTo2MuKSequence = cms.Sequence( BTo2MuKSequence + BTo2MuKTable )
    return process

def nanoAOD_customizeBTo2MuP(process):
    process.nanoBTo2MuPSequence = cms.Sequence( BTo2MuPSequence + BTo2MuPTable )
    return process

def nanoAOD_customizeBTo2Mu3P(process):
    process.nanoBTo2Mu3PSequence = cms.Sequence( BTo2Mu3PSequence + BTo2Mu3PTable )
    return process

from FWCore.ParameterSet.MassReplace import massSearchReplaceAnyInputTag
def nanoAOD_customizeMC(process):
    for name, path in process.paths.iteritems():
        # replace all the non-match embedded inputs with the matched ones
        massSearchReplaceAnyInputTag(path, 'muonTrgSelector:SelectedMuons', 'selectedMuonsMCMatchEmbedded')
        massSearchReplaceAnyInputTag(path, 'tracksBPark:SelectedTracks', 'tracksBParkMCMatchEmbedded')

        # modify the path to include mc-specific info
        path.insert(0, nanoSequenceMC)
        path.replace(process.muonBParkSequence, process.muonBParkMC)
        path.replace(process.tracksBParkSequence, process.tracksBParkMC)

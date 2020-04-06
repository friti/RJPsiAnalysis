import FWCore.ParameterSet.Config as cms
from PhysicsTools.BParkingNano.common_cff import *

electronPairsForKee = cms.EDProducer(
    'DiElectronBuilder',
    src = cms.InputTag('electronsForAnalysis', 'SelectedElectrons'),
    transientTracksSrc = cms.InputTag('electronsForAnalysis', 'SelectedTransientElectrons'),
    lep1Selection = cms.string('pt > 1.5 && userFloat("unBiased") >= 3'),
    lep2Selection = cms.string(''),
    preVtxSelection = cms.string(
        'abs(userCand("l1").vz - userCand("l2").vz) <= 1. && mass() < 5 '
        '&& mass() > 0 && charge() == 0 && userFloat("lep_deltaR") > 0.03'
    ),
    postVtxSelection = cms.string('userFloat("sv_chi2") < 998 && userFloat("sv_prob") > 1.e-5'),
)

BToKee = cms.EDProducer(
    'BToKLLBuilder',
    dileptons = cms.InputTag('electronPairsForKee'),
    leptonTransientTracks = electronPairsForKee.transientTracksSrc,
    kaons = cms.InputTag('tracksBPark', 'SelectedTracks'),
    kaonsTransientTracks = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    tracks = cms.InputTag("packedPFCandidates"),
    lostTracks = cms.InputTag("lostTracks"),
    kaonSelection = cms.string(''),
    isoTracksSelection = cms.string('pt > 0.7 && abs(eta)<2.5'),
    preVtxSelection = cms.string(
        'pt > 3. && userFloat("min_dr") > 0.03 '
        '&& mass < 7. && mass > 4.'
        ),
    postVtxSelection = cms.string(
        'userInt("sv_OK") == 1 && userFloat("sv_prob") > 0.001 '
        '&& userFloat("fitted_cos_theta_2D") >= 0 '
        '&& userFloat("fitted_mass") > 4.5 && userFloat("fitted_mass") < 6.'
    )
)

muonPairsFormmm = cms.EDProducer(
    'DiMuonBuilder',
    src = cms.InputTag('muonTrgSelector', 'SelectedMuons'),
    transientTracksSrc = cms.InputTag('muonTrgSelector', 'SelectedTransientMuons'),
    lep1Selection = cms.string('pt > 1.5'),
    lep2Selection = cms.string(''),
    preVtxSelection = cms.string('abs(userCand("l1").vz - userCand("l2").vz) <= 1. && mass() < 5 '
                                 '&& mass() > 0 && charge() == 0 && userFloat("lep_deltaR") > 0.03'),
    postVtxSelection = electronPairsForKee.postVtxSelection,
)

BTommm = cms.EDProducer(
    'BTommmBuilder',
    dileptons = cms.InputTag('muonPairsFormmm'),
    leptonTransientTracks = muonPairsFormmm.transientTracksSrc,
    kaons = cms.InputTag('muonTrgSelector', 'SelectedMuons'),
    kaonsTransientTracks = cms.InputTag('muonTrgSelector', 'SelectedTransientMuons'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    tracks = cms.InputTag("packedPFCandidates"),
    lostTracks = cms.InputTag("lostTracks"),
    kaonSelection = cms.string(''),
    isoTracksSelection = BToKee.isoTracksSelection,
    # This in principle can be different between electrons and muons
    preVtxSelection = cms.string(
        'pt > 3. && userFloat("min_dr") > 0.03'
        #'&& mass < 7. && mass > 4.'
        ),
    postVtxSelection = cms.string(
        'userInt("sv_OK") == 1 && userFloat("sv_prob") > 0.001 '
        '&& userFloat("fitted_cos_theta_2D") >= 0'
        #'&& userFloat("fitted_mass") > 4.5 && userFloat("fitted_mass") < 6.'
    ),
    #TRIGGER                                        
     
    bits=cms.InputTag("TriggerResults","","HLT"),               
    objects=cms.InputTag("slimmedPatTrigger"), 

)

BToKeeTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("BToKee"),
    cut = cms.string(""),
    name = cms.string("BToKEE"),
    doc = cms.string("BToKEE Variable"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        # pre-fit quantities
        CandVars,
        l1Idx = uint('l1_idx'),
        l2Idx = uint('l2_idx'),
        kIdx = uint('k_idx'),
        minDR = ufloat('min_dr'),
        maxDR = ufloat('max_dr'),
        # fit and vtx info
        #chi2 = ufloat('sv_chi2'),
        svprob = ufloat('sv_prob'),
        l_xy = ufloat('l_xy'),
        l_xy_unc = ufloat('l_xy_unc'),
        vtx_x = ufloat('vtx_x'),
        vtx_y = ufloat('vtx_y'),
        vtx_z = ufloat('vtx_z'),
        vtx_ex = ufloat('vtx_ex'), ## only saving diagonal elements of the cov matrix
        vtx_ey = ufloat('vtx_ey'),
        vtx_ez = ufloat('vtx_ez'),
        # Mll
        mll_raw = Var('userCand("dilepton").mass()', float),
        mll_llfit = Var('userCand("dilepton").userFloat("fitted_mass")', float), # this might not work
        mllErr_llfit = Var('userCand("dilepton").userFloat("fitted_massErr")', float), # this might not work
        mll_fullfit = ufloat('fitted_mll'),
        # Cos(theta)
        cos2D = ufloat('cos_theta_2D'),
        fit_cos2D = ufloat('fitted_cos_theta_2D'),
        # post-fit momentum
        fit_mass = ufloat('fitted_mass'),
        fit_massErr = ufloat('fitted_massErr'),
        fit_pt = ufloat('fitted_pt'),
        fit_eta = ufloat('fitted_eta'),
        fit_phi = ufloat('fitted_phi'),
        fit_l1_pt = ufloat('fitted_l1_pt'),
        fit_l1_eta = ufloat('fitted_l1_eta'),
        fit_l1_phi = ufloat('fitted_l1_phi'),
        fit_l2_pt = ufloat('fitted_l2_pt'),
        fit_l2_eta = ufloat('fitted_l2_eta'),
        fit_l2_phi = ufloat('fitted_l2_phi'),
        fit_k_pt = ufloat('fitted_k_pt'),
        fit_k_eta = ufloat('fitted_k_eta'),
        fit_k_phi = ufloat('fitted_k_phi'),
        l1_iso03 = ufloat('l1_iso03'),
        l1_iso04 = ufloat('l1_iso04'),
        l2_iso03 = ufloat('l2_iso03'),
        l2_iso04 = ufloat('l2_iso04'),
        k_iso03  = ufloat('k_iso03'),
        k_iso04  = ufloat('k_iso04'),
        b_iso03  = ufloat('b_iso03'),
        b_iso04  = ufloat('b_iso04'),
        n_k_used = uint('n_k_used'),
        n_l1_used = uint('n_l1_used'),
        n_l2_used = uint('n_l2_used'),
    )
)

BTommmTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("BTommm"),
    cut = cms.string(""),
    name = cms.string("BTommm"),
    doc = cms.string("BTommm Variable"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        # pre-fit quantities                                                      
        CandVars,
        #nome branch= nome variabile del .cc
        l1Idx = uint('l1_idx'),
        l2Idx = uint('l2_idx'),
        kIdx = uint('k_idx'),
        minDR = ufloat('min_dr'),
        maxDR = ufloat('max_dr'),
        # fit and vtx info                                                                                                    
        #chi2 = ufloat('sv_chi2'),
                         
        svprob = ufloat('sv_prob'),
        l_xy = ufloat('l_xy'),
        l_xy_unc = ufloat('l_xy_unc'),
        vtx_x = ufloat('vtx_x'),
        vtx_y = ufloat('vtx_y'),
        vtx_z = ufloat('vtx_z'),
        vtx_ex = ufloat('vtx_ex'), ## only saving diagonal elements of the cov matrix                                         
        vtx_ey = ufloat('vtx_ey'),
        vtx_ez = ufloat('vtx_ez'),
        # Mll                                                                                                                 
        mll_raw = Var('userCand("dilepton").mass()', float),
        mll_llfit = Var('userCand("dilepton").userFloat("fitted_mass")', float), # this might not work                        
        mllErr_llfit = Var('userCand("dilepton").userFloat("fitted_massErr")', float), # this might not work                  
        mll_fullfit = ufloat('fitted_mll'),
        # Cos(theta)                                                                                                          
        cos2D = ufloat('cos_theta_2D'),
        fit_cos2D = ufloat('fitted_cos_theta_2D'),
        # post-fit momentum                                                                                                   
        fit_mass = ufloat('fitted_mass'),
        fit_massErr = ufloat('fitted_massErr'),
        fit_pt = ufloat('fitted_pt'),
        fit_eta = ufloat('fitted_eta'),
        fit_phi = ufloat('fitted_phi'),
        fit_l1_pt = ufloat('fitted_l1_pt'),
        fit_l1_eta = ufloat('fitted_l1_eta'),
        fit_l1_phi = ufloat('fitted_l1_phi'),
        fit_l2_pt = ufloat('fitted_l2_pt'),
        fit_l2_eta = ufloat('fitted_l2_eta'),
        fit_l2_phi = ufloat('fitted_l2_phi'),
        fit_k_pt = ufloat('fitted_k_pt'),
        fit_k_eta = ufloat('fitted_k_eta'),
        fit_k_phi = ufloat('fitted_k_phi'),
        l1_iso03 = ufloat('l1_iso03'),
        l1_iso04 = ufloat('l1_iso04'),
        l2_iso03 = ufloat('l2_iso03'),
        l2_iso04 = ufloat('l2_iso04'),
        k_iso03  = ufloat('k_iso03'),
        k_iso04  = ufloat('k_iso04'),
        b_iso03  = ufloat('b_iso03'),
        b_iso04  = ufloat('b_iso04'),
        n_k_used = uint('n_k_used'),
        n_l1_used = uint('n_l1_used'),
        n_l2_used = uint('n_l2_used'),
        #my variables
        pass_3mu=uint('pass_3mu'),
        m_miss_sq=ufloat('m_miss_2'),
        Q_sq=ufloat('Q_2'),
        pt_miss=ufloat('pt_miss'),
        pt_miss_vec=ufloat('pt_miss_vec'),
        pt_var=ufloat('pt_var'),
        DR=ufloat('DR'),
        E_mu_star=ufloat('E_mu_star'),
        E_mu_canc=ufloat('E_mu_#'),
        m_jpsi=ufloat('m_jpsi'),
        #jPsi_mass_online=ufloat('jPsi_mass_online')                                                         
    )
)

CountBToKee = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("BToKee")
)    
CountBTommm = CountBToKee.clone(
    minNumber = cms.uint32(0),
    src = cms.InputTag("BTommm")
)


BTommmSequence = cms.Sequence(
    (muonPairsFormmm * BTommm)
)
BToKEESequence = cms.Sequence(
    (electronPairsForKee * BToKee)
)

#BToKLLSequence = cms.Sequence(
#    (electronPairsForKee * BToKee) +
#    (muonPairsFormmm * BTommm)
#)
BToKLLTables = cms.Sequence(BTommmTable)


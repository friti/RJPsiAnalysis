import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var, CandVars
from PhysicsTools.RJPsiNano.common_cff import RJpsiCandVars, ufloat, uint, ubool

JpsiMuonPairs = cms.EDProducer(
    'DiMuonBuilder',
    src                    = cms.InputTag('muonTrgSelector', 'SelectedMuons'),
    transientTracksSrc     = cms.InputTag('muonTrgSelector', 'SelectedTransientMuons'),
    vertexCollection       = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muon1Selection         = cms.string('pt > 2.5'),
    muon2Selection         = cms.string(''),
    preVtxSelection        = cms.string(' && '.join([
        'abs(userCand("mu1").bestTrack.dz - userCand("mu2").bestTrack.dz) <= 0.4 ',
        'mass() > 2',
        'mass() < 4',
        'userFloat("muons12_deltaR") > 0.01',
        #'mass() > 0.0',
        #'mass() < 5.0',
        #'userFloat("muons12_deltaR") > 0.03',
        ])
    ),
    postVtxSelection   = cms.string(
        'userFloat("sv_prob") > 1.e-5 '
        '&& pt > 3 '
        #        '&& userFloat("sv_chi2") < 998 ' 
    ),
    beamSpot              = cms.InputTag("offlineBeamSpot"),
)

BuilderDefaultCfg = cms.PSet(
    dimuons             = cms.InputTag('JpsiMuonPairs','muonPairsForB'),
    #pvSelected = cms.InputTag('pvSelector', 'bestVertex'),
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    #muons            = cms.InputTag('muonTrgSelector', 'SelectedMuons'),
    muonsTransientTracks = JpsiMuonPairs.transientTracksSrc,
    #kaons                 = cms.InputTag('tracksBPark', 'SelectedTracks'),
    #kaonsTransientTracks  = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    beamSpot              = cms.InputTag("offlineBeamSpot"),
    tracks                = cms.InputTag("packedPFCandidates"),
    lostTracks            = cms.InputTag("lostTracks"),
    #kaonSelection         = cms.string(''),
    isoTracksSelection    = cms.string('pt > 0.5 && abs(eta)<2.5'),
    preVtxSelection       = cms.string(''),
    postVtxSelection      = cms.string(' && '.join([
        'userFloat("fitted_cos_theta_2D") >= 0',
        'mass < 8.',
        'userInt("sv_OK") == 1',
        'userFloat("sv_prob") > 1e-8',
        ])
    ),
    bits                  = cms.InputTag("TriggerResults","","HLT"),               
    objects               = cms.InputTag("slimmedPatTrigger"), 
)

TableDefault = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src       = cms.InputTag("your_cand_name"),
    cut       = cms.string(""),
    name      = cms.string("your_cand_name"),
    doc       = cms.string("your_cand_name Variable"),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    variables = cms.PSet(),
)

TableDefaultVariables = cms.PSet(
    RJpsiCandVars,
    mu1Idx = uint('mu1_idx'),
    mu2Idx = uint('mu2_idx'),
    
    minDR = ufloat('min_dr'),
    maxDR = ufloat('max_dr'),
    

    
    #all final particles vertex
    bodies3_chi2     = ufloat('sv_chi2'),
    bodies3_chi2_2     = ufloat('vtx_chi2'),
    bodies3_svprob   = ufloat('sv_prob'),
    bodies3_l_xy     = ufloat('l_xy'),
    bodies3_l_xy_unc = ufloat('l_xy_unc'),
    bodies3_vtx_x    = ufloat('vtx_x'),
    bodies3_vtx_y    = ufloat('vtx_y'),
    bodies3_vtx_z    = ufloat('vtx_z'),
    bodies3_vtx_ex   = ufloat('vtx_ex'), 
    bodies3_vtx_ey   = ufloat('vtx_ey'),
    bodies3_vtx_ez   = ufloat('vtx_ez'),
    bodies3_cos2D     = ufloat('cos_theta_2D'),
    bodies3_svOK     = uint('sv_OK'),
    
    #post fit all particles vertex
    bodies3_fit_mass    = ufloat('fitted_mass'),
    bodies3_fit_massErr = ufloat('fitted_massErr'),
    bodies3_fit_pt      = ufloat('fitted_pt'),
    bodies3_fit_eta     = ufloat('fitted_eta'),
    bodies3_fit_phi     = ufloat('fitted_phi'),
    bodies3_fit_mu1_pt   = ufloat('fitted_mu1_pt'),
    bodies3_fit_mu1_eta  = ufloat('fitted_mu1_eta'),
    bodies3_fit_mu1_phi  = ufloat('fitted_mu1_phi'),
    bodies3_fit_mu2_pt   = ufloat('fitted_mu2_pt'),
    bodies3_fit_mu2_eta  = ufloat('fitted_mu2_eta'),
    bodies3_fit_mu2_phi  = ufloat('fitted_mu2_phi'),
    bodies3_fit_cos2D = ufloat('fitted_cos_theta_2D'),

    #2 particles vertex
    jpsivtx_chi2 = Var('userCand("dilepton").userFloat("sv_chi2")', float),
    jpsivtx_svprob = Var('userCand("dilepton").userFloat("sv_prob")', float),
    jpsivtx_l_xy = Var('userCand("dilepton").userFloat("l_xy")', float),
    jpsivtx_l_xy_unc = Var('userCand("dilepton").userFloat("l_xy_unc")', float),
    jpsivtx_vtx_x = Var('userCand("dilepton").userFloat("vtx_x")', float),
    jpsivtx_vtx_y = Var('userCand("dilepton").userFloat("vtx_y")', float),
    jpsivtx_vtx_z = Var('userCand("dilepton").userFloat("vtx_z")', float),
    jpsivtx_vtx_ex = Var('userCand("dilepton").userFloat("vtx_ex")', float),
    jpsivtx_vtx_ey = Var('userCand("dilepton").userFloat("vtx_ey")', float),
    jpsivtx_vtx_ez = Var('userCand("dilepton").userFloat("vtx_ez")', float),
    jpsivtx_cos2D = Var('userCand("dilepton").userFloat("cos_theta_2D")', float),

    #post fit 2 particles vertex
    jpsivtx_fit_mass    = Var('userCand("dilepton").userFloat("fitted_mass")', float),
    jpsivtx_fit_massErr = Var('userCand("dilepton").userFloat("fitted_massErr")', float),
    jpsivtx_fit_pt    = Var('userCand("dilepton").userFloat("fitted_pt")', float),
    jpsivtx_fit_eta    = Var('userCand("dilepton").userFloat("fitted_eta")', float),
    jpsivtx_fit_phi   = Var('userCand("dilepton").userFloat("fitted_phi")', float),
    jpsivtx_fit_mu1_pt   = Var('userCand("dilepton").userFloat("fitted_mu1_pt")', float),
    jpsivtx_fit_mu1_eta   = Var('userCand("dilepton").userFloat("fitted_mu1_eta")', float),
    jpsivtx_fit_mu1_phi   = Var('userCand("dilepton").userFloat("fitted_mu1_phi")', float),
    jpsivtx_fit_mu2_pt   = Var('userCand("dilepton").userFloat("fitted_mu2_pt")', float),
    jpsivtx_fit_mu2_eta   = Var('userCand("dilepton").userFloat("fitted_mu2_eta")', float),
    jpsivtx_fit_mu2_phi   = Var('userCand("dilepton").userFloat("fitted_mu2_phi")', float),
    jpsivtx_fit_cos2D = Var('userCand("dilepton").userFloat("fitted_cos_theta_2D")', float),

    #isolation
    mu1_iso03 = ufloat('mu1_iso03'),
    mu1_iso04 = ufloat('mu1_iso04'),
    mu2_iso03 = ufloat('mu2_iso03'),
    mu2_iso04 = ufloat('mu2_iso04'),
    b_iso03  = ufloat('b_iso03'),
    b_iso04  = ufloat('b_iso04'),

    #beamspot
    beamspot_x     = ufloat('beamspot_x'),
    beamspot_y     = ufloat('beamspot_y'),
    beamspot_z     = ufloat('beamspot_z'),    


    #new variables
    m_miss_sq   = ufloat('m_miss_2'),
    Q_sq        = ufloat('Q_2'),
    pt_miss     = ufloat('pt_miss'),
    pt_miss_vec = ufloat('pt_miss_vec'),
    pt_var      = ufloat('pt_var'),
    DR          = ufloat('DR'),
    m_jpsi      = ufloat('m_jpsi'),

    #PV vertex
    pv_idx = uint('pv_idx'),
    nPrimaryVertices = uint('nPrimaryVertices'),
    pv_x = ufloat('pv_x'),
    pv_y = ufloat('pv_y'),
    pv_z = ufloat('pv_z'),
    pv_ex = ufloat('pv_ex'),
    pv_ey = ufloat('pv_ey'),
    pv_ez = ufloat('pv_ez'),
    pv_exy = ufloat('pv_exz'),
    pv_eyz = ufloat('pv_eyz'),
    pv_exz = ufloat('pv_exz'),
    pv_chi2 = ufloat('pv_chi2'),

    # Mll (do we need this???)                                                                                               
    mll_raw = Var('userCand("dimuon").mass()', float),
    mll_llfit = Var('userCand("dimuon").userFloat("fitted_mass")', float), # this might not work       
    mllErr_llfit = Var('userCand("dimuon").userFloat("fitted_massErr")', float), # this might not work 
    mll_fullfit = ufloat('fitted_mll'),
    mll_vtxex= Var('userCand("dimuon").userFloat("vtx_ex")',float),

    #number of muons used
    n_mu1_used = uint('n_mu1_used'),
    n_mu2_used = uint('n_mu2_used'),

    mu1_dxy = ufloat('mu1_dxy'),
    mu1_dz = ufloat('mu1_dz'),
    mu2_dxy = ufloat('mu2_dxy'),
    mu2_dz = ufloat('mu2_dz'),

    mu1_dxyErr = ufloat('mu1_dxyErr'),
    mu1_dzErr = ufloat('mu1_dzErr'),
    mu2_dxyErr = ufloat('mu2_dxyErr'),
    mu2_dzErr = ufloat('mu2_dzErr'),
)

#builder for final states with 3 particles
Final3PartTableVariables = TableDefaultVariables.clone(
    kIdx     = uint('k_idx'),
    bodies3_fit_k_pt    = ufloat('fitted_k_pt'),
    bodies3_fit_k_eta   = ufloat('fitted_k_eta'),
    bodies3_fit_k_phi   = ufloat('fitted_k_phi'),
    k_iso03     = ufloat('k_iso03'),
    k_iso04     = ufloat('k_iso04'),
    E_mu_star   = ufloat('E_mu_star'),
    E_mu_canc   = ufloat('E_mu_#'),
    n_k_used = uint('n_k_used'),
    ip3D = ufloat('ip3D'),
    ip3D_e = ufloat('ip3D'),

    #dz and dxy  for muon particle w.r.t. best pv.

    k_dxy = ufloat('k_dxy'),
    k_dz = ufloat('k_dz'),
    k_dxyErr = ufloat('k_dxyErr'),
    k_dzErr = ufloat('k_dzErr'),
)

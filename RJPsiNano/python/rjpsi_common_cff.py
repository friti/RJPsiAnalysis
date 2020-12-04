import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var, CandVars
from PhysicsTools.RJPsiNano.common_cff import RJpsiCandVars, ufloat, uint, ubool

JpsiMuonPairs = cms.EDProducer(
    'DiMuonBuilder',
    src                = cms.InputTag('muonTrgSelector', 'SelectedMuons'),
    transientTracksSrc = cms.InputTag('muonTrgSelector', 'SelectedTransientMuons'),
    muon1Selection      = cms.string('pt > 2.5'),
    muon2Selection      = cms.string(''),
    preVtxSelection    = cms.string(' && '.join([
#         'abs(userCand("l1").dz - userCand("l2").dz) <= 0.4 ',
        'abs(userCand("l1").bestTrack.dz - userCand("l2").bestTrack.dz) <= 0.4 ',
        'mass() > 2',
        'mass() < 4',
        'userFloat("muons12_deltaR") > 0.01',
        ])
    ),
    postVtxSelection   = cms.string(
        'pt > 3 '
#         '&& userFloat("sv_chi2") < 998 ' 
        '&& userFloat("sv_prob") > 1.e-5 '
    ),
)

BuilderDefaultCfg = cms.PSet(
    dimuons             = cms.InputTag('JpsiMuonPairs','muonPairsForBTo3Mu'),
    pvSelected = cms.InputTag('pvSelector', 'bestVertex'),
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
        'userInt("sv_OK") == 1',
        'userFloat("sv_prob") > 1e-8',
        'userFloat("fitted_cos_theta_2D") >= 0',
#         'userFloat("fitted_mass") > 4.5',
#         'userFloat("fitted_mass") < 8.',
        'mass > 4.5',
        'mass < 8.',
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
  # pre-fit quantities                                                      
  RJpsiCandVars,
  #nome branch= nome variabile del .cc
  mu1Idx = uint('mu1_idx'),
  mu2Idx = uint('mu2_idx'),
  kIdx = uint('k_idx'),
  minDR = ufloat('min_dr'),
  maxDR = ufloat('max_dr'),
  # fit and vtx info                                                                                                    
  #chi2 = ufloat('sv_chi2'),
  ip3D = ufloat('ip3D'),
  ip3D_e = ufloat('ip3D'),
                   
  svprob = ufloat('sv_prob'),
  l_xy = ufloat('l_xy'),
  l_xy_unc = ufloat('l_xy_unc'),
  vtx_x = ufloat('vtx_x'),
  vtx_y = ufloat('vtx_y'),
  vtx_z = ufloat('vtx_z'),
  vtx_ex = ufloat('vtx_ex'), ## only saving diagonal elements of the cov matrix                                         
  vtx_ey = ufloat('vtx_ey'),
  vtx_ez = ufloat('vtx_ez'),
  vtx_chi2 = ufloat('vtx_chi2'),

  jpsi_vtx_x = ufloat('jpsi_vtx_x'),
  jpsi_vtx_y = ufloat('jpsi_vtx_y'),
  jpsi_vtx_z = ufloat('jpsi_vtx_z'),
  jpsi_vtx_ex = ufloat('jpsi_vtx_ex'),
  jpsi_vtx_ey = ufloat('jpsi_vtx_ey'),
  jpsi_vtx_ez = ufloat('jpsi_vtx_ez'),
  jpsi_vtx_chi2 = ufloat('jpsi_vtx_chi2'),

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
  # Mll                                                                                                                 
  mll_raw = Var('userCand("dimuon").mass()', float),
  mll_llfit = Var('userCand("dimuon").userFloat("fitted_mass")', float), # this might not work                        
  mllErr_llfit = Var('userCand("dimuon").userFloat("fitted_massErr")', float), # this might not work                  
  mll_fullfit = ufloat('fitted_mll'),
  mll_vtxex= Var('userCand("dimuon").userFloat("vtx_ex")',float),
#  mll_vtxx= Var('userCand("dimuon").userFloat("vtx_x")',float),
#  mll_vtxy= Var('userCand("dimuon").userFloat("vtx_y")',float),
#  mll_vtxz= Var('userCand("dimuon").userFloat("vtx_z")',float),
#  mll_vtxex= Var('userCand("dimuon").userFloat("vtx_ex")',float),
#  mll_vtxey= Var('userCand("dimuon").userFloat("vtx_ey")',float),
#  mll_vtxez= Var('userCand("dimuon").userFloat("vtx_ez")',float),

  # Cos(theta)                                                                                                          
  cos2D = ufloat('cos_theta_2D'),
  fit_cos2D = ufloat('fitted_cos_theta_2D'),
  # post-fit momentum                                                                                                   
  fit_mass = ufloat('fitted_mass'),
  fit_massErr = ufloat('fitted_massErr'),
  fit_pt = ufloat('fitted_pt'),
  fit_eta = ufloat('fitted_eta'),
  fit_phi = ufloat('fitted_phi'),
  fit_mu1_pt = ufloat('fitted_mu1_pt'),
  fit_mu1_eta = ufloat('fitted_mu1_eta'),
  fit_mu1_phi = ufloat('fitted_mu1_phi'),
  fit_mu2_pt = ufloat('fitted_mu2_pt'),
  fit_mu2_eta = ufloat('fitted_mu2_eta'),
  fit_mu2_phi = ufloat('fitted_mu2_phi'),
  fit_k_pt = ufloat('fitted_k_pt'),
  fit_k_eta = ufloat('fitted_k_eta'),
  fit_k_phi = ufloat('fitted_k_phi'),
  mu1_iso03 = ufloat('mu1_iso03'),
  mu1_iso04 = ufloat('mu1_iso04'),
  mu2_iso03 = ufloat('mu2_iso03'),
  mu2_iso04 = ufloat('mu2_iso04'),
  k_iso03  = ufloat('k_iso03'),
  k_iso04  = ufloat('k_iso04'),
  b_iso03  = ufloat('b_iso03'),
  b_iso04  = ufloat('b_iso04'),
  n_k_used = uint('n_k_used'),
  n_mu1_used = uint('n_mu1_used'),
  n_mu2_used = uint('n_mu2_used'),
  #my variables
  #pass_3mu=uint('pass_3mu'),
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

  #Gen Variables
  is_jpsi_mu=uint("is_jpsi_mu"),
  is_psi2s_mu=uint("is_psi2s_mu"),
  is_chic0_mu=uint("is_chic0_mu"),
  is_chic1_mu=uint("is_chic1_mu"),
  is_chic2_mu=uint("is_chic2_mu"),
  is_hc_mu=uint("is_hc_mu"),
  is_jpsi_tau=uint("is_jpsi_tau"),
  is_psi2s_tau=uint("is_psi2s_tau"),
  is_jpsi_pi=uint("is_jpsi_pi"),
  is_jpsi_3pi=uint("is_jpsi_3pi"),
  is_jpsi_hc=uint("is_jpsi_hc"),
  is_error=uint("is_error"),
  weightGen= ufloat("weightGen")
)
Final3PartTableVariables = TableDefaultVariables.clone(
    kIdx     = uint('k_idx'),
    bodies3_fit_k_pt    = ufloat('fitted_k_pt'),
    bodies3_fit_k_eta   = ufloat('fitted_k_eta'),
    bodies3_fit_k_phi   = ufloat('fitted_k_phi'),
    k_iso03     = ufloat('k_iso03'),
    k_iso04     = ufloat('k_iso04'),
    E_mu_star   = ufloat('E_mu_star'),
    E_mu_canc   = ufloat('E_mu_#'),
)

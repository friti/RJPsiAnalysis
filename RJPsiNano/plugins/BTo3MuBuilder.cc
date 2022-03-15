#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "TLorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "RecoTauTag/ImpactParameter/interface/ImpactParameterAlgorithm.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include <RecoBTag/BTagTools/interface/SignedImpactParameter3D.h>
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

#include <vector>
#include <memory>
#include <map>
#include <string>
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "helper.h"
#include <limits>
#include <algorithm>
#include "KinVtxFitter.h"

constexpr bool debugGen = false;
constexpr bool debug = false;


class BTo3MuBuilder : public edm::global::EDProducer<> {

  // perhaps we need better structure here (begin run etc)
public:
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;
  typedef std::vector<reco::GenParticle> GenParticleCollection;
  //typedef std::vector<KinVtxFitter> KinVtxFitterCollection;

  explicit BTo3MuBuilder(const edm::ParameterSet &cfg):
    k_selection_{cfg.getParameter<std::string>("muonSelection")},
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    dimuons_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("dimuons") )},
    primaryVertices_{consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("primaryVertices"))},
    muons_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("muonsTransientTracks") )},
    muons_{consumes<pat::MuonCollection>( cfg.getParameter<edm::InputTag>("muons") )},
    isotracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
    isolostTracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))),

    isotrk_selection_{cfg.getParameter<std::string>("isoTracksSelection")},
    beamspot_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") )} 
  {
      produces<pat::CompositeCandidateCollection>();
  }

  ~BTo3MuBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
  int  getPVIdx(const reco::VertexCollection*,const reco::TransientTrack&) const;
  Measurement1D getIP(edm::Ptr<pat::CompositeCandidate> ll_ptr, reco::Vertex pv, reco::TransientTrack transientTrackMu) const;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
  private:
  const StringCutObjectSelector<pat::Muon> k_selection_; 
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; // cut on the di-muon before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; // cut on the di-muon after the SV fit

  const edm::EDGetTokenT<pat::CompositeCandidateCollection> dimuons_;
  const edm::EDGetTokenT<reco::VertexCollection> primaryVertices_;
  const edm::EDGetTokenT<TransientTrackCollection> muons_ttracks_;
  const edm::EDGetTokenT<pat::MuonCollection> muons_;
  //const edm::EDGetTokenT<TransientTrackCollection> muons_ttracks_;

  const edm::EDGetTokenT<pat::PackedCandidateCollection> isotracksToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> isolostTracksToken_;

  const StringCutObjectSelector<pat::PackedCandidate> isotrk_selection_; 

  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;  
};

void BTo3MuBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  //input
  edm::Handle<pat::CompositeCandidateCollection> dimuons;
  evt.getByToken(dimuons_, dimuons);

  edm::Handle<reco::VertexCollection> primaryVertices;
  evt.getByToken(primaryVertices_, primaryVertices);
  
  edm::Handle<TransientTrackCollection> muons_ttracks;
  evt.getByToken(muons_ttracks_, muons_ttracks);

  edm::Handle<pat::MuonCollection> muons;
  evt.getByToken(muons_, muons);
  
  edm::Handle<reco::BeamSpot> beamspot;
  evt.getByToken(beamspot_, beamspot);  

  //for isolation
  edm::Handle<pat::PackedCandidateCollection> iso_tracks;
  evt.getByToken(isotracksToken_, iso_tracks);
  edm::Handle<pat::PackedCandidateCollection> iso_lostTracks;
  evt.getByToken(isolostTracksToken_, iso_lostTracks);

//////

  unsigned int nTracks     = iso_tracks->size();
  unsigned int totalTracks = nTracks + iso_lostTracks->size();

  std::vector<int> used_muon1_id, used_muon2_id, used_trk_id;

  std::unique_ptr<pat::CompositeCandidateCollection> ret_val(new pat::CompositeCandidateCollection());

  const reco::VertexCollection* vertices = primaryVertices.product();
  int nPrimaryVertices = vertices->size();
  // output
  if (debug) std::cout <<"number of dimuons"<<dimuons->size()<< std::endl;
  for(size_t ll_idx = 0; ll_idx < dimuons->size(); ++ll_idx) 
  {
    //std::cout << "PV " << ll_idx << ": " << vertices->at(ll_idx).position() << std::endl;
    edm::Ptr<pat::CompositeCandidate> ll_ptr(dimuons, ll_idx);
    edm::Ptr<reco::Candidate> mu1_ptr = ll_ptr->userCand("mu1");
    edm::Ptr<reco::Candidate> mu2_ptr = ll_ptr->userCand("mu2");
    size_t mu1_idx = abs(ll_ptr->userInt("mu1_idx"));
    size_t mu2_idx = abs(ll_ptr->userInt("mu2_idx"));

    // PV shortest dz from jpsi
    int pvIdx = ll_ptr->userInt("pvIdx");
    reco::Vertex pv_jpsi = vertices->at(pvIdx);


    double mu1_pvjpsi_dxy = mu1_ptr->bestTrack()->dxy(pv_jpsi.position());
    double mu2_pvjpsi_dxy = mu2_ptr->bestTrack()->dxy(pv_jpsi.position());
    double mu1_pvjpsi_dz = mu1_ptr->bestTrack()->dz(pv_jpsi.position());
    double mu2_pvjpsi_dz = mu2_ptr->bestTrack()->dz(pv_jpsi.position());
    double mu1_pvjpsi_dxyErr = mu1_ptr->bestTrack()->dxyError(pv_jpsi.position(),pv_jpsi.covariance());
    double mu2_pvjpsi_dxyErr = mu2_ptr->bestTrack()->dxyError(pv_jpsi.position(),pv_jpsi.covariance());
    double mu1_pvjpsi_dzErr = mu1_ptr->bestTrack()->dzError();
    double mu2_pvjpsi_dzErr = mu2_ptr->bestTrack()->dzError();

    // first pv of the collection
    reco::Vertex pv_first = vertices->at(0);
    double mu1_pvfirst_dxy = mu1_ptr->bestTrack()->dxy(pv_first.position());
    double mu2_pvfirst_dxy = mu2_ptr->bestTrack()->dxy(pv_first.position());
    double mu1_pvfirst_dz = mu1_ptr->bestTrack()->dz(pv_first.position());
    double mu2_pvfirst_dz = mu2_ptr->bestTrack()->dz(pv_first.position());
    double mu1_pvfirst_dxyErr = mu1_ptr->bestTrack()->dxyError(pv_first.position(),pv_first.covariance());
    double mu2_pvfirst_dxyErr = mu2_ptr->bestTrack()->dxyError(pv_first.position(),pv_first.covariance());
    double mu1_pvfirst_dzErr = mu1_ptr->bestTrack()->dzError();
    double mu2_pvfirst_dzErr = mu2_ptr->bestTrack()->dzError();

    //if(debug) std::cout << "mu1_dxy" << mu1_dxy << std::endl;
    //if(debug) std::cout << "mu1_dz" << mu1_dz << std::endl;
    //if(debug) std::cout << "mu2_dxy" << mu2_dxy << std::endl;
    //if(debug) std::cout << "mu2_dz" << mu2_dz << std::endl;
    
    size_t isDimuon_dimuon0Trg = abs(ll_ptr->userInt("muonpair_fromdimuon0"));
    size_t isDimuon_jpsiTrkTrg = abs(ll_ptr->userInt("muonpair_fromjpsitrk"));
    size_t isDimuon_jpsiTrk_PsiPrimeTrg = abs(ll_ptr->userInt("muonpair_fromjpsitrk_PsiPrime"));
    size_t isDimuon_jpsiTrk_NonResonantTrg = abs(ll_ptr->userInt("muonpair_fromjpsitrk_NonResonant"));
    //size_t isDimuon_jpsiTrkTrg = abs(ll_ptr->userInt("isJpsiTrkTrg"));
    //size_t isDimuon_dimuon0Trg = abs(ll_ptr->userInt("isDimuon0Trg"));
    if(!isDimuon_dimuon0Trg && !isDimuon_jpsiTrk_PsiPrimeTrg && !isDimuon_jpsiTrkTrg && !isDimuon_jpsiTrk_NonResonantTrg) {
      if(debug) std::cout<<"Not dimuon0 trigger couple"<<std::endl;
      continue;
    }
    
    //Loop  on displaced muons    
    if(debug) std::cout <<"Number of muons"<<std::endl;
    for(size_t k_idx = 0; k_idx < muons->size(); ++k_idx) {
      edm::Ptr<pat::Muon> k_ptr(muons, k_idx);
      if((mu1_idx == k_idx) || (mu2_idx == k_idx)) continue;
      if( !k_selection_(*k_ptr) ) continue;
      double k_pvjpsi_dxy = k_ptr->bestTrack()->dxy(pv_jpsi.position());
      double k_pvjpsi_dz = k_ptr->bestTrack()->dz(pv_jpsi.position());
      double k_pvjpsi_dxyErr = k_ptr->bestTrack()->dxyError(pv_jpsi.position(),pv_jpsi.covariance());
      double k_pvjpsi_dzErr = k_ptr->bestTrack()->dzError();

      double k_pvfirst_dxy = k_ptr->bestTrack()->dxy(pv_first.position());
      double k_pvfirst_dz = k_ptr->bestTrack()->dz(pv_first.position());
      double k_pvfirst_dxyErr = k_ptr->bestTrack()->dxyError(pv_first.position(),pv_first.covariance());
      double k_pvfirst_dzErr = k_ptr->bestTrack()->dzError();
  
      //std::cout << "here1" << std::endl;
      //ha trovato il mu displaced
      bool isDimuon0Trg = k_ptr->userInt("isDimuon0Trg");
      bool isJpsiTrkTrg = k_ptr->userInt("isJpsiTrkTrg");
      bool isJpsiTrk_PsiPrimeTrg = k_ptr->userInt("isJpsiTrk_PsiPrimeTrg");
      bool isJpsiTrk_NonResonantTrg = k_ptr->userInt("isJpsiTrk_NonResonantTrg");
      
      bool isMuonFromJpsi_dimuon0Trg = k_ptr->userInt("isMuonFromJpsi_dimuon0Trg");
      bool isMuonFromJpsi_jpsiTrkTrg = k_ptr->userInt("isMuonFromJpsi_jpsiTrkTrg");
      bool isMuonFromJpsi_jpsiTrk_PsiPrimeTrg = k_ptr->userInt("isMuonFromJpsi_jpsiTrk_PsiPrimeTrg");
      bool isMuonFromJpsi_jpsiTrk_NonResonantTrg = k_ptr->userInt("isMuonFromJpsi_jpsiTrk_NonResonantTrg");

      bool isUnpairedMuon_dimuon0 = (isDimuon0Trg && !isMuonFromJpsi_dimuon0Trg) && isDimuon_dimuon0Trg;
      bool isUnpairedMuon_jpsiTrk = (isJpsiTrkTrg && !isMuonFromJpsi_jpsiTrkTrg) && isDimuon_jpsiTrkTrg;
      bool isUnpairedMuon_jpsiTrk_PsiPrime = (isJpsiTrk_PsiPrimeTrg && !isMuonFromJpsi_jpsiTrk_PsiPrimeTrg) && isDimuon_jpsiTrk_PsiPrimeTrg;
      bool isUnpairedMuon_jpsiTrk_NonResonant = (isJpsiTrk_NonResonantTrg && !isMuonFromJpsi_jpsiTrk_NonResonantTrg) && isDimuon_jpsiTrk_NonResonantTrg;
      if(debug) std::cout<<"displaced muon:"<<k_ptr->pt()<<" isDImuon0Trg "<<isDimuon0Trg<<" isMuonFromJpsi_dimuon0Trg "<<isMuonFromJpsi_dimuon0Trg<<" isUnpairedMuon "<<isUnpairedMuon_dimuon0<<std::endl;
      //      if(!(isDimuon_jpsiTrkTrg || isUnpairedMuon)) continue;
      if(!(isUnpairedMuon_dimuon0 && isDimuon_dimuon0Trg) && 
         !(isUnpairedMuon_jpsiTrk && isDimuon_jpsiTrkTrg) && 
         !(isUnpairedMuon_jpsiTrk_PsiPrime && isDimuon_jpsiTrk_PsiPrimeTrg) && 
         !(isUnpairedMuon_jpsiTrk_NonResonant && isDimuon_jpsiTrk_NonResonantTrg)) continue;
      
      math::PtEtaPhiMLorentzVector k_p4(
                k_ptr->pt(), 
                k_ptr->eta(),
                k_ptr->phi(),
                k_ptr->mass()
                );
  
      //std::cout << "Here" << std::endl;
      // Use UserCands as they should not use memory but keep the Ptr itself
      // Put the muon passing the corresponding selection


      pat::CompositeCandidate cand;
      cand.setP4(ll_ptr->p4() + k_p4);
      cand.setCharge(ll_ptr->charge() + k_ptr->charge());

      // pv info

      cand.addUserInt("pvjpsi_idx", pvIdx);
      cand.addUserInt("pvfirst_idx", 0);
      cand.addUserInt("nPrimaryVertices", nPrimaryVertices);

      // tracks info
      cand.addUserCand("mu1", mu1_ptr);
      cand.addUserCand("mu2", mu2_ptr);
      cand.addUserCand("k", k_ptr);
      cand.addUserCand("dimuon", ll_ptr);
      
      cand.addUserInt("mu1_idx", mu1_idx);
      cand.addUserInt("mu2_idx", mu2_idx);
      cand.addUserInt("k_idx", k_idx);

      cand.addUserFloat("mu1_pvjpsi_dxy", mu1_pvjpsi_dxy);
      cand.addUserFloat("mu1_pvjpsi_dz", mu1_pvjpsi_dz);
      cand.addUserFloat("mu2_pvjpsi_dxy", mu2_pvjpsi_dxy);
      cand.addUserFloat("mu2_pvjpsi_dz", mu2_pvjpsi_dz);
      cand.addUserFloat("k_pvjpsi_dxy", k_pvjpsi_dxy);
      cand.addUserFloat("k_pvjpsi_dz", k_pvjpsi_dz);
      cand.addUserFloat("mu1_pvjpsi_dxyErr", mu1_pvjpsi_dxyErr);
      cand.addUserFloat("mu1_pvjpsi_dzErr", mu1_pvjpsi_dzErr);
      cand.addUserFloat("mu2_pvjpsi_dxyErr", mu2_pvjpsi_dxyErr);
      cand.addUserFloat("mu2_pvjpsi_dzErr", mu2_pvjpsi_dzErr);
      cand.addUserFloat("k_pvjpsi_dxyErr", k_pvjpsi_dxyErr);
      cand.addUserFloat("k_pvjpsi_dzErr", k_pvjpsi_dzErr);

      cand.addUserFloat("mu1_pvfirst_dxy", mu1_pvfirst_dxy);
      cand.addUserFloat("mu1_pvfirst_dz", mu1_pvfirst_dz);
      cand.addUserFloat("mu2_pvfirst_dxy", mu2_pvfirst_dxy);
      cand.addUserFloat("mu2_pvfirst_dz", mu2_pvfirst_dz);
      cand.addUserFloat("k_pvfirst_dxy", k_pvfirst_dxy);
      cand.addUserFloat("k_pvfirst_dz", k_pvfirst_dz);
      cand.addUserFloat("mu1_pvfirst_dxyErr", mu1_pvfirst_dxyErr);
      cand.addUserFloat("mu1_pvfirst_dzErr", mu1_pvfirst_dzErr);
      cand.addUserFloat("mu2_pvfirst_dxyErr", mu2_pvfirst_dxyErr);
      cand.addUserFloat("mu2_pvfirst_dzErr", mu2_pvfirst_dzErr);
      cand.addUserFloat("k_pvfirst_dxyErr", k_pvfirst_dxyErr);
      cand.addUserFloat("k_pvfirst_dzErr", k_pvfirst_dzErr);


      if(debug) std::cout<<"cand pt "<<cand.pt()<<std::endl;
      if(debug) std::cout<<"displ mu "<<k_ptr->pt()<<std::endl;
      if(debug) std::cout<<"displ m1 "<<mu1_ptr->pt()<<std::endl;
      if(debug) std::cout<<"displ m2 "<<mu2_ptr->pt()<<std::endl;

      auto dr_info = min_max_dr({mu1_ptr, mu2_ptr, k_ptr});

      cand.addUserFloat("min_dr", dr_info.first);
      cand.addUserFloat("max_dr", dr_info.second);
      // TODO add meaningful variables
      
      if( !pre_vtx_selection_(cand) ) continue;
      //if(debug) std::cout << "here2" << std::endl;
      
      //        if(debug) std::cout<<"PRIMA"<<std::endl;
      KinVtxFitter fitter(
        {muons_ttracks->at(mu1_idx), muons_ttracks->at(mu2_idx), muons_ttracks->at(k_idx)},
        {mu1_ptr->mass(), mu2_ptr->mass(), k_ptr->mass()},
        {LEP_SIGMA, LEP_SIGMA, LEP_SIGMA} //some small sigma for the muon mass
        );
      used_muon1_id.emplace_back(mu1_idx);
      used_muon2_id.emplace_back(mu2_idx);
      used_trk_id.emplace_back(k_idx);
      if(fitter.success()) {
      
	cand.setVertex( 
		       reco::Candidate::Point( 
					      fitter.fitted_vtx().x(),
					      fitter.fitted_vtx().y(),
					      fitter.fitted_vtx().z()
             )  
          );

	// pv shortest dz from jpsi + mu
	const reco::TransientTrack& threemuonTT = fitter.fitted_candidate_ttrk();
	int pvbIdx = getPVIdx(vertices, threemuonTT);
	if(debug) std::cout<<pvbIdx<<" "<<pvIdx<<std::endl;
	reco::Vertex pv_b = vertices->at(pvbIdx);
	double mu1_pvb_dxy = mu1_ptr->bestTrack()->dxy(pv_b.position());
	double mu2_pvb_dxy = mu2_ptr->bestTrack()->dxy(pv_b.position());
	double mu1_pvb_dz = mu1_ptr->bestTrack()->dz(pv_b.position());
	double mu2_pvb_dz = mu2_ptr->bestTrack()->dz(pv_b.position());
	double mu1_pvb_dxyErr = mu1_ptr->bestTrack()->dxyError(pv_b.position(),pv_b.covariance());
	double mu2_pvb_dxyErr = mu2_ptr->bestTrack()->dxyError(pv_b.position(),pv_b.covariance());
	double mu1_pvb_dzErr = mu1_ptr->bestTrack()->dzError();
	double mu2_pvb_dzErr = mu2_ptr->bestTrack()->dzError();
	double k_pvb_dxy = k_ptr->bestTrack()->dxy(pv_b.position());
	double k_pvb_dz = k_ptr->bestTrack()->dz(pv_b.position());
	double k_pvb_dxyErr = k_ptr->bestTrack()->dxyError(pv_b.position(),pv_b.covariance());
	double k_pvb_dzErr = k_ptr->bestTrack()->dzError();


	Measurement1D ip3D_pvjpsi = getIP(ll_ptr, pv_jpsi, muons_ttracks->at(k_idx));
	Measurement1D ip3D_pvb = getIP(ll_ptr, pv_b, muons_ttracks->at(k_idx));
	Measurement1D ip3D_pvfirst = getIP(ll_ptr, pv_first, muons_ttracks->at(k_idx));


	cand.addUserInt("pvb_idx", pvbIdx);
	cand.addUserFloat("mu1_pvb_dxy", mu1_pvb_dxy);
	cand.addUserFloat("mu1_pvb_dz", mu1_pvb_dz);
	cand.addUserFloat("mu2_pvb_dxy", mu2_pvb_dxy);
	cand.addUserFloat("mu2_pvb_dz", mu2_pvb_dz);
	cand.addUserFloat("k_pvb_dxy", k_pvb_dxy);
	cand.addUserFloat("k_pvb_dz", k_pvb_dz);
	cand.addUserFloat("mu1_pvb_dxyErr", mu1_pvb_dxyErr);
	cand.addUserFloat("mu1_pvb_dzErr", mu1_pvb_dzErr);
	cand.addUserFloat("mu2_pvb_dxyErr", mu2_pvb_dxyErr);
	cand.addUserFloat("mu2_pvb_dzErr", mu2_pvb_dzErr);
	cand.addUserFloat("k_pvb_dxyErr", k_pvb_dxyErr);
	cand.addUserFloat("k_pvb_dzErr", k_pvb_dzErr);

        auto lxy = l_xy(fitter, *beamspot);
        cand.addUserFloat("l_xy", lxy.value());
        cand.addUserFloat("l_xy_unc", lxy.error());
        cand.addUserFloat("ip3D_pvjpsi", ip3D_pvjpsi.value());
        cand.addUserFloat("ip3D_pvjpsi_e", ip3D_pvjpsi.error());
        cand.addUserFloat("ip3D_pvb", ip3D_pvb.value());
        cand.addUserFloat("ip3D_pvb_e", ip3D_pvb.error());
        cand.addUserFloat("ip3D_pvfirst", ip3D_pvfirst.value());
        cand.addUserFloat("ip3D_pvfirst_e", ip3D_pvfirst.error());
        cand.addUserInt("sv_OK" , fitter.success());
        cand.addUserFloat("sv_chi2", fitter.chi2());
        cand.addUserFloat("sv_ndof", fitter.dof()); // float??
        cand.addUserFloat("sv_prob", fitter.prob());
        cand.addUserFloat("fitted_mll" , (fitter.daughter_p4(0) + fitter.daughter_p4(1)).mass());
        auto fit_p4 = fitter.fitted_p4();
        cand.addUserFloat("fitted_pt"  , fit_p4.pt()); 
        cand.addUserFloat("fitted_eta" , fit_p4.eta());
        cand.addUserFloat("fitted_phi" , fit_p4.phi());
        cand.addUserFloat("fitted_mass", fitter.fitted_candidate().mass());      
        cand.addUserFloat("fitted_massErr", sqrt(fitter.fitted_candidate().kinematicParametersError().matrix()(6,6)));      
        cand.addUserFloat(
            "cos_theta_2D", 
            cos_theta_2D(fitter, *beamspot, cand.p4())
            );
        cand.addUserFloat(
            "fitted_cos_theta_2D", 
            cos_theta_2D(fitter, *beamspot, fit_p4)
            );

        cand.addUserFloat("vtx_ex", sqrt(fitter.fitted_vtx_uncertainty().cxx()));
        cand.addUserFloat("vtx_ey", sqrt(fitter.fitted_vtx_uncertainty().cyy()));
        cand.addUserFloat("vtx_ez", sqrt(fitter.fitted_vtx_uncertainty().czz()));
        cand.addUserFloat("vtx_chi2", ChiSquaredProbability(fitter.chi2(), fitter.dof()));
        cand.addUserFloat("fitted_mu1_pt" , fitter.daughter_p4(0).pt()); 
        cand.addUserFloat("fitted_mu1_eta", fitter.daughter_p4(0).eta());
        cand.addUserFloat("fitted_mu1_phi", fitter.daughter_p4(0).phi());
        cand.addUserFloat("fitted_mu2_pt" , fitter.daughter_p4(1).pt()); 
        cand.addUserFloat("fitted_mu2_eta", fitter.daughter_p4(1).eta());
        cand.addUserFloat("fitted_mu2_phi", fitter.daughter_p4(1).phi());
        cand.addUserFloat("fitted_k_pt"  , fitter.daughter_p4(2).pt()); 
        cand.addUserFloat("fitted_k_eta" , fitter.daughter_p4(2).eta());
        cand.addUserFloat("fitted_k_phi" , fitter.daughter_p4(2).phi());

	cand.addUserFloat("pvb_x", pv_b.position().x());
	cand.addUserFloat("pvb_y", pv_b.position().y());
	cand.addUserFloat("pvb_z", pv_b.position().z());
	cand.addUserFloat("pvb_ex", pv_b.covariance(0,0));
	cand.addUserFloat("pvb_ey", pv_b.covariance(1,1));
	cand.addUserFloat("pvb_ez", pv_b.covariance(2,2));
	cand.addUserFloat("pvb_exy", pv_b.covariance(0,1));
	cand.addUserFloat("pvb_eyz", pv_b.covariance(0,2));
	cand.addUserFloat("pvb_exz", pv_b.covariance(1,2));
	cand.addUserFloat("pvb_chi2", ChiSquaredProbability(pv_b.chi2(), pv_b.ndof()));

      }
      else
      {
        cand.setVertex(reco::Candidate::Point(0.,0.,0.));
        cand.addUserFloat("l_xy", -99.);
        cand.addUserFloat("l_xy_unc", -99.);
        cand.addUserFloat("ip3D_pvjpsi", -99.);
        cand.addUserFloat("ip3D_pvjpsi_e", -99.);
        cand.addUserFloat("ip3D_pvb", -99.);
        cand.addUserFloat("ip3D_pvb_e", -99.);
        cand.addUserFloat("ip3D_pvfirst", -99.);
        cand.addUserFloat("ip3D_pvfirst_e", -99.);
        cand.addUserInt("sv_OK" , fitter.success());
        cand.addUserFloat("sv_chi2", -99.);
        cand.addUserFloat("sv_ndof", -99.); // float??
        cand.addUserFloat("sv_prob", -99.);
        cand.addUserFloat("fitted_mll" , -99.);
        cand.addUserFloat("fitted_pt"  , -99.); 
        cand.addUserFloat("fitted_eta" , -99.);
        cand.addUserFloat("fitted_phi" , -99.);
        cand.addUserFloat("fitted_mass", -99.);      
        cand.addUserFloat("fitted_massErr", -99.);      
        cand.addUserFloat("cos_theta_2D", -99.);
        cand.addUserFloat("fitted_cos_theta_2D", -99.);
        cand.addUserFloat("vtx_ex", -99.);
        cand.addUserFloat("vtx_ey", -99.);
        cand.addUserFloat("vtx_ez", -99.);
        cand.addUserFloat("vtx_chi2", -99.);
        cand.addUserFloat("fitted_mu1_pt" , -99.); 
        cand.addUserFloat("fitted_mu1_eta", -99.);
        cand.addUserFloat("fitted_mu1_phi", -99.);
        cand.addUserFloat("fitted_mu2_pt" , -99.); 
        cand.addUserFloat("fitted_mu2_eta", -99.);
        cand.addUserFloat("fitted_mu2_phi", -99.);
        cand.addUserFloat("fitted_k_pt"  , -99.); 
        cand.addUserFloat("fitted_k_eta" , -99.);
        cand.addUserFloat("fitted_k_phi" , -99.);

	cand.addUserInt("pvb_idx", -99.);
	cand.addUserFloat("mu1_pvb_dxy", -99.);
	cand.addUserFloat("mu1_pvb_dz", -99.);
	cand.addUserFloat("mu2_pvb_dxy", -99.);
	cand.addUserFloat("mu2_pvb_dz", -99.);
	cand.addUserFloat("k_pvb_dxy", -99.);
	cand.addUserFloat("k_pvb_dz", -99.);
	cand.addUserFloat("mu1_pvb_dxyErr", -99.);
	cand.addUserFloat("mu1_pvb_dzErr", -99.);
	cand.addUserFloat("mu2_pvb_dxyErr", -99.);
	cand.addUserFloat("mu2_pvb_dzErr", -99.);
	cand.addUserFloat("k_pvb_dxyErr", -99.);
	cand.addUserFloat("k_pvb_dzErr", -99.);
	cand.addUserFloat("pvb_x", -99.);
	cand.addUserFloat("pvb_y", -99.);
	cand.addUserFloat("pvb_z", -99.);
	cand.addUserFloat("pvb_ex", -99.);
	cand.addUserFloat("pvb_ey", -99.);
	cand.addUserFloat("pvb_ez", -99.);
	cand.addUserFloat("pvb_exy", -99.);
	cand.addUserFloat("pvb_eyz", -99.);
	cand.addUserFloat("pvb_exz", -99.);
	cand.addUserFloat("pvb_chi2", -99.);

      }
      cand.addUserFloat("vtx_x", cand.vx());
      cand.addUserFloat("vtx_y", cand.vy());
      cand.addUserFloat("vtx_z", cand.vz());

      //alredy defined in dimuon builder
      /*
      cand.addUserFloat("jpsi_vtx_x", ll_ptr->userFloat("vtx_x"));
      cand.addUserFloat("jpsi_vtx_y", ll_ptr->userFloat("vtx_y"));
      cand.addUserFloat("jpsi_vtx_z", ll_ptr->userFloat("vtx_z"));
      cand.addUserFloat("jpsi_vtx_ex", ll_ptr->userFloat("vtx_ex"));
      cand.addUserFloat("jpsi_vtx_ey", ll_ptr->userFloat("vtx_ey"));
      cand.addUserFloat("jpsi_vtx_ez", ll_ptr->userFloat("vtx_ez"));
      cand.addUserFloat("jpsi_vtx_chi2", ll_ptr->userFloat("vtx_chi2"));
      */
      cand.addUserFloat("pvjpsi_x", pv_jpsi.position().x());
      cand.addUserFloat("pvjpsi_y", pv_jpsi.position().y());
      cand.addUserFloat("pvjpsi_z", pv_jpsi.position().z());
      cand.addUserFloat("pvjpsi_ex", pv_jpsi.covariance(0,0));
      cand.addUserFloat("pvjpsi_ey", pv_jpsi.covariance(1,1));
      cand.addUserFloat("pvjpsi_ez", pv_jpsi.covariance(2,2));
      cand.addUserFloat("pvjpsi_exy", pv_jpsi.covariance(0,1));
      cand.addUserFloat("pvjpsi_eyz", pv_jpsi.covariance(0,2));
      cand.addUserFloat("pvjpsi_exz", pv_jpsi.covariance(1,2));
      cand.addUserFloat("pvjpsi_chi2", ChiSquaredProbability(pv_jpsi.chi2(), pv_jpsi.ndof()));

      cand.addUserFloat("pvfirst_x", pv_first.position().x());
      cand.addUserFloat("pvfirst_y", pv_first.position().y());
      cand.addUserFloat("pvfirst_z", pv_first.position().z());
      cand.addUserFloat("pvfirst_ex", pv_first.covariance(0,0));
      cand.addUserFloat("pvfirst_ey", pv_first.covariance(1,1));
      cand.addUserFloat("pvfirst_ez", pv_first.covariance(2,2));
      cand.addUserFloat("pvfirst_exy", pv_first.covariance(0,1));
      cand.addUserFloat("pvfirst_eyz", pv_first.covariance(0,2));
      cand.addUserFloat("pvfirst_exz", pv_first.covariance(1,2));
      cand.addUserFloat("pvfirst_chi2", ChiSquaredProbability(pv_first.chi2(), pv_first.ndof()));

      const reco::BeamSpot &bm = *beamspot;

      cand.addUserFloat("beamspot_x", bm.x0());
      cand.addUserFloat("beamspot_y", bm.y0());
      cand.addUserFloat("beamspot_z", bm.z0());

      //mie variabili                                                                                
      //conto quanti mu totali aveva l'evento
      //cand.addUserInt("pass_3mu",pass_3mu.size());                                 
      
      float B_pt=(Bc_MASS/cand.mass())*cand.pt();
      
      TLorentzVector P_b;
      P_b.SetPtEtaPhiM(B_pt,cand.eta(),cand.phi(),Bc_MASS);
      
      TLorentzVector P_k;
      P_k.SetPtEtaPhiM(k_ptr->pt(),k_ptr->eta(),k_ptr->phi(),k_ptr->mass());
      
      TLorentzVector P_mu1;
      P_mu1.SetPtEtaPhiM(mu1_ptr->pt(),mu1_ptr->eta(),mu1_ptr->phi(),mu1_ptr->mass());
      
      TLorentzVector P_mu2;
      P_mu2.SetPtEtaPhiM(mu2_ptr->pt(),mu2_ptr->eta(),mu2_ptr->phi(),mu2_ptr->mass());
      
      float m_miss_2=(P_b-P_k-P_mu1-P_mu2)*(P_b-P_k-P_mu1-P_mu2);
      float Q_2=(P_b-P_mu1-P_mu2)*(P_b-P_mu1-P_mu2);
      float pt_miss=(P_b.Pt()-P_k.Pt()-P_mu1.Pt()-P_mu2.Pt());
      float pt_miss_vec=((P_b-P_k-P_mu1-P_mu2).Pt());
      float pt_var=((P_mu1+P_mu2).Pt()-P_k.Pt());
      float DR=deltaR(P_mu1.Eta(),P_mu1.Phi(),P_mu2.Eta(),P_mu2.Phi());
      //float deta = P_mu1.Eta() - P_mu2.Eta();
      //float dphi = P_mu1.Phi() - P_mu2.Phi();
      //float DR_2=std::sqrt(deta * deta + dphi * dphi);        
      float m_jpsi=sqrt((P_mu1+P_mu2)*(P_mu1+P_mu2));
      cand.addUserFloat("m_miss_2", m_miss_2);
      cand.addUserFloat("Q_2",Q_2);
      cand.addUserFloat("pt_miss",pt_miss);
      cand.addUserFloat("pt_miss_vec",pt_miss_vec);
      cand.addUserFloat("pt_var",pt_var);
      cand.addUserFloat("DR",DR);
      cand.addUserFloat("m_jpsi",m_jpsi);

      //energia del mu unpaired in diversi sistemi di riferimento                  
      TLorentzVector P_mu=P_k;        
      TVector3 mu_beta_lab=P_b.BoostVector();

      P_mu.Boost(-mu_beta_lab);
      cand.addUserFloat("E_mu_star",P_mu.E());
      P_mu=P_k;        
      TLorentzVector jpsi=P_mu1+P_mu2;
      TVector3 jpsi_beta_lab=jpsi.BoostVector();
      P_mu.Boost(-jpsi_beta_lab);
      cand.addUserFloat("E_mu_#",P_mu.E());
      //if(debug) std::cout<<"post vertex features "<<std::cout << "postvtx" <<fitter.success()<<" "<<fitter.prob()<<" "<<cos_theta_2D(fitter, *beamspot, fit_p4)<<" "<<cand.mass()<< std::endl;
      if( !post_vtx_selection_(cand) ) {
	if(debug) std::cout<<"post vertx sel dies"<<std::endl;
	continue;        
      }
      if(debug) std::cout<<"pass post vertx sel"<<std::endl;

      //compute isolation
      float mu1_iso03 = 0;
      float mu1_iso04 = 0;
      float mu2_iso03 = 0;
      float mu2_iso04 = 0;
      float k_iso03  = 0;
      float k_iso04  = 0;
      float b_iso03  = 0;
      float b_iso04  = 0;
      
      for( unsigned int iTrk=0; iTrk<totalTracks; ++iTrk ) 
      {
        const pat::PackedCandidate & trk = (iTrk < nTracks) ? (*iso_tracks)[iTrk] : (*iso_lostTracks)[iTrk-nTracks];
        // define selections for iso tracks (pT, eta, ...)
        if( !isotrk_selection_(trk) ) continue;
        // check if the track is the muon

	// only consider tracks originating close to the three bodies
	if ( !mu1_ptr->bestTrack() || fabs(trk.dz() - mu1_ptr->bestTrack()->dz()) > 0.4 ) continue;
	if ( !mu2_ptr->bestTrack() || fabs(trk.dz() - mu2_ptr->bestTrack()->dz()) > 0.4 ) continue;
	if ( !k_ptr ->bestTrack() || fabs(trk.dz() - k_ptr ->bestTrack()->dz()) > 0.4 ) continue;

        if (k_ptr->userCand("cand") ==  edm::Ptr<reco::Candidate> ( iso_tracks, iTrk ) ) {
          
          if(debug) std::cout<<"old"<<std::endl;
          continue;
        }
        
        if(track_to_muon_match(k_ptr, iso_tracks.id(), iTrk)) 
        {
          // if(debug) std::cout<<"new"<<std::endl;
          continue;
        }
        // check if the track is one of the two muons
        if (track_to_muon_match(mu1_ptr, iso_tracks.id(), iTrk) || 
            track_to_muon_match(mu2_ptr, iso_tracks.id(), iTrk) ) continue;
        
        // add to final particle iso if dR < cone
        float dr_to_mu1 = deltaR(cand.userFloat("fitted_mu1_eta"), cand.userFloat("fitted_mu1_phi"), trk.eta(), trk.phi());
        float dr_to_mu2 = deltaR(cand.userFloat("fitted_mu2_eta"), cand.userFloat("fitted_mu2_phi"), trk.eta(), trk.phi());
        float dr_to_k  = deltaR(cand.userFloat("fitted_k_eta") , cand.userFloat("fitted_k_phi") , trk.eta(), trk.phi());
        float dr_to_b  = deltaR(cand.userFloat("fitted_eta")   , cand.userFloat("fitted_phi") , trk.eta(), trk.phi());
        
	if (dr_to_mu1 < 0.4 && dr_to_mu1>0.01){
	  mu1_iso04 += trk.pt();
	  if ( dr_to_mu1 < 0.3) mu1_iso03 += trk.pt();
	}
	if (dr_to_mu2 < 0.4 && dr_to_mu2>0.01){
	  mu2_iso04 += trk.pt();
	  if (dr_to_mu2 < 0.3)  mu2_iso03 += trk.pt();
	}
	if (dr_to_k < 0.4 && dr_to_k>0.01){
	  k_iso04 += trk.pt();
	  if (dr_to_k < 0.3) k_iso03 += trk.pt();
	}
	if (dr_to_b < 0.4){
	  b_iso04 += trk.pt();
	  if (dr_to_b < 0.3) b_iso03 += trk.pt();
	}      }

      cand.addUserFloat("mu1_iso03", mu1_iso03);
      cand.addUserFloat("mu1_iso04", mu1_iso04);
      cand.addUserFloat("mu2_iso03", mu2_iso03);
      cand.addUserFloat("mu2_iso04", mu2_iso04);
      cand.addUserFloat("k_iso03" , k_iso03 );
      cand.addUserFloat("k_iso04" , k_iso04 );
      cand.addUserFloat("b_iso03" , b_iso03 );
      cand.addUserFloat("b_iso04" , b_iso04 );

      //salvo il candidato
      ret_val->push_back(cand);
    }//for(size_t ll_idx = 0; ll_idx < dimuons->size(); ++ll_idx) 
  }//for(size_t k_idx = 0; k_idx < muons->size(); ++k_idx)
  
  for (auto & cand: *ret_val){
    cand.addUserInt("n_k_used", std::count(used_trk_id.begin(),used_trk_id.end(),cand.userInt("k_idx")));
    cand.addUserInt("n_mu1_used", std::count(used_muon1_id.begin(),used_muon1_id.end(),cand.userInt("mu1_idx"))+std::count(used_muon2_id.begin(),used_muon2_id.end(),cand.userInt("mu1_idx")));
    cand.addUserInt("n_mu2_used", std::count(used_muon1_id.begin(),used_muon1_id.end(),cand.userInt("mu2_idx"))+std::count(used_muon2_id.begin(),used_muon2_id.end(),cand.userInt("mu2_idx")));
  }
  evt.put(std::move(ret_val));
}//produce  

Measurement1D BTo3MuBuilder::getIP(edm::Ptr<pat::CompositeCandidate> ll_ptr, reco::Vertex pv, reco::TransientTrack transientTrackMu) const
{
  // computing 3d impact parameter
  //GlobalVector jpsiDirection(fitter.fitted_p4().x(), fitter.fitted_p4().y(), fitter.fitted_p4().z());
  //reco::Vertex::Point jpsiVertexPosition(fitter.fitted_vtx().x(),fitter.fitted_vtx().y(),fitter.fitted_vtx().z());
  reco::Vertex::Point jpsiVertexPosition(ll_ptr->userFloat("vtx_x"),ll_ptr->userFloat("vtx_y"),ll_ptr->userFloat("vtx_z"));

  /*const double err00 = fitter.fitted_candidate().kinematicParametersError().matrix()(0,0);
  const double err11 = fitter.fitted_candidate().kinematicParametersError().matrix()(1,1);
  const double err22 = fitter.fitted_candidate().kinematicParametersError().matrix()(2,2);
  const double err01 = fitter.fitted_candidate().kinematicParametersError().matrix()(0,1);
  const double err02 = fitter.fitted_candidate().kinematicParametersError().matrix()(0,2);
  const double err12 = fitter.fitted_candidate().kinematicParametersError().matrix()(1,2);
  */
  
  reco::Vertex::Error jpsiVertexError;

  jpsiVertexError(0,0) = ll_ptr->userFloat("jpsi_err00");
  jpsiVertexError(0,1) = ll_ptr->userFloat("jpsi_err01");
  jpsiVertexError(0,2) = ll_ptr->userFloat("jpsi_err02");
  jpsiVertexError(1,0) = ll_ptr->userFloat("jpsi_err01");
  jpsiVertexError(1,1) = ll_ptr->userFloat("jpsi_err11");
  jpsiVertexError(1,2) = ll_ptr->userFloat("jpsi_err12");
  jpsiVertexError(2,0) = ll_ptr->userFloat("jpsi_err02");
  jpsiVertexError(2,1) = ll_ptr->userFloat("jpsi_err12");
  jpsiVertexError(2,2) = ll_ptr->userFloat("jpsi_err22");

  //jpsi vtx - pv
  //  GlobalVector jpsiGlobalVector(fitter.fitted_p4().x() - pv.position().x(), fitter.fitted_p4().y() - pv.position().y(), fitter.fitted_p4().z() - pv.position().z());
  GlobalVector jpsiGlobalVector(ll_ptr->userFloat("vtx_x") - pv.position().x(), ll_ptr->userFloat("vtx_y") - pv.position().y(), ll_ptr->userFloat("vtx_z") - pv.position().z());
  //const reco::Vertex jpsiVertex(jpsiVertexPosition, jpsiVertexError, fitter.chi2(), fitter.dof(), 2);
  const reco::Vertex jpsiVertex(jpsiVertexPosition, jpsiVertexError, ll_ptr->userFloat("sv_chi2"), ll_ptr->userFloat("sv_ndof"), 2);

  SignedImpactParameter3D signed_ip3D;
  Measurement1D ip3D = signed_ip3D.apply(transientTrackMu,jpsiGlobalVector,jpsiVertex).second;
  return ip3D;
}

int BTo3MuBuilder::getPVIdx(const reco::VertexCollection* vertices,const reco::TransientTrack& dimuonTT) const
{
    double dzMin = 1000000.;
    reco::Vertex bestVertex;
    int pvIdx = 0;
    //const reco::VertexCollection* vertices = thePrimaryVerticesHandle.product();
    for(size_t i = 0; i < vertices->size() ; i++)
    {
      reco::Vertex primVertex = vertices->at(i);
      //std::cout << "prim vertex z: " << primVertex->z() << std::endl;
      if (abs(dzMin) > abs(dimuonTT.track().dz(primVertex.position())))
      {
        bestVertex = primVertex;
        pvIdx = i;
        //bestVertex = primVertex;
        dzMin = dimuonTT.track().dz(primVertex.position());
      }
    }
    if(debug) std::cout<< "Best vertex x: " << bestVertex.x() << std::endl;
    if(debug) std::cout<< "Best vertex id: " << pvIdx << std::endl;
  return pvIdx;
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(BTo3MuBuilder);

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


class BTo2Mu3PiBuilder : public edm::global::EDProducer<> {

  // perhaps we need better structure here (begin run etc)
public:
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;
  typedef std::vector<reco::GenParticle> GenParticleCollection;

  explicit BTo2Mu3PiBuilder(const edm::ParameterSet &cfg):
    particle_selection_{cfg.getParameter<std::string>("particleSelection")},
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    dimuons_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("dimuons") )},
    primaryVertices_{consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("primaryVertices"))},
    particles_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("particles") )},
    particles_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("particlesTransientTracks") )},
    muons_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("muonsTransientTracks") )},
    //kaons_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("kaonsTransientTracks") )},
    isotracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
    isolostTracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))),
    isotrk_selection_{cfg.getParameter<std::string>("isoTracksSelection")},
    beamspot_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") )} 
  {
      produces<pat::CompositeCandidateCollection>();
  }

  ~BTo2Mu3PiBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
  int  getPVIdx(const reco::VertexCollection*,const reco::TransientTrack&) const;
  Measurement1D getIP(KinVtxFitter fitter, reco::TransientTrack transientTrackMu) const;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
  private:
  const StringCutObjectSelector<pat::CompositeCandidate> particle_selection_; 
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; // cut on the di-muon before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; // cut on the di-muon after the SV fit

  const edm::EDGetTokenT<pat::CompositeCandidateCollection> dimuons_;
  const edm::EDGetTokenT<reco::VertexCollection> primaryVertices_;
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> particles_;
  const edm::EDGetTokenT<TransientTrackCollection> particles_ttracks_;
  const edm::EDGetTokenT<TransientTrackCollection> muons_ttracks_;

  const edm::EDGetTokenT<pat::PackedCandidateCollection> isotracksToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> isolostTracksToken_;

  const StringCutObjectSelector<pat::PackedCandidate> isotrk_selection_; 


  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;  
};

void BTo2Mu3PiBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  //input
  edm::Handle<pat::CompositeCandidateCollection> dimuons;
  evt.getByToken(dimuons_, dimuons);

  edm::Handle<reco::VertexCollection> primaryVertices;
  evt.getByToken(primaryVertices_, primaryVertices);
  
  edm::Handle<TransientTrackCollection> particles_ttracks;
  evt.getByToken(particles_ttracks_, particles_ttracks);

  edm::Handle<pat::CompositeCandidateCollection> particles;
  evt.getByToken(particles_, particles);
  
  edm::Handle<TransientTrackCollection> muons_ttracks;
  evt.getByToken(muons_ttracks_, muons_ttracks);  

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

  std::vector<int> used_muon1_id, used_muon2_id, used_pi1_id,used_pi2_id,used_pi3_id;

  std::unique_ptr<pat::CompositeCandidateCollection> ret_val(new pat::CompositeCandidateCollection());

  const reco::VertexCollection* vertices = primaryVertices.product();
  int nPrimaryVertices = vertices->size();
  // output
  for(size_t ll_idx = 0; ll_idx < dimuons->size(); ++ll_idx) 
  {
    if(debug) std::cout<<"Begining of the dimuon loop"<<particles->size()<<std::endl;
    //std::cout << "PV " << ll_idx << ": " << vertices->at(ll_idx).position() << std::endl;
    edm::Ptr<pat::CompositeCandidate> ll_ptr(dimuons, ll_idx);
    edm::Ptr<reco::Candidate> mu1_ptr = ll_ptr->userCand("mu1");
    edm::Ptr<reco::Candidate> mu2_ptr = ll_ptr->userCand("mu2");
    size_t mu1_idx = abs(ll_ptr->userInt("mu1_idx"));
    size_t mu2_idx = abs(ll_ptr->userInt("mu2_idx"));

    int pvIdx = ll_ptr->userInt("pvIdx");
    reco::Vertex pv_jpsi = vertices->at(pvIdx);
    double mu1_pvjpsi_dxy = mu1_ptr->bestTrack()->dxy(pv_jpsi.position());
    double mu1_pvjpsi_dz = mu1_ptr->bestTrack()->dz(pv_jpsi.position());
    double mu2_pvjpsi_dxy = mu2_ptr->bestTrack()->dxy(pv_jpsi.position());
    double mu2_pvjpsi_dz = mu2_ptr->bestTrack()->dz(pv_jpsi.position());
    double mu1_pvjpsi_dxyErr = mu1_ptr->bestTrack()->dxyError(pv_jpsi.position(),pv_jpsi.covariance());
    double mu2_pvjpsi_dxyErr = mu2_ptr->bestTrack()->dxyError(pv_jpsi.position(),pv_jpsi.covariance());
    double mu1_pvjpsi_dzErr = mu1_ptr->bestTrack()->dzError();
    double mu2_pvjpsi_dzErr = mu2_ptr->bestTrack()->dzError();

    size_t isDimuon_dimuon0Trg = abs(ll_ptr->userInt("muonpair_fromdimuon0"));
    size_t isDimuon_jpsiTrkTrg = abs(ll_ptr->userInt("muonpair_fromjpsitrk"));
    //size_t isDimuon_jpsiTrkTrg = abs(ll_ptr->userInt("isJpsiTrkTrg"));
    //size_t isDimuon_dimuon0Trg = abs(ll_ptr->userInt("isDimuon0Trg"));
    if(!(isDimuon_jpsiTrkTrg)) continue;

    // first loop on pion- this one with trigger matching
    if(debug) std::cout<<"paerticles size "<<particles->size()<<std::endl;
    for(size_t pi1_idx = 0; pi1_idx < particles->size(); ++pi1_idx) {
      edm::Ptr<pat::CompositeCandidate> pi1_ptr(particles, pi1_idx);
      if( !particle_selection_(*pi1_ptr) ) continue;
      //dz requirement
      if ( fabs(particles_ttracks->at(pi1_idx).track().dz() - mu1_ptr->bestTrack()->dz()) > 0.4 || fabs(particles_ttracks->at(pi1_idx).track().dz() - mu2_ptr->bestTrack()->dz()) > 0.4) continue;
      
      //dR requirement
      if(deltaR(muons_ttracks->at(mu2_idx).track().eta(),muons_ttracks->at(mu2_idx).track().phi(),particles_ttracks->at(pi1_idx).track().eta(),particles_ttracks->at(pi1_idx).track().phi()) < 0.1 ) continue;
      if(deltaR(muons_ttracks->at(mu1_idx).track().eta(),muons_ttracks->at(mu1_idx).track().phi(),particles_ttracks->at(pi1_idx).track().eta(),particles_ttracks->at(pi1_idx).track().phi()) < 0.1 ) continue;
      
      bool isPartTrg = pi1_ptr->userInt("isTriggering");
      //ha trovato il mu displaced
      if(!(isPartTrg)) {
	//if(debug) std::cout<<"is NOT track triggered "<<k_ptr->pt()<<std::endl;
	continue;
      }
      

      math::PtEtaPhiMLorentzVector pi1_p4(
					  pi1_ptr->pt(), 
					  pi1_ptr->eta(),
					  pi1_ptr->phi(),
					  PI_MASS
					  );
  
      //loop pion 2
      for(size_t pi2_idx = 0; pi2_idx < particles->size(); ++pi2_idx) {
	edm::Ptr<pat::CompositeCandidate> pi2_ptr(particles, pi2_idx);
	if( !particle_selection_(*pi2_ptr) ) continue;
	if(pi2_idx == pi1_idx) continue;
	// dz between track and leptons
	if ( fabs(particles_ttracks->at(pi2_idx).track().dz() - mu1_ptr->bestTrack()->dz()) > 0.4 ||  fabs(particles_ttracks->at(pi2_idx).track().dz() - mu2_ptr->bestTrack()->dz()) > 0.4) continue;
	//DR between tracks and leptons
	if(deltaR(muons_ttracks->at(mu2_idx).track().eta(),muons_ttracks->at(mu2_idx).track().phi(),particles_ttracks->at(pi2_idx).track().eta(),particles_ttracks->at(pi2_idx).track().phi()) < 0.1 ) continue;
	if(deltaR(muons_ttracks->at(mu1_idx).track().eta(),muons_ttracks->at(mu1_idx).track().phi(),particles_ttracks->at(pi2_idx).track().eta(),particles_ttracks->at(pi2_idx).track().phi()) < 0.1 ) continue;

	math::PtEtaPhiMLorentzVector pi2_p4(
					    pi2_ptr->pt(),
					    pi2_ptr->eta(),
					    pi2_ptr->phi(),
                                                  PI_MASS
					    );
      // Use UserCands as they should not use memory but keep the Ptr itself
      // Put the muon passing the corresponding selection

	  //loop on the third pion
	  for(size_t pi3_idx = 0; pi3_idx < particles->size(); ++pi3_idx) {
	    edm::Ptr<pat::CompositeCandidate> pi3_ptr(particles, pi3_idx);
	    if( !particle_selection_(*pi3_ptr) ) continue;
	    if(pi3_idx == pi1_idx or pi3_idx == pi2_idx) continue;
	    //dz requirement
	    if ( fabs(particles_ttracks->at(pi3_idx).track().dz() - mu1_ptr->bestTrack()->dz()) > 0.4 ||  fabs(particles_ttracks->at(pi3_idx).track().dz() - mu2_ptr->bestTrack()->dz()) > 0.4) continue;
	    //DR requirements
	    if(deltaR(muons_ttracks->at(mu2_idx).track().eta(),muons_ttracks->at(mu2_idx).track().phi(),particles_ttracks->at(pi3_idx).track().eta(),particles_ttracks->at(pi3_idx).track().phi()) < 0.1 ) continue;
	    if(deltaR(muons_ttracks->at(mu1_idx).track().eta(),muons_ttracks->at(mu1_idx).track().phi(),particles_ttracks->at(pi3_idx).track().eta(),particles_ttracks->at(pi3_idx).track().phi()) < 0.1 ) continue;

	    if(debug) std::cout<<"before dz "<<std::endl;

      double pi1_dxy = particles_ttracks->at(pi1_idx).track().dxy(pv_jpsi.position());
      double pi1_dz = particles_ttracks->at(pi1_idx).track().dz(pv_jpsi.position());
      double pi2_dxy = particles_ttracks->at(pi2_idx).track().dxy(pv_jpsi.position());
      double pi2_dz = particles_ttracks->at(pi2_idx).track().dz(pv_jpsi.position());
      double pi3_dxy = particles_ttracks->at(pi3_idx).track().dxy(pv_jpsi.position());
      double pi3_dz = particles_ttracks->at(pi3_idx).track().dz(pv_jpsi.position());
      double pi1_dxyErr = particles_ttracks->at(pi1_idx).track().dxyError(pv_jpsi.position(),pv_jpsi.covariance());
      double pi1_dzErr = particles_ttracks->at(pi1_idx).track().dzError();
      double pi2_dxyErr = particles_ttracks->at(pi2_idx).track().dxyError(pv_jpsi.position(),pv_jpsi.covariance());
      double pi2_dzErr = particles_ttracks->at(pi2_idx).track().dzError();
      double pi3_dxyErr = particles_ttracks->at(pi3_idx).track().dxyError(pv_jpsi.position(),pv_jpsi.covariance());
      double pi3_dzErr = particles_ttracks->at(pi3_idx).track().dzError();

	    if(debug) std::cout<<"p1 dxy "<<pi1_dxy<<std::endl;
	    if(debug) std::cout<<"p1 dz "<<pi1_dz<<std::endl;

	    if(debug) std::cout<<"after dz "<<std::endl;
	    math::PtEtaPhiMLorentzVector pi3_p4(
						pi3_ptr->pt(),
						pi3_ptr->eta(),
						pi3_ptr->phi(),
                                                    PI_MASS
						);

	    //the other two pions don't need the trigger matching
	    pat::CompositeCandidate cand;
	    cand.setP4(ll_ptr->p4() + + pi1_p4 + pi2_p4 + pi3_p4);
	    cand.setCharge(ll_ptr->charge() + pi1_ptr->charge() + pi2_ptr->charge() + pi3_ptr->charge());
	    if(debug) std::cout<<"cand pt "<<cand.pt()<<std::endl;
	    if(debug) std::cout<<"displ p1 "<<pi1_ptr->pt()<<std::endl;
	    if(debug) std::cout<<"displ p2 "<<pi2_ptr->pt()<<std::endl;
	    if(debug) std::cout<<"displ p3 "<<pi3_ptr->pt()<<std::endl;
	    if(debug) std::cout<<"displ m1 "<<mu1_ptr->pt()<<std::endl;
	    if(debug) std::cout<<"displ m2 "<<mu2_ptr->pt()<<std::endl;

      // pv info

      cand.addUserInt("pvjpsi_idx", pvIdx);
      cand.addUserInt("nPrimaryVertices", nPrimaryVertices);

      // tracks info
	    
	    cand.addUserCand("mu1", mu1_ptr);
	    cand.addUserCand("mu2", mu2_ptr);
	    cand.addUserCand("pi1", pi1_ptr);
	    cand.addUserCand("pi2", pi2_ptr);
	    cand.addUserCand("pi3", pi3_ptr);

	    cand.addUserCand("dimuon", ll_ptr);
	    
	    cand.addUserInt("mu1_idx", mu1_idx);
	    cand.addUserInt("mu2_idx", mu2_idx);
	    cand.addUserInt("pi1_idx", pi1_idx);
	    cand.addUserInt("pi2_idx", pi2_idx);
	    cand.addUserInt("pi3_idx", pi3_idx);

      cand.addUserFloat("mu1_pvjpsi_dxy", mu1_pvjpsi_dxy);
      cand.addUserFloat("mu1_pvjpsi_dz", mu1_pvjpsi_dz);
      cand.addUserFloat("mu2_pvjpsi_dxy", mu2_pvjpsi_dxy);
      cand.addUserFloat("mu2_pvjpsi_dz", mu2_pvjpsi_dz);
      cand.addUserFloat("pi1_dxy", pi1_dxy);
      cand.addUserFloat("pi1_dz", pi1_dz);
      cand.addUserFloat("pi2_dxy", pi2_dxy);
      cand.addUserFloat("pi2_dz", pi2_dz);
      cand.addUserFloat("pi3_dxy", pi3_dxy);
      cand.addUserFloat("pi3_dz", pi3_dz);

      cand.addUserFloat("mu1_pvjpsi_dxyErr", mu1_pvjpsi_dxyErr);
      cand.addUserFloat("mu1_pvjpsi_dzErr", mu1_pvjpsi_dzErr);
      cand.addUserFloat("mu2_pvjpsi_dxyErr", mu2_pvjpsi_dxyErr);
      cand.addUserFloat("mu2_pvjpsi_dzErr", mu2_pvjpsi_dzErr);
      cand.addUserFloat("pi1_dxyErr", pi1_dxyErr);
      cand.addUserFloat("pi1_dzErr", pi1_dzErr);
      cand.addUserFloat("pi2_dxyErr", pi2_dxyErr);
      cand.addUserFloat("pi2_dzErr", pi2_dzErr);
      cand.addUserFloat("pi3_dxyErr", pi3_dxyErr);
      cand.addUserFloat("pi3_dzErr", pi3_dzErr);

	    auto dr_info = min_max_dr({mu1_ptr, mu2_ptr, pi1_ptr, pi2_ptr, pi3_ptr});
	    
	    cand.addUserFloat("min_dr", dr_info.first);
	    cand.addUserFloat("max_dr", dr_info.second);
	    // TODO add meaningful variables
	    
	    if( !pre_vtx_selection_(cand) ) continue;
	    //std::cout << "here2" << std::endl;
	    
	    //        std::cout<<"PRIMA"<<std::endl;
	    KinVtxFitter fitter(
				{muons_ttracks->at(mu1_idx), muons_ttracks->at(mu2_idx), particles_ttracks->at(pi1_idx),particles_ttracks->at(pi2_idx),particles_ttracks->at(pi3_idx)},
				{mu1_ptr->mass(), mu2_ptr->mass(), PI_MASS, PI_MASS,PI_MASS},
				{LEP_SIGMA, LEP_SIGMA, PI_SIGMA,PI_SIGMA,PI_SIGMA,} //some small sigma for the muon mass
				);
	    //std::cout<<"DOPO"<<std::endl;
	    if(!fitter.success()) continue; // hardcoded, but do we need otherwise?
	    //std::cout << "here3" << std::endl;
	    cand.setVertex( 
			   reco::Candidate::Point( 
						  fitter.fitted_vtx().x(),
						  fitter.fitted_vtx().y(),
						  fitter.fitted_vtx().z()
						   )  
			    );
	    /*Measurement1D ip3D = getIP(fitter, particles_ttracks->at(k_idx));
	    cand.addUserFloat("ip3D", ip3D.value());
	    cand.addUserFloat("ip3D_e", ip3D.error());
	    */
	    // pv shortest dz from jpsi + mu
	    const reco::TransientTrack& threemuonTT = fitter.fitted_candidate_ttrk();
	    
	    used_muon1_id.emplace_back(mu1_idx);
	    used_muon2_id.emplace_back(mu2_idx);
	    used_pi1_id.emplace_back(pi1_idx);
	    used_pi2_id.emplace_back(pi2_idx);
	    used_pi3_id.emplace_back(pi3_idx);
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
	    auto lxy = l_xy(fitter, *beamspot);
	    cand.addUserFloat("l_xy", lxy.value());
	    cand.addUserFloat("l_xy_unc", lxy.error());
	    cand.addUserFloat("vtx_x", cand.vx());
	    cand.addUserFloat("vtx_y", cand.vy());
	    cand.addUserFloat("vtx_z", cand.vz());
	    cand.addUserFloat("vtx_ex", sqrt(fitter.fitted_vtx_uncertainty().cxx()));
	    cand.addUserFloat("vtx_ey", sqrt(fitter.fitted_vtx_uncertainty().cyy()));
	    cand.addUserFloat("vtx_ez", sqrt(fitter.fitted_vtx_uncertainty().czz()));
	    cand.addUserFloat("vtx_chi2", ChiSquaredProbability(fitter.chi2(), fitter.dof()));

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
	    


	    cand.addUserFloat("fitted_mu1_pt" , fitter.daughter_p4(0).pt()); 
	    cand.addUserFloat("fitted_mu1_eta", fitter.daughter_p4(0).eta());
	    cand.addUserFloat("fitted_mu1_phi", fitter.daughter_p4(0).phi());
	    cand.addUserFloat("fitted_mu2_pt" , fitter.daughter_p4(1).pt()); 
	    cand.addUserFloat("fitted_mu2_eta", fitter.daughter_p4(1).eta());
	    cand.addUserFloat("fitted_mu2_phi", fitter.daughter_p4(1).phi());
	    cand.addUserFloat("fitted_pi1_pt"  , fitter.daughter_p4(2).pt()); 
	    cand.addUserFloat("fitted_pi1_eta" , fitter.daughter_p4(2).eta());
	    cand.addUserFloat("fitted_pi1_phi" , fitter.daughter_p4(2).phi());
	    cand.addUserFloat("fitted_pi2_eta" , fitter.daughter_p4(3).eta());
	    cand.addUserFloat("fitted_pi2_phi" , fitter.daughter_p4(3).phi());
	    cand.addUserFloat("fitted_pi2_pt" , fitter.daughter_p4(3).pt());
	    cand.addUserFloat("fitted_pi3_pt"  , fitter.daughter_p4(4).pt());
	    cand.addUserFloat("fitted_pi3_eta" , fitter.daughter_p4(4).eta());
	    cand.addUserFloat("fitted_pi3_phi" , fitter.daughter_p4(4).phi());
	    
	    
	    const reco::BeamSpot &bm = *beamspot;

	    cand.addUserFloat("beamspot_x", bm.x0());
	    cand.addUserFloat("beamspot_y", bm.y0());
	    cand.addUserFloat("beamspot_z", bm.z0());
	    //	    float B_pt=(Bc_MASS/cand.mass())*cand.pt();
	    TLorentzVector P_b;
	    P_b.SetPtEtaPhiM(cand.pt(),cand.eta(),cand.phi(),cand.mass());
	    TLorentzVector P_pi1;
	    P_pi1.SetPtEtaPhiM(pi1_ptr->pt(),pi1_ptr->eta(),pi1_ptr->phi(),pi1_ptr->mass());

	    TLorentzVector P_pi2;
	    P_pi2.SetPtEtaPhiM(pi2_ptr->pt(),pi2_ptr->eta(),pi2_ptr->phi(),pi2_ptr->mass());

	    TLorentzVector P_pi3;
	    P_pi3.SetPtEtaPhiM(pi3_ptr->pt(),pi3_ptr->eta(),pi3_ptr->phi(),pi3_ptr->mass());
    
	    TLorentzVector P_mu1;
	    P_mu1.SetPtEtaPhiM(mu1_ptr->pt(),mu1_ptr->eta(),mu1_ptr->phi(),mu1_ptr->mass());
	    
	    TLorentzVector P_mu2;
	    P_mu2.SetPtEtaPhiM(mu2_ptr->pt(),mu2_ptr->eta(),mu2_ptr->phi(),mu2_ptr->mass());
      
	    float m_miss_2=(P_b-P_pi1 - P_pi2 - P_pi3 -P_mu1-P_mu2)*(P_b-P_pi1 - P_pi2 - P_pi3-P_mu1-P_mu2);
	    float Q_2=(P_b-P_mu1-P_mu2)*(P_b-P_mu1-P_mu2);
	    float pt_miss=(P_b.Pt()-P_pi1.Pt() - P_pi2.Pt() - P_pi3.Pt()-P_mu1.Pt()-P_mu2.Pt());
	    float pt_miss_vec=((P_b-P_pi1 - P_pi2 - P_pi3-P_mu1-P_mu2).Pt());
	    float pt_var=((P_mu1+P_mu2).Pt()-(P_pi1 - P_pi2 - P_pi3).Pt());
	    float DR=deltaR(P_mu1.Eta(),P_mu1.Phi(),P_mu2.Eta(),P_mu2.Phi());

	    float m_jpsi=sqrt((P_mu1+P_mu2)*(P_mu1+P_mu2));
	    cand.addUserFloat("m_miss_2", m_miss_2);
	    cand.addUserFloat("Q_2",Q_2);
	    cand.addUserFloat("pt_miss",pt_miss);
	    cand.addUserFloat("pt_miss_vec",pt_miss_vec);
	    cand.addUserFloat("pt_var",pt_var);
	    cand.addUserFloat("DR",DR);
	    cand.addUserFloat("m_jpsi",m_jpsi);
	    
	    //energia del mu unpaired in diversi sistemi di riferimento                  
	    /*
	    TLorentzVector P_mu=P_k;        
	    TVector3 mu_beta_lab=P_b.BoostVector();
	    
	    P_mu.Boost(-mu_beta_lab);
	    cand.addUserFloat("E_mu_star",P_mu.E());
	    P_mu=P_k;        
	    TLorentzVector jpsi=P_mu1+P_mu2;
	    TVector3 jpsi_beta_lab=jpsi.BoostVector();
	    P_mu.Boost(-jpsi_beta_lab);
	    cand.addUserFloat("E_mu_#",P_mu.E());
	    */
	    if( !post_vtx_selection_(cand) ) {
	      if(debug) std::cout<<"post vrxt dies "<<std::endl;
	      continue;        
	    }
	    //std::cout << "here4" << std::endl;
	    if(debug) std::cout<<"post vrxt Survives!! "<<std::endl;
	    //compute isolation
	    float mu1_iso03 = 0;
	    float mu1_iso04 = 0;
	    float mu2_iso03 = 0;
	    float mu2_iso04 = 0;
	    float pi1_iso03  = 0;
	    float pi1_iso04  = 0;
	    float pi2_iso03  = 0;
	    float pi2_iso04  = 0;
	    float pi3_iso03  = 0;
	    float pi3_iso04  = 0;
	    float b_iso03  = 0;
	    float b_iso04  = 0;
      
	    for( unsigned int iTrk=0; iTrk<totalTracks; ++iTrk ) 
	      {
		const pat::PackedCandidate & trk = (iTrk < nTracks) ? (*iso_tracks)[iTrk] : (*iso_lostTracks)[iTrk-nTracks];
		// define selections for iso tracks (pT, eta, ...)
		if( !isotrk_selection_(trk) ) continue;


		// only consider tracks originating close to the three bodies
		if ( !mu1_ptr->bestTrack() || fabs(trk.dz() - mu1_ptr->bestTrack()->dz()) > 0.4 ) continue;
		if ( !mu2_ptr->bestTrack() || fabs(trk.dz() - mu2_ptr->bestTrack()->dz()) > 0.4 ) continue;
		if ( fabs(trk.dz() - particles_ttracks->at(pi1_idx).track().dz()) > 0.4 ) continue;
		if ( fabs(trk.dz() - particles_ttracks->at(pi2_idx).track().dz()) > 0.4 ) continue;
		if (  fabs(trk.dz() - particles_ttracks->at(pi3_idx).track().dz()) > 0.4 ) continue;

		if(track_to_muon_match(pi1_ptr, iso_tracks.id(), iTrk)  )  continue;
		if(track_to_muon_match(pi2_ptr, iso_tracks.id(), iTrk)  )  continue;
		if(track_to_muon_match(pi3_ptr, iso_tracks.id(), iTrk)  )  continue;

		// check if the track is one of the two muons
		if (track_to_muon_match(mu1_ptr, iso_tracks.id(), iTrk) ||
		    track_to_muon_match(mu2_ptr, iso_tracks.id(), iTrk) ) continue;
		
		// add to final particle iso if dR < cone
		float dr_to_mu1 = deltaR(cand.userFloat("fitted_mu1_eta"), cand.userFloat("fitted_mu1_phi"), trk.eta(), trk.phi());
		float dr_to_mu2 = deltaR(cand.userFloat("fitted_mu2_eta"), cand.userFloat("fitted_mu2_phi"), trk.eta(), trk.phi());
		float dr_to_pi1  = deltaR(cand.userFloat("fitted_pi1_eta") , cand.userFloat("fitted_pi1_phi") , trk.eta(), trk.phi());
		float dr_to_pi2  = deltaR(cand.userFloat("fitted_pi2_eta") , cand.userFloat("fitted_pi2_phi") , trk.eta(), trk.phi());
		float dr_to_pi3  = deltaR(cand.userFloat("fitted_pi3_eta") , cand.userFloat("fitted_pi3_phi") , trk.eta(), trk.phi());

		float dr_to_b  = deltaR(cand.userFloat("fitted_eta")   , cand.userFloat("fitted_phi") , trk.eta(), trk.phi());
        
		if (dr_to_mu1 < 0.4 && dr_to_mu1>0.01){
		  mu1_iso04 += trk.pt();
		  if ( dr_to_mu1 < 0.3) mu1_iso03 += trk.pt();
		}
		if (dr_to_mu2 < 0.4 && dr_to_mu2>0.01){
		  mu2_iso04 += trk.pt();
		  if (dr_to_mu2 < 0.3)  mu2_iso03 += trk.pt();
		}
		if (dr_to_pi1 < 0.4 && dr_to_pi1>0.01){
		  pi1_iso04 += trk.pt();
		  if (dr_to_pi1 < 0.3) pi1_iso03 += trk.pt();
		}
		if (dr_to_pi2 < 0.4 && dr_to_pi2>0.01){
		  pi2_iso04 += trk.pt();
		  if (dr_to_pi2 < 0.3) pi2_iso03 += trk.pt();
		}
		if (dr_to_pi3 < 0.4 && dr_to_pi3>0.01){
		  pi3_iso04 += trk.pt();
		  if (dr_to_pi3 < 0.3) pi3_iso03 += trk.pt();
		}
		if (dr_to_b < 0.4){
		  b_iso04 += trk.pt();
		  if (dr_to_b < 0.3) b_iso03 += trk.pt();
		}      
	      }

	    cand.addUserFloat("mu1_iso03", mu1_iso03);
	    cand.addUserFloat("mu1_iso04", mu1_iso04);
	    cand.addUserFloat("mu2_iso03", mu2_iso03);
	    cand.addUserFloat("mu2_iso04", mu2_iso04);
	    cand.addUserFloat("pi1_iso03" , pi1_iso03 );
	    cand.addUserFloat("pi1_iso04" , pi1_iso04 );
	    cand.addUserFloat("pi2_iso03" , pi2_iso03 );
	    cand.addUserFloat("pi2_iso04" , pi2_iso04 );
	    cand.addUserFloat("pi3_iso03" , pi3_iso03 );
	    cand.addUserFloat("pi3_iso04" , pi3_iso04 );      
	    cand.addUserFloat("b_iso03" , b_iso03 );
	    cand.addUserFloat("b_iso04" , b_iso04 );

	    //save the candidate
	    ret_val->push_back(cand);
	  }//for(size_t ll_idx = 0; ll_idx < dimuons->size(); ++ll_idx) 
      }//loop pion3
    }//loop pion2
  }//loop pion1
  for (auto & cand: *ret_val){
    cand.addUserInt("n_pi1_used", std::count(used_pi1_id.begin(),used_pi1_id.end(),cand.userInt("pi1_idx")));
    cand.addUserInt("n_pi2_used", std::count(used_pi2_id.begin(),used_pi2_id.end(),cand.userInt("pi2_idx")));
    cand.addUserInt("n_pi3_used", std::count(used_pi3_id.begin(),used_pi3_id.end(),cand.userInt("pi3_idx")));
    cand.addUserInt("n_mu1_used", std::count(used_muon1_id.begin(),used_muon1_id.end(),cand.userInt("mu1_idx"))+std::count(used_muon2_id.begin(),used_muon2_id.end(),cand.userInt("mu1_idx")));
    cand.addUserInt("n_mu2_used", std::count(used_muon1_id.begin(),used_muon1_id.end(),cand.userInt("mu2_idx"))+std::count(used_muon2_id.begin(),used_muon2_id.end(),cand.userInt("mu2_idx")));
  }
  evt.put(std::move(ret_val));
}//produce  

Measurement1D BTo2Mu3PiBuilder::getIP(KinVtxFitter fitter,reco::TransientTrack transientTrackMu) const
{
  // computing 3d impact parameter
  GlobalVector jpsiDirection(fitter.fitted_p4().x(), fitter.fitted_p4().y(), fitter.fitted_p4().z());
  reco::Vertex::Point jpsiVertexPosition(fitter.fitted_vtx().x(),fitter.fitted_vtx().y(),fitter.fitted_vtx().z());

  const double err00 = fitter.fitted_candidate().kinematicParametersError().matrix()(0,0);
  const double err11 = fitter.fitted_candidate().kinematicParametersError().matrix()(1,1);
  const double err22 = fitter.fitted_candidate().kinematicParametersError().matrix()(2,2);
  const double err01 = fitter.fitted_candidate().kinematicParametersError().matrix()(0,1);
  const double err02 = fitter.fitted_candidate().kinematicParametersError().matrix()(0,2);
  const double err12 = fitter.fitted_candidate().kinematicParametersError().matrix()(1,2);
  reco::Vertex::Error jpsiVertexError;

  jpsiVertexError(0,0) = err00;
  jpsiVertexError(0,1) = err01;
  jpsiVertexError(0,2) = err02;
  jpsiVertexError(1,0) = err01;
  jpsiVertexError(1,1) = err11;
  jpsiVertexError(1,2) = err12;
  jpsiVertexError(2,0) = err02;
  jpsiVertexError(2,1) = err12;
  jpsiVertexError(2,2) = err22;

  GlobalVector jpsiGlobalVector(fitter.fitted_p4().x(), fitter.fitted_p4().y(), fitter.fitted_p4().z());
  const reco::Vertex jpsiVertex(jpsiVertexPosition, jpsiVertexError, fitter.chi2(), fitter.dof(), 2);

  SignedImpactParameter3D signed_ip3D;
  Measurement1D ip3D = signed_ip3D.apply(transientTrackMu,jpsiGlobalVector,jpsiVertex).second;
  return ip3D;
}

int BTo2Mu3PiBuilder::getPVIdx(const reco::VertexCollection* vertices,const reco::TransientTrack& dimuonTT) const
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
DEFINE_FWK_MODULE(BTo2Mu3PiBuilder);

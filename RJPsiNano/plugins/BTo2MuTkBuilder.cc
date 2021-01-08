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


class BTo2MuTkBuilder : public edm::global::EDProducer<> {

  // perhaps we need better structure here (begin run etc)
public:
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;
  typedef std::vector<reco::GenParticle> GenParticleCollection;

  explicit BTo2MuTkBuilder(const edm::ParameterSet &cfg):
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
    particle_mass{cfg.getParameter<double>("particle_mass")},
    beamspot_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") )} 
  {
      produces<pat::CompositeCandidateCollection>();
  }

  ~BTo2MuTkBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
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
  const double particle_mass;


  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;  
};

void BTo2MuTkBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

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

  std::vector<int> used_muon1_id, used_muon2_id, used_trk_id;

  std::unique_ptr<pat::CompositeCandidateCollection> ret_val(new pat::CompositeCandidateCollection());

  const reco::VertexCollection* vertices = primaryVertices.product();
  // output

  if(debug) std::cout<<"dimuons->size() "<<dimuons->size()<<std::endl;
  for(size_t ll_idx = 0; ll_idx < dimuons->size(); ++ll_idx) 
  {
    //std::cout << "PV " << ll_idx << ": " << vertices->at(ll_idx).position() << std::endl;
    edm::Ptr<pat::CompositeCandidate> ll_ptr(dimuons, ll_idx);
    edm::Ptr<reco::Candidate> mu1_ptr = ll_ptr->userCand("mu1");
    edm::Ptr<reco::Candidate> mu2_ptr = ll_ptr->userCand("mu2");
    size_t mu1_idx = abs(ll_ptr->userInt("mu1_idx"));
    size_t mu2_idx = abs(ll_ptr->userInt("mu2_idx"));

    int pvIdx = ll_ptr->userInt("pvIdx");
    reco::Vertex bestVertex = vertices->at(pvIdx);
    double mu1_dxy = mu1_ptr->bestTrack()->dxy(bestVertex.position());
    double mu2_dxy = mu2_ptr->bestTrack()->dxy(bestVertex.position());
    double mu1_dz = mu1_ptr->bestTrack()->dz(bestVertex.position());
    double mu2_dz = mu2_ptr->bestTrack()->dz(bestVertex.position());

    size_t isDimuon_dimuon0Trg = abs(ll_ptr->userInt("muonpair_fromdimuon0"));
    size_t isDimuon_jpsiTrkTrg = abs(ll_ptr->userInt("muonpair_fromjpsitrk"));
    //size_t isDimuon_jpsiTrkTrg = abs(ll_ptr->userInt("isJpsiTrkTrg"));
    //size_t isDimuon_dimuon0Trg = abs(ll_ptr->userInt("isDimuon0Trg"));
    if(debug) std::cout<<"isDimuon_jpsiTrkTrg  "<<isDimuon_jpsiTrkTrg<<std::endl;
    if(!(isDimuon_jpsiTrkTrg)) continue;

    //Loop  on displaced muons    
    if(debug) std::cout<<"paerticles size "<<particles->size()<<std::endl;
    for(size_t k_idx = 0; k_idx < particles->size(); ++k_idx) {
      edm::Ptr<pat::CompositeCandidate> k_ptr(particles, k_idx);
      if( !particle_selection_(*k_ptr) ) continue;

      //double k_dxy = k_ptr->bestTrack()->dxy(bestVertex.position());
      //double k_dz = k_ptr->bestTrack()->dz(bestVertex.position());
      double k_dxy = particles_ttracks->at(k_idx).track().dxy(bestVertex.position());
      double k_dz = particles_ttracks->at(k_idx).track().dz(bestVertex.position());

      bool isPartTrg = k_ptr->userInt("isTriggering");
      if(debug) std::cout<<"isTriggering "<<isPartTrg<<std::endl;
      //ha trovato il mu displaced
      if(!(isPartTrg)) {
	//if(debug) std::cout<<"is NOT track triggered "<<k_ptr->pt()<<std::endl;
	continue;
      }
      

      math::PtEtaPhiMLorentzVector k_p4(
                k_ptr->pt(), 
                k_ptr->eta(),
                k_ptr->phi(),
                particle_mass 
                );
  
      // Use UserCands as they should not use memory but keep the Ptr itself
      // Put the muon passing the corresponding selection

      pat::CompositeCandidate cand;
      cand.setP4(ll_ptr->p4() + k_p4);
      cand.setCharge(ll_ptr->charge() + k_ptr->charge());

      // pv info

      cand.addUserInt("pv_idx", pvIdx);

      // tracks info
      if(debug) std::cout<<"cand pt "<<cand.pt()<<std::endl;
      if(debug) std::cout<<"displ mu "<<k_ptr->pt()<<std::endl;
      if(debug) std::cout<<"displ m1 "<<mu1_ptr->pt()<<std::endl;
      if(debug) std::cout<<"displ m2 "<<mu2_ptr->pt()<<std::endl;

      cand.addUserCand("mu1", mu1_ptr);
      cand.addUserCand("mu2", mu2_ptr);
      cand.addUserCand("k", k_ptr);
      cand.addUserCand("dimuon", ll_ptr);
      
      cand.addUserInt("mu1_idx", mu1_idx);
      cand.addUserInt("mu2_idx", mu2_idx);
      cand.addUserInt("k_idx", k_idx);

      cand.addUserFloat("mu1_dxy", mu1_dxy);
      cand.addUserFloat("mu1_dz", mu1_dz);
      cand.addUserFloat("mu2_dxy", mu2_dxy);
      cand.addUserFloat("mu2_dz", mu2_dz);
      cand.addUserFloat("k_dxy", k_dxy);
      cand.addUserFloat("k_dz", k_dz);
      
      auto dr_info = min_max_dr({mu1_ptr, mu2_ptr, k_ptr});

      cand.addUserFloat("min_dr", dr_info.first);
      cand.addUserFloat("max_dr", dr_info.second);
      // TODO add meaningful variables
      
      if( !pre_vtx_selection_(cand) ) continue;
      //std::cout << "here2" << std::endl;
      
      //        std::cout<<"PRIMA"<<std::endl;
      KinVtxFitter fitter(
        {muons_ttracks->at(mu1_idx), muons_ttracks->at(mu2_idx), particles_ttracks->at(k_idx)},
        {mu1_ptr->mass(), mu2_ptr->mass(), particle_mass},
        {LEP_SIGMA, LEP_SIGMA, LEP_SIGMA} //some small sigma for the muon mass
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
      Measurement1D ip3D = getIP(fitter, particles_ttracks->at(k_idx));
      cand.addUserFloat("ip3D", ip3D.value());
      cand.addUserFloat("ip3D_e", ip3D.error());

      used_muon1_id.emplace_back(mu1_idx);
      used_muon2_id.emplace_back(mu2_idx);
      used_trk_id.emplace_back(k_idx);
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

      cand.addUserFloat("pv_x", bestVertex.position().x());
      cand.addUserFloat("pv_y", bestVertex.position().y());
      cand.addUserFloat("pv_z", bestVertex.position().z());
      cand.addUserFloat("pv_ex", bestVertex.covariance(0,0));
      cand.addUserFloat("pv_ey", bestVertex.covariance(1,1));
      cand.addUserFloat("pv_ez", bestVertex.covariance(2,2));
      cand.addUserFloat("pv_exy", bestVertex.covariance(0,1));
      cand.addUserFloat("pv_eyz", bestVertex.covariance(0,2));
      cand.addUserFloat("pv_exz", bestVertex.covariance(1,2));
      cand.addUserFloat("pv_chi2", ChiSquaredProbability(bestVertex.chi2(), bestVertex.ndof()));
      

      cand.addUserFloat("fitted_mu1_pt" , fitter.daughter_p4(0).pt()); 
      cand.addUserFloat("fitted_mu1_eta", fitter.daughter_p4(0).eta());
      cand.addUserFloat("fitted_mu1_phi", fitter.daughter_p4(0).phi());
      cand.addUserFloat("fitted_mu2_pt" , fitter.daughter_p4(1).pt()); 
      cand.addUserFloat("fitted_mu2_eta", fitter.daughter_p4(1).eta());
      cand.addUserFloat("fitted_mu2_phi", fitter.daughter_p4(1).phi());
      cand.addUserFloat("fitted_k_pt"  , fitter.daughter_p4(2).pt()); 
      cand.addUserFloat("fitted_k_eta" , fitter.daughter_p4(2).eta());
      cand.addUserFloat("fitted_k_phi" , fitter.daughter_p4(2).phi());


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
	if ( fabs(trk.dz() - particles_ttracks->at(k_idx).track().dz()) > 0.4 ) continue;

        if (k_ptr->userCand("cand") ==  edm::Ptr<reco::Candidate> ( iso_tracks, iTrk ) ) {
          
          std::cout<<"old"<<std::endl;
          continue;
        }
        
        if(track_to_muon_match(k_ptr, iso_tracks.id(), iTrk)) 
        {
          // std::cout<<"new"<<std::endl;
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
	}
      }

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

Measurement1D BTo2MuTkBuilder::getIP(KinVtxFitter fitter,reco::TransientTrack transientTrackMu) const
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
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(BTo2MuTkBuilder);

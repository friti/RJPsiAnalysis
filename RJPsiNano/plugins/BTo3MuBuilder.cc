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


class BTo3MuBuilder : public edm::global::EDProducer<> {

  // perhaps we need better structure here (begin run etc)
public:
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;
  typedef std::vector<reco::GenParticle> GenParticleCollection;

  explicit BTo3MuBuilder(const edm::ParameterSet &cfg):
    mu3_selection_{cfg.getParameter<std::string>("muonSelection")},
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    dimuons_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("dimuons") )},
    pvSelected_{consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("pvSelected"))},
    muons_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("muonTransientTracks") )},
    muons_{consumes<pat::MuonCollection>( cfg.getParameter<edm::InputTag>("muons") )},
    //kaons_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("kaonsTransientTracks") )},
    isotracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
    isolostTracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))),

    isotrk_selection_{cfg.getParameter<std::string>("isoTracksSelection")},
   //GEN
    srcToken_(consumes<GenParticleCollection>(cfg.getParameter<edm::InputTag>("srcGen"))), 
    beamspot_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") )} 
  {
      produces<pat::CompositeCandidateCollection>();
  }

  ~BTo3MuBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
  Measurement1D getIP(KinVtxFitter fitter, reco::TransientTrack transientTrackMu) const;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
  private:
  const StringCutObjectSelector<pat::Muon> mu3_selection_; 
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; // cut on the di-muon before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; // cut on the di-muon after the SV fit

  const edm::EDGetTokenT<pat::CompositeCandidateCollection> dimuons_;
  const edm::EDGetTokenT<reco::VertexCollection> pvSelected_;
  const edm::EDGetTokenT<TransientTrackCollection> muons_ttracks_;
  const edm::EDGetTokenT<pat::MuonCollection> muons_;
  //const edm::EDGetTokenT<TransientTrackCollection> muons_ttracks_;

  const edm::EDGetTokenT<pat::PackedCandidateCollection> isotracksToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> isolostTracksToken_;

  const StringCutObjectSelector<pat::PackedCandidate> isotrk_selection_; 

  //GEN
  edm::EDGetTokenT<reco::GenParticleCollection> srcToken_;

  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;  
};

void BTo3MuBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  //input
  edm::Handle<pat::CompositeCandidateCollection> dimuons;
  evt.getByToken(dimuons_, dimuons);

  edm::Handle<reco::VertexCollection> pvSelected;
  evt.getByToken(pvSelected_, pvSelected);
  
  edm::Handle<TransientTrackCollection> muons_ttracks;
  evt.getByToken(muons_ttracks_, muons_ttracks);

  edm::Handle<pat::MuonCollection> muons;
  evt.getByToken(muons_, muons);
  
  //edm::Handle<TransientTrackCollection> kaons_ttracks;
  //evt.getByToken(kaons_ttracks_, kaons_ttracks);  

  edm::Handle<reco::BeamSpot> beamspot;
  evt.getByToken(beamspot_, beamspot);  

  //for isolation
  edm::Handle<pat::PackedCandidateCollection> iso_tracks;
  evt.getByToken(isotracksToken_, iso_tracks);
  edm::Handle<pat::PackedCandidateCollection> iso_lostTracks;
  evt.getByToken(isolostTracksToken_, iso_lostTracks);

  //GEN
  int flag_jpsi_mu = -99. ;
  int flag_psi2s_mu = -99.;
  int flag_chic0_mu = -99.;
  int flag_chic1_mu = -99.;
  int flag_chic2_mu = -99.;
  int flag_hc_mu = -99.;
  int flag_jpsi_tau = -99.;
  int flag_psi2s_tau = -99.;
  int flag_jpsi_pi = -99.;
  int flag_jpsi_3pi = -99.;
  int flag_jpsi_hc = -99.;
  int flag_error = -99.;
  
  if(evt.eventAuxiliary().run() == 1){ //only if it is a MC sample
    //GEN
    if(debugGen) std::cout<<"In the flag code..."<<std::endl;
    edm::Handle<GenParticleCollection> src;
    evt.getByToken(srcToken_, src);
    const size_t n = src->size();
    std::vector<int> final_daus;

    if(debugGen) std::cout<<"Number of gen particles: "<<n<<std::endl;
    
    for(unsigned int  i = 0; i < n; ++i) {  //loop on gen particles
      const reco::GenParticle & gen = (*src)[i];
      const reco::Candidate* daughter; 
      const reco::Candidate* the_b;
      int is_doublemu = 0;
      int is_b = 0;
      final_daus.clear();
    
      if(abs(gen.pdgId()) == 443){  // looking for jpsi      

      if(debugGen) std::cout<<"There is a jpsi and she has "<<gen.numberOfDaughters()<<" daughters"<<std::endl;
      for(unsigned int dau = 0; dau < gen.numberOfDaughters(); dau++){  //loop of jpsi daughters
        if(debugGen) std::cout<<"Jpsi daughter: "<<gen.daughter(dau)->pdgId()<<std::endl;
        is_b = 0;
        if (abs(gen.daughter(dau)->pdgId())==13){
          is_doublemu += 1;
        }
      } //end loop on daughters
      if(is_doublemu>=2){  // jpsi -> mu mu
        if(debugGen) std::cout<<"The daughters are muons"<<std::endl;
        the_b = gen.mother(0); // jpsi mother
        if(abs(the_b->pdgId()) == 541){ //Bc->jpsi
          if(debugGen) std::cout<<"The direct mother is a Bc"<<std::endl;
          is_b = 1;
        }  
        else if(the_b->numberOfMothers() > 0)
        {
          the_b = gen.mother(0)->mother(0); // Bc->X->jpsi
          if(abs(the_b->pdgId()) == 541 ){
            if(debugGen) std::cout<<"The non direct mother is a Bc"<<std::endl;
            is_b = 1;
          }
        }
        if(is_b == 1){
          if(debugGen) std::cout<<"The Bc has "<<the_b->numberOfDaughters()<<"daughters"<<std::endl;

          for(unsigned int bdau=0; bdau < the_b->numberOfDaughters(); bdau ++){
            daughter = the_b->daughter(bdau);
            if(abs(daughter->pdgId())!= 541 and abs(daughter->pdgId())!= 22){    //not gamma
              final_daus.push_back(abs(daughter->pdgId()));
              //        cout<<daughter->pdgId()<<endl;
              if(debugGen) std::cout<<"The Bc daughters are "<< daughter->pdgId()<<std::endl;
            }
          }
          
          std::sort(final_daus.begin(), final_daus.end());  //sort the pdgIds of the daughters
          /*
          for(unsigned int item=0; item< final_daus.size(); item ++){
            if(debugGen) std::cout<<final_daus[item]<<std::endl;
            if(item == final_daus.size() -1) if(debugGen) std::cout<<" "<<std::endl;
          }
          */

          flag_jpsi_mu = 0;
          flag_psi2s_mu = 0;
          flag_chic0_mu = 0;
          flag_chic1_mu = 0;
          flag_chic2_mu = 0;
          flag_hc_mu = 0;
          flag_jpsi_tau = 0;
          flag_psi2s_tau = 0;
          flag_jpsi_pi = 0;
          flag_jpsi_3pi = 0;
          flag_jpsi_hc = 0;
          flag_error = 0;
          
          if(final_daus[0] == 13){  //muon
            if(final_daus[1] == 14){
              if(final_daus[2] == 443)  flag_jpsi_mu=1;
              else if (final_daus[2] == 100443) flag_psi2s_mu = 1;
              else if (final_daus[2] == 10441) flag_chic0_mu = 1;
              else if (final_daus[2] == 20443) flag_chic1_mu = 1;
              else if (final_daus[2] == 445) flag_chic2_mu = 1;
              else if (final_daus[2] == 10443) flag_hc_mu = 1;
            }
          }
          else if(final_daus[0] == 15){ //tau
            if (final_daus[1] == 16){
              if(final_daus[2] == 443) flag_jpsi_tau = 1;
              else if(final_daus[2] == 100443) flag_psi2s_tau = 1;
            }
          }
          else if(final_daus[0] == 211){
            if (final_daus[1] == 443) flag_jpsi_pi = 1;
            if (final_daus[1] == 211 && final_daus[2] ==211 && final_daus[1] == 443) flag_jpsi_3pi = 1;
          }
          else if ((final_daus[0] == 431 && final_daus[1] == 443) || (final_daus[0] == 433 && final_daus[1] == 443)) flag_jpsi_hc = 1;  
          else{
            flag_error = 1;
          }
        } // if(is_b == 1)
      } //if(is_doublemu>=2)
    }//if(abs(gen.pdgId()) == 443)
  }//for(unsigned int  i = 0; i < n; ++i)
}//if(evt.eventAuxiliary().run() == 1)

//////

  unsigned int nTracks     = iso_tracks->size();
  unsigned int totalTracks = nTracks + iso_lostTracks->size();

  std::vector<int> used_muon1_id, used_muon2_id, used_trk_id;

  std::unique_ptr<pat::CompositeCandidateCollection> ret_val(new pat::CompositeCandidateCollection());

  const reco::VertexCollection* vertices = pvSelected.product();
  // output
  for(size_t ll_idx = 0; ll_idx < dimuons->size(); ++ll_idx) 
  {
    //std::cout << "PV " << ll_idx << ": " << vertices->at(ll_idx).position() << std::endl;
    reco::Vertex bestVertex = vertices->at(ll_idx);
    edm::Ptr<pat::CompositeCandidate> ll_prt(dimuons, ll_idx);
    edm::Ptr<reco::Candidate> mu1_ptr = ll_prt->userCand("mu1");
    edm::Ptr<reco::Candidate> mu2_ptr = ll_prt->userCand("mu2");
    size_t mu1_idx = abs(ll_prt->userInt("mu1_idx"));
    size_t mu2_idx = abs(ll_prt->userInt("mu2_idx"));
    size_t isDimuon_dimuon0Trg = abs(ll_prt->userInt("isDimuon0Trg"));
    size_t isDimuon_jpsiTrkTrg = abs(ll_prt->userInt("isJpsiTrkTrg"));


    //Loop  on displaced muons    
    for(size_t mu3_idx = 0; mu3_idx < muons->size(); ++mu3_idx) {
      edm::Ptr<pat::Muon> mu3_ptr(muons, mu3_idx);
      if((mu1_idx == mu3_idx) || (mu2_idx == mu3_idx)) continue;
      if( !mu3_selection_(*mu3_ptr) ) continue;
  
      //ha trovato il mu displaced
      bool isDimuon0Trg = mu3_ptr->userInt("isDimuon0Trg");
      bool isJpsiTrkTrg = mu3_ptr->userInt("isJpsiTrkTrg");
      bool isJpsiMuon = mu3_ptr->userInt("isJpsiMuon");
      bool isUnpairedMuon = isDimuon0Trg && !isJpsiMuon;
      if(!(isDimuon_jpsiTrkTrg || isUnpairedMuon)) continue;
      
      math::PtEtaPhiMLorentzVector mu3_p4(
                mu3_ptr->pt(), 
                mu3_ptr->eta(),
                mu3_ptr->phi(),
                mu3_ptr->mass()
                );
  
      // Use UserCands as they should not use memory but keep the Ptr itself
      // Put the muon passing the corresponding selection

      pat::CompositeCandidate cand;
      cand.setP4(ll_prt->p4() + mu3_p4);
      cand.setCharge(ll_prt->charge() + mu3_ptr->charge());

      cand.addUserCand("mu1", mu1_ptr);
      cand.addUserCand("mu2", mu2_ptr);
      cand.addUserCand("mu3", mu3_ptr);
      cand.addUserCand("dimuon", ll_prt);
      
      cand.addUserInt("mu1_idx", mu1_idx);
      cand.addUserInt("mu2_idx", mu2_idx);
      cand.addUserInt("mu3_idx", mu3_idx);
      
      auto dr_info = min_max_dr({mu1_ptr, mu2_ptr, mu3_ptr});

      cand.addUserFloat("min_dr", dr_info.first);
      cand.addUserFloat("max_dr", dr_info.second);
      // TODO add meaningful variables
      //GEN variables
      cand.addUserInt("is_jpsi_mu", flag_jpsi_mu);
      cand.addUserInt("is_psi2s_mu", flag_psi2s_mu);
      cand.addUserInt("is_chic0_mu", flag_chic0_mu);
      cand.addUserInt("is_chic1_mu", flag_chic1_mu);
      cand.addUserInt("is_chic2_mu", flag_chic2_mu);
      cand.addUserInt("is_hc_mu", flag_hc_mu);
      cand.addUserInt("is_jpsi_tau", flag_jpsi_tau);
      cand.addUserInt("is_psi2s_tau", flag_psi2s_tau);
      cand.addUserInt("is_jpsi_pi", flag_jpsi_pi);
      cand.addUserInt("is_jpsi_3pi", flag_jpsi_3pi);
      cand.addUserInt("is_jpsi_hc", flag_jpsi_hc);
      cand.addUserInt("is_error", flag_error);
      int weight;
      
      if(evt.eventAuxiliary().run() == 1){
        if(flag_jpsi_mu == 1) weight = 1.;
        else if(flag_psi2s_mu == 1) weight = 0.5474;
        else if(flag_chic0_mu == 1) weight = 0.0116;
        else if(flag_chic1_mu == 1) weight = 0.3440;
        else if(flag_chic2_mu == 1) weight = 0.1950;
        else if(flag_hc_mu == 1) weight = 0.01;
        else if(flag_jpsi_tau == 1) weight = 1.;
        else if(flag_psi2s_tau == 1) weight = 0.5474;
        else if(flag_jpsi_pi == 1) weight = 1.;
        else if(flag_jpsi_3pi == 1) weight = 1.;
        else if(flag_jpsi_hc == 1) weight = 1.;
        else weight = -1.;
      }
      else weight = -99.;
      cand.addUserFloat("weightGen", weight);
      
      if( !pre_vtx_selection_(cand) ) continue;
      //std::cout << "here2" << std::endl;
      
      //        std::cout<<"PRIMA"<<std::endl;
      KinVtxFitter fitter(
        {muons_ttracks->at(mu1_idx), muons_ttracks->at(mu2_idx), muons_ttracks->at(mu3_idx)},
        {mu1_ptr->mass(), mu2_ptr->mass(), mu3_ptr->mass()},
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
      Measurement1D ip3D = getIP(fitter, muons_ttracks->at(mu3_idx));
      cand.addUserFloat("ip3D", ip3D.value());
      cand.addUserFloat("ip3D_e", ip3D.error());

      used_muon1_id.emplace_back(mu1_idx);
      used_muon2_id.emplace_back(mu2_idx);
      used_trk_id.emplace_back(mu3_idx);
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

      cand.addUserFloat("jpsi_vtx_x", ll_prt->userFloat("vtx_x"));
      cand.addUserFloat("jpsi_vtx_y", ll_prt->userFloat("vtx_y"));
      cand.addUserFloat("jpsi_vtx_z", ll_prt->userFloat("vtx_z"));
      cand.addUserFloat("jpsi_vtx_ex", ll_prt->userFloat("vtx_ex"));
      cand.addUserFloat("jpsi_vtx_ey", ll_prt->userFloat("vtx_ey"));
      cand.addUserFloat("jpsi_vtx_ez", ll_prt->userFloat("vtx_ez"));
      cand.addUserFloat("jpsi_vtx_chi2", ll_prt->userFloat("vtx_chi2"));

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
      cand.addUserFloat("fitted_mu3_pt"  , fitter.daughter_p4(2).pt()); 
      cand.addUserFloat("fitted_mu3_eta" , fitter.daughter_p4(2).eta());
      cand.addUserFloat("fitted_mu3_phi" , fitter.daughter_p4(2).phi());
      
      //mie variabili                                                                                
      //conto quanti mu totali aveva l'evento
      //cand.addUserInt("pass_3mu",pass_3mu.size());                                 
      
      float B_pt=(Bc_MASS/cand.mass())*cand.pt();
      
      TLorentzVector P_b;
      P_b.SetPtEtaPhiM(B_pt,cand.eta(),cand.phi(),Bc_MASS);
      
      TLorentzVector P_k;
      P_k.SetPtEtaPhiM(mu3_ptr->pt(),mu3_ptr->eta(),mu3_ptr->phi(),mu3_ptr->mass());
      
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
      if( !post_vtx_selection_(cand) ) continue;        
      //std::cout << "here4" << std::endl;

      //compute isolation
      float mu1_iso03 = 0;
      float mu1_iso04 = 0;
      float mu2_iso03 = 0;
      float mu2_iso04 = 0;
      float mu3_iso03  = 0;
      float mu3_iso04  = 0;
      float b_iso03  = 0;
      float b_iso04  = 0;
      
      for( unsigned int iTrk=0; iTrk<totalTracks; ++iTrk ) 
      {
        const pat::PackedCandidate & trk = (iTrk < nTracks) ? (*iso_tracks)[iTrk] : (*iso_lostTracks)[iTrk-nTracks];
        // define selections for iso tracks (pT, eta, ...)
        if( !isotrk_selection_(trk) ) continue;
        // check if the track is the muon
        if (mu3_ptr->userCand("cand") ==  edm::Ptr<reco::Candidate> ( iso_tracks, iTrk ) ) {
          
          std::cout<<"old"<<std::endl;
          continue;
        }
        
        if(track_to_muon_match(mu3_ptr, iso_tracks.id(), iTrk)) 
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
        float dr_to_k  = deltaR(cand.userFloat("fitted_mu3_eta") , cand.userFloat("fitted_mu3_phi") , trk.eta(), trk.phi());
        float dr_to_b  = deltaR(cand.userFloat("fitted_eta")   , cand.userFloat("fitted_phi") , trk.eta(), trk.phi());
        
        if (dr_to_mu1 < 0.4)
        {
          mu1_iso04 += trk.pt();
          if ( dr_to_mu1 < 0.3) mu1_iso03 += trk.pt();
        }
        if (dr_to_mu2 < 0.4)
        {
          mu2_iso04 += trk.pt();
          if (dr_to_mu2 < 0.3)  mu2_iso03 += trk.pt();
        }
        if (dr_to_k < 0.4)
        {
          mu3_iso04 += trk.pt();
          if (dr_to_k < 0.3) mu3_iso03 += trk.pt();
        }
        if (dr_to_b < 0.4)
        {
          b_iso04 += trk.pt();
          if (dr_to_b < 0.3) b_iso03 += trk.pt();
        }
      }

      cand.addUserFloat("mu1_iso03", mu1_iso03);
      cand.addUserFloat("mu1_iso04", mu1_iso04);
      cand.addUserFloat("mu2_iso03", mu2_iso03);
      cand.addUserFloat("mu2_iso04", mu2_iso04);
      cand.addUserFloat("mu3_iso03" , mu3_iso03 );
      cand.addUserFloat("mu3_iso04" , mu3_iso04 );
      cand.addUserFloat("b_iso03" , b_iso03 );
      cand.addUserFloat("b_iso04" , b_iso04 );

      //salvo il candidato
      ret_val->push_back(cand);
    }//for(size_t ll_idx = 0; ll_idx < dimuons->size(); ++ll_idx) 
  }//for(size_t mu3_idx = 0; mu3_idx < muons->size(); ++mu3_idx)
  
  for (auto & cand: *ret_val){
    cand.addUserInt("n_mu3_used", std::count(used_trk_id.begin(),used_trk_id.end(),cand.userInt("mu3_idx")));
    cand.addUserInt("n_mu1_used", std::count(used_muon1_id.begin(),used_muon1_id.end(),cand.userInt("mu1_idx"))+std::count(used_muon2_id.begin(),used_muon2_id.end(),cand.userInt("mu1_idx")));
    cand.addUserInt("n_mu2_used", std::count(used_muon1_id.begin(),used_muon1_id.end(),cand.userInt("mu2_idx"))+std::count(used_muon2_id.begin(),used_muon2_id.end(),cand.userInt("mu2_idx")));
  }
  evt.put(std::move(ret_val));
}//produce  

Measurement1D BTo3MuBuilder::getIP(KinVtxFitter fitter,reco::TransientTrack transientTrackMu) const
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
DEFINE_FWK_MODULE(BTo3MuBuilder);

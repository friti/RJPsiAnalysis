#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/Math/interface/deltaR.h" 
#include <algorithm>

#include <vector>
#include <memory>
#include <map>
#include <string>
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "helper.h"
#include <limits>
#include <algorithm>
#include "KinVtxFitter.h"
#include "TLorentzVector.h"

class BTommmBuilder : public edm::global::EDProducer<> {

  // perhaps we need better structure here (begin run etc)
public:
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;

  explicit BTommmBuilder(const edm::ParameterSet &cfg):
    k_selection_{cfg.getParameter<std::string>("kaonSelection")},
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    dileptons_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("dileptons") )},
    leptons_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("leptonTransientTracks") )},
    kaons_{consumes<pat::MuonCollection>( cfg.getParameter<edm::InputTag>("kaons") )},
    kaons_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("kaonsTransientTracks") )},   
    isotracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
    isolostTracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))),

    //TRIGGER
    triggerBits_(consumes<edm::TriggerResults>(cfg.getParameter<edm::InputTag>("bits"))),
    triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone>>(cfg.getParameter<edm::InputTag>("objects"))),
    


    //GEN
    //src_(consumes<reco::CandidateView>(cfg.getParameter<edm::InputTag>("src"))),



    isotrk_selection_{cfg.getParameter<std::string>("isoTracksSelection")},
    beamspot_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") )} {
      produces<pat::CompositeCandidateCollection>();
    }

    
    
    
  ~BTommmBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const StringCutObjectSelector<pat::Muon> k_selection_; 
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; // cut on the di-lepton before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; // cut on the di-lepton after the SV fit

  const edm::EDGetTokenT<pat::CompositeCandidateCollection> dileptons_;
  const edm::EDGetTokenT<TransientTrackCollection> leptons_ttracks_;
  const edm::EDGetTokenT<pat::MuonCollection> kaons_;
  const edm::EDGetTokenT<TransientTrackCollection> kaons_ttracks_;

  const edm::EDGetTokenT<pat::PackedCandidateCollection> isotracksToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> isolostTracksToken_;
  
  //TRIGGER
  const edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  const edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggerObjects_;

  //GEN
  //  const edm::EDGetTokenT<reco::CandidateView> src_;

  const StringCutObjectSelector<pat::PackedCandidate> isotrk_selection_; 

  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;  

  
   
};

void BTommmBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  //lui fa da solo evento per evento
  //input
  //vuole accedere ai dati

  edm::Handle<pat::CompositeCandidateCollection> dileptons;
  evt.getByToken(dileptons_, dileptons);
  
  edm::Handle<TransientTrackCollection> leptons_ttracks;
  evt.getByToken(leptons_ttracks_, leptons_ttracks);

  //edm::Handle<pat::CompositeCandidateCollection> kaons;
  //evt.getByToken(kaons_, kaons);
  
  edm::Handle<pat::MuonCollection> kaons;
  evt.getByToken(kaons_,kaons);

  edm::Handle<TransientTrackCollection> kaons_ttracks;
  evt.getByToken(kaons_ttracks_, kaons_ttracks);  

  edm::Handle<reco::BeamSpot> beamspot;
  evt.getByToken(beamspot_, beamspot);  

  //for isolation
  edm::Handle<pat::PackedCandidateCollection> iso_tracks;
  evt.getByToken(isotracksToken_, iso_tracks);
  edm::Handle<pat::PackedCandidateCollection> iso_lostTracks;
  evt.getByToken(isolostTracksToken_, iso_lostTracks);
  
  //TRIGGER
  edm::Handle<edm::TriggerResults> triggerBits;
  evt.getByToken(triggerBits_, triggerBits);
  const edm::TriggerNames &names = evt.triggerNames(*triggerBits);
  std::vector<pat::TriggerObjectStandAlone> triggeringMuons;
  
  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerObjects;
  evt.getByToken(triggerObjects_, triggerObjects);
   

  /*std::cout << "\n === TRIGGER PATHS === " << std::endl;
  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
    std::cout << "Trigger " << names.triggerName(i) << std::endl;
    }*/

  
  //std::cout<<names.triggerNames()<<std::endl;
  

  //GEN
  /*  edm::Handle<reco::CandidateView> cands;
  evt.getByToken(src_, cands);
  unsigned int ncand = cands->size();
  //ma che sono questi candidati? hanno stessa dimensione di kaons, cioè i muoni displaced.
  std::cout<<"numero di candidati: "<<ncand<<", kaons size: "<<kaons->size()<<std::endl;
  //std::cout<<"cand[0]"<<cands->ptrAt(0)->pt()<<std::endl;
  */


  unsigned int nTracks     = iso_tracks->size();
  unsigned int totalTracks = nTracks + iso_lostTracks->size();

  std::vector<int> used_lep1_id, used_lep2_id, used_trk_id;

  
  // output
   std::unique_ptr<pat::CompositeCandidateCollection> ret_val(new pat::CompositeCandidateCollection());
  

  //TRIGGER
  //Controllo se ha passato il PATH giusto di trigger
  unsigned int index = names.triggerIndex("HLT_Dimuon0_Jpsi3p5_Muon2_v5");
  //std::cout<<index<<", "<<names.size()<<std::endl;    
  if(index==triggerBits->size()){
    evt.put(std::move(ret_val));
  }

  else{
  bool pass_hlt=triggerBits->accept(index); 
  
  //li salvo in questi array
  std::vector<pat::TriggerObjectStandAlone> pass_jpsi;
  std::vector<pat::TriggerObjectStandAlone> pass_3mu;

  if(pass_hlt){
   
    //loop sugli oggetti triggerati
    for (pat::TriggerObjectStandAlone obj : *triggerObjects){
      obj.unpackFilterLabels(evt, *triggerBits);
      obj.unpackPathNames(names);
    
      //controlla che sia un muone (necessario?)
      /*bool isTriggerMuon = false;
      for (unsigned h = 0; h < obj.filterIds().size(); ++h){
	if(obj.filterIds()[h] == 83){ 
	  isTriggerMuon = true; 
	  break;
	} 
	if(!isTriggerMuon) continue; 
      */
      /*      for (unsigned h = 0; h < obj.filterLabels().size(); ++h){
	std::string filterName = obj.filterLabels()[h];
	std::cout<<filterName<<std::endl;
	}*/


      //muoni della jpsi
      if(obj.hasFilterLabel("hltVertexmumuFilterJpsiMuon3p5")){
	pass_jpsi.push_back(obj);
      }
      
      if (obj.hasFilterLabel("hltTripleMuL3PreFiltered222")){
	pass_3mu.push_back(obj);
      }
      

    }
  }


  
  /*if(pass_jpsi.size()>0){
    std::cout<<"EHIIIIIIIIIIIIIIII"<<pass_jpsi.size()<<std::endl;
  } 
  if(pass_3mu.size()>0){
    std::cout<<"EHIIIIIIIIIIIIIIII"<<pass_3mu.size()<<std::endl;
    } */
  /*for(unsigned int i=0; i< pass_jpsi.size(); i++){
    std::cout<<pass_jpsi[i].pt()<<std::endl;
    }*/

  /*
      //    for (unsigned h = 0; h < obj.filterLabels().size(); ++h){
      //std::string filterName = obj.filterLabels()[h];
	//std::cout<<filterName<<std::endl;
        //if(filterName=="hltL3") != std::string::npos  && filterName.find("Park") != std::string::npos){
	//isTriggerMuon = true;
	//if(debug) std::cout << "\t   Filters:   " << filterName; 
	//break;
	}
	//else{ isTriggerMuon = false; }
	
  }
*/
  /*
    if(!isTriggerMuon) continue;
    triggeringMuons.push_back(obj);
    if(debug){ std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
    // Print trigger object collection and type
    std::cout << "\t   Collection: " << obj.collection() << std::endl;
    }
    }//trigger objects
    
    if(debug){
    std::cout << "\n total n of triggering muons = " << triggeringMuons.size() << std::endl;
      for(auto ij : triggeringMuons){
	std::cout << " >>> components (pt, eta, phi) = " << ij.pt() << " " << ij.eta() << " " << ij.phi() << std::endl;
      }
    */  
  
  
  /*//match online-offline di muoni che formano la jpsi
  std::vector<TLorentzVector> trig_obj_jpsi;
  std::vector<TLorentzVector> muon_jpsi;
  TLorentzVector online;
  TLorentzVector offline;

  for(pat::TriggerObjectStandAlone obj: pass_jpsi){
    for(pat::Composite muon: *kaons){
      //std::cout<<deltaR(obj,muon)<<std::endl;
      if(deltaR(obj,muon)<0.05){	
	online.SetPtEtaPhiM(obj.pt(),obj.eta(),obj.phi(),obj.mass());
	trig_obj_jpsi.push_back(online);
	offline.SetPtEtaPhiM(muon.pt(),muon.eta(),muon.phi(),muon.mass());
	muon_jpsi.push_back(offline);
      }
    }
    }*/


  //il match dei muoni che formano la jpsi e fatto dopo!!

  //vanno riconosciuti i muoni displaced
  std::vector<pat::TriggerObjectStandAlone> displaced;
  
  //std::vector<edm::Handle<pat::MuonCollection>> displaced;
  //  std::vector<unsigned int> i_jp;
  int flag=0;
  for(unsigned int i=0;i<pass_3mu.size() ;i++){
    flag=0;
    for(unsigned int j=0;j<pass_jpsi.size();j++){
      if(pass_3mu[i].pt()==pass_jpsi[j].pt()){
	//i_jp.push_back(i);
	//break
	//std::cout<<"pt1= "<<pass_3mu[i].pt()<<" ; pt2= "<<pass_jpsi[j].pt()<<std::endl;
	flag=1;
      }
    }
    if(flag==0){
      displaced.push_back(pass_3mu[i]);
    }
  }
  
  //matching online-offline del mu displaced
  std::vector<pat::TriggerObjectStandAlone> trig_obj_displ;
  std::vector<pat::Muon> muon_displ;
  for(pat::TriggerObjectStandAlone obj: displaced){
    for(pat::Muon muon: *kaons){
      if(deltaR(obj,muon)<0.05){	
	//online.SetPtEtaPhiM(obj.pt(),obj.eta(),obj.phi(),obj.mass());
	trig_obj_displ.push_back(obj);
	//offline.SetPtEtaPhiM(muon.pt(),muon.eta(),muon.phi(),muon.mass());
	muon_displ.push_back(muon);
      }
    }
  }
  
  //voglio che muon_displ sia il nuovo kaons.
  //muon_displ contiene i muoni displaced matchati con quelli triggerati online

  //std::cout << __PRETTY_FUNCTION__ << "]\t" << __LINE__ << std::endl;

  //Ciclo su tutti i muoni displaced 
  for(size_t k_idx = 0; k_idx < kaons->size(); ++k_idx) {
    
    edm::Ptr<pat::Muon> k_ptr(kaons, k_idx);
    if( !k_selection_(*k_ptr) ) continue;
    
    //controllo che questo kaon matchi con uno di quelli triggerati online
    //altrimenti lo skippa
    flag=0;
    for(pat::Muon muon: muon_displ){
      if(k_ptr->pt()==muon.pt()){
	  flag=1;
	  //std::cout<<k_ptr->pt()<<", "<<muon.pt()<<std::endl;
      }
    }
    if(flag==0){
      //std::cout<<"ADIOS"<<std::endl;
      continue;
      
    }
    
    
    math::PtEtaPhiMLorentzVector k_p4(
				      k_ptr->pt(), 
				      k_ptr->eta(),
				      k_ptr->phi(),
				      0.105658
				      );
    
    //math::PtEtaPhiMLorentzVector k_p4(muon_displ[k_idx].pt(),muon_displ[k_idx].eta(),muon_displ[k_idx].phi(),muon_displ[k_idx].mass());
    
    //ciclo sui leptoni che fanno la jpsi
    for(size_t ll_idx = 0; ll_idx < dileptons->size(); ++ll_idx) {
      edm::Ptr<pat::CompositeCandidate> ll_prt(dileptons, ll_idx);
      edm::Ptr<reco::Candidate> l1_ptr = ll_prt->userCand("l1");
      edm::Ptr<reco::Candidate> l2_ptr = ll_prt->userCand("l2");
      int l1_idx = ll_prt->userInt("l1_idx");
      int l2_idx = ll_prt->userInt("l2_idx");
    

      //matching online-offline di jpsi 
      std::vector<pat::TriggerObjectStandAlone> trig_obj_jpsi;
      std::vector<edm::Ptr<reco::Candidate>> muon_jpsi;
      flag=0;
      for(pat::TriggerObjectStandAlone obj: pass_jpsi){
	for(pat::TriggerObjectStandAlone obj2: pass_jpsi){
	  //std::cout<<deltaR(obj,muon)<<std::endl;
	  if((deltaR(obj,*l1_ptr)<0.05) & (deltaR(obj2,*l2_ptr)<0.05)){	
	    
	    trig_obj_jpsi.push_back(obj);
	    trig_obj_jpsi.push_back(obj2);
	    muon_jpsi.push_back(l1_ptr);
	    muon_jpsi.push_back(l2_ptr);
	    flag=1;
	    // std::cout<<obj.pt()<<", "<<l1_ptr->pt()<<", "<<obj2.pt()<<", "<<l2_ptr->pt()<<std::endl;
	  }
	}
      }
      //se le coppie non matchano, continua con la prossima coppia di leptoni
      if(flag==0){
	
	//std::cout<<"ADIOS2"<<std::endl;
	continue;
      }


      //si crea lo stato finale cand
      pat::CompositeCandidate cand;
      cand.setP4(ll_prt->p4() + k_p4);
      cand.setCharge(ll_prt->charge() + k_ptr->charge());
      // Use UserCands as they should not use memory but keep the Ptr itself
      // Put the lepton passing the corresponding selection
      cand.addUserCand("l1", l1_ptr);
      cand.addUserCand("l2", l2_ptr);
      cand.addUserCand("K",  k_ptr);
      cand.addUserCand("dilepton", ll_prt);

      cand.addUserInt("l1_idx", l1_idx);
      cand.addUserInt("l2_idx", l2_idx);
      cand.addUserInt("k_idx", k_idx);
    
      auto dr_info = min_max_dr({l1_ptr, l2_ptr,  k_ptr});
      //cand.addUserFloat("eMu",(k_ptr->pt()*k_ptr->pt()));
      cand.addUserFloat("min_dr", dr_info.first);
      cand.addUserFloat("max_dr", dr_info.second);
      // TODO add meaningful variables
      
      

      if( !pre_vtx_selection_(cand) ) continue;
    
      KinVtxFitter fitter(
        {leptons_ttracks->at(l1_idx), leptons_ttracks->at(l2_idx), kaons_ttracks->at(k_idx)},
        {l1_ptr->mass(), l2_ptr->mass(), 0.105658},
        {LEP_SIGMA, LEP_SIGMA, K_SIGMA} //some small sigma for the lepton mass
        );
      if(!fitter.success()) continue; // hardcoded, but do we need otherwise?
      cand.setVertex( 
        reco::Candidate::Point( 
          fitter.fitted_vtx().x(),
          fitter.fitted_vtx().y(),
          fitter.fitted_vtx().z()
          )  
        );
      used_lep1_id.emplace_back(l1_idx);
      used_lep2_id.emplace_back(l2_idx);
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

      cand.addUserFloat("fitted_l1_pt" , fitter.daughter_p4(0).pt()); 
      cand.addUserFloat("fitted_l1_eta", fitter.daughter_p4(0).eta());
      cand.addUserFloat("fitted_l1_phi", fitter.daughter_p4(0).phi());
      cand.addUserFloat("fitted_l2_pt" , fitter.daughter_p4(1).pt()); 
      cand.addUserFloat("fitted_l2_eta", fitter.daughter_p4(1).eta());
      cand.addUserFloat("fitted_l2_phi", fitter.daughter_p4(1).phi());
      cand.addUserFloat("fitted_k_pt"  , fitter.daughter_p4(2).pt()); 
      cand.addUserFloat("fitted_k_eta" , fitter.daughter_p4(2).eta());
      cand.addUserFloat("fitted_k_phi" , fitter.daughter_p4(2).phi());
 
      //mie variabili

      TLorentzVector P_b;
      P_b.SetPtEtaPhiM(fit_p4.pt(),fit_p4.eta(),fit_p4.phi(),fitter.fitted_candidate().mass());

      TLorentzVector P_k;
      P_k.SetPtEtaPhiM(fitter.daughter_p4(2).pt(),fitter.daughter_p4(2).eta(),fitter.daughter_p4(2).phi(),MUON_MASS);

      TLorentzVector P_l1;
      P_l1.SetPtEtaPhiM(fitter.daughter_p4(0).pt(),fitter.daughter_p4(0).eta(),fitter.daughter_p4(0).phi(),MUON_MASS);

      TLorentzVector P_l2;
      P_l2.SetPtEtaPhiM(fitter.daughter_p4(1).pt(),fitter.daughter_p4(1).eta(),fitter.daughter_p4(1).phi(),MUON_MASS);


      float m_miss_2=(P_b-P_k-P_l1-P_l2)*(P_b-P_k-P_l1-P_l2);
      
      float Q_2=(P_b-P_l1-P_l2)*(P_b-P_l1-P_l2);

      float pt_miss=(P_b.Pt()-P_k.Pt()-P_l1.Pt()-P_l2.Pt());

      float pt_var=((P_l1.Pt()+P_l2.Pt())-P_k.Pt());

      float dr_l1l2=sqrt((P_l1.Phi()-P_l2.Phi())*(P_l1.Phi()-P_l2.Phi())+(P_l1.Eta()-P_l2.Eta())*(P_l1.Eta()-P_l2.Eta()));

      cand.addUserFloat("m_miss_2", m_miss_2);
      cand.addUserFloat("Q_2",Q_2);
      cand.addUserFloat("pt_miss",pt_miss);
      cand.addUserFloat("pt_var",pt_var);
      cand.addUserFloat("dr_l1l2",dr_l1l2);
    
      //energia del mu unpaired in diversi sistemi di riferimento

      TLorentzVector P_mu;
      P_mu.SetPtEtaPhiM(fit_p4.pt(),fit_p4.eta(),fit_p4.phi(),MUON_MASS);

      TVector3 mu_beta_lab=P_mu.BoostVector();

      P_mu.Boost(-mu_beta_lab);
      

      cand.addUserFloat("E_mu_star",P_mu.E());

      //Risetto P_mu al valore originario perche altrimenti è boostato
      P_mu.SetPtEtaPhiM(fit_p4.pt(),fit_p4.eta(),fit_p4.phi(),MUON_MASS);
      
      //quadrivettori dei due muoni della jpsi
      //TLorentzVector P_mu_1;
      //TLorentzVector P_mu_2;

      //P_mu_1.SetPtEtaPhiM(fitter.daughter_p4(0).pt(), fitter.daughter_p4(0).eta(), fitter.daughter_p4(0).phi(),MUON_MASS);
      //P_mu_2.SetPtEtaPhiM(fitter.daughter_p4(1).pt(), fitter.daughter_p4(1).eta(), fitter.daughter_p4(1).phi(),MUON_MASS);

      TLorentzVector jpsi=P_l1+P_l1;

      TVector3 jpsi_beta_lab=jpsi.BoostVector();

      P_mu.Boost(jpsi_beta_lab);
      
      cand.addUserFloat("E_mu_#",P_mu.E());


      

      if( !post_vtx_selection_(cand) ) continue;        

      //compute isolation
      float l1_iso03 = 0;
      float l1_iso04 = 0;
      float l2_iso03 = 0;
      float l2_iso04 = 0;
      float k_iso03  = 0;
      float k_iso04  = 0;
      float b_iso03  = 0;
      float b_iso04  = 0;

      for( unsigned int iTrk=0; iTrk<totalTracks; ++iTrk ) {
      
        const pat::PackedCandidate & trk = (iTrk < nTracks) ? (*iso_tracks)[iTrk] : (*iso_lostTracks)[iTrk-nTracks];
        // define selections for iso tracks (pT, eta, ...)
        if( !isotrk_selection_(trk) ) continue;
        // check if the track is the kaon
        if ( k_ptr->userCand("cand") ==  edm::Ptr<reco::Candidate> ( iso_tracks, iTrk ) ) continue;
        // check if the track is one of the two leptons
        if (track_to_lepton_match(l1_ptr, iso_tracks.id(), iTrk) || 
            track_to_lepton_match(l2_ptr, iso_tracks.id(), iTrk) ) continue;

        // add to final particle iso if dR < cone
        float dr_to_l1 = deltaR(cand.userFloat("fitted_l1_eta"), cand.userFloat("fitted_l1_phi"), trk.eta(), trk.phi());
        float dr_to_l2 = deltaR(cand.userFloat("fitted_l2_eta"), cand.userFloat("fitted_l2_phi"), trk.eta(), trk.phi());
        float dr_to_k  = deltaR(cand.userFloat("fitted_k_eta") , cand.userFloat("fitted_k_phi") , trk.eta(), trk.phi());
        float dr_to_b  = deltaR(cand.userFloat("fitted_eta")   , cand.userFloat("fitted_phi") , trk.eta(), trk.phi());

        if (dr_to_l1 < 0.4){
          l1_iso04 += trk.pt();
          if ( dr_to_l1 < 0.3) l1_iso03 += trk.pt();
        }
        if (dr_to_l2 < 0.4){
          l2_iso04 += trk.pt();
          if (dr_to_l2 < 0.3)  l2_iso03 += trk.pt();
        }
        if (dr_to_k < 0.4){
          k_iso04 += trk.pt();
          if (dr_to_k < 0.3) k_iso03 += trk.pt();
        }
        if (dr_to_b < 0.4){
          b_iso04 += trk.pt();
          if (dr_to_b < 0.3) b_iso03 += trk.pt();
        }
      }
      cand.addUserFloat("l1_iso03", l1_iso03);
      cand.addUserFloat("l1_iso04", l1_iso04);
      cand.addUserFloat("l2_iso03", l2_iso03);
      cand.addUserFloat("l2_iso04", l2_iso04);
      cand.addUserFloat("k_iso03" , k_iso03 );
      cand.addUserFloat("k_iso04" , k_iso04 );
      cand.addUserFloat("b_iso03" , b_iso03 );
      cand.addUserFloat("b_iso04" , b_iso04 );

      ret_val->push_back(cand);}
    // for(size_t ll_idx = 0; ll_idx < dileptons->size(); ++ll_idx) {
  

    
  } // for(size_t k_idx = 0; k_idx < kaons->size(); ++k_idx)
 
  /*
  //controllo della massa della jpsi
  pat::CompositeCandidate jPsi;
  for (pat::TriggerObjectStandAlone obj: pass_jpsi){
      for (pat::TriggerObjectStandAlone obj2: pass_jpsi){
	 math::PtEtaPhiMLorentzVector l1(obj.pt(),obj.eta(),obj.phi(),obj.mass());
	 math::PtEtaPhiMLorentzVector l2(obj2.pt(),obj2.eta(),obj2.phi(),obj2.mass());
	 jPsi.setP4(l1+l2);
	//std::cout<<jPsi.M()<<std::endl;
	 jPsi.addUserFloat("jPsi_mass_online",jPsi.mass());
      }
  }
  
  if(pass_jpsi.size()==0 || pass_hlt!=1){
    jPsi.addUserFloat("jPsi_mass_online",-99);
    }*/
  

  for (auto & cand: *ret_val){
    cand.addUserInt("n_k_used", std::count(used_trk_id.begin(),used_trk_id.end(),cand.userInt("k_idx")));
    cand.addUserInt("n_l1_used", std::count(used_lep1_id.begin(),used_lep1_id.end(),cand.userInt("l1_idx"))+std::count(used_lep2_id.begin(),used_lep2_id.end(),cand.userInt("l1_idx")));
    cand.addUserInt("n_l2_used", std::count(used_lep1_id.begin(),used_lep1_id.end(),cand.userInt("l2_idx"))+std::count(used_lep2_id.begin(),used_lep2_id.end(),cand.userInt("l2_idx")));
  }

  evt.put(std::move(ret_val));}
}  
    

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(BTommmBuilder);

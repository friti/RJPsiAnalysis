#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <vector>
#include <memory>
#include <map>
#include <string>
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "helper.h"
#include <limits>
#include <algorithm>
#include "KinVtxFitter.h"

class DiMuonBuilder : public edm::global::EDProducer<> {

  // perhaps we need better structure here (begin run etc)
public:
  typedef std::vector<pat::Muon> MuonCollection;
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;

  explicit DiMuonBuilder(const edm::ParameterSet &cfg):
    mu1_selection_{cfg.getParameter<std::string>("muon1Selection")},
    mu2_selection_{cfg.getParameter<std::string>("muon2Selection")},
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    src_{consumes<MuonCollection>( cfg.getParameter<edm::InputTag>("src") )},
    ttracks_src_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("transientTracksSrc") )} {
       produces<pat::CompositeCandidateCollection>();
    }

  ~DiMuonBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const StringCutObjectSelector<pat::Muon> mu1_selection_; // cut on leading muon
  const StringCutObjectSelector<pat::Muon> mu2_selection_; // cut on sub-leading muon
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; // cut on the di-muon before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; // cut on the di-muon after the SV fit
  const edm::EDGetTokenT<MuonCollection> src_;
  const edm::EDGetTokenT<TransientTrackCollection> ttracks_src_;
};

void DiMuonBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  //input
  edm::Handle<MuonCollection> muons;
  evt.getByToken(src_, muons);
  
  edm::Handle<TransientTrackCollection> ttracks;
  evt.getByToken(ttracks_src_, ttracks);

  // output
  std::unique_ptr<pat::CompositeCandidateCollection> ret_value(new pat::CompositeCandidateCollection());
  
  for(size_t mu1_idx = 0; mu1_idx < muons->size(); ++mu1_idx) {
    edm::Ptr<pat::Muon> mu1_ptr(muons, mu1_idx);
    if(!mu1_selection_(*mu1_ptr)) continue; 
    int isJpsiMuon1 = mu1_ptr->userInt("isJpsiMuon");
    //int isDimuon0Trg1 = mu1_ptr->userInt("isDimuon0Trg");
    int isJpsiTrkTrg1 = mu1_ptr->userInt("isJpsiTrkTrg");
    
    for(size_t mu2_idx = mu1_idx + 1; mu2_idx < muons->size(); ++mu2_idx) {
      edm::Ptr<pat::Muon> mu2_ptr(muons, mu2_idx);
      if(!mu2_selection_(*mu2_ptr)) continue;
      // Form pairs only with triggered muons
      int isJpsiMuon2 = mu2_ptr->userInt("isJpsiMuon");
      //int isDimuon0Trg2 = mu2_ptr->userInt("isDimuon0Trg");
      int isJpsiTrkTrg2 = mu2_ptr->userInt("isJpsiTrkTrg");
      bool trg1 = (isJpsiTrkTrg1 && isJpsiTrkTrg2);
      bool trg2 = (isJpsiMuon1 && isJpsiMuon2);

      //if(!trg1 && !trg2) continue;
      if(!trg2) continue;
      //std::cout << isDimuon0Trg1 

      pat::CompositeCandidate muon_pair;
      muon_pair.setP4(mu1_ptr->p4() + mu2_ptr->p4());
      muon_pair.setCharge(mu1_ptr->charge() + mu2_ptr->charge());
      muon_pair.addUserFloat("muons12_deltaR", reco::deltaR(*mu1_ptr, *mu2_ptr));
      // Put the muon passing the corresponding selection
      muon_pair.addUserInt("mu1_idx", mu1_idx );
      muon_pair.addUserInt("mu2_idx", mu2_idx );
      // Use UserCands as they should not use memory but keep the Ptr itself
      muon_pair.addUserCand("mu1", mu1_ptr );
      muon_pair.addUserCand("mu2", mu2_ptr );
           if( !pre_vtx_selection_(muon_pair) ) continue; // before making the SV, cut on the info we have

      KinVtxFitter fitter(
        {ttracks->at(mu1_idx), ttracks->at(mu2_idx)},
        {mu1_ptr->mass(), mu2_ptr->mass()},
        {LEP_SIGMA, LEP_SIGMA} //some small sigma for the particle mass
        );

      muon_pair.addUserFloat("sv_chi2", fitter.chi2());
      muon_pair.addUserFloat("sv_ndof", fitter.dof()); // float??
      muon_pair.addUserFloat("sv_prob", fitter.prob());
      muon_pair.addUserFloat("fitted_mass", fitter.success() ? fitter.fitted_candidate().mass() : -1);
      muon_pair.addUserFloat("fitted_massErr", fitter.success() ? sqrt(fitter.fitted_candidate().kinematicParametersError().matrix()(6,6)) : -1);
      muon_pair.addUserFloat("vtx_ex",1.0);
      /* muon_pair.addUserFloat("vtx_x",muon_pair.vx());
      muon_pair.addUserFloat("vtx_y",muon_pair.vy());
      muon_pair.addUserFloat("vtx_z",muon_pair.vz());
      muon_pair.addUserFloat("vtx_ex",sqrt(fitter.fitted_vtx_uncertainty().cxx()));
      muon_pair.addUserFloat("vtx_ey",sqrt(fitter.fitted_vtx_uncertainty().cyy()));
      muon_pair.addUserFloat("vtx_ez",sqrt(fitter.fitted_vtx_uncertainty().czz()));
      */
     
      // if needed, add here more stuff

      // cut on the SV info
      if( !post_vtx_selection_(muon_pair) ) continue;
      ret_value->push_back(muon_pair);
    }
  }
  
  evt.put(std::move(ret_value));
}

DEFINE_FWK_MODULE(DiMuonBuilder);

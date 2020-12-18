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
constexpr bool debug = false;

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
       produces<pat::CompositeCandidateCollection>("muonPairsForBTo3Mu");
       produces<TransientTrackCollection>("dimuonTransientTracks");
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
  std::unique_ptr<TransientTrackCollection> dimuon_tt(new TransientTrackCollection);
  
  for(size_t mu1_idx = 0; mu1_idx < muons->size(); ++mu1_idx) {
    edm::Ptr<pat::Muon> mu1_ptr(muons, mu1_idx);
    if(!mu1_selection_(*mu1_ptr)) continue; 
    int isJpsiMuon1 = mu1_ptr->userInt("isJpsiMuon");
    int isDimuon0Trg1 = mu1_ptr->userInt("isDimuon0Trg");
    int isJpsiTrkTrg1 = mu1_ptr->userInt("isJpsiTrkTrg");
    for(size_t mu2_idx = mu1_idx + 1; mu2_idx < muons->size(); ++mu2_idx) {
      edm::Ptr<pat::Muon> mu2_ptr(muons, mu2_idx);
      if(!mu2_selection_(*mu2_ptr)) continue;
      // Form pairs only with triggered muons
      int isJpsiMuon2 = mu2_ptr->userInt("isJpsiMuon");
      int isDimuon0Trg2 = mu2_ptr->userInt("isDimuon0Trg");
      int isJpsiTrkTrg2 = mu2_ptr->userInt("isJpsiTrkTrg");
      bool trg1 = ((isJpsiTrkTrg1 && isJpsiMuon1) && (isJpsiTrkTrg2 && isJpsiMuon2));
      bool trg2 = ((isDimuon0Trg1 && isJpsiMuon1) && (isDimuon0Trg2 && isJpsiMuon2));

      int dimuon0_trigger = 0;
      int jpsitrk_trigger = 0;
      if(!trg1 && !trg2) continue;
      //if(!trg2) continue;
      if(trg1) jpsitrk_trigger = 1;
      if(trg2) dimuon0_trigger = 1;
      
      pat::CompositeCandidate muon_pair;
      muon_pair.setP4(mu1_ptr->p4() + mu2_ptr->p4());
      muon_pair.setCharge(mu1_ptr->charge() + mu2_ptr->charge());
      muon_pair.addUserFloat("muons12_deltaR", reco::deltaR(*mu1_ptr, *mu2_ptr));
      // Put the muon passing the corresponding selection
      muon_pair.addUserInt("mu1_idx", mu1_idx );
      muon_pair.addUserInt("mu2_idx", mu2_idx );
      muon_pair.addUserInt("isJpsiTrkTrg", trg1);
      muon_pair.addUserInt("isDimuon0Trg", trg2);
      // Use UserCands as they should not use memory but keep the Ptr itself
      muon_pair.addUserCand("mu1", mu1_ptr );
      muon_pair.addUserCand("mu2", mu2_ptr );
      if(debug) std::cout<<"l1 "<<mu1_ptr->pt()<<" l2 "<<mu2_ptr->pt()<<" mass "<<muon_pair.mass()<<" deltaR"<<reco::deltaR(*mu1_ptr, *mu2_ptr)<<" dz "<<mu1_ptr->bestTrack()->dz()-mu2_ptr->bestTrack()->dz()<<std::endl;
      if( !pre_vtx_selection_(muon_pair) ) {
	if(debug) std::cout<<"pre vtx selection dies"<<std::endl;

	continue; // before making the SV, cut on the info we have
      }
      KinVtxFitter fitter(
        {ttracks->at(mu1_idx), ttracks->at(mu2_idx)},
        {mu1_ptr->mass(), mu2_ptr->mass()},
        {LEP_SIGMA, LEP_SIGMA} //some small sigma for the particle mass
        );

      muon_pair.addUserFloat("sv_chi2", fitter.chi2());
      //if(debug) std::cout << "vx_: " << fitter.fitted_candidate() << std::endl;
      muon_pair.addUserFloat("sv_position", fitter.fitted_vtx().x()); // float??
      muon_pair.addUserFloat("sv_ndof", fitter.dof()); // float??
      muon_pair.addUserFloat("sv_prob", fitter.prob());
      muon_pair.addUserFloat("fitted_mass", fitter.success() ? fitter.fitted_candidate().mass() : -1);
      muon_pair.addUserFloat("fitted_massErr", fitter.success() ? sqrt(fitter.fitted_candidate().kinematicParametersError().matrix()(6,6)) : -1);
      muon_pair.addUserFloat("vtx_x",muon_pair.vx());
      muon_pair.addUserFloat("vtx_y",muon_pair.vy());
      muon_pair.addUserFloat("vtx_z",muon_pair.vz());
      muon_pair.addUserFloat("vtx_ex",sqrt(fitter.fitted_vtx_uncertainty().cxx()));
      muon_pair.addUserFloat("vtx_ey",sqrt(fitter.fitted_vtx_uncertainty().cyy()));
      muon_pair.addUserFloat("vtx_ez",sqrt(fitter.fitted_vtx_uncertainty().czz()));
      muon_pair.addUserFloat("vtx_chi2", fitter.chi2());
     
      // if needed, add here more stuff

      // cut on the SV info
      // const reco::TransientTrack& fitted_candidate_ttrk()
      if( !post_vtx_selection_(muon_pair) ) continue;
      if(!fitter.fitted_candidate_ttrk().isValid()) continue;
      dimuon_tt->emplace_back(fitter.fitted_candidate_ttrk());
      //ret_value->push_back(muon_pair);
      ret_value->emplace_back(muon_pair);
      ret_value->back().addUserInt("muonpair_fromdimuon0", dimuon0_trigger);
      ret_value->back().addUserInt("muonpair_fromjpsitrk", jpsitrk_trigger);
    }
  }
  
  evt.put(std::move(ret_value), "muonPairsForBTo3Mu");
  evt.put(std::move(dimuon_tt), "dimuonTransientTracks");
}

DEFINE_FWK_MODULE(DiMuonBuilder);

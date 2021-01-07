// Class to select the PV from the collection
// selecting the one closest in z to 
// every secondary vertex


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"


#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "helper.h"

using namespace std;

constexpr bool debug = false;

class PrimaryVertexSelector : public edm::EDProducer {
public:
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;
  
  explicit PrimaryVertexSelector(const edm::ParameterSet &iConfig);
  ~PrimaryVertexSelector() override {};

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);
  const edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> dimuonsSrc_;
  const edm::EDGetTokenT<TransientTrackCollection> dimuonsTTSrc_;
};

PrimaryVertexSelector::PrimaryVertexSelector(const edm::ParameterSet &iConfig):
  vertexSrc_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"))),
  dimuonsSrc_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("dimuonCollection"))),
  dimuonsTTSrc_(consumes<TransientTrackCollection>(iConfig.getParameter<edm::InputTag>("dimuonTTCollection")))
{
  produces<reco::VertexCollection>("bestVertex");
}

void PrimaryVertexSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // Inputs
  edm::Handle<pat::CompositeCandidateCollection> dimuons;
  iEvent.getByToken(dimuonsSrc_, dimuons);

  edm::Handle<TransientTrackCollection> dimuons_tt;
  iEvent.getByToken(dimuonsTTSrc_,dimuons_tt);

  edm::Handle<reco::VertexCollection> thePrimaryVerticesHandle;
  iEvent.getByToken(vertexSrc_, thePrimaryVerticesHandle);
  // Chosing the closest PV in Z direction to the JPsi trajectory projection.

  std::unique_ptr<reco::VertexCollection> bestVertex_out ( new reco::VertexCollection );

  for(size_t tt_idx = 0; tt_idx < dimuons_tt->size(); ++tt_idx)
  {
    const reco::TransientTrack& dimuonTT = dimuons_tt->at(tt_idx);
    double dzMin = 1000000.;
    reco::Vertex bestVertex;
    const reco::VertexCollection* vertices = thePrimaryVerticesHandle.product();
    //for(reco::VertexCollection::const_iterator  primVertex = vertices->begin(); primVertex!= vertices->end(); primVertex++)
    for(size_t i = 0; i < vertices->size() ; i++)
    {
      reco::Vertex primVertex = vertices->at(i);
      //std::cout << "prim vertex z: " << primVertex->z() << std::endl;
      if (abs(dzMin) > abs(dimuonTT.track().dz(primVertex.position())))
      {
        bestVertex = primVertex;
        //bestVertex = primVertex;
        dzMin = dimuonTT.track().dz(primVertex.position());
      }
    }   
    //cout<< "Best vertex: " << bestVertex.x() << endl;
    bestVertex_out->emplace_back(bestVertex);
  }
  iEvent.put(std::move(bestVertex_out), "bestVertex");
}
DEFINE_FWK_MODULE(PrimaryVertexSelector);

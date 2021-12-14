// Class to include the flags for the mother of jpsi -> for the sample Hb -> jpsi X 
// and the the truth match

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h" 
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "helper.h"

using namespace std;

constexpr bool debug = false;

bool isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle)
{
  //particle is already the ancestor
  if(ancestor == particle ) return true;

  //otherwise loop on mothers, if any and return true if the ancestor is found
  for(size_t i=0;i< particle->numberOfMothers();i++)
    {
      //         std::cout << "checking if ancestor: " << particle->mother(i)->pdgId() << "   " << ancestor->pdgId() << std::endl;
      if(isAncestor(ancestor,particle->mother(i))) return true;
    }
  //if we did not return yet, then particle and ancestor are not relatives
  return false;
}

int diquarks[] = { 1103,2101,2103,2203,3101,3103,3201,3203,3303,4101,4103,4201,4203,4301,4303,4403,5101,5103,5201,5203,5301,5303,5401, 5403,5503};


const reco::Candidate* checkMom(const reco::Candidate * candMom){

  if (candMom == nullptr) return nullptr;
  
  if (candMom->mother(0) == nullptr) {
    return candMom;
  }  
  int * p = std::find (diquarks, diquarks+25, candMom->mother(0)->pdgId());
  if (abs(candMom->mother(0)->pdgId()) < 8  || \
      abs(candMom->mother(0)->pdgId())== 21 || \
      abs(candMom->mother(0)->pdgId())== 2212 || \
      (p != (diquarks+25))
      ){ 
    //   else if (abs(candMom->mother(0)->pdgId()) < 9 || abs(candMom->mother(0)->pdgId())==21 || abs(candMom->mother(0)->pdgId())==2212){ 
    //     std::cout << "check mom - input particle mother is prime: " << candMom->mother(0)->pdgId() << std::endl;
    return candMom;
  }
  else {
    candMom = checkMom(candMom->mother(0));
    return candMom;
  }  
}



class jpsiMotherFlagsTableProducer : public edm::EDProducer {
public:
  typedef std::vector<reco::GenParticle> GenParticleCollection;
  typedef std::vector<pat::PackedGenParticle> PackedGenParticleCollection;

  explicit jpsiMotherFlagsTableProducer(const edm::ParameterSet &iConfig):
    prunedToken_{consumes<GenParticleCollection>(iConfig.getParameter<edm::InputTag>("prunedGen"))},
    packedToken_{consumes<PackedGenParticleCollection>(iConfig.getParameter<edm::InputTag>("packedGen"))}
  {
    produces<nanoaod::FlatTable>();
  }
  ~jpsiMotherFlagsTableProducer() override {};

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);
  //GEN
  const edm::EDGetTokenT<reco::GenParticleCollection> prunedToken_; 
  const edm::EDGetTokenT<pat::PackedGenParticleCollection> packedToken_; 
};

//jpsiMotherFlagsTableProducer::jpsiMotherFlagsTableProducer(const edm::ParameterSet &iConfig):
//  prunedToken_(consumes<GenParticleCollection>(iConfig.getParameter<edm::InputTag>("prunedGen")))
//packedToken_(consumes<PackedGenParticleCollection>(iConfig.getParameter<edm::InputTag>("packedGen")))
//{
//  produces<nanoaod::FlatTable>();
//}
void jpsiMotherFlagsTableProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    
  if(debug) std::cout<<"In the jpsi mother flag code..."<<std::endl;
  edm::Handle<GenParticleCollection> pruned;
  iEvent.getByToken(prunedToken_, pruned);
  edm::Handle<PackedGenParticleCollection> packed;
  iEvent.getByToken(packedToken_, packed);


  bool foundJpsiMup = false;
  bool foundJpsiMum = false;
    
  std::vector<int> broPdgId;
  std::vector<double> broPt, broEta, broPhi;
  std::vector<const reco::Candidate*> anc_vec;
  //const reco::Candidate *theMom = nullptr; 
  int theMom_pdgId; 

  //Looking for jpsi daughters -> muons
  for (const reco::GenParticle &jpsiMeson : *pruned) {
    if (abs(jpsiMeson.pdgId()) == 443 ){
      if(debug) std::cout<<"There is a jpsi"<<std::endl;
      foundJpsiMup = false; // plus muon
      foundJpsiMum = false; //minus muon

      for (const pat::PackedGenParticle &dau : *packed) {
	const reco::Candidate * motherInPrunedCollection = dau.mother(0) ;
	if(motherInPrunedCollection != nullptr && isAncestor( &jpsiMeson , motherInPrunedCollection)){
	  if (dau.pdgId() == -13)     {
	    foundJpsiMum = true; 

	  }
	  else if (dau.pdgId() == 13) {
	    foundJpsiMup = true; 

	  }  
	}
      }
      if(debug) std::cout<<"Is there a plus and a minus mu?"<<foundJpsiMup<<" "<<foundJpsiMum<<std::endl;
      if (!(foundJpsiMup && foundJpsiMum)) continue;    
      if(debug) std::cout<<"YES"<<std::endl;

    
      // find jpsi original ancestor -> To Save!
      if(debug) std::cout<<"Number of mothers of the jpsi "<<jpsiMeson.numberOfMothers()<<std::endl;
      if ( jpsiMeson.numberOfMothers()>1 && debug) std::cout << "number of jpsi mom: " << jpsiMeson.numberOfMothers() << std::endl;
            
      for (size_t imom = 0; imom < jpsiMeson.numberOfMothers(); imom++){
	if(debug) std::cout<<"jpsimother"<<jpsiMeson.mother(imom)->pdgId()<<std::endl;
	const reco::Candidate * candMom = checkMom( jpsiMeson.mother(imom)) ;
	if(debug) std::cout<<"candMom"<<candMom->pdgId()<<std::endl;
	
	if (candMom == nullptr) continue;
	if (candMom != nullptr && isAncestor(candMom, &jpsiMeson )){
	  if ( candMom->mother(0) != nullptr) {

	    //theMom = candMom;
	    theMom_pdgId = candMom->pdgId();
	    if(debug) std::cout<<"We found a mom "<<theMom_pdgId<<std::endl;
	    break;
	  }
	}
      }
    }
  }        


  //GEN
  if(debug) std::cout<<"After the loops the pdg of the mom is  "<<theMom_pdgId<<std::endl;
  
  uint8_t flag_bzero = 0;
  uint8_t flag_bplus = 0;
  uint8_t flag_bzero_s = 0;
  uint8_t flag_bplus_c = 0;
  uint8_t flag_sigmaminus_b = 0;
  uint8_t flag_lambdazero_b = 0;
  uint8_t flag_ximinus_b = 0;
  uint8_t flag_sigmazero_b = 0;
  uint8_t flag_xizero_b = 0;
  uint8_t flag_other = 0;

  if(abs(theMom_pdgId) == 511 || abs(theMom_pdgId) == 513) flag_bzero = 1.;
  else if(abs(theMom_pdgId) == 521 || abs(theMom_pdgId) == 523) flag_bplus = 1;
  else if(abs(theMom_pdgId) == 531 || abs(theMom_pdgId) == 533) flag_bzero_s = 1;
  else if(abs(theMom_pdgId) == 541 || abs(theMom_pdgId) == 543) flag_bplus_c = 1;
  else if(abs(theMom_pdgId) == 5112 || abs(theMom_pdgId) == 5114) flag_sigmaminus_b = 1;
  else if(abs(theMom_pdgId) == 5122) flag_lambdazero_b = 1;
  else if(abs(theMom_pdgId) == 5132 || abs(theMom_pdgId) == 5314) flag_ximinus_b = 1;
  else if(abs(theMom_pdgId) == 5212 || abs(theMom_pdgId) == 5214) flag_sigmazero_b = 1;
  else if(abs(theMom_pdgId) == 5232 || abs(theMom_pdgId) == 5324) flag_xizero_b = 1;
  else flag_other = 1;

  auto tab = std::make_unique<nanoaod::FlatTable>(1,"",true);

  tab->addColumnValue<uint8_t>("JpsiMotherFlag_bplus", flag_bplus, "decay flag", nanoaod::FlatTable::UInt8Column);
  tab->addColumnValue<uint8_t>("JpsiMotherFlag_bzero", flag_bzero, "decay flag", nanoaod::FlatTable::UInt8Column);
  tab->addColumnValue<uint8_t>("JpsiMotherFlag_bzero_s", flag_bzero_s, "decay flag", nanoaod::FlatTable::UInt8Column);
  tab->addColumnValue<uint8_t>("JpsiMotherFlag_bplus_c", flag_bplus_c, "decay flag", nanoaod::FlatTable::UInt8Column);
  tab->addColumnValue<uint8_t>("JpsiMotherFlag_sigmaminus_b", flag_sigmaminus_b, "decay flag", nanoaod::FlatTable::UInt8Column);
  tab->addColumnValue<uint8_t>("JpsiMotherFlag_lambdazero_b", flag_lambdazero_b, "decay flag", nanoaod::FlatTable::UInt8Column);
  tab->addColumnValue<uint8_t>("JpsiMotherFlag_ximinus_b", flag_ximinus_b, "decay flag", nanoaod::FlatTable::UInt8Column);
  tab->addColumnValue<uint8_t>("JpsiMotherFlag_sigmazero_b", flag_sigmazero_b, "decay flag", nanoaod::FlatTable::UInt8Column);
  tab->addColumnValue<uint8_t>("JpsiMotherFlag_xizero_b", flag_xizero_b, "decay flag", nanoaod::FlatTable::UInt8Column);
  tab->addColumnValue<uint8_t>("JpsiMotherFlag_other", flag_other, "decay flag", nanoaod::FlatTable::UInt8Column);

  iEvent.put(std::move(tab));
}

DEFINE_FWK_MODULE(jpsiMotherFlagsTableProducer);

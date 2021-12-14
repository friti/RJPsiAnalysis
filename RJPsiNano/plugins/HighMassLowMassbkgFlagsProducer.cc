// Class to include the flags for the mother of jpsi -> for the sample Hb -> jpsi X 
// and the the truth match

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
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

bool Ancestor(const reco::Candidate* ancestor, const reco::Candidate * particle)
{
  //particle is already the ancestor
  if(ancestor == particle ) return true;

  //otherwise loop on mothers, if any and return true if the ancestor is found
  for(size_t i=0;i< particle->numberOfMothers();i++)
    {
      //         std::cout << "checking if ancestor: " << particle->mother(i)->pdgId() << "   " << ancestor->pdgId() << std::endl;
      if(Ancestor(ancestor,particle->mother(i))) return true;
    }
  //if we did not return yet, then particle and ancestor are not relatives
  return false;
}

int diquark[] = { 1103,2101,2103,2203,3101,3103,3201,3203,3303,4101,4103,4201,4203,4301,4303,4403,5101,5103,5201,5203,5301,5303,5401, 5403,5503};

template <class T>
unsigned int findIdx(T cand, edm::Handle<std::vector<reco::GenParticle>> pruned){
  //unsigned int findIdx(const reco::Candidate * cand, edm::Handle<std::vector<reco::GenParticle>> pruned){

  unsigned int idx;
  for (unsigned int ipart = 0; ipart < pruned->size(); ipart++) {
    const reco::GenParticle particle = pruned->at(ipart);
    if (cand->pt() == particle.pt() && cand->pdgId() == particle.pdgId() && cand->eta() == particle.eta() && cand->phi() == particle.phi()){
      idx = ipart;
      break;
    }
  }
  return idx;
}

const reco::Candidate* checkmom(const reco::Candidate * candMom){

  if (candMom == nullptr) return nullptr;
  
  if (candMom->mother(0) == nullptr) {
    return candMom;
  }  
  int * p = std::find (diquark, diquark+25, candMom->mother(0)->pdgId());
  if (abs(candMom->mother(0)->pdgId()) < 8  || \
      abs(candMom->mother(0)->pdgId())== 21 || \
      abs(candMom->mother(0)->pdgId())== 2212 || \
      (p != (diquark+25))
      ){ 
    //   else if (abs(candMom->mother(0)->pdgId()) < 9 || abs(candMom->mother(0)->pdgId())==21 || abs(candMom->mother(0)->pdgId())==2212){ 
    //     std::cout << "check mom - input particle mother is prime: " << candMom->mother(0)->pdgId() << std::endl;
    return candMom;
  }
  else {
    candMom = checkmom(candMom->mother(0));
    return candMom;
  }  
}


class HighMassLowMassbkgFlagsProducer : public edm::EDProducer {
public:
  typedef std::vector<reco::GenParticle> GenParticleCollection;
  typedef std::vector<pat::PackedGenParticle> PackedGenParticleCollection;

  explicit HighMassLowMassbkgFlagsProducer(const edm::ParameterSet &iConfig):
    prunedToken_{consumes<GenParticleCollection>(iConfig.getParameter<edm::InputTag>("prunedGen"))},
    packedToken_{consumes<PackedGenParticleCollection>(iConfig.getParameter<edm::InputTag>("packedGen"))}
  {
    produces<pat::CompositeCandidateCollection>("AncestorFlags");
  }
  ~HighMassLowMassbkgFlagsProducer() override {};

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);
  //GEN
  const edm::EDGetTokenT<reco::GenParticleCollection> prunedToken_; 
  const edm::EDGetTokenT<pat::PackedGenParticleCollection> packedToken_; 
};

//HighMassLowMassbkgFlagsProducer::HighMassLowMassbkgFlagsProducer(const edm::ParameterSet &iConfig):
//  prunedToken_(consumes<GenParticleCollection>(iConfig.getParameter<edm::InputTag>("prunedGen")))
//packedToken_(consumes<PackedGenParticleCollection>(iConfig.getParameter<edm::InputTag>("packedGen")))
//{
//  produces<nanoaod::FlatTable>();
//}
void HighMassLowMassbkgFlagsProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    
  if(debug) std::cout<<"In the bkg hb flag code..."<<std::endl;
  edm::Handle<GenParticleCollection> pruned;
  iEvent.getByToken(prunedToken_, pruned);
  edm::Handle<PackedGenParticleCollection> packed;
  iEvent.getByToken(packedToken_, packed);

  //Define the vectors for the final Table
  std::unique_ptr<pat::CompositeCandidateCollection> jpsi_out( new pat::CompositeCandidateCollection() );

  bool foundJpsiMup = false;
  bool foundJpsiMum = false;
    
  std::vector<int> broPdgId;
  std::vector<double> broPt, broEta, broPhi;
  std::vector<const reco::Candidate*> anc_vec;

  int theMom_pdgId;
  unsigned int mu1_idx = -1;
  unsigned int mu2_idx = -1;
  unsigned int jpsi_ancestor_idx = -1;

  // JPSI
  //Looking for jpsi daughters -> muons
  for (unsigned int ipart = 0; ipart < pruned->size(); ipart++) {
    //for (const reco::GenParticle &jpsiMeson : *pruned) {
    const reco::GenParticle &jpsiMeson = pruned->at(ipart);
    if (abs(jpsiMeson.pdgId()) == 443 ){
      if(debug) std::cout<<"There is a jpsi"<<std::endl;
      if(debug) std::cout<<"The idx is "<<ipart<<std::endl;
      foundJpsiMup = false; // plus muon
      foundJpsiMum = false; //minus muon
      
      //for (const pat::PackedGenParticle* dau : &packed) {
      for (auto dau : *packed) {
	const reco::Candidate * motherInPrunedCollection = dau.mother(0) ;
	if(motherInPrunedCollection != nullptr && Ancestor( &jpsiMeson , motherInPrunedCollection)){
	  if (dau.pdgId() == -13)     {
	    foundJpsiMum = true; 
	    //mu1_idx = findIdx<const pat::PackedGenParticle* >(&dau, pruned);
	    mu1_idx = findIdx<const reco::Candidate * >(jpsiMeson.daughter(0), pruned);
	    if(debug) std::cout<<"mu1 indx "<<mu1_idx<<std::endl;
	  }
	  else if (dau.pdgId() == 13) {
	    foundJpsiMup = true; 
	    //mu2_idx = findIdx<const pat::PackedGenParticle* >(&dau, pruned);
	    mu2_idx = findIdx<const reco::Candidate * >(jpsiMeson.daughter(1), pruned);
	    if(debug) std::cout<<"mu2 indx "<<mu2_idx<<std::endl;
	  }  
	}
      }
      if(debug) std::cout<<"Is there a plus and a minus mu?"<<foundJpsiMup<<" "<<foundJpsiMum<<std::endl;
      if (!(foundJpsiMup && foundJpsiMum)) continue;    
      if(debug) std::cout<<"YES"<<std::endl;

      // Not sure if it's the best way to do this
      // define a compositecandidate jpsiCand, such that ic an create my table easier (I don't know how to create a table from gen)
      

      // find jpsi original ancestor -> To Save!
      if(debug) std::cout<<"Number of mothers of the jpsi "<<jpsiMeson.numberOfMothers()<<std::endl;
      if ( jpsiMeson.numberOfMothers()>1 && debug) std::cout << "number of jpsi mom: " << jpsiMeson.numberOfMothers() << std::endl;
            
      for (size_t imom = 0; imom < jpsiMeson.numberOfMothers(); imom++){
	if(debug) std::cout<<"jpsimother"<<jpsiMeson.mother(imom)->pdgId()<<std::endl;
	const reco::Candidate * candMom = checkmom( jpsiMeson.mother(imom)) ;
	if(debug) std::cout<<"candMom"<<candMom->pdgId()<<std::endl;
	
	if (candMom == nullptr) continue;
	if (candMom != nullptr && Ancestor(candMom, &jpsiMeson )){
	  if ( candMom->mother(0) != nullptr) {

	    if(debug) std::cout<<"We found a mom "<<candMom->pdgId()<<std::endl;

	    jpsi_ancestor_idx = findIdx<const reco::Candidate*>(candMom, pruned);
	    if(debug) std::cout<<"The mom has idx "<<jpsi_ancestor_idx<<std::endl;
	    break;
	  }// if ( candMom->mother(0) != nullptr)
	}//if (candMom != nullptr && Ancestor(candMom, &jpsiMeson ))
      }//for (size_t imom = 0; imom < jpsiMeson.numberOfMothers(); imom++)
      
      std::vector<unsigned int> mu3_idx;
      std::vector<unsigned int> mu3_ancestor_idx;
      //find other Muon
      unsigned int number_of_other_muons = 0;
      for (unsigned int imuon = 0; imuon < pruned->size(); imuon++) {
	const reco::GenParticle &muon3 = pruned->at(imuon);
	if (abs(muon3.pdgId()) == 13 ){
	  if (!Ancestor( &jpsiMeson , &muon3)){
	    number_of_other_muons +=1;
	    if(debug) std::cout<<"Muon NOT from jpsi "<<muon3.pdgId()<<std::endl;
	    // save index for third muon
	    mu3_idx.push_back(imuon);
	    //find muon ancestor
	    if(debug) std::cout<<"Number of mothers of the muon "<<muon3.numberOfMothers()<<std::endl;
	    if ( muon3.numberOfMothers()>1 && debug) std::cout << "number of muon3 moms: " << muon3.numberOfMothers() << std::endl;
            
	    for (size_t imom = 0; imom < muon3.numberOfMothers(); imom++){
	      if(debug) std::cout<<"muon3mother"<<muon3.mother(imom)->pdgId()<<std::endl;
	      const reco::Candidate * candMom = checkmom( muon3.mother(imom)) ;
	      if(debug) std::cout<<"candMom"<<candMom->pdgId()<<std::endl;
	
	      if (candMom == nullptr) continue;
	      if (candMom != nullptr && Ancestor(candMom, &muon3 )){
		if ( candMom->mother(0) != nullptr) {

		  if(debug) std::cout<<"We found a mom "<<candMom->pdgId()<<std::endl;

		  mu3_ancestor_idx.push_back(findIdx<const reco::Candidate*>(candMom, pruned));
		  if(debug) std::cout<<"The mom has idx, "<<findIdx<const reco::Candidate*>(candMom, pruned)<<std::endl;
		  break;
		}// if ( candMom->mother(0) != nullptr)
	      }//if (candMom != nullptr && Ancestor(candMom, &jpsiMeson ))
	    }//for (size_t imom = 0; imom < muon3.numberOfMothers(); imom++)

	  } //(!Ancestor( &jpsiMeson , &muon3))
	}
      } // for muons in pruned
      if (debug) std::cout<<"number of other muons "<<number_of_other_muons<<std::endl;
      // for each third muon I save a different candidate, to be checked later
      for(unsigned int n=0;n< number_of_other_muons;n++){
	      pat::CompositeCandidate jpsiCand;
	      jpsiCand.setP4(jpsiMeson.p4());
	      jpsiCand.setCharge(jpsiMeson.charge());
	      // So i can find back which jpsi I am using as a candidate
	      jpsiCand.addUserInt("mu1_idx", mu1_idx);      
	      jpsiCand.addUserInt("mu2_idx", mu2_idx);
	      jpsiCand.addUserInt("jpsi_idx", ipart);
	      jpsiCand.addUserInt("jpsi_ancestor_idx", jpsi_ancestor_idx);
	      jpsiCand.addUserInt("mu3_idx", mu3_idx[n]);
	      jpsiCand.addUserInt("mu3_ancestor_idx", mu3_ancestor_idx[n]);
	      jpsi_out->push_back(jpsiCand);
	      if(debug) std::cout<<"IDX "<<mu1_idx<<std::endl;
	      if(debug) std::cout<<"IDX "<<jpsi_ancestor_idx<<std::endl;
	      if(debug) std::cout<<"IDX "<<mu3_idx[n]<<std::endl;
	      if(debug) std::cout<<"IDX "<<mu3_ancestor_idx[n]<<std::endl;
	      if(debug) std::cout<<jpsiCand.pt()<<std::endl;
      }

    }// if jpsi == 443
  }        
  if(debug) std::cout<<"number of events "<<jpsi_out->size()<<std::endl;
  /*
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
  */
  //export the jpsi_out vector
  if(debug) std::cout<<jpsi_out->size()<<std::endl;
  iEvent.put(std::move(jpsi_out), "AncestorFlags");
}

DEFINE_FWK_MODULE(HighMassLowMassbkgFlagsProducer);

// Class to include the gen infos of Bc and the daughters, to have a proper computation of hammer weights:
// otherwise if ht ecandidate is recostructed badly I don't have the right hammer weight

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "helper.h"

using namespace std;

constexpr bool debug = false;

class BcGenInfoTableProducer : public edm::EDProducer {
public:
  typedef std::vector<reco::GenParticle> GenParticleCollection;

  explicit BcGenInfoTableProducer(const edm::ParameterSet &iConfig);
  ~BcGenInfoTableProducer() override {};

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);
  //GEN
  const edm::EDGetTokenT<reco::GenParticleCollection> srcToken_; 
};

BcGenInfoTableProducer::BcGenInfoTableProducer(const edm::ParameterSet &iConfig):
  srcToken_(consumes<GenParticleCollection>(iConfig.getParameter<edm::InputTag>("srcGen")))
{
  produces<nanoaod::FlatTable>();
}
void BcGenInfoTableProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // Inputs


  //GEN  
  float bc_gen_pt = 0.;
  float bc_gen_eta = 0.;
  float bc_gen_phi = 0.;
  float bc_gen_mass = 0.;
  float bc_gen_vx = 0.;
  float bc_gen_vy = 0.;
  float bc_gen_vz = 0.;
  float mu1_gen_pt = 0.;
  float mu1_gen_eta = 0.;
  float mu1_gen_phi = 0.;
  float mu1_gen_mass = 0.;
  float mu2_gen_pt = 0.;
  float mu2_gen_eta = 0.;
  float mu2_gen_phi = 0.;
  float mu2_gen_mass = 0.;
  float mu3_gen_pt = 0.;
  float mu3_gen_eta = 0.;
  float mu3_gen_phi = 0.;
  float mu3_gen_mass = 0.;
  float tau_gen_pt = 0.;
  float tau_gen_eta = 0.;
  float tau_gen_phi = 0.;
  float tau_gen_mass = 0.;
  float jpsi_gen_pt = 0.;
  float jpsi_gen_eta = 0.;
  float jpsi_gen_phi = 0.;
  float jpsi_gen_mass = 0.;
  float jpsi_gen_vx = 0.;
  float jpsi_gen_vy = 0.;
  float jpsi_gen_vz = 0.;

  //GEN
  if(debug) std::cout<<"In the bc gen code..."<<std::endl;
  edm::Handle<GenParticleCollection> src;
  iEvent.getByToken(srcToken_, src);
  const size_t n = src->size();
  std::vector<const reco::Candidate*> final_daus;
  //std::vector<reco::Candidate> final_daus_part;

  if(debug) std::cout<<"Number of gen particles: "<<n<<std::endl;
  
  for(unsigned int  i = 0; i < n; ++i) {  //loop on gen particles
    const reco::GenParticle & gen = (*src)[i];
    const reco::Candidate* daughter; 
    const reco::Candidate* the_b;
    int is_doublemu = 0;
    int is_b = 0;
    final_daus.clear();
  
    if(abs(gen.pdgId()) == 443){  // looking for jpsi      
      if(debug) std::cout<<"There is a jpsi and she has "<<gen.numberOfDaughters()<<" daughters"<<std::endl;
      for(unsigned int dau = 0; dau < gen.numberOfDaughters(); dau++){  //loop of jpsi daughters
        if(debug) std::cout<<"Jpsi daughter: "<<gen.daughter(dau)->pdgId()<<std::endl;
        is_b = 0;
        if (abs(gen.daughter(dau)->pdgId())==13){
          is_doublemu += 1;
        }
      } //end loop on daughters
      if(is_doublemu>=2){  // jpsi -> mu mu
        if(debug) std::cout<<"The daughters are muons"<<std::endl;
        the_b = gen.mother(0); // jpsi mother
        if(abs(the_b->pdgId()) == 541){ //Bc->jpsi
          if(debug) std::cout<<"The direct mother is a Bc"<<std::endl;
          is_b = 1;
        }  
        else if(the_b->numberOfMothers() > 0)
        {
          the_b = gen.mother(0)->mother(0); // Bc->X->jpsi
          if(abs(the_b->pdgId()) == 541 ){
            if(debug) std::cout<<"The non direct mother is a Bc"<<std::endl;
            is_b = 1;
          }
        }
        if(is_b == 1){
          if(debug) std::cout<<"The Bc has "<<the_b->numberOfDaughters()<<"daughters"<<std::endl;

          for(unsigned int bdau=0; bdau < the_b->numberOfDaughters(); bdau ++){
            daughter = the_b->daughter(bdau);
            if(abs(daughter->pdgId())!= 541 and abs(daughter->pdgId())!= 22){    //not gamma
              final_daus.push_back(daughter);
              //final_daus_part.push_back(daughter);
              //        cout<<daughter->pdgId()<<endl;
              if(debug) std::cout<<"The Bc daughters are "<< daughter->pdgId()<<std::endl;
            }
          }
          
          std::sort(final_daus.begin(), final_daus.end(), [](const reco::Candidate* a, const reco::Candidate* b){ return abs(a->pdgId()) < abs(b->pdgId()); });  //sort the pdgIds of the daughters
          /*
          for(unsigned int item=0; item< final_daus.size(); item ++){
            if(debug) std::cout<<final_daus[item]<<std::endl;
            if(item == final_daus.size() -1) if(debug) std::cout<<" "<<std::endl;
          }
          */

	  mu3_gen_pt = 0.;
	  mu3_gen_eta = 0.;
	  mu3_gen_phi = 0.;
	  mu3_gen_mass = 0.;
	  tau_gen_pt = 0.;
	  tau_gen_eta = 0.;
	  tau_gen_phi = 0.;
	  tau_gen_mass = 0.;
          
	  bc_gen_pt = the_b->pt();
	  bc_gen_eta = the_b->eta();
	  bc_gen_phi = the_b->phi();
	  bc_gen_mass = the_b->mass();

	  bc_gen_vx = the_b->vx();
	  bc_gen_vy = the_b->vy();
	  bc_gen_vz = the_b->vz();

	  if(debug) std::cout<<"The Bc production vertex is "<<the_b->vx()<<std::endl;
	  jpsi_gen_pt = gen.pt();
	  jpsi_gen_eta = gen.eta();
	  jpsi_gen_phi = gen.phi();
	  jpsi_gen_mass = gen.mass();

	  jpsi_gen_vx = gen.vx();
	  jpsi_gen_vy = gen.vy();
	  jpsi_gen_vz = gen.vz();


	  // Bc. jpsi+mu
          if(abs(final_daus[0]->pdgId()) == 13){  //muon
            if(abs(final_daus[1]->pdgId()) == 14){ //muonic neutrino
              if(abs(final_daus[2]->pdgId()) == 443){  //jpsi
		mu3_gen_pt = final_daus[0]->pt();
		mu3_gen_eta = final_daus[0]->eta();
		mu3_gen_phi = final_daus[0]->phi();
		mu3_gen_mass = final_daus[0]->mass();

	      }
	    }
          }
	  // Bc-> jpsi+tau
          else if(abs(final_daus[0]->pdgId()) == 15){ //tau
            if (abs(final_daus[1]->pdgId()) == 16){ //tauonic neutrino
              if(abs(final_daus[2]->pdgId()) == 443) {//jpsi
		tau_gen_pt = final_daus[0]->pt();
		tau_gen_eta = final_daus[0]->eta();
		tau_gen_phi = final_daus[0]->phi();
		tau_gen_mass = final_daus[0]->mass();
	      }
            }
          }
          else{
	    if(debug) std::cout<<"neither tau or mu, somthing else"<<std::endl;
          }
        } // if(is_b == 1)
      } //if(is_doublemu>=2)
    }//if(abs(gen.pdgId()) == 443)
  }//for(unsigned int  i = 0; i < n; ++i)


  auto tab = std::make_unique<nanoaod::FlatTable>(1,"",true);
  tab->addColumnValue<float>("BcGenInfo_bc_gen_pt", bc_gen_pt, "bc gen information", nanoaod::FlatTable::FloatColumn);
  tab->addColumnValue<float>("BcGenInfo_bc_gen_eta", bc_gen_eta, "bc gen information", nanoaod::FlatTable::FloatColumn);
  tab->addColumnValue<float>("BcGenInfo_bc_gen_phi", bc_gen_phi, "bc gen information", nanoaod::FlatTable::FloatColumn);
  tab->addColumnValue<float>("BcGenInfo_bc_gen_mass", bc_gen_mass, "bc gen information", nanoaod::FlatTable::FloatColumn);
  tab->addColumnValue<float>("BcGenInfo_mu1_gen_pt", mu1_gen_pt, "bc gen information", nanoaod::FlatTable::FloatColumn);
  tab->addColumnValue<float>("BcGenInfo_mu1_gen_eta", mu1_gen_eta, "bc gen information", nanoaod::FlatTable::FloatColumn);
  tab->addColumnValue<float>("BcGenInfo_mu1_gen_phi", mu1_gen_phi, "bc gen information", nanoaod::FlatTable::FloatColumn);
  tab->addColumnValue<float>("BcGenInfo_mu1_gen_mass", mu1_gen_mass, "bc gen information", nanoaod::FlatTable::FloatColumn);
  tab->addColumnValue<float>("BcGenInfo_mu2_gen_pt", mu2_gen_pt, "bc gen information", nanoaod::FlatTable::FloatColumn);
  tab->addColumnValue<float>("BcGenInfo_mu2_gen_eta", mu2_gen_eta, "bc gen information", nanoaod::FlatTable::FloatColumn);
  tab->addColumnValue<float>("BcGenInfo_mu2_gen_phi", mu2_gen_phi, "bc gen information", nanoaod::FlatTable::FloatColumn);
  tab->addColumnValue<float>("BcGenInfo_mu2_gen_mass", mu2_gen_mass, "bc gen information", nanoaod::FlatTable::FloatColumn);
  tab->addColumnValue<float>("BcGenInfo_mu3_gen_pt", mu3_gen_pt, "bc gen information", nanoaod::FlatTable::FloatColumn);
  tab->addColumnValue<float>("BcGenInfo_mu3_gen_eta", mu3_gen_eta, "bc gen information", nanoaod::FlatTable::FloatColumn);
  tab->addColumnValue<float>("BcGenInfo_mu3_gen_phi", mu3_gen_phi, "bc gen information", nanoaod::FlatTable::FloatColumn);
  tab->addColumnValue<float>("BcGenInfo_mu3_gen_mass", mu3_gen_mass, "bc gen information", nanoaod::FlatTable::FloatColumn);
  tab->addColumnValue<float>("BcGenInfo_tau_gen_pt", tau_gen_pt, "bc gen information", nanoaod::FlatTable::FloatColumn);
  tab->addColumnValue<float>("BcGenInfo_tau_gen_eta", tau_gen_eta, "bc gen information", nanoaod::FlatTable::FloatColumn);
  tab->addColumnValue<float>("BcGenInfo_tau_gen_phi", tau_gen_phi, "bc gen information", nanoaod::FlatTable::FloatColumn);
  tab->addColumnValue<float>("BcGenInfo_tau_gen_mass", tau_gen_mass, "bc gen information", nanoaod::FlatTable::FloatColumn);
  tab->addColumnValue<float>("BcGenInfo_jpsi_gen_pt", jpsi_gen_pt, "bc gen information", nanoaod::FlatTable::FloatColumn);
  tab->addColumnValue<float>("BcGenInfo_jpsi_gen_eta", jpsi_gen_eta, "bc gen information", nanoaod::FlatTable::FloatColumn);
  tab->addColumnValue<float>("BcGenInfo_jpsi_gen_phi", jpsi_gen_phi, "bc gen information", nanoaod::FlatTable::FloatColumn);
  tab->addColumnValue<float>("BcGenInfo_jpsi_gen_mass", jpsi_gen_mass, "bc gen information", nanoaod::FlatTable::FloatColumn);
  tab->addColumnValue<float>("BcGenInfo_bc_gen_vx", bc_gen_vx, "bc gen information", nanoaod::FlatTable::FloatColumn);
  tab->addColumnValue<float>("BcGenInfo_bc_gen_vy", bc_gen_vy, "bc gen information", nanoaod::FlatTable::FloatColumn);
  tab->addColumnValue<float>("BcGenInfo_bc_gen_vz", bc_gen_vz, "bc gen information", nanoaod::FlatTable::FloatColumn);
  tab->addColumnValue<float>("BcGenInfo_jpsi_gen_vx", jpsi_gen_vx, "bc gen information", nanoaod::FlatTable::FloatColumn);
  tab->addColumnValue<float>("BcGenInfo_jpsi_gen_vy", jpsi_gen_vy, "bc gen information", nanoaod::FlatTable::FloatColumn);
  tab->addColumnValue<float>("BcGenInfo_jpsi_gen_vz", jpsi_gen_vz, "bc gen information", nanoaod::FlatTable::FloatColumn);
  iEvent.put(std::move(tab));
}

DEFINE_FWK_MODULE(BcGenInfoTableProducer);

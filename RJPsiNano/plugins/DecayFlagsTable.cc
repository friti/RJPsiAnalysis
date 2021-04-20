// Class to include the flags of the different decays
// and the the truth match

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

class DecayFlagsTableProducer : public edm::EDProducer {
public:
  typedef std::vector<reco::GenParticle> GenParticleCollection;

  explicit DecayFlagsTableProducer(const edm::ParameterSet &iConfig);
  ~DecayFlagsTableProducer() override {};

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);
  //GEN
  const edm::EDGetTokenT<reco::GenParticleCollection> srcToken_; 
};

DecayFlagsTableProducer::DecayFlagsTableProducer(const edm::ParameterSet &iConfig):
  srcToken_(consumes<GenParticleCollection>(iConfig.getParameter<edm::InputTag>("srcGen")))
{
  produces<nanoaod::FlatTable>();
}
void DecayFlagsTableProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // Inputs


  //GEN
  
  uint8_t flag_jpsi_mu = 0;
  uint8_t flag_psi2s_mu = 0;
  uint8_t flag_chic0_mu = 0;
  uint8_t flag_chic1_mu = 0;
  uint8_t flag_chic2_mu = 0;
  uint8_t flag_hc_mu = 0;
  uint8_t flag_jpsi_tau = 0;
  uint8_t flag_psi2s_tau = 0;
  uint8_t flag_jpsi_pi = 0;
  uint8_t flag_jpsi_3pi = 0;
  uint8_t flag_jpsi_hc = 0;
  uint8_t flag_error = 0;
  //GEN
  if(debug) std::cout<<"In the flag code..."<<std::endl;
  edm::Handle<GenParticleCollection> src;
  iEvent.getByToken(srcToken_, src);
  const size_t n = src->size();
  std::vector<int> final_daus;

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
              final_daus.push_back(abs(daughter->pdgId()));
              //        cout<<daughter->pdgId()<<endl;
              if(debug) std::cout<<"The Bc daughters are "<< daughter->pdgId()<<std::endl;
            }
          }
          
          std::sort(final_daus.begin(), final_daus.end());  //sort the pdgIds of the daughters
          /*
          for(unsigned int item=0; item< final_daus.size(); item ++){
            if(debug) std::cout<<final_daus[item]<<std::endl;
            if(item == final_daus.size() -1) if(debug) std::cout<<" "<<std::endl;
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
            if (final_daus[1] == 211 && final_daus[2] ==211 && final_daus[3] == 443) flag_jpsi_3pi = 1;
          }
          else if ((final_daus[0] == 431 && final_daus[1] == 443) || (final_daus[0] == 433 && final_daus[1] == 443)) flag_jpsi_hc = 1;  
          else{
            flag_error = 1;
          }
        } // if(is_b == 1)
      } //if(is_doublemu>=2)
    }//if(abs(gen.pdgId()) == 443)
  }//for(unsigned int  i = 0; i < n; ++i)

  float weight=0;
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

  auto tab = std::make_unique<nanoaod::FlatTable>(1,"",true);

  tab->addColumnValue<uint8_t>("DecayFlag_is_jpsi_mu", flag_jpsi_mu, "decay flag", nanoaod::FlatTable::UInt8Column);
  tab->addColumnValue<uint8_t>("DecayFlag_is_psi2s_mu", flag_psi2s_mu, "decay flag", nanoaod::FlatTable::UInt8Column);
  tab->addColumnValue<uint8_t>("DecayFlag_is_chic0_mu", flag_chic0_mu, "decay flag", nanoaod::FlatTable::UInt8Column);
  tab->addColumnValue<uint8_t>("DecayFlag_is_chic1_mu", flag_chic1_mu, "decay flag", nanoaod::FlatTable::UInt8Column);
  tab->addColumnValue<uint8_t>("DecayFlag_is_chic2_mu", flag_chic2_mu, "decay flag", nanoaod::FlatTable::UInt8Column);
  tab->addColumnValue<uint8_t>("DecayFlag_is_hc_mu", flag_hc_mu, "decay flag", nanoaod::FlatTable::UInt8Column);
  tab->addColumnValue<uint8_t>("DecayFlag_is_jpsi_tau", flag_jpsi_tau, "decay flag", nanoaod::FlatTable::UInt8Column);
  tab->addColumnValue<uint8_t>("DecayFlag_is_psi2s_tau", flag_psi2s_tau, "decay flag", nanoaod::FlatTable::UInt8Column);
  tab->addColumnValue<uint8_t>("DecayFlag_is_jpsi_pi", flag_jpsi_pi, "decay flag", nanoaod::FlatTable::UInt8Column);
  tab->addColumnValue<uint8_t>("DecayFlag_is_jpsi_3pi", flag_jpsi_3pi, "decay flag", nanoaod::FlatTable::UInt8Column);
  tab->addColumnValue<uint8_t>("DecayFlag_is_jpsi_hc", flag_jpsi_hc, "decay flag", nanoaod::FlatTable::UInt8Column);
  tab->addColumnValue<uint8_t>("DecayFlag_is_error", flag_error, "decay flag", nanoaod::FlatTable::UInt8Column);
  tab->addColumnValue<float>("DecayFlag_weight", weight, "decay weight", nanoaod::FlatTable::FloatColumn, 10);
  iEvent.put(std::move(tab));
}

DEFINE_FWK_MODULE(DecayFlagsTableProducer);

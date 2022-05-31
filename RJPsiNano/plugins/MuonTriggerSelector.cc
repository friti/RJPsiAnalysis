// class to produce 2 pat::MuonCollections
// one matched to the Park triggers
// another fitered wrt the Park triggers


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

#include "DataFormats/PatCandidates/interface/Muon.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerAlgorithm.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include <TLorentzVector.h>
#include "helper.h"

using namespace std;

constexpr bool debug = false;
constexpr bool debugTrg = false;

class MuonTriggerSelector : public edm::EDProducer {
    
public:
    
    explicit MuonTriggerSelector(const edm::ParameterSet &iConfig);
    
    ~MuonTriggerSelector() override {};
    
    
private:

    virtual void produce(edm::Event&, const edm::EventSetup&);

    edm::EDGetTokenT<std::vector<pat::Muon>> muonSrc_;
    edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
    edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggerObjects_;
    edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
    edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;

    //for trigger match
    const double maxdR_;

    //for filter wrt trigger
    const double dzTrg_cleaning_; // selects primary vertex

    const double ptMin_;          // min pT in all muons for B candidates
    const double absEtaMax_;      //max eta ""
    const bool softMuonsOnly_;    //cuts muons without soft ID
};


MuonTriggerSelector::MuonTriggerSelector(const edm::ParameterSet &iConfig):
  muonSrc_( consumes<std::vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "muonCollection" ) ) ),
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("objects"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  vertexSrc_( consumes<reco::VertexCollection> ( iConfig.getParameter<edm::InputTag>( "vertexCollection" ) ) ), 
  maxdR_(iConfig.getParameter<double>("maxdR_matching")),
  dzTrg_cleaning_(iConfig.getParameter<double>("dzForCleaning_wrtTrgMuon")),
  ptMin_(iConfig.getParameter<double>("ptMin")),
  absEtaMax_(iConfig.getParameter<double>("absEtaMax")),
  softMuonsOnly_(iConfig.getParameter<bool>("softMuonsOnly"))
{
  // produce 2 collections: trgMuons (tags) and SelectedMuons (probes & tags if survive preselection cuts)
    produces<pat::MuonCollection>("trgMuons"); 
    produces<pat::MuonCollection>("SelectedMuons");
    produces<TransientTrackCollection>("SelectedTransientMuons");  
}



void MuonTriggerSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    
  edm::ESHandle<MagneticField> bFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

  edm::Handle<reco::VertexCollection> vertexHandle;
  iEvent.getByToken(vertexSrc_, vertexHandle);
  //const reco::Vertex & PV = vertexHandle->front();

  if(debug) std::cout << " MuonTriggerSelector::produce " << std::endl;

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

  std::vector<pat::TriggerObjectStandAlone> triggeringMuons;

  //taken from https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2016#Trigger
  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);
  if(debug) std::cout << "\n TRIGGER OBJECTS " << std::endl;

  edm::Handle<std::vector<pat::Muon>> muons;
  iEvent.getByToken(muonSrc_, muons);

  std::unique_ptr<pat::MuonCollection>      trgmuons_out   ( new pat::MuonCollection );
  std::unique_ptr<pat::MuonCollection>      muons_out      ( new pat::MuonCollection );
  std::unique_ptr<TransientTrackCollection> trans_muons_out( new TransientTrackCollection );

  // Getting the indexes of the HLT paths
  unsigned int index_dimuon01 = names.triggerIndex("HLT_Dimuon0_Jpsi3p5_Muon2_v5");
  unsigned int index_dimuon02 = names.triggerIndex("HLT_Dimuon0_Jpsi3p5_Muon2_v6");
  unsigned int index_jpsiTrk1 = names.triggerIndex("HLT_DoubleMu4_JpsiTrk_Displaced_v14");
  unsigned int index_jpsiTrk2 = names.triggerIndex("HLT_DoubleMu4_JpsiTrk_Displaced_v15");
  unsigned int index_jpsiTrk3 = names.triggerIndex("HLT_DoubleMu4_PsiPrimeTrk_Displaced_v14");
  unsigned int index_jpsiTrk4 = names.triggerIndex("HLT_DoubleMu4_PsiPrimeTrk_Displaced_v15");
  unsigned int index_jpsiTrk5 = names.triggerIndex("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v14");
  unsigned int index_jpsiTrk6 = names.triggerIndex("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v15");

  
  
  bool pass_dimuon01_path= false;
  bool pass_dimuon02_path= false;
  bool pass_jpsiTrk1_path = false;
  bool pass_jpsiTrk2_path = false;
  bool pass_jpsiTrk3_path = false;
  bool pass_jpsiTrk4_path = false;
  bool pass_jpsiTrk5_path = false;
  bool pass_jpsiTrk6_path = false;
  if(index_dimuon01 != triggerBits->size()) 
    pass_dimuon01_path = triggerBits->accept(index_dimuon01);
  if(index_dimuon02 != triggerBits->size()) 
    pass_dimuon02_path = triggerBits->accept(index_dimuon02);
  if(index_jpsiTrk1 != triggerBits->size()) 
    pass_jpsiTrk1_path = triggerBits->accept(index_jpsiTrk1);
  if(index_jpsiTrk2 != triggerBits->size()) 
    pass_jpsiTrk2_path = triggerBits->accept(index_jpsiTrk2);
  if(index_jpsiTrk3 != triggerBits->size()) 
    pass_jpsiTrk3_path = triggerBits->accept(index_jpsiTrk3);
  if(index_jpsiTrk4 != triggerBits->size()) 
    pass_jpsiTrk4_path = triggerBits->accept(index_jpsiTrk4);
  if(index_jpsiTrk5 != triggerBits->size()) 
    pass_jpsiTrk5_path = triggerBits->accept(index_jpsiTrk5);
  if(index_jpsiTrk6 != triggerBits->size()) 
    pass_jpsiTrk6_path = triggerBits->accept(index_jpsiTrk6);

  //  if(debug) std::cout << "pass_dimuuon0 " << pass_dimuon0_path<<" pass_trk "<<pass_jpsiTrk_path<<std::endl;
  //if(pass_dimuon0_path) std::cout << "pass_dimuuon0" << std::endl;
  std::vector<bool> jpsiFromMuon_fromDimuon0_flags;
  std::vector<bool> jpsiFromMuon_fromJpsiTrk_flags;
  std::vector<bool> jpsiFromMuon_fromJpsiTrk_PsiPrime_flags;
  std::vector<bool> jpsiFromMuon_fromJpsiTrk_NonResonant_flags;
  std::vector<bool> dimuon0Flags;
  std::vector<bool> jpsiTrkFlags;
  std::vector<bool> jpsiTrk_PsiPrimeFlags;
  std::vector<bool> jpsiTrk_NonResonantFlags;

  if(pass_dimuon01_path || pass_dimuon02_path || pass_jpsiTrk1_path || pass_jpsiTrk2_path || pass_jpsiTrk3_path || pass_jpsiTrk4_path || pass_jpsiTrk5_path || pass_jpsiTrk6_path) {  
    //if(pass_dimuon01_path || pass_jpsiTrk1_path || pass_jpsiTrk2_path || pass_jpsiTrk3_path || pass_jpsiTrk4_path || pass_jpsiTrk5_path || pass_jpsiTrk6_path) {  
    for (pat::TriggerObjectStandAlone obj : *triggerObjects) 
    { // note: not "const &" since we want to call unpackPathNames
      obj.unpackFilterLabels(iEvent, *triggerBits);
      obj.unpackPathNames(names);

      bool muonFromJpsi_fromDimuon0Path = false;
      bool muonFromJpsi_fromJpsiTrkPath = false;
      bool muonFromJpsi_fromJpsiTrk_PsiPrimePath = false;
      bool muonFromJpsi_fromJpsiTrk_NonResonantPath = false;
      bool dimuon0_seed = false;
      bool jpsitrk_seed = false;
      bool jpsitrk_PsiPrime_seed = false;
      bool jpsitrk_NonResonant_seed = false;

      //if(pass_jpsiTrk5_path || pass_jpsiTrk6_path) 
      //  if(debugTrg) std::cout << "pass_nonREsonant path" << std::endl;

      if(obj.hasFilterLabel("hltVertexmumuFilterJpsiMuon3p5") )
        muonFromJpsi_fromDimuon0Path = true;
      if(obj.hasFilterLabel("hltTripleMuL3PreFiltered222")){
        dimuon0_seed = true;
      }

      if(obj.hasFilterLabel("hltDisplacedmumuFilterDoubleMu4Jpsi") )
        muonFromJpsi_fromJpsiTrkPath = true;
      if(obj.hasFilterLabel("hltJpsiTkVertexFilter")){
        jpsitrk_seed = true;
      }

      if(obj.hasFilterLabel("hltDisplacedmumuFilterDoubleMu4PsiPrime") )
      {
        muonFromJpsi_fromJpsiTrk_PsiPrimePath = true;
      }
      if(obj.hasFilterLabel("hltPsiPrimeTkVertexFilter")){
        jpsitrk_PsiPrime_seed = true;
      }

      if(obj.hasFilterLabel("hltDisplacedmumuFilterDoubleMu4LowMassNonResonant") )
      {
        muonFromJpsi_fromJpsiTrk_NonResonantPath = true;
        //if(debugTrg) std::cout << "hltDisplacedmumuFilterDoubleMu4NonResonant" << std::endl;
      }
      if(obj.hasFilterLabel("hltLowMassNonResonantTkVertexFilter")){
        jpsitrk_NonResonant_seed = true;
      }
      
      //for each triggered muon I know which trigger it passes
      dimuon0Flags.push_back(dimuon0_seed); //unpaired muon passes the dimuon0
      jpsiTrkFlags.push_back(jpsitrk_seed); //unpaired muon passes the jpsitrk
      jpsiTrk_PsiPrimeFlags.push_back(jpsitrk_PsiPrime_seed); //unpaired muon passes the psiprimetrk
      jpsiTrk_NonResonantFlags.push_back(jpsitrk_NonResonant_seed); //unpaired muon passes the nonResonantrk

      jpsiFromMuon_fromDimuon0_flags.push_back(muonFromJpsi_fromDimuon0Path); //the muon is from the jpsi
      jpsiFromMuon_fromJpsiTrk_flags.push_back(muonFromJpsi_fromJpsiTrkPath); //the muon is from the jpsi
      jpsiFromMuon_fromJpsiTrk_PsiPrime_flags.push_back(muonFromJpsi_fromJpsiTrk_PsiPrimePath); //the muon is from the jpsi
      jpsiFromMuon_fromJpsiTrk_NonResonant_flags.push_back(muonFromJpsi_fromJpsiTrk_NonResonantPath); //the muon is from the jpsi
      triggeringMuons.push_back(obj);
      if(debug){ 

        std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
      	// Print trigger object collection and type
      	std::cout << "\t   Collection: " << obj.collection() << std::endl;
      	std::cout<< "\t muonFromJpsi_fromDimuon0Path "<<muonFromJpsi_fromDimuon0Path<<std::endl;
      	std::cout<<"\t muonFromJpsi_fromJpsiTrkPath "<<muonFromJpsi_fromJpsiTrkPath<<std::endl;
      	std::cout<<"\t pass_dimuon0_path"<<pass_dimuon01_path<<std::endl;
      	std::cout<<"\t dimuon0_seed"<<dimuon0_seed<<std::endl;
      	std::cout<<"\t jpsitrk_seed"<<jpsitrk_seed<<std::endl;
      	std::cout<<"\t jpsitrk_PsiPrime_seed"<<jpsitrk_PsiPrime_seed<<std::endl;
      	std::cout<<"\t jpsitrk_NonResonant_seed"<<jpsitrk_NonResonant_seed<<std::endl;
	
      }
      
    }//trigger objects
  }

  if(debug)
  {
    std::cout << "\n total n of triggering muons = " << triggeringMuons.size() << std::endl;
    for(auto ij : triggeringMuons)
    {
	    std::cout << " >>> components (pt, eta, phi) = " << ij.pt() << " " << ij.eta() << " " << ij.phi() << std::endl;
    }
  }
  //now check for reco muons matched to triggering muons
  std::vector<int> muonIsTrigger(muons->size(), 0);
  std::vector<int> muonIsFromJpsi_dimuon0Path(muons->size(), 0);
  std::vector<int> muonIsFromJpsi_jpsiTrkPath(muons->size(), 0);
  std::vector<int> muonIsFromJpsi_jpsiTrk_PsiPrimePath(muons->size(), 0);
  std::vector<int> muonIsFromJpsi_jpsiTrk_NonResonantPath(muons->size(), 0);
  std::vector<int> muonIsDimuon0Trg(muons->size(), 0);
  std::vector<int> muonIsJpsiTrkTrg(muons->size(), 0);
  std::vector<int> muonIsJpsiTrk_PsiPrimeTrg(muons->size(), 0);
  std::vector<int> muonIsJpsiTrk_NonResonantTrg(muons->size(), 0);

  if(debugTrg) std::cout << "Event -------------------- " << std::endl;
  
  for(const pat::Muon & muon : *muons)
  {
    //this is for triggering muon not really need to be configurable
    unsigned int iMuo(&muon - &(muons->at(0)) );
    //if(!(muon.isLooseMuon() && muon.isSoftMuon(PV))) continue;
    bool isMuonMatchedToDimuon0Path = !(muon.triggerObjectMatchByPath("HLT_Dimuon0_Jpsi3p5_Muon2_v5")==nullptr) || !(muon.triggerObjectMatchByPath("HLT_Dimuon0_Jpsi3p5_Muon2_v6")==nullptr);
    bool isMuonMatchedToJpsiTrkPath = !(muon.triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrk_Displaced_v14")==nullptr) || !(muon.triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrk_Displaced_v15")==nullptr);
    bool isMuonMatchedToJpsiTrk_PsiPrimePath = !(muon.triggerObjectMatchByPath("HLT_DoubleMu4_PsiPrimeTrk_Displaced_v14")==nullptr) || !(muon.triggerObjectMatchByPath("HLT_DoubleMu4_PsiPrimeTrk_Displaced_v15")==nullptr);
    bool isMuonMatchedToJpsiTrk_NonResonantPath = !(muon.triggerObjectMatchByPath("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v14")==nullptr) || ! (muon.triggerObjectMatchByPath("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v15")==nullptr);

    //if( isMuonMatchedToJpsiTrk_NonResonantPath ) std::cout << " isMuonMatchedToJpsiTrk_NonResonantPath-------------------- " << std::endl;


    float dRMuonMatchingDimuon0 = -1.;
    int recoMuonMatchingDimuon0_index = -1;
    int trgMuonMatchingDimuon0_index = -1;

    float dRMuonMatchingJpsiTrk = -1.;
    int recoMuonMatchingJpsiTrk_index = -1;
    int trgMuonMatchingJpsiTrk_index = -1;

    float dRMuonMatchingJpsiTrk_PsiPrime = -1.;
    int recoMuonMatchingJpsiTrk_PsiPrime_index = -1;
    int trgMuonMatchingJpsiTrk_PsiPrime_index = -1;

    float dRMuonMatchingJpsiTrk_NonResonant = -1.;
    int recoMuonMatchingJpsiTrk_NonResonant_index = -1;
    int trgMuonMatchingJpsiTrk_NonResonant_index = -1;

    for(unsigned int iTrg=0; iTrg<triggeringMuons.size(); ++iTrg)
    {
      if(!dimuon0Flags[iTrg] && !jpsiTrkFlags[iTrg] && !jpsiTrk_PsiPrimeFlags[iTrg] && !jpsiTrk_NonResonantFlags[iTrg]) continue;
      //it passes the dimuon0 trigger
      

      float dR = reco::deltaR(triggeringMuons[iTrg], muon);
      if(isMuonMatchedToDimuon0Path && dimuon0Flags[iTrg] && (dR < dRMuonMatchingDimuon0 || dRMuonMatchingDimuon0 == -1)  && dR < maxdR_)
      {
	      dRMuonMatchingDimuon0 = dR;
	      recoMuonMatchingDimuon0_index = iMuo;
	      trgMuonMatchingDimuon0_index = iTrg;
	      if(debug) std::cout<<"trigger dimuon0 "<< dimuon0Flags[iTrg]<<std::endl;
	      if(debug) std::cout << " dR = " << dR <<std::endl;
	      if(debug) std::cout << " HLT = " << triggeringMuons[iTrg].pt() << " " << triggeringMuons[iTrg].eta() << " " << triggeringMuons[iTrg].phi()          << std::endl;
	      if(debug) std::cout << " reco = " << muon.pt() << " " << muon.eta() << " " << muon.phi()          << std::endl;
      }
      if(isMuonMatchedToJpsiTrkPath && jpsiTrkFlags[iTrg] && (dR < dRMuonMatchingJpsiTrk || dRMuonMatchingJpsiTrk == -1)  && dR < maxdR_)
      {
	      dRMuonMatchingJpsiTrk = dR;
	      recoMuonMatchingJpsiTrk_index = iMuo;
	      trgMuonMatchingJpsiTrk_index = iTrg;
	      if(debug) std::cout<<"trigger jpsitrk "<< jpsiTrkFlags[iTrg]<<std::endl;
	      if(debug) std::cout << " dR = " << dR <<std::endl;
	      if(debug) std::cout << " HLT = " << triggeringMuons[iTrg].pt() << " " << triggeringMuons[iTrg].eta() << " " << triggeringMuons[iTrg].phi()          << std::endl;
	      if(debug) std::cout << " reco = " << muon.pt() << " " << muon.eta() << " " << muon.phi()          << std::endl;
      }
      if(isMuonMatchedToJpsiTrk_PsiPrimePath && jpsiTrk_PsiPrimeFlags[iTrg] && (dR < dRMuonMatchingJpsiTrk_PsiPrime || dRMuonMatchingJpsiTrk_PsiPrime == -1)  && dR < maxdR_)
      {
	      dRMuonMatchingJpsiTrk_PsiPrime = dR;
	      recoMuonMatchingJpsiTrk_PsiPrime_index = iMuo;
	      trgMuonMatchingJpsiTrk_PsiPrime_index = iTrg;
	      if(debug) std::cout<<"trigger jpsitrk "<< jpsiTrk_PsiPrimeFlags[iTrg]<<std::endl;
	      if(debug) std::cout << " dR = " << dR <<std::endl;
	      if(debug) std::cout << " HLT = " << triggeringMuons[iTrg].pt() << " " << triggeringMuons[iTrg].eta() << " " << triggeringMuons[iTrg].phi()          << std::endl;
	      if(debug) std::cout << " reco = " << muon.pt() << " " << muon.eta() << " " << muon.phi()          << std::endl;
      }
      if(isMuonMatchedToJpsiTrk_NonResonantPath && jpsiTrk_NonResonantFlags[iTrg] && (dR < dRMuonMatchingJpsiTrk_NonResonant || dRMuonMatchingJpsiTrk_NonResonant == -1)  && dR < maxdR_)
      {
	      dRMuonMatchingJpsiTrk_NonResonant = dR;
	      recoMuonMatchingJpsiTrk_NonResonant_index = iMuo;
	      trgMuonMatchingJpsiTrk_NonResonant_index = iTrg;
	      if(debug) std::cout<<"trigger jpsitrk "<< jpsiTrk_NonResonantFlags[iTrg]<<std::endl;
	      if(debug) std::cout << " dR = " << dR <<std::endl;
	      if(debug) std::cout << " HLT = " << triggeringMuons[iTrg].pt() << " " << triggeringMuons[iTrg].eta() << " " << triggeringMuons[iTrg].phi()          << std::endl;
	      if(debug) std::cout << " reco = " << muon.pt() << " " << muon.eta() << " " << muon.phi()          << std::endl;
      }

    }
    //save reco muon 
    if(recoMuonMatchingDimuon0_index != -1 || recoMuonMatchingJpsiTrk_index !=-1 || recoMuonMatchingJpsiTrk_PsiPrime_index !=-1 || recoMuonMatchingJpsiTrk_NonResonant_index !=-1)
    {
	    muonIsTrigger[iMuo] = 1;
	    pat::Muon recoTriggerMuonCand(muon);
	    recoTriggerMuonCand.addUserInt("trgMuonDimuon0_index", trgMuonMatchingDimuon0_index);
	    recoTriggerMuonCand.addUserInt("trgMuonJpsiTrk_index", trgMuonMatchingJpsiTrk_index);
	    recoTriggerMuonCand.addUserInt("trgMuonJpsiTrk_PsiPrime_index", trgMuonMatchingJpsiTrk_PsiPrime_index);
	    recoTriggerMuonCand.addUserInt("trgMuonJpsiTrk_NonResonant_index", trgMuonMatchingJpsiTrk_NonResonant_index);
	    trgmuons_out->emplace_back(recoTriggerMuonCand);
	    //keep track of original muon index for SelectedMuons collection

      if(recoMuonMatchingDimuon0_index != -1)
      {
	      muonIsFromJpsi_dimuon0Path[iMuo] = jpsiFromMuon_fromDimuon0_flags[trgMuonMatchingDimuon0_index];
	      muonIsDimuon0Trg[iMuo] = dimuon0Flags[trgMuonMatchingDimuon0_index];
      }
      else
      {
	      muonIsFromJpsi_dimuon0Path[iMuo] = 0;//jpsiFromMuon_fromDimuon0_flags[trgMuonMatchingDimuon0_index];
	      muonIsDimuon0Trg[iMuo] = 0;//dimuon0Flags[trgMuonMatchingDimuon0_index];
      }

      if(recoMuonMatchingJpsiTrk_index != -1)
      {
	      if(debug) std::cout<<"here it fills the vector for "<<muon.pt()<<" flag "<<jpsiFromMuon_fromJpsiTrk_flags[trgMuonMatchingJpsiTrk_index]<<std::endl;
	      muonIsFromJpsi_jpsiTrkPath[iMuo] = jpsiFromMuon_fromJpsiTrk_flags[trgMuonMatchingJpsiTrk_index];
	      muonIsJpsiTrkTrg[iMuo] = jpsiTrkFlags[trgMuonMatchingJpsiTrk_index];
      }
      else
      {
	      muonIsFromJpsi_jpsiTrkPath[iMuo] = 0;//jpsiFromMuon_fromJpsiTrk_flags[trgMuonMatchingJpsiTrk_index];
	      muonIsJpsiTrkTrg[iMuo] = 0;//jpsiTrkFlags[trgMuonMatchingJpsiTrk_index];
      }
      if(recoMuonMatchingJpsiTrk_PsiPrime_index != -1)
      {
	      muonIsFromJpsi_jpsiTrk_PsiPrimePath[iMuo] = jpsiFromMuon_fromJpsiTrk_PsiPrime_flags[trgMuonMatchingJpsiTrk_PsiPrime_index];
	      muonIsJpsiTrk_PsiPrimeTrg[iMuo] = jpsiTrk_PsiPrimeFlags[trgMuonMatchingJpsiTrk_PsiPrime_index];
      }
      else
      {
	      muonIsFromJpsi_jpsiTrk_PsiPrimePath[iMuo] = 0;//jpsiFromMuon_fromJpsiTrk_flags[trgMuonMatchingJpsiTrk_index];
	      muonIsJpsiTrk_PsiPrimeTrg[iMuo] = 0;//jpsiTrkFlags[trgMuonMatchingJpsiTrk_index];
      }

      if(debugTrg) if( isMuonMatchedToJpsiTrk_NonResonantPath ) std::cout << "recoMuonMatchingJpsiTrk_NonResonant_flag=====================  " << jpsiFromMuon_fromJpsiTrk_NonResonant_flags[recoMuonMatchingJpsiTrk_NonResonant_index] << std::endl;
      if(recoMuonMatchingJpsiTrk_NonResonant_index != -1)
      {
	      muonIsFromJpsi_jpsiTrk_NonResonantPath[iMuo] = jpsiFromMuon_fromJpsiTrk_NonResonant_flags[trgMuonMatchingJpsiTrk_NonResonant_index];
	      muonIsJpsiTrk_NonResonantTrg[iMuo] = jpsiTrk_NonResonantFlags[trgMuonMatchingJpsiTrk_NonResonant_index];
      }
      else
      {
	      muonIsFromJpsi_jpsiTrk_NonResonantPath[iMuo] = 0;//jpsiFromMuon_fromJpsiTrk_flags[trgMuonMatchingJpsiTrk_index];
	      muonIsJpsiTrk_NonResonantTrg[iMuo] = 0;//jpsiTrkFlags[trgMuonMatchingJpsiTrk_index];
      }
      /*
      if(debug){      
        std::cout << "HERE" << std::endl;
        std::cout << "trgMuonMatching_index: " << trgMuonMatching_index << std::endl;
        std::cout << "jpsiFromMuon_fromDimuon0_flags.push_back(muonFromJpsi_fromDimuon0Path)"<< jpsiFromMuon_fromDimuon0_flags[trgMuonMatchingDimuon0_index] << std::endl;
        std::cout << "jpsiFromMuon_fromJpsiTrk_flags.push_back(muonFromJpsi_fromDimuon0Path)"<< jpsiFromMuon_fromJpsiTrk_flags[trgMuonMatching_index] << std::endl;
        std::cout << "dimuon0Flags.push_back(pass_dimuon0_path)"<< dimuon0Flags[trgMuonMatchingDimuon0_index] << std::endl;
        std::cout << "jpsiTrkFlags.push_back(pass_jpsiTrk_path)"<< jpsiTrkFlags[trgMuonMatching_index] << std::endl;
           std::cout  << "----- reco = " << muon.pt() << " " << muon.eta() << " " << muon.phi() << " " 
		      << " HLT Dimuon0 = " << triggeringMuons[trgMuonMatchingDimuon0_index].pt() << " " << triggeringMuons[trgMuonMatchingDimuon0_index].eta() << " " << triggeringMuons[trgMuonMatchingDimuon0_index].phi()
		      << " HLT JpsiTrk = " << triggeringMuons[trgMuonMatchingJpsiTrk_index].pt() << " " << triggeringMuons[trgMuonMatchingJpsiTrk_index].eta() << " " << triggeringMuons[trgMuonMatchingJpsiTrk_index].phi()
		      << std::endl;
      }
      */

    }


  }


 
  // now produce output for analysis (code simplified loop of trg inside)
  // trigger muon + all compatible in dz with any tag
  for(unsigned int muIdx=0; muIdx<muons->size(); ++muIdx) 
  {
    const pat::Muon& mu = (*muons)[muIdx];
    //selection cuts
    if (mu.pt() < ptMin_) continue;
    if (fabs(mu.eta()) > absEtaMax_) continue;
    //following ID is needed for trigger muons not here
    // anyway it is off in the configuration
    //G: if (softMuonsOnly_ && !mu.isSoftMuon(PV)) continue;

    /* // same PV as the tag muon, both tag and probe only dz selection
    bool SkipMuon=true;
    for (const pat::Muon & trgmu : *trgmuons_out) {
	    if( fabs(mu.vz()-trgmu.vz()) > dzTrg_cleaning_ && dzTrg_cleaning_ >0 )
	      continue;
	    SkipMuon=false;
    } 
    // needs decission: what about events without trg muon? now we SKIP them
    if (SkipMuon)  continue;
    */

    // build transient track
    const reco::TransientTrack muonTT((*(mu.bestTrack())), &(*bFieldHandle)); //sara: check, why not using inner track for muons? 
    if (!muonTT.isValid()) continue;
    

    muons_out->emplace_back(mu);
    muons_out->back().addUserInt("isTriggering", muonIsTrigger[muIdx]);
    muons_out->back().addUserInt("isMuonFromJpsi_dimuon0Trg", muonIsFromJpsi_dimuon0Path[muIdx]);
    muons_out->back().addUserInt("isMuonFromJpsi_jpsiTrkTrg", muonIsFromJpsi_jpsiTrkPath[muIdx]);
    muons_out->back().addUserInt("isMuonFromJpsi_jpsiTrk_PsiPrimeTrg", muonIsFromJpsi_jpsiTrk_PsiPrimePath[muIdx]);
    muons_out->back().addUserInt("isMuonFromJpsi_jpsiTrk_NonResonantTrg", muonIsFromJpsi_jpsiTrk_NonResonantPath[muIdx]);
    muons_out->back().addUserInt("isDimuon0Trg", muonIsDimuon0Trg[muIdx]);
    muons_out->back().addUserInt("isJpsiTrkTrg", muonIsJpsiTrkTrg[muIdx]);
    muons_out->back().addUserInt("isJpsiTrk_PsiPrimeTrg", muonIsJpsiTrk_PsiPrimeTrg[muIdx]);
    muons_out->back().addUserInt("isJpsiTrk_NonResonantTrg", muonIsJpsiTrk_NonResonantTrg[muIdx]);

    trans_muons_out->emplace_back(muonTT);
  }

  iEvent.put(std::move(trgmuons_out),    "trgMuons");
  iEvent.put(std::move(muons_out),       "SelectedMuons");
  iEvent.put(std::move(trans_muons_out), "SelectedTransientMuons");
}

DEFINE_FWK_MODULE(MuonTriggerSelector);


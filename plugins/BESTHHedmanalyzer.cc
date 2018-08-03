// -*- C++ -*-
//===============================================================================================================
// Package:    BESTHHedmanalyzer/BESTHHedmanalyzer --------------------------------------------------------------
// Class:      BESTHHedmanalyzer                   --------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------
/**\class BESTHHedmanalyzer BESTHHedmanalyzer.cc BESTHHedmanalyzer/BESTHHedmanalyzer/plugins/BESTHHedmanalyzer.cc
-----------------------------------------------------------------------------------------------------------------
 Description: This class uses the BEST standalone function to analyze HH events                               ---
 ----------------------------------------------------------------------------------------------------------------
 Implementation:                                                                                              ---
     This EDAnalyzer is meant to be used with CMSSW_9_4_8 and utilized BEST and lwtnn                         ---
*/
//===============================================================================================================
// Original Author:  Brendan Regnery               --------------------------------------------------------------
//         Created:  Fri, 29 Jun 2018 22:32:28 GMT --------------------------------------------------------------
//===============================================================================================================
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

// ROOT include files
#include "TTree.h"
#include "TFile.h" 

// include the BEST files
#include "BESTAnalysis/BoostedEventShapeTagger/interface/BoostedEventShapeTagger.h"

/////////////////////////////////////////////////////////////////////
// Class declaration ------------------------------------------------
/////////////////////////////////////////////////////////////////////

class BESTHHedmanalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit BESTHHedmanalyzer(const edm::ParameterSet&);
      ~BESTHHedmanalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      
      //==================================================
      // User functions ----------------------------------
      //==================================================

      std::vector<pat::Jet> * sortJets(std::vector<pat::Jet> *jets);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //==================================================
      // ----------member data ---------------------------
      //==================================================

      // BEST Variables
      BoostedEventShapeTagger *m_BEST; 

      // Tokens
      edm::EDGetTokenT<std::vector<pat::Jet> > ak8JetsToken_;
      edm::EDGetTokenT<std::vector<reco::GenParticle> > genPartToken_;

      // tree variables
      TTree *bestTree;
      std::map<std::string, float> treeVars;
      std::vector<std::string> listOfVars;
};

//////////////////////////////////////////////////////////////////////
// constants, enums and typedefs -------------------------------------
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// static data member definitions ------------------------------------
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Constructors ------------------------------------------------------
//////////////////////////////////////////////////////////////////////

BESTHHedmanalyzer::BESTHHedmanalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   usesResource("TFileService");

   // pass the name of the BEST configuration file
   m_BEST = new BoostedEventShapeTagger("/afs/cern.ch/work/b/bregnery/public/HHwwwwMCgenerator/BEST_HH/CMSSW_9_4_8/src/BESTAnalysis/BoostedEventShapeTagger/data/config.txt");

   // Use TFile service to create a tree to store histogram variables
   edm::Service<TFileService> fs;
   bestTree = fs->make<TTree>("bestTree","bestTree");

   //----------------------------------------------------------------
   // Create tree variables and branches ----------------------------
   //----------------------------------------------------------------

   // Variabls for BEST results
   listOfVars.push_back("jet1_particleType");
   listOfVars.push_back("jet1_dnn_w");
   listOfVars.push_back("jet1_dnn_z");
   listOfVars.push_back("jet1_dnn_h");
   listOfVars.push_back("jet1_dnn_top");
   listOfVars.push_back("jet1_dnn_qcd");

   listOfVars.push_back("jet2_particleType");
   listOfVars.push_back("jet2_dnn_w");
   listOfVars.push_back("jet2_dnn_z");
   listOfVars.push_back("jet2_dnn_h");
   listOfVars.push_back("jet2_dnn_top");
   listOfVars.push_back("jet2_dnn_qcd");

   listOfVars.push_back("jet3_particleType");
   listOfVars.push_back("jet3_dnn_w");
   listOfVars.push_back("jet3_dnn_z");
   listOfVars.push_back("jet3_dnn_h");
   listOfVars.push_back("jet3_dnn_top");
   listOfVars.push_back("jet3_dnn_qcd");

   // AK8 jet variables
   listOfVars.push_back("nJets");

   listOfVars.push_back("jet1AK8_phi");
   listOfVars.push_back("jet1AK8_eta");

   listOfVars.push_back("jet2AK8_phi");
   listOfVars.push_back("jet2AK8_eta");

   listOfVars.push_back("jet3AK8_phi");
   listOfVars.push_back("jet3AK8_eta");

   // Generator Particle variables
   listOfVars.push_back("genHiggs1_phi");
   listOfVars.push_back("genHiggs1_eta");

   listOfVars.push_back("genHiggs2_phi");
   listOfVars.push_back("genHiggs2_eta");

   listOfVars.push_back("genW1_phi");
   listOfVars.push_back("genW1_eta");

   listOfVars.push_back("genW2_phi");
   listOfVars.push_back("genW2_eta");

   listOfVars.push_back("genW3_phi");
   listOfVars.push_back("genW3_eta");

   listOfVars.push_back("genW4_phi");
   listOfVars.push_back("genW4_eta");

   listOfVars.push_back("genVirW1_phi");
   listOfVars.push_back("genVirW1_eta");

   listOfVars.push_back("genVirW2_phi");
   listOfVars.push_back("genVirW2_eta");

   // Make Branches for each variable
   for (unsigned i = 0; i < listOfVars.size(); i++){
      treeVars[ listOfVars[i] ] = -999.99;
      bestTree->Branch( (listOfVars[i]).c_str() , &(treeVars[ listOfVars[i] ]), (listOfVars[i]+"/F").c_str() );
   }

   //----------------------------------------------------------------
   // Define input tags ---------------------------------------------
   //----------------------------------------------------------------

   // AK8 jets
   edm::InputTag ak8JetsTag_;
   ak8JetsTag_ = edm::InputTag("slimmedJetsAK8", "", "PAT");
   ak8JetsToken_ = consumes<std::vector<pat::Jet> >(ak8JetsTag_);
  
   // Gen particles
   edm::InputTag genPartTag_;
   genPartTag_ = edm::InputTag("prunedGenParticles", "", "PAT"); 
   genPartToken_ = consumes<std::vector<reco::GenParticle> >(genPartTag_);
}

//////////////////////////////////////////////////////////////////////
// Destructor --------------------------------------------------------
//////////////////////////////////////////////////////////////////////

BESTHHedmanalyzer::~BESTHHedmanalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

   // delete the BEST variable
   delete m_BEST;
}


/////////////////////////////////////////////////////////////////////
// member functions -------------------------------------------------
/////////////////////////////////////////////////////////////////////

//========================================================
// ------------ method called for each event  ------------
//========================================================
void
BESTHHedmanalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace fastjet;
   using namespace std;

   //---------------------------------------------------------------
   // Creat miniAOD object collections ----------------------------- 
   //---------------------------------------------------------------

   // Find objects corresponding to the token and link to the handle
   Handle< std::vector<pat::Jet> > ak8JetsCollection;
   iEvent.getByToken(ak8JetsToken_, ak8JetsCollection);
   std::vector<pat::Jet> ak8JetsVec = *ak8JetsCollection.product();
 
   Handle< std::vector<reco::GenParticle> > genPartCollection;
   iEvent.getByToken(genPartToken_, genPartCollection);
   std::vector<reco::GenParticle> genPart = *genPartCollection.product();

   //----------------------------------------------------------------
   // AK8 Jet Loop --------------------------------------------------
   //----------------------------------------------------------------

   std::vector<pat::Jet> *ak8Jets = new std::vector<pat::Jet>;
   ak8Jets = sortJets(&ak8JetsVec);
   int nJets = 0;
   for (std::vector<pat::Jet>::const_iterator jetBegin = ak8Jets->begin(), jetEnd = ak8Jets->end(), ijet = jetBegin; ijet != jetEnd; ++ijet){

      // AK8 Jets must meet BEST criteria
      if(ijet->numberOfDaughters() >= 2 && ijet->pt() >= 500 && ijet->userFloat("ak8PFJetsCHSSoftDropMass") > 40 ){

         nJets++;

         // Now sent jets to BEST for classification
         std::map<std::string,double> NNresults = m_BEST->execute(*ijet);  // ijet is a pat::Jet
         int particleType = m_BEST->getParticleID();                      // automatically calculate the particle classification

         // Link BEST and Jet results to the tree
         if(nJets == 1){
            // BEST
            treeVars["jet1_particleType"] = particleType;
            treeVars["jet1_dnn_w"] = NNresults["dnn_w"];
            treeVars["jet1_dnn_z"] = NNresults["dnn_z"];
            treeVars["jet1_dnn_h"] = NNresults["dnn_higgs"];
            treeVars["jet1_dnn_top"] = NNresults["dnn_top"];
            treeVars["jet1_dnn_qcd"] = NNresults["dnn_qcd"];

            // AK8 Jets
            treeVars["jet1AK8_phi"] = ijet->phi();
            treeVars["jet1AK8_eta"] = ijet->eta();
         }
         if(nJets == 2){
            // BEST
            treeVars["jet2_particleType"] = particleType;
            treeVars["jet2_dnn_w"] = NNresults["dnn_w"];
            treeVars["jet2_dnn_z"] = NNresults["dnn_z"];
            treeVars["jet2_dnn_h"] = NNresults["dnn_higgs"];
            treeVars["jet2_dnn_top"] = NNresults["dnn_top"];
            treeVars["jet2_dnn_qcd"] = NNresults["dnn_qcd"];

            // AK8 Jets
            treeVars["jet2AK8_phi"] = ijet->phi();
            treeVars["jet2AK8_eta"] = ijet->eta();
         }
         if(nJets == 3){
            // BEST
            treeVars["jet3_particleType"] = particleType;
            treeVars["jet3_dnn_w"] = NNresults["dnn_w"];
            treeVars["jet3_dnn_z"] = NNresults["dnn_z"];
            treeVars["jet3_dnn_h"] = NNresults["dnn_higgs"];
            treeVars["jet3_dnn_top"] = NNresults["dnn_top"];
            treeVars["jet3_dnn_qcd"] = NNresults["dnn_qcd"];

            // AK8 Jets
            treeVars["jet3AK8_phi"] = ijet->phi();
            treeVars["jet3AK8_eta"] = ijet->eta();
         }
      }
   }
   treeVars["nJets"] = nJets;

   //----------------------------------------------------------------
   // Gen Particles Loop --------------------------------------------
   //----------------------------------------------------------------

   int nHiggs = 0;
   int nWbosons = 0;
   for (vector<reco::GenParticle>::const_iterator genBegin = genPart.begin(), genEnd = genPart.end(), ipart = genBegin; ipart != genEnd; ++ipart){

      // Higgs
      if(abs(ipart->pdgId() ) == 25) nHiggs++;
      if(nHiggs == 1 && abs(ipart->pdgId() ) == 25 && ipart->status() < 50 && ipart->status() > 20 ){
         treeVars["genHiggs1_phi"] = ipart->phi();
         treeVars["genHiggs1_eta"] = ipart->eta();
      }
      if(nHiggs == 2 && abs(ipart->pdgId() ) == 25 && ipart->status() < 50 && ipart->status() > 20 ){
         treeVars["genHiggs2_phi"] = ipart->phi();
         treeVars["genHiggs2_eta"] = ipart->eta();
      }

      // W bosons
      if(abs(ipart->pdgId() ) == 24) nWbosons++;
      if(nWbosons == 1 && abs(ipart->pdgId() ) == 24 && ipart->status() < 50 && ipart->status() > 20 ){
         treeVars["genW1_phi"] = ipart->phi();
         treeVars["genW1_eta"] = ipart->eta();
      }
      if(nWbosons == 2 && abs(ipart->pdgId() ) == 24 && ipart->status() < 50 && ipart->status() > 20 ){
         treeVars["genW2_phi"] = ipart->phi();
         treeVars["genW2_eta"] = ipart->eta();
      }
      if(nWbosons == 3 && abs(ipart->pdgId() ) == 24 && ipart->status() < 50 && ipart->status() > 20 ){
         treeVars["genW3_phi"] = ipart->phi();
         treeVars["genW3_eta"] = ipart->eta();
      }
      if(nWbosons == 4 && abs(ipart->pdgId() ) == 24 && ipart->status() < 50 && ipart->status() > 20 ){
         treeVars["genW4_phi"] = ipart->phi();
         treeVars["genW4_eta"] = ipart->eta();
      }

      // Virtual W bosons
      if(nHiggs == 1 && abs(ipart->pdgId() ) == 25 && ipart->status() < 50 && ipart->status() > 20 && ipart->numberOfDaughters() == 3){
         if(abs(ipart->daughter(0)->pdgId() ) == 24){
            TLorentzVector daughter1 = TLorentzVector(ipart->daughter(1)->px(), ipart->daughter(1)->py(), ipart->daughter(1)->pz(), ipart->daughter(1)->energy());
            TLorentzVector daughter2 = TLorentzVector(ipart->daughter(2)->px(), ipart->daughter(2)->py(), ipart->daughter(2)->pz(), ipart->daughter(2)->energy());
            TLorentzVector virW1 = daughter1 + daughter2;
            treeVars["genVirW1_phi"] = virW1.Phi();
            treeVars["genVirW1_eta"] = virW1.Eta();
         }
         if(abs(ipart->daughter(2)->pdgId() ) == 24){
            TLorentzVector daughter1 = TLorentzVector(ipart->daughter(0)->px(), ipart->daughter(0)->py(), ipart->daughter(0)->pz(), ipart->daughter(0)->energy());
            TLorentzVector daughter2 = TLorentzVector(ipart->daughter(1)->px(), ipart->daughter(1)->py(), ipart->daughter(1)->pz(), ipart->daughter(1)->energy());
            TLorentzVector virW1 = daughter1 + daughter2;
            treeVars["genVirW1_phi"] = virW1.Phi();
            treeVars["genVirW1_eta"] = virW1.Eta();
         }
      }
      if(nHiggs == 2 && abs(ipart->pdgId() ) == 25 && ipart->status() < 50 && ipart->status() > 20 && ipart->numberOfDaughters() == 3){
         if(abs(ipart->daughter(0)->pdgId() ) == 24){
            TLorentzVector daughter1 = TLorentzVector(ipart->daughter(1)->px(), ipart->daughter(1)->py(), ipart->daughter(1)->pz(), ipart->daughter(1)->energy());
            TLorentzVector daughter2 = TLorentzVector(ipart->daughter(2)->px(), ipart->daughter(2)->py(), ipart->daughter(2)->pz(), ipart->daughter(2)->energy());
            TLorentzVector virW2 = daughter1 + daughter2;
            treeVars["genVirW2_phi"] = virW2.Phi();
            treeVars["genVirW2_eta"] = virW2.Eta();
         }
         if(abs(ipart->daughter(2)->pdgId() ) == 24){
            TLorentzVector daughter1 = TLorentzVector(ipart->daughter(0)->px(), ipart->daughter(0)->py(), ipart->daughter(0)->pz(), ipart->daughter(0)->energy());
            TLorentzVector daughter2 = TLorentzVector(ipart->daughter(1)->px(), ipart->daughter(1)->py(), ipart->daughter(1)->pz(), ipart->daughter(1)->energy());
            TLorentzVector virW2 = daughter1 + daughter2;
            treeVars["genVirW2_phi"] = virW2.Phi();
            treeVars["genVirW2_eta"] = virW2.Eta();
         }
      }
   }

   // Fill the tree that stores the NNresults
   bestTree->Fill();

}

//===================================================================
// Sort Jets --------------------------------------------------------
//-------------------------------------------------------------------
// This function takes in a Jet handle and returns a handle with  ---
// the top four jets sorted by pT                                 ---
//-------------------------------------------------------------------
 
std::vector<pat::Jet> *
BESTHHedmanalyzer::sortJets(std::vector<pat::Jet> *jets)
{
   std::vector<pat::Jet> *sortedJets = new std::vector<pat::Jet>;
   pat::Jet *jet1 = new pat::Jet;
   pat::Jet *jet2 = new pat::Jet;
   pat::Jet *jet3 = new pat::Jet;
   pat::Jet *jet4 = new pat::Jet;
   
   // Get leading jet
   if(jets->size() > 0){
      for (std::vector<pat::Jet>::const_iterator jetBegin = jets->begin(), jetEnd = jets->end(), ijet = jetBegin; ijet != jetEnd; ++ijet){
          if(!jet1){
             *jet1 = *ijet;
          }
          if(jet1 && ijet->pt() > jet1->pt() ){
             *jet1 = *ijet;
          }
       }
       sortedJets->push_back(*jet1);
       //std::cout << "Jet 1 pT: " << jet1->pt() << std::endl;
   }
   
   // Get subleading jet
   if(jets->size() > 1){
      for (std::vector<pat::Jet>::const_iterator jetBegin = jets->begin(), jetEnd = jets->end(), ijet = jetBegin; ijet != jetEnd; ++ijet){
          if(!jet2 && jet1->pt() > ijet->pt() ){
             *jet2 = *ijet;
          }
          if(jet2 && jet1->pt() > ijet->pt() && ijet->pt() > jet2->pt() ){
             *jet2 = *ijet;
          }
       }
       sortedJets->push_back(*jet2);
       //std::cout << "Jet 2 pT: " << jet2->pt() << std::endl;
   } 

   // Get third leading jet
   if(jets->size() > 2){
      for (std::vector<pat::Jet>::const_iterator jetBegin = jets->begin(), jetEnd = jets->end(), ijet = jetBegin; ijet != jetEnd; ++ijet){
          if(!jet3 && jet2->pt() > ijet->pt() ){
             *jet3 = *ijet;
          }
          if(jet3 && jet2->pt() > ijet->pt() && ijet->pt() > jet3->pt() ){
             *jet3 = *ijet;
          }
       }
       sortedJets->push_back(*jet3);
       //std::cout << "Jet 3 pT: " << jet3->pt() << std::endl;
   } 

   // Get fourth leading jet
   if(jets->size() > 3){
      for (std::vector<pat::Jet>::const_iterator jetBegin = jets->begin(), jetEnd = jets->end(), ijet = jetBegin; ijet != jetEnd; ++ijet){
          if(!jet4 && jet3->pt() > ijet->pt() ){
             *jet4 = *ijet;
          }
          if(jet4 && jet3->pt() > ijet->pt() && ijet->pt() > jet4->pt() ){
             *jet4 = *ijet;
          }
       }
       sortedJets->push_back(*jet4);
       //std::cout << "Jet 4 pT: " << jet4->pt() << std::endl;
   }

   return sortedJets;
} 

//=======================================================================================
// ------------ method called once each job just before starting event loop  ------------
//=======================================================================================

void 
BESTHHedmanalyzer::beginJob()
{
}

//========================================================================================
// ------------ method called once each job just after ending the event loop  ------------
//========================================================================================

void 
BESTHHedmanalyzer::endJob() 
{
}

//==================================================================================================
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
//==================================================================================================

void
BESTHHedmanalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BESTHHedmanalyzer);

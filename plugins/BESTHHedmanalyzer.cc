// -*- C++ -*-
//
// Package:    BESTHHedmanalyzer/BESTHHedmanalyzer
// Class:      BESTHHedmanalyzer
// 
/**\class BESTHHedmanalyzer BESTHHedmanalyzer.cc BESTHHedmanalyzer/BESTHHedmanalyzer/plugins/BESTHHedmanalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Brendan Regnery
//         Created:  Fri, 29 Jun 2018 22:32:28 GMT
//
//


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

// ROOT include files
#include "TTree.h"
#include "TFile.h" 

// include the BEST files
#include "BESTAnalysis/BoostedEventShapeTagger/interface/BoostedEventShapeTagger.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class BESTHHedmanalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit BESTHHedmanalyzer(const edm::ParameterSet&);
      ~BESTHHedmanalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      BoostedEventShapeTagger *m_BEST; 
      TTree *bestTree;
      edm::EDGetTokenT<std::vector<pat::Jet> > ak8JetsToken_;
      std::map<std::string, float> treeVars;
      std::vector<std::string> listOfVars;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
BESTHHedmanalyzer::BESTHHedmanalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   usesResource("TFileService");

   // pass the name of the BEST configuration file
   //m_BEST = new BoostedEventShapeTagger("BESTAnalysis/BoostedEventShapeTagger/data/config.txt");
   m_BEST = new BoostedEventShapeTagger("/afs/cern.ch/work/b/bregnery/public/HHwwwwMCgenerator/BEST_HH/CMSSW_9_4_8/src/BESTAnalysis/BoostedEventShapeTagger/data/config.txt");

   // Use TFile service to create a tree to store histogram variables
   edm::Service<TFileService> fs;
   bestTree = fs->make<TTree>("bestTree","bestTree");

   // Create tree variables and branches
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
   listOfVars.push_back("nJets");
   for (unsigned i = 0; i < listOfVars.size(); i++){
      treeVars[ listOfVars[i] ] = -999.99;
      bestTree->Branch( (listOfVars[i]).c_str() , &(treeVars[ listOfVars[i] ]), (listOfVars[i]+"/F").c_str() );
   }
 
   // Define input tags
   edm::InputTag ak8JetsTag_;
   ak8JetsTag_ = edm::InputTag("slimmedJetsAK8", "", "PAT");
   ak8JetsToken_ = consumes<std::vector<pat::Jet> >(ak8JetsTag_);

}


BESTHHedmanalyzer::~BESTHHedmanalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

   // delete the BEST variable
   delete m_BEST;
}


//
// member functions
//

// ------------ method called for each event  ------------
void
BESTHHedmanalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace fastjet;
   using namespace std;

   // Find objects corresponding to the token and link to the handle
   Handle< std::vector<pat::Jet> > ak8Jets;
   iEvent.getByToken(ak8JetsToken_, ak8Jets);

   // loop over the jets
   int nJets = 0;
   for (std::vector<pat::Jet>::const_iterator jetBegin = ak8Jets->begin(), jetEnd = ak8Jets->end(), ijet = jetBegin; ijet != jetEnd; ++ijet){
      if(ijet->numberOfDaughters() >= 2 && ijet->pt() >= 500 && ijet->userFloat("ak8PFJetsCHSSoftDropMass") > 40 ){
         std::map<std::string,double> NNresults = m_BEST->execute(*ijet);  // ijet is a pat::Jet
         int particleType = m_BEST->getParticleID();                      // automatically calculate the particle classification
         if(nJets == 0){
            treeVars["jet1_particleType"] = particleType;
            treeVars["jet1_dnn_w"] = NNresults["dnn_w"];
            treeVars["jet1_dnn_z"] = NNresults["dnn_z"];
            treeVars["jet1_dnn_h"] = NNresults["dnn_higgs"];
            treeVars["jet1_dnn_top"] = NNresults["dnn_top"];
            treeVars["jet1_dnn_qcd"] = NNresults["dnn_qcd"];
         }
         if(nJets == 1){
            treeVars["jet2_particleType"] = particleType;
            treeVars["jet2_dnn_w"] = NNresults["dnn_w"];
            treeVars["jet2_dnn_z"] = NNresults["dnn_z"];
            treeVars["jet2_dnn_h"] = NNresults["dnn_higgs"];
            treeVars["jet2_dnn_top"] = NNresults["dnn_top"];
            treeVars["jet2_dnn_qcd"] = NNresults["dnn_qcd"];
         }
         if(nJets == 2){
            treeVars["jet3_particleType"] = particleType;
            treeVars["jet3_dnn_w"] = NNresults["dnn_w"];
            treeVars["jet3_dnn_z"] = NNresults["dnn_z"];
            treeVars["jet3_dnn_h"] = NNresults["dnn_higgs"];
            treeVars["jet3_dnn_top"] = NNresults["dnn_top"];
            treeVars["jet3_dnn_qcd"] = NNresults["dnn_qcd"];
         }
         nJets++;
      }
   }
   treeVars["nJets"] = nJets;

  // Fill the tree that stores the NNresults
  bestTree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
BESTHHedmanalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
BESTHHedmanalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
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

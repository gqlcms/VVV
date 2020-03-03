// -*- C++ -*-
//
// Package:    hlt/hlt
// Class:      hlt
// 
/**\class hlt hlt.cc hlt/hlt/plugins/hlt.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Qilong Guo
//         Created:  Thu, 27 Feb 2020 13:17:52 GMT
//
//


#include <iostream>
// system include files
#include <memory>

// user include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//jet
#include "DataFormats/PatCandidates/interface/Jet.h"


#include "TTree.h"

using namespace std;

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class hlt : public edm::EDAnalyzer  {
   public:
      explicit hlt(const edm::ParameterSet&);
      ~hlt();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      virtual bool tightJetIDpuppi( const pat::Jet& j);
      void setDummyValues();


      // ----------member data ---------------------------
    TTree* outTree_;

    edm::Handle<pat::JetCollection> puppijets_;
    edm::EDGetTokenT<pat::JetCollection> puppijetInputToken_;
    double jetAK8puppi_dnnDecorrW,jetAK8puppi_dnnDecorrW_2,jetAK8puppi_dnnDecorrW_3;
    double jetAK8puppi_dnnDecorrZ,jetAK8puppi_dnnDecorrZ_2,jetAK8puppi_dnnDecorrZ_3;
    double jetAK8puppi_ptJEC, jetAK8puppi_eta, jetAK8puppi_phi, jetAK8puppi_sd;
    double jetAK8puppi_ptJEC_2, jetAK8puppi_eta_2, jetAK8puppi_phi_2, jetAK8puppi_sd_2;
    double jetAK8puppi_ptJEC_3, jetAK8puppi_eta_3, jetAK8puppi_phi_3, jetAK8puppi_sd_3;



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
hlt::hlt(const edm::ParameterSet& iConfig) 
{
   //now do what ever initialization is needed
   //usesResource("TFileService");

    edm::Service<TFileService> fs;

    outTree_ = fs->make<TTree>("EDBRCandidates","EDBR Candidates");

    outTree_->Branch("jetAK8puppi_dnnDecorrW"           ,&jetAK8puppi_dnnDecorrW         ,"jetAK8puppi_dnnDecorrW/D"           );
    outTree_->Branch("jetAK8puppi_dnnDecorrW_2"           ,&jetAK8puppi_dnnDecorrW_2         ,"jetAK8puppi_dnnDecorrW_2/D"           );
    outTree_->Branch("jetAK8puppi_dnnDecorrW_3"           ,&jetAK8puppi_dnnDecorrW_3         ,"jetAK8puppi_dnnDecorrW_3/D"           );
    outTree_->Branch("jetAK8puppi_dnnDecorrZ"           ,&jetAK8puppi_dnnDecorrZ         ,"jetAK8puppi_dnnDecorrZ/D"           );
    outTree_->Branch("jetAK8puppi_dnnDecorrZ_2"           ,&jetAK8puppi_dnnDecorrZ_2         ,"jetAK8puppi_dnnDecorrZ_2/D"           );
    outTree_->Branch("jetAK8puppi_dnnDecorrZ_3"           ,&jetAK8puppi_dnnDecorrZ_3         ,"jetAK8puppi_dnnDecorrZ_3/D"           );
    outTree_->Branch("jetAK8puppi_ptJEC"          ,&jetAK8puppi_ptJEC         ,"jetAK8puppi_ptJEC/D"         );
    outTree_->Branch("jetAK8puppi_ptJEC_2"          ,&jetAK8puppi_ptJEC_2         ,"jetAK8puppi_ptJEC_2/D"         );
    outTree_->Branch("jetAK8puppi_ptJEC_3"          ,&jetAK8puppi_ptJEC_3         ,"jetAK8puppi_ptJEC_3/D"         );
    outTree_->Branch("jetAK8puppi_eta"          ,&jetAK8puppi_eta         ,"jetAK8puppi_eta/D"         );
    outTree_->Branch("jetAK8puppi_eta_2"          ,&jetAK8puppi_eta_2         ,"jetAK8puppi_eta_2/D"         );
    outTree_->Branch("jetAK8puppi_eta_3"          ,&jetAK8puppi_eta_3         ,"jetAK8puppi_eta_3/D"         );
    outTree_->Branch("jetAK8puppi_phi"          ,&jetAK8puppi_phi         ,"jetAK8puppi_phi/D"         );
    outTree_->Branch("jetAK8puppi_phi_2"          ,&jetAK8puppi_phi_2         ,"jetAK8puppi_phi_2/D"         );
    outTree_->Branch("jetAK8puppi_phi_3"          ,&jetAK8puppi_phi_3         ,"jetAK8puppi_phi_3/D"         );
    outTree_->Branch("jetAK8puppi_sd"          ,&jetAK8puppi_sd         ,"jetAK8puppi_sd/D"         );
    outTree_->Branch("jetAK8puppi_sd_2"          ,&jetAK8puppi_sd_2         ,"jetAK8puppi_sd_2/D"         );
    outTree_->Branch("jetAK8puppi_sd_3"          ,&jetAK8puppi_sd_3         ,"jetAK8puppi_sd_3/D"         );

    puppijetInputToken_ = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("puppijets"));

}


hlt::~hlt()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
hlt::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

    setDummyValues(); //Initalize variables with dummy values

//        bool doPuppi  = iEvent.getByToken(puppijetInputToken_, puppijets_ );
        iEvent.getByToken(puppijetInputToken_, puppijets_ );
        bool doPuppi = true;
        if( doPuppi ){//1

//leading pt jet:
            int usenumber3 = -1; double pt_larger=0;
            int numvhad = puppijets_->size();
            for( int inum = 0; inum< numvhad; inum++){
                const pat::Jet& Vpuppi = puppijets_->at(inum);
                if(tightJetIDpuppi(Vpuppi)<1) continue;
                if(Vpuppi.pt() > pt_larger && fabs(Vpuppi.eta())<2.4 && inum<4) {pt_larger = Vpuppi.pt(); usenumber3 = inum; continue;}
            }
            if (usenumber3>-1) {//2
                const pat::Jet& hadronicVpuppi = puppijets_->at(usenumber3);
                jetAK8puppi_dnnDecorrW         = hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:WvsQCD");
                jetAK8puppi_dnnDecorrZ         = hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZvsQCD");
                jetAK8puppi_ptJEC       = hadronicVpuppi.pt();
                jetAK8puppi_eta     = hadronicVpuppi.eta();
                jetAK8puppi_phi     = hadronicVpuppi.phi();
                jetAK8puppi_sd       =  hadronicVpuppi.userFloat("ak8PFJetsPuppiSoftDropMass"); 
}
//2
            int usenumber2 = -1; double pt_larger2=0;
            for( int inum = 0; inum< numvhad; inum++){
                const pat::Jet& Vpuppi = puppijets_->at(inum);
                if(tightJetIDpuppi(Vpuppi)<1) continue;
                if(Vpuppi.pt() > pt_larger2 && fabs(Vpuppi.eta())<2.4 && inum != usenumber3 && inum<4) {pt_larger2 = Vpuppi.pt(); usenumber2 = inum; continue;}
            }

            if (usenumber2>-1) {//2
                const pat::Jet& hadronicVpuppi = puppijets_->at(usenumber2);
                jetAK8puppi_dnnDecorrW_2         = hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:WvsQCD");
                jetAK8puppi_dnnDecorrZ_2         = hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZvsQCD");
                jetAK8puppi_ptJEC_2       = hadronicVpuppi.pt();
                jetAK8puppi_eta_2     = hadronicVpuppi.eta();
                jetAK8puppi_phi_2     = hadronicVpuppi.phi();
                jetAK8puppi_sd_2       =  hadronicVpuppi.userFloat("ak8PFJetsPuppiSoftDropMass"); 
}
//3
            int usenumber1 = -1; double pt_larger1=0;
            for( int inum = 0; inum< numvhad; inum++){
                const pat::Jet& Vpuppi = puppijets_->at(inum);
                if(tightJetIDpuppi(Vpuppi)<1) continue;
                if(Vpuppi.pt() > pt_larger1 && fabs(Vpuppi.eta())<2.4 && inum != usenumber3 && inum != usenumber2 && inum<4) {pt_larger1 = Vpuppi.pt(); usenumber1 = inum; continue;}
            }

            if (usenumber1>-1) {//2
                const pat::Jet& hadronicVpuppi = puppijets_->at(usenumber1);
                jetAK8puppi_dnnDecorrW_3         = hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:WvsQCD");
                jetAK8puppi_dnnDecorrZ_3         = hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZvsQCD");
                jetAK8puppi_ptJEC_3       = hadronicVpuppi.pt();
                jetAK8puppi_eta_3     = hadronicVpuppi.eta();
                jetAK8puppi_phi_3     = hadronicVpuppi.phi();
                jetAK8puppi_sd_3       =  hadronicVpuppi.userFloat("ak8PFJetsPuppiSoftDropMass"); 
}




outTree_->Fill();
        }

}


// ------------ method called once each job just before starting event loop  ------------
void 
hlt::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
hlt::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
hlt::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void hlt::setDummyValues() {

    jetAK8puppi_dnnDecorrW          = -99;
    jetAK8puppi_dnnDecorrW_2          = -99;
    jetAK8puppi_dnnDecorrW_3          = -99;
    jetAK8puppi_dnnDecorrZ          = -99;
    jetAK8puppi_dnnDecorrZ_2          = -99;
    jetAK8puppi_dnnDecorrZ_3          = -99;
    jetAK8puppi_ptJEC         = -99;
    jetAK8puppi_eta         = -99;
    jetAK8puppi_phi         = -99;
    jetAK8puppi_ptJEC_2         = -99;
    jetAK8puppi_eta_2         = -99;
    jetAK8puppi_phi_2         = -99;
    jetAK8puppi_ptJEC_3         = -99;
    jetAK8puppi_eta_3         = -99;
    jetAK8puppi_phi_3         = -99;
    jetAK8puppi_sd         = -99;
    jetAK8puppi_sd_2         = -99;
    jetAK8puppi_sd_3         = -99;


}

bool
hlt::tightJetIDpuppi( const pat::Jet& j ) {
    // refer to https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_data
    if(j.pt()>0.){
    double NHF = j.neutralHadronEnergyFraction();
    double NEMF = j.neutralEmEnergyFraction();
    double CHF = j.chargedHadronEnergyFraction();
    //double MUF = j.muonEnergyFraction();
    int NumConst = j.chargedMultiplicity()+j.neutralMultiplicity();
    int NumNeutralParticle =j.neutralMultiplicity();
    int CHM = j.chargedMultiplicity();
    double eta = j.eta();
    return ((  (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 ) || (abs(eta)>2.4 && abs(eta)<=2.7) )) || (NHF<0.99 && abs(eta)>2.7 && abs(eta)<=3.0 ) || (NEMF<0.90 && NHF>0.02 &&NumNeutralParticle>2 && NumNeutralParticle<15 && abs(eta)>3.0) ) ;
}
else{
return (0);
    }
}




//define this as a plug-in
DEFINE_FWK_MODULE(hlt);

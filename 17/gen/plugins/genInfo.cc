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

#include "DataFormats/Candidate/interface/Candidate.h"

//Gen 
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


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
      void setDummyValues();
    virtual const reco::Candidate* findLastW(const reco::Candidate *particle,int IDpdg);
    virtual const reco::Candidate* findLasttau(const reco::Candidate *particle,int IDpdg);
    virtual const reco::Candidate* findFirstW(const reco::Candidate *particle,int IDpdg);




      // ----------member data ---------------------------
    TTree* outTree_;





      //-----------------------------  GenInfo-----------------------------
    bool RunOnMC_;
    edm::EDGetTokenT<edm::View<reco::GenParticle>> genSrc_;

//----------------------------------  top info ----------------------------

//top,
    double gentop_pt, gentop_eta, gentop_phi, gentop_mass;
//anti top 
    double genantitop_pt, genantitop_eta, genantitop_phi, genantitop_mass;
//topw,
    double gent_w_pt,gent_w_phi,gent_w_eta,gent_w_mass;
    double gent_w_tag;
    //top w daughter
    double gent_w_q1_pt,gent_w_q1_phi,gent_w_q1_eta,gent_w_q1_e,gent_w_q1_pdg;
    double gent_w_q2_pt,gent_w_q2_phi,gent_w_q2_eta,gent_w_q2_e,gent_w_q2_pdg;
//top b
    double gent_b_pt,gent_b_phi,gent_b_eta,gent_b_mass;
//anti top w
    double genantit_w_pt,genantit_w_phi,genantit_w_eta,genantit_w_mass;
    double genantit_w_tag;
    // anti top w daughter
    double genantit_w_q1_pt,genantit_w_q1_phi,genantit_w_q1_eta,genantit_w_q1_e,genantit_w_q1_pdg;
    double genantit_w_q2_pt,genantit_w_q2_phi,genantit_w_q2_eta,genantit_w_q2_e,genantit_w_q2_pdg;
//anti top b
    double genantit_b_pt,genantit_b_phi,genantit_b_eta,genantit_b_mass;

// ================================ top ===================================

// ---------------------------------- w ----------------------------------

// last w
    double ptgenwl[5],etagenwl[5],phigenwl[5],massgenwl[5],taggenwl[5],taggenwmother[5];
// first w 
    double ptgenwf[5],etagenwf[5],phigenwf[5],massgenwf[5];
// last w daughter
    double genw_q1_pt[5],genw_q1_eta[5],genw_q1_phi[5],genw_q1_e[5],genw_q1_pdg[5];
    double genw_q2_pt[5],genw_q2_eta[5],genw_q2_phi[5],genw_q2_e[5],genw_q2_pdg[5];
// ================================== w ==================================

// --------------------------------- z ----------------------------------
//first z
    double ptgenzf[5],etagenzf[5],phigenzf[5],massgenzf[5];
//last z
    double ptgenzl[5],etagenzl[5],phigenzl[5],massgenzl[5],taggenzl[5];

// ================================= z ==================================
     //#############################  GenInfo  #############################


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
//top
    outTree_->Branch("gentop_pt"        ,&gentop_pt       ,"gentop_pt/D"       );
    outTree_->Branch("gentop_eta"        ,&gentop_eta       ,"gentop_eta/D"       );
    outTree_->Branch("gentop_phi"        ,&gentop_phi       ,"gentop_phi/D"       );
    outTree_->Branch("gentop_mass"        ,&gentop_mass       ,"gentop_mass/D"       );
//anti top
    outTree_->Branch("genantitop_pt"        ,&genantitop_pt       ,"genantitop_pt/D"       );
    outTree_->Branch("genantitop_eta"        ,&genantitop_eta       ,"genantitop_eta/D"       );
    outTree_->Branch("genantitop_phi"        ,&genantitop_phi       ,"genantitop_phi/D"       );
    outTree_->Branch("genantitop_mass"        ,&genantitop_mass       ,"genantitop_mass/D"       );

//topw
        outTree_->Branch("gent_w_pt"           ,&gent_w_pt         ,"gent_w_pt/D"          );
        outTree_->Branch("gent_w_eta"           ,&gent_w_eta         ,"gent_w_eta/D"          );
        outTree_->Branch("gent_w_phi"           ,&gent_w_phi         ,"gent_w_phi/D"          );
        outTree_->Branch("gent_w_mass"           ,&gent_w_mass         ,"gent_w_mass/D"          );
        //decay method 
        outTree_->Branch("gent_w_tag"           ,&gent_w_tag         ,"gent_w_tag/D"          );
        //w daughter
        outTree_->Branch("gent_w_q1_pt"           ,&gent_w_q1_pt         ,"gent_w_q1_pt/D"          );
        outTree_->Branch("gent_w_q1_eta"           ,&gent_w_q1_eta         ,"gent_w_q1_eta/D"          );
        outTree_->Branch("gent_w_q1_phi"           ,&gent_w_q1_phi         ,"gent_w_q1_phi/D"          );
        outTree_->Branch("gent_w_q1_e"           ,&gent_w_q1_e         ,"gent_w_q1_e/D"          );
        outTree_->Branch("gent_w_q1_pdg"           ,&gent_w_q1_pdg         ,"gent_w_q1_pdg/D"          );
        outTree_->Branch("gent_w_q2_pt"           ,&gent_w_q2_pt         ,"gent_w_q2_pt/D"          );
        outTree_->Branch("gent_w_q2_eta"           ,&gent_w_q2_eta         ,"gent_w_q2_eta/D"          );
        outTree_->Branch("gent_w_q2_phi"           ,&gent_w_q2_phi         ,"gent_w_q2_phi/D"          );
        outTree_->Branch("gent_w_q2_e"           ,&gent_w_q2_e         ,"gent_w_q2_e/D"          );
        outTree_->Branch("gent_w_q2_pdg"           ,&gent_w_q2_pdg         ,"gent_w_q2_pdg/D"          );
//top b
        outTree_->Branch("gent_b_pt"           ,&gent_b_pt         ,"gent_b_pt/D"          );
        outTree_->Branch("gent_b_eta"           ,&gent_b_eta         ,"gent_b_eta/D"          );
        outTree_->Branch("gent_b_phi"           ,&gent_b_phi         ,"gent_b_phi/D"          );
        outTree_->Branch("gent_b_mass"           ,&gent_b_mass         ,"gent_b_mass/D"          );
//anti top w
        outTree_->Branch("genantit_w_pt"           ,&genantit_w_pt         ,"genantit_w_pt/D"          );
        outTree_->Branch("genantit_w_eta"           ,&genantit_w_eta         ,"genantit_w_eta/D"          );
        outTree_->Branch("genantit_w_phi"           ,&genantit_w_phi         ,"genantit_w_phi/D"          );
        outTree_->Branch("genantit_w_mass"           ,&genantit_w_mass         ,"genantit_w_mass/D"          );

        outTree_->Branch("genantit_w_tag"           ,&genantit_w_tag         ,"genantit_w_tag/D"          );
        //anti top daughter
        outTree_->Branch("genantit_w_q1_pt"           ,&genantit_w_q1_pt         ,"genantit_w_q1_pt/D"          );
        outTree_->Branch("genantit_w_q1_eta"           ,&genantit_w_q1_eta         ,"genantit_w_q1_eta/D"          );
        outTree_->Branch("genantit_w_q1_phi"           ,&genantit_w_q1_phi         ,"genantit_w_q1_phi/D"          );
        outTree_->Branch("genantit_w_q1_e"           ,&genantit_w_q1_e         ,"genantit_w_q1_e/D"          );
        outTree_->Branch("genantit_w_q1_pdg"           ,&genantit_w_q1_pdg         ,"genantit_w_q1_pdg/D"          );
        outTree_->Branch("genantit_w_q2_pt"           ,&genantit_w_q2_pt         ,"genantit_w_q2_pt/D"          );
        outTree_->Branch("genantit_w_q2_eta"           ,&genantit_w_q2_eta         ,"genantit_w_q2_eta/D"          );
        outTree_->Branch("genantit_w_q2_phi"           ,&genantit_w_q2_phi         ,"genantit_w_q2_phi/D"          );
        outTree_->Branch("genantit_w_q2_e"           ,&genantit_w_q2_e         ,"gent_w_q2_e/D"          );
        outTree_->Branch("genantit_w_q2_pdg"           ,&genantit_w_q2_pdg         ,"genantit_w_q2_pdg/D"          );

//anti top b
        outTree_->Branch("genantit_b_pt"           ,&genantit_b_pt         ,"genantit_b_pt/D"          );
        outTree_->Branch("genantit_b_eta"           ,&genantit_b_eta         ,"genantit_b_eta/D"          );
        outTree_->Branch("genantit_b_phi"           ,&genantit_b_phi         ,"genantit_b_phi/D"          );
        outTree_->Branch("genantit_b_mass"           ,&genantit_b_mass         ,"genantit_b_mass/D"          );


// --------------------------------------------------------- w ------------------------------------------------------
// last w
        outTree_->Branch("ptgenwl"           ,ptgenwl         ,"ptgenwl[5]/D"          );
        outTree_->Branch("etagenwl"           ,etagenwl         ,"etagenwl[5]/D"          );
        outTree_->Branch("phigenwl"           ,phigenwl       ,"phigenwl[5]/D"          );
        outTree_->Branch("massgenwl"           ,massgenwl         ,"massgenwl[5]/D"          );
        outTree_->Branch("taggenwl"           ,taggenwl         ,"taggenwl[5]/D"          );
        outTree_->Branch("taggenwmother"           ,taggenwmother         ,"taggenwmother[5]/D"          );
// first w
        outTree_->Branch("ptgenwf"           ,ptgenwf         ,"ptgenwf[5]/D"          );
        outTree_->Branch("etagenwf"           ,etagenwf         ,"etagenwf[5]/D"          );
        outTree_->Branch("phigenwf"           ,phigenwf       ,"phigenwf[5]/D"          );
        outTree_->Branch("massgenwf"           ,massgenwf         ,"massgenwf[5]/D"          );
// last w daughter
        outTree_->Branch("genw_q1_pt"           ,genw_q1_pt         ,"genw_q1_pt[5]/D"          );
        outTree_->Branch("genw_q1_phi"           ,genw_q1_phi         ,"genw_q1_phi[5]/D"          );
        outTree_->Branch("genw_q1_eta"           ,genw_q1_eta         ,"genw_q1_eta[5]/D"          );
        outTree_->Branch("genw_q1_e"           ,genw_q1_e         ,"genw_q1_e[5]/D"          );
        outTree_->Branch("genw_q1_pdg"           ,genw_q1_pdg         ,"genw_q1_pdg[5]/D"          );
        outTree_->Branch("genw_q2_pt"           ,genw_q2_pt         ,"genw_q2_pt[5]/D"          );
        outTree_->Branch("genw_q2_phi"           ,genw_q2_phi         ,"genw_q2_phi[5]/D"          );
        outTree_->Branch("genw_q2_eta"           ,genw_q2_eta         ,"genw_q2_eta[5]/D"          );
        outTree_->Branch("genw_q2_e"           ,genw_q2_e         ,"genw_q2_e[5]/D"          );
        outTree_->Branch("genw_q2_pdg"           ,genw_q2_pdg         ,"genw_q2_pdg[5]/D"          );


// ========================================================= w =======================================================

//---------------------------------------------------------- z -------------------------------------------------------

        outTree_->Branch("ptgenzf"           ,ptgenzf         ,"ptgenzf[5]/D"          );
        outTree_->Branch("etagenzf"           ,etagenzf         ,"etagenzf[5]/D"          );
        outTree_->Branch("phigenzf"           ,phigenzf       ,"phigenzf[5]/D"          );
        outTree_->Branch("massgenzf"           ,massgenzf         ,"massgenzf[5]/D"          );
        outTree_->Branch("taggenzl"           ,taggenzl         ,"taggenzl[5]/D"          );
        outTree_->Branch("ptgenzf"           ,ptgenzf         ,"ptgenzf[5]/D"          );
        outTree_->Branch("etagenzf"           ,etagenzf         ,"etagenzf[5]/D"          );
        outTree_->Branch("phigenzf"           ,phigenzf       ,"phigenzf[5]/D"          );
        outTree_->Branch("massgenzf"           ,massgenzf         ,"massgenzf[5]/D"          );


// ========================================================= z ========================================================

      //-----------------------------  GenInfo -----------------------------
    RunOnMC_        = iConfig.getParameter<bool>("RunOnMC");
    genSrc_      = consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>( "genSrc") ) ;

     //#############################  GenInfo  #############################

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


    edm::Handle<edm::View<reco::GenParticle> > genParticles;
    iEvent.getByToken(genSrc_, genParticles);

    if(RunOnMC_)
    {//MC Info
        for(size_t ik=0; ik<genParticles->size();ik++)
        {// loop on gen
            //for top
            const reco::Candidate* ptop0 = &(*genParticles)[ik];
            //findLasttau(if it's top.return last top, else,return itself)
            const reco::Candidate* ptop=findLasttau(ptop0,6);
                if(ptop0->pdgId()== 6 && gentop_pt==-99) {
                    gentop_pt = ptop->pt();
                    gentop_eta = ptop->eta();
                    gentop_phi = ptop->phi();
                    gentop_mass = ptop->mass();
                    for(int i=0;ptop->daughter(i)!=NULL;i++){
                    //top->w
                        if(abs(ptop->daughter(i)->pdgId())==24){
                            gent_w_pt=ptop->daughter(i)->pt();
                            gent_w_eta=ptop->daughter(i)->eta();
                            gent_w_phi=ptop->daughter(i)->phi();
                            gent_w_mass=ptop->daughter(i)->mass();
                            const reco::Candidate* ptw0 = ptop->daughter(i);
                            const reco::Candidate* ptw= findLastW(ptw0,24);
                            if(ptw->daughter(0)!=NULL)
                            {

                                if( abs(ptw->daughter(0)->pdgId())<=6 ){
                                    gent_w_tag=4;
                                    gent_w_q1_pt=ptw->daughter(0)->pt();
                                    gent_w_q1_eta=ptw->daughter(0)->eta();
                                    gent_w_q1_phi=ptw->daughter(0)->phi();
                                    gent_w_q1_e=ptw->daughter(0)->energy();
                                    gent_w_q1_pdg=ptw->daughter(0)->pdgId();
                                    gent_w_q2_pt=ptw->daughter(1)->pt();
                                    gent_w_q2_eta=ptw->daughter(1)->eta();
                                    gent_w_q2_phi=ptw->daughter(1)->phi();
                                    gent_w_q2_e=ptw->daughter(1)->energy();
                                    gent_w_q2_pdg=ptw->daughter(1)->pdgId();
                                }
                                //lepton w
                                if( abs(ptw->daughter(0)->pdgId())==11 ||abs(ptw->daughter(0)->pdgId())==12 ) gent_w_tag=1;
                                if( abs(ptw->daughter(0)->pdgId())==12 ||abs(ptw->daughter(0)->pdgId())==13 ) gent_w_tag=2;
                                if( abs(ptw->daughter(0)->pdgId())==14 ||abs(ptw->daughter(0)->pdgId())==15 ) gent_w_tag=3;
                                if( abs(ptw->daughter(0)->pdgId())==11 ||abs(ptw->daughter(0)->pdgId())==12||abs(ptw->daughter(0)->pdgId())==13 ||abs(ptw->daughter(0)->pdgId())==14||abs(ptw->daughter(0)->pdgId())==15 ||abs(ptw->daughter(0)->pdgId())==16)
                                {
                                    gent_w_q1_pt=ptw->daughter(0)->pt();
                                    gent_w_q1_eta=ptw->daughter(0)->eta();
                                    gent_w_q1_phi=ptw->daughter(0)->phi();
                                    gent_w_q1_e=ptw->daughter(0)->energy();
                                    gent_w_q1_pdg=ptw->daughter(0)->pdgId();
                                    gent_w_q2_pt=ptw->daughter(1)->pt();
                                    gent_w_q2_eta=ptw->daughter(1)->eta();
                                    gent_w_q2_phi=ptw->daughter(1)->phi();
                                    gent_w_q2_e=ptw->daughter(1)->energy();
                                    gent_w_q2_pdg=ptw->daughter(1)->pdgId();
                                }
                        }
                        }
                        // top -> b
                        if(abs(ptop->daughter(i)->pdgId())==5){
                            gent_b_pt=ptop->daughter(i)->pt();
                            gent_b_eta=ptop->daughter(i)->eta();
                            gent_b_phi=ptop->daughter(i)->phi();
                            gent_b_mass=ptop->daughter(i)->mass();
                        }
                }
                }
                // anti top
                if(ptop0->pdgId()== -6 && genantitop_pt==-99) {
                    genantitop_pt = ptop->pt();
                    genantitop_eta = ptop->eta();
                    genantitop_phi = ptop->phi();
                    genantitop_mass = ptop->mass();
                    for(int i=0;ptop->daughter(i)!=NULL;i++){
                        if(abs(ptop->daughter(i)->pdgId())==24){
                            genantit_w_pt=ptop->daughter(i)->pt();
                            genantit_w_eta=ptop->daughter(i)->eta();
                            genantit_w_phi=ptop->daughter(i)->phi();
                            genantit_w_mass=ptop->daughter(i)->mass();
                            const reco::Candidate* ptw0 = ptop->daughter(i);
                            const reco::Candidate* ptw= findLastW(ptw0,24);
                            if(ptw->daughter(0)!=NULL)
                            {
                                if( abs(ptw->daughter(0)->pdgId())<=6 ){
                                    genantit_w_tag=4;
                                    genantit_w_q1_pt=ptw->daughter(0)->pt();
                                    genantit_w_q1_eta=ptw->daughter(0)->eta();
                                    genantit_w_q1_phi=ptw->daughter(0)->phi();
                                    genantit_w_q1_e=ptw->daughter(0)->energy();
                                    genantit_w_q1_pdg=ptw->daughter(0)->pdgId();
                                    genantit_w_q2_pt=ptw->daughter(1)->pt();
                                    genantit_w_q2_eta=ptw->daughter(1)->eta();
                                    genantit_w_q2_phi=ptw->daughter(1)->phi();
                                    genantit_w_q2_e=ptw->daughter(1)->energy();
                                    genantit_w_q2_pdg=ptw->daughter(1)->pdgId();
                                }
                                if( abs(ptw->daughter(0)->pdgId())==11 ||abs(ptw->daughter(0)->pdgId())==12 ) genantit_w_tag=1;
                                if( abs(ptw->daughter(0)->pdgId())==12 ||abs(ptw->daughter(0)->pdgId())==13 ) genantit_w_tag=2;
                                if( abs(ptw->daughter(0)->pdgId())==14 ||abs(ptw->daughter(0)->pdgId())==15 ) genantit_w_tag=3;
                                if( abs(ptw->daughter(0)->pdgId())==11 ||abs(ptw->daughter(0)->pdgId())==12||abs(ptw->daughter(0)->pdgId())==13 ||abs(ptw->daughter(0)->pdgId())==14||abs(ptw->daughter(0)->pdgId())==15 ||abs(ptw->daughter(0)->pdgId())==16)
                                {
                                    genantit_w_q1_pt=ptw->daughter(0)->pt();
                                    genantit_w_q1_eta=ptw->daughter(0)->eta();
                                    genantit_w_q1_phi=ptw->daughter(0)->phi();
                                    genantit_w_q1_e=ptw->daughter(0)->energy();
                                    genantit_w_q1_pdg=ptw->daughter(0)->pdgId();
                                    genantit_w_q2_pt=ptw->daughter(1)->pt();
                                    genantit_w_q2_eta=ptw->daughter(1)->eta();
                                    genantit_w_q2_phi=ptw->daughter(1)->phi();
                                    genantit_w_q2_e=ptw->daughter(1)->energy();
                                    genantit_w_q2_pdg=ptw->daughter(1)->pdgId();
                                }
                            }
                        }
                        if(abs(ptop->daughter(i)->pdgId())==5){
                            genantit_b_pt=ptop->daughter(i)->pt();
                            genantit_b_eta=ptop->daughter(i)->eta();
                            genantit_b_phi=ptop->daughter(i)->phi();
                            genantit_b_mass=ptop->daughter(i)->mass();
                        }
                    }
                }
                
    }//end of loop on gen

        
// --------------------------------------------- w -----------------------------------------------------
    //w and top info      
        int igenw=0;
        int sizew=5;
        for(size_t ik=0; ik<genParticles->size();ik++)
        {
            if(abs((*genParticles)[ik].pdgId())==24)
                {
                    const reco::Candidate* pwtmp1 = &(*genParticles)[ik];
                    const reco::Candidate* pwtmp=findLastW(pwtmp1,24);
                    int woverlap=0;
                    for (int ia=0;ia<igenw;ia++){
                        if(pwtmp->pt()==ptgenwl[ia]) woverlap=1;
                    }
                    if(pwtmp->pt()>50&&igenw<sizew&&woverlap==0){
                    ptgenwl[igenw] = pwtmp->pt();
                    etagenwl[igenw] = pwtmp->eta();
                    phigenwl[igenw] = pwtmp->phi();
                    massgenwl[igenw] = pwtmp->mass();
                    const reco::Candidate* pwtmp2=findFirstW(pwtmp1,24);
                    ptgenwf[igenw] = pwtmp2->pt();
                    etagenwf[igenw] = pwtmp2->eta();
                    phigenwf[igenw] = pwtmp2->phi();
                    massgenwf[igenw] = pwtmp2->mass();
                        taggenwmother[igenw]=pwtmp2->mother(0)->pdgId();

                    if(pwtmp->daughter(0)!=NULL)//loop on w daughter
                    {
                        const reco::Candidate* pltmp = pwtmp->daughter(0);
                         if( (abs(pltmp->pdgId())==11) || (abs(pltmp->pdgId())==12) ){
                                taggenwl[igenw]=1;                    }//end of w daugter loop
                        if( (abs(pltmp->pdgId())==13) || (abs(pltmp->pdgId())==14) ){
                            taggenwl[igenw]=2;                    }//end of w daugter loop
                        if( (abs(pltmp->pdgId())==15) || (abs(pltmp->pdgId())==16) ){
                            taggenwl[igenw]=3;                    }//end of w daugter loop
                        if( (abs(pltmp->pdgId())<6) ) {
                            taggenwl[igenw]=4;                    }
                        if( (abs(pltmp->pdgId())<6) || (abs(pltmp->pdgId())==11) || (abs(pltmp->pdgId())==12) ||(abs(pltmp->pdgId())==13) || (abs(pltmp->pdgId())==14) ||(abs(pltmp->pdgId())==15) || (abs(pltmp->pdgId())==16)){
                            genw_q1_pt[igenw]=pwtmp->daughter(0)->pt();
                            genw_q1_eta[igenw]=pwtmp->daughter(0)->eta();
                            genw_q1_phi[igenw]=pwtmp->daughter(0)->phi();
                            genw_q1_e[igenw]=pwtmp->daughter(0)->energy();
                            genw_q1_pdg[igenw]=pwtmp->daughter(0)->pdgId();
                            genw_q2_pt[igenw]=pwtmp->daughter(1)->pt();
                            genw_q2_eta[igenw]=pwtmp->daughter(1)->eta();
                            genw_q2_phi[igenw]=pwtmp->daughter(1)->phi();
                            genw_q2_e[igenw]=pwtmp->daughter(1)->energy();
                            genw_q2_pdg[igenw]=pwtmp->daughter(1)->pdgId();

                        }
                    }
                    igenw+=1;
                    }
                }//end of if w
        }

// --------------------------------------------- z ---------------------------------------

        int igenz=0;
        for(size_t ik=0; ik<genParticles->size();ik++)
        {
            if(abs((*genParticles)[ik].pdgId())==23)
            {
                const reco::Candidate* pztmp1 = &(*genParticles)[ik];
                const reco::Candidate* pztmp=findLasttau(pztmp1,23);
                int zoverlap=0;
                for (int ia=0;ia<igenz;ia++){
                    if(pztmp->pt()==ptgenzl[ia]) zoverlap=1;}
                if(pztmp->pt()>50&&igenz<sizew&&zoverlap==0){
                    ptgenzl[igenz] = pztmp->pt();
                    etagenzl[igenz] = pztmp->eta();
                    phigenzl[igenz] = pztmp->phi();
                    massgenzl[igenz] = pztmp->mass();
                    const reco::Candidate* pztmp2=findFirstW(pztmp1,23);
                    ptgenzf[igenz] = pztmp2->pt();
                    etagenzf[igenz] = pztmp2->eta();
                    phigenzf[igenz] = pztmp2->phi();
                    massgenzf[igenz] = pztmp2->mass();
                    //for(int i=0;pz->daughter(i)!=NULL;i++)//loop on w daughter
                    if(pztmp->daughter(0)!=NULL)//loop on w daughter
                    {
                        const reco::Candidate* pltmp = pztmp->daughter(0);
                        //std::cout<< "pl pdgId" << pl->pdgId() << std::endl;
                        if( (abs(pltmp->pdgId())==11) || (abs(pltmp->pdgId())==12) ){
                            taggenzl[igenz]=1;                    }//end of w daugter loop
                        if( (abs(pltmp->pdgId())==13) || (abs(pltmp->pdgId())==14) ){
                            taggenzl[igenz]=2;                    }//end of w daugter loop
                        if( (abs(pltmp->pdgId())==15) || (abs(pltmp->pdgId())==16) ){
                            taggenzl[igenz]=3;                    }//end of w daugter loop
                        if(abs(pltmp->pdgId())<6 ) {
                            taggenzl[igenz]=4;}
                    }
                    igenz+=1;
                }
            }//end of if w
        }
// ====================================== z ==========================================


    }//end of MC Info

outTree_->Fill();

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

//top
    gentop_pt  = -99;
    gentop_eta  = -99;
    gentop_phi  = -99;
    gentop_mass  = -99;
//anti top
    genantitop_pt  = -99;
    genantitop_eta  = -99;
    genantitop_phi  = -99;
    genantitop_mass  = -99;


//topw
    gent_w_pt=-99;gent_w_phi=-99;gent_w_eta=-99;gent_w_mass=-99;
    gent_w_tag=-99;
//anti top w
    genantit_w_pt=-99;genantit_w_phi=-99;genantit_w_eta=-99;genantit_w_mass=-99;
    genantit_w_tag=-99;

//top w daughter
    gent_w_q1_pt=-99;gent_w_q1_phi=-99;gent_w_q1_eta=-99;gent_w_q1_e=-99;gent_w_q1_pdg=-99;
    gent_w_q2_pt=-99;gent_w_q2_phi=-99;gent_w_q2_eta=-99;gent_w_q2_e=-99;gent_w_q2_pdg=-99;
// anti top w daughter
    genantit_w_q1_pt=-99;genantit_w_q1_phi=-99;genantit_w_q1_eta=-99;genantit_w_q1_e=-99;genantit_w_q1_pdg=-99;
    genantit_w_q2_pt=-99;genantit_w_q2_phi=-99;genantit_w_q2_eta=-99;genantit_w_q2_e=-99;genantit_w_q2_pdg=-99;


//top b
    gent_b_pt=-99;gent_b_phi=-99;gent_b_eta=-99;gent_b_mass=-99;
//anti top b
    genantit_b_pt=-99;genantit_b_phi=-99;genantit_b_eta=-99;genantit_b_mass=-99;

// ---------------------------------------- w ----------------------------------------------------
    for(int i=0;i<5;i++){
        ptgenwl[i]=-99;etagenwl[i]=-99;phigenwl[i]=-99;massgenwl[i]=-99;taggenwl[i]=-99;taggenwmother[i]=-99;
        ptgenwf[i]=-99;etagenwf[i]=-99;phigenwf[i]=-99;massgenwf[i]=-99;
        genw_q1_pt[i]=-99;genw_q1_phi[i]=-99;genw_q1_eta[i]=-99;genw_q1_e[i]=-99;genw_q1_pdg[i]=-99;
        genw_q2_pt[i]=-99;genw_q2_phi[i]=-99;genw_q2_eta[i]=-99;genw_q2_e[i]=-99;genw_q2_pdg[i]=-99;
    //  z
        ptgenzl[i]=-99;etagenzl[i]=-99;phigenzl[i]=-99;massgenzl[i]=-99;taggenzl[i]=-99;
        ptgenzf[i]=-99;etagenzf[i]=-99;phigenzf[i]=-99;massgenzf[i]=-99;

}
// ======================================== w ====================================================

}

const reco::Candidate*  hlt::findLasttau(const reco::Candidate *particle,int IDpdg){
    int iw=0;
    int pidw=0;
    const reco::Candidate* pw=particle;
    //cout<<"check 1 "<<pw->pdgId()<<"    "<<pw->status()<<"   "<<endl;
    for(int ii=0;particle->daughter(ii)!=NULL;ii++){
        if(abs(particle->daughter(ii)->pdgId())== IDpdg) {
            iw=ii;
            pidw=abs(particle->daughter(ii)->pdgId());
            //cout<<"check 2 "<<iw<<"    "<<pidw<<"   "<<endl;
        }
    }
    if( abs(pidw) == IDpdg ){
        pw = particle->daughter(iw);
        //cout<<"check 5 "<<pw->pdgId()<<"    "<<pw->status()<<"   "<<endl;

        return (findLasttau(pw,IDpdg));
    }
    //cout<<"check 3 "<<pw->pdgId()<<"    "<<pw->status()<<"   "<<pw->daughter(0)->pdgId()<<"    "<<endl;
    return pw;
}

const reco::Candidate*  hlt::findLastW(const reco::Candidate *particle,int IDpdg){
    int iw=0;
    int pidw=0;
    const reco::Candidate* pw=particle;
    //cout<<"check 1 "<<pw->pdgId()<<"    "<<pw->status()<<"   "<<endl;
    for(int ii=0;particle->daughter(ii)!=NULL;ii++){
        if(abs(particle->daughter(ii)->pdgId())>pidw) {
            iw=ii;
            pidw=abs(particle->daughter(ii)->pdgId());
            //cout<<"check 2 "<<iw<<"    "<<pidw<<"   "<<endl;
        }
    }
    if( abs(pidw) == IDpdg ){
        pw = particle->daughter(iw);
        //cout<<"check 5 "<<pw->pdgId()<<"    "<<pw->status()<<"   "<<endl;
        return (findLastW(pw,IDpdg));
    }
    //cout<<"check 3 "<<pw->pdgId()<<"    "<<pw->status()<<"   "<<pw->daughter(0)->pdgId()<<"    "<<endl;
    return pw;
}

const reco::Candidate*  hlt::findFirstW(const reco::Candidate *particle,int IDpdg){
    if (particle->mother(0)!=NULL){
        if(abs(particle->mother(0)->pdgId()) == IDpdg )
        return (findFirstW(particle->mother(0),IDpdg));
    }
    //cout<<"check 3 "<<pw->pdgId()<<"    "<<pw->status()<<"   "<<pw->daughter(0)->pdgId()<<"    "<<endl;
    return particle;
}



//define this as a plug-in
DEFINE_FWK_MODULE(hlt);

// -*- C++ -*-
//
// Package:    RPCInclusionAnalyzer
// Class:      RPCInclusionAnalyzer
//
/**\class RPCInclusionAnalyzer RPCInclusionAnalyzer.cc L1TCSCTFUpgrade/RPCInclusionAnalyzer/plugins/RPCInclusionAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Renjie Wang
//         Created:  Thu, 04 Jun 2015 16:50:38 GMT
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// includes to fetch all reguired data products from the edm::Event
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
#include "DataFormats/L1CSCTrackFinder/interface/L1CSCTrackCollection.h"
#include "DataFormats/L1CSCTrackFinder/interface/L1CSCStatusDigiCollection.h"
#include "L1Trigger/CSCCommonTrigger/interface/CSCTriggerGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "DataFormats/L1CSCTrackFinder/interface/CSCTriggerContainer.h"
#include "DataFormats/L1CSCTrackFinder/interface/TrackStub.h"

#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuRegionalCand.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTCand.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTExtendedCand.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTReadoutCollection.h"

// Sector Receiver LUT class to transform wire/strip numbers to eta/phi observables
#include "L1Trigger/CSCTrackFinder/interface/CSCSectorReceiverLUT.h"
#include "CondFormats/L1TObjects/interface/L1MuTriggerScales.h"
#include "CondFormats/DataRecord/interface/L1MuTriggerScalesRcd.h"
#include "CondFormats/L1TObjects/interface/L1MuTriggerPtScale.h"
#include "CondFormats/DataRecord/interface/L1MuTriggerPtScaleRcd.h"

#include "DataFormats/L1CSCTrackFinder/interface/L1CSCStatusDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
#include "DataFormats/L1CSCTrackFinder/interface/L1CSCTrackCollection.h"
#include "DataFormats/L1CSCTrackFinder/interface/CSCTriggerContainer.h"
#include "DataFormats/L1CSCTrackFinder/interface/TrackStub.h"

#include "TMath.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <unistd.h>

#include "L1TCSCTFUpgrade/RPCInclusionAnalyzer/interface/DataEvtSummaryHandler.h"
//
// class declaration
//

class RPCInclusionAnalyzer : public edm::EDAnalyzer {
public:
    explicit RPCInclusionAnalyzer(const edm::ParameterSet&);
    ~RPCInclusionAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;


    bool isMC_;
    DataEvtSummaryHandler summaryHandler_;

    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

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
RPCInclusionAnalyzer::RPCInclusionAnalyzer(const edm::ParameterSet& iConfig):
	isMC_(              iConfig.getParameter<bool>("isMC"))

{
    //now do what ever initialization is needed
    edm::Service<TFileService> fs;
    summaryHandler_.initTree(  fs->make<TTree>("data","Event Summary") );
    TFileDirectory baseDir=fs->mkdir(iConfig.getParameter<std::string>("dtag"));


}


RPCInclusionAnalyzer::~RPCInclusionAnalyzer()
{

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
RPCInclusionAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace std;

    summaryHandler_.resetStruct();
    //event summary to be filled
    DataEvtSummary_t &ev=summaryHandler_.getEvent();

    //event header
    ev.run    = event.id().run();
    ev.lumi   = event.luminosityBlock();
    ev.event  = event.id().event();

    summaryHandler_.fillTree();
}


// ------------ method called once each job just before starting event loop  ------------
void
RPCInclusionAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
RPCInclusionAnalyzer::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
/*
void
RPCInclusionAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
RPCInclusionAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
RPCInclusionAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
RPCInclusionAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
RPCInclusionAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RPCInclusionAnalyzer);

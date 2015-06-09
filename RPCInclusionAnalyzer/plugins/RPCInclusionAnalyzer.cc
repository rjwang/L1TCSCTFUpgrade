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


#include "L1TriggerDPGUpgrade/DataFormats/interface/L1TMuonTriggerPrimitive.h"
#include "L1TriggerDPGUpgrade/DataFormats/interface/L1TMuonTriggerPrimitiveFwd.h"
#include "L1TriggerDPGUpgrade/DataFormats/interface/L1TMuonInternalTrack.h"
#include "L1TriggerDPGUpgrade/DataFormats/interface/L1TMuonInternalTrackFwd.h"

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

//CSCTF
#include "DataFormats/L1CSCTrackFinder/interface/L1CSCTrackCollection.h"
#include "L1Trigger/CSCTrackFinder/interface/CSCTFPtLUT.h"
#include "L1Trigger/CSCCommonTrigger/interface/CSCPatternLUT.h"
#include "DataFormats/MuonDetId/interface/CSCTriggerNumbering.h"
#include "L1Trigger/CSCTrackFinder/interface/CSCSectorReceiverLUT.h"
#include "CondFormats/L1TObjects/interface/L1MuTriggerScales.h"
#include "CondFormats/L1TObjects/interface/L1MuTriggerPtScale.h"
#include "CondFormats/DataRecord/interface/L1MuTriggerScalesRcd.h"
#include "CondFormats/DataRecord/interface/L1MuTriggerPtScaleRcd.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"


#include "DataFormats/MuonDetId/interface/RPCDetId.h"

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

    int calcPhiBits(float GblPhi, int sector);
    int calcEtaBits(float GblEta);
    int sectorRPC2CSC(float rpcphi);

    bool isMC_;
    edm::InputTag RPCTPTag_, CSCTFTag_;
    DataEvtSummaryHandler summaryHandler_;

    // Needed for CSCTF LUTs
    CSCSectorReceiverLUT* srLUTs_[5][2];
    const L1MuTriggerScales  *scale;
    const L1MuTriggerPtScale *ptScale;

    // ptLUT. Input: packed pt/eta bit from track converter
    const float ptscale[33] = {
        -1.,   0.0,   1.5,   2.0,   2.5,   3.0,   3.5,   4.0,
        4.5,   5.0,   6.0,   7.0,   8.0,  10.0,  12.0,  14.0,
        16.0,  18.0,  20.0,  25.0,  30.0,  35.0,  40.0,  45.0,
        50.0,  60.0,  70.0,  80.0,  90.0, 100.0, 120.0, 140.0, 1.E6
    };

    const float etabins[16] = {
        0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
        1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4
    };

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

    RPCTPTag_      = iConfig.getParameter<edm::InputTag>("RPCTPTag");
    CSCTFTag_      = iConfig.getParameter<edm::InputTag>("CSCTFTag");


    bzero(srLUTs_ , sizeof(srLUTs_));
    int sector=1;    // assume SR LUTs are all same for every sector
    bool TMB07=true; // specific TMB firmware
    edm::ParameterSet srLUTset;
    srLUTset.addUntrackedParameter<bool>("ReadLUTs", false);
    srLUTset.addUntrackedParameter<bool>("Binary",   false);
    srLUTset.addUntrackedParameter<std::string>("LUTPath", "./");


    // positive endcap
    int endcap = 1;
    for(int station=1,fpga=0; station<=4 && fpga<5; station++) {
        if(station==1)
            for(int subSector=0; subSector<2 && fpga<5; subSector++)
                srLUTs_[fpga++][1] = new CSCSectorReceiverLUT(endcap,sector,subSector+1,
                        station, srLUTset, TMB07);
        else
            srLUTs_[fpga++][1]   = new CSCSectorReceiverLUT(endcap,  sector,   0,
                    station, srLUTset, TMB07);
    }

    // negative endcap
    endcap = 2;
    for(int station=1,fpga=0; station<=4 && fpga<5; station++) {
        if(station==1)
            for(int subSector=0; subSector<2 && fpga<5; subSector++)
                srLUTs_[fpga++][0] = new CSCSectorReceiverLUT(endcap,sector,subSector+1,
                        station, srLUTset, TMB07);
        else
            srLUTs_[fpga++][0]   = new CSCSectorReceiverLUT(endcap,  sector,   0,
                    station, srLUTset, TMB07);
    }


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


    Handle< vector<L1TMuon::TriggerPrimitive> > rpcs;
    event.getByLabel(RPCTPTag_, rpcs);


    Handle< vector<L1TMuon::InternalTrack> > CSCTFTracks;
    event.getByLabel(CSCTFTag_, CSCTFTracks);


    // Initialize CSCTF pT LUTs
    ESHandle< L1MuTriggerScales > scales;
    iSetup.get< L1MuTriggerScalesRcd >().get(scales);
    scale = scales.product();

    ESHandle< L1MuTriggerPtScale > ptscales;
    iSetup.get< L1MuTriggerPtScaleRcd >().get(ptscales);
    ptScale = ptscales.product();

//    cout << "\n\n=============== NEW EVENTS ============" << endl;

    vector<int> cluster_types;
    vector<int> cluster_regions,cluster_rings,cluster_stations,cluster_sectors,cluster_layers,cluster_subsectors,cluster_rolls,cluster_nhits;
    vector<float> cluster_phis, cluster_etas;

    for(auto rpc = rpcs->cbegin(); rpc < rpcs->cend(); rpc++) {

        RPCDetId rpc_id = rpc -> detId<RPCDetId>();

        float phi  = rpc->getCMSGlobalPhi();
        float eta  = rpc->getCMSGlobalEta();

        int region     = rpc_id.region();
        int ring       = rpc_id.ring();
        int station    = rpc_id.station();
        int sector     = rpc_id.sector();
        int layer 	= rpc->getRPCData().layer;
        int subsector  = rpc_id.subsector();
        int roll       = rpc_id.roll();

        int type=0;
        type += roll*1E0;
        type += layer*1E1;
        type += (ring+2)*1E2;
        type += subsector*1E3;
        type += sector*1E4;
        type += station*1E6;
        type += (region+1)*1E7;

        auto it = std::find(cluster_types.begin(), cluster_types.end(), type);
        if (it == cluster_types.end()) {
            //not in the list yet
            cluster_types.push_back(type);

            cluster_regions.push_back(region);
            cluster_rings.push_back(ring);
            cluster_stations.push_back(station);
            cluster_sectors.push_back(sector);
            cluster_layers.push_back(layer);
            cluster_subsectors.push_back(subsector);
            cluster_rolls.push_back(roll);
            cluster_nhits.push_back(1);

            cluster_phis.push_back(phi);
            cluster_etas.push_back(eta);

        } else {
            auto index = std::distance(cluster_types.begin(), it);
            cluster_phis[index] += phi;
            cluster_etas[index] += eta;
            cluster_nhits[index]++;
        }
        /*
                    cout << "Phi: " << phi
                         << " Eta: " << eta
                         << " region: " << (region+1)
                         << " Station: " << station
                         << " sector: " << sector
                         << " subsector: " << subsector
                         << " ring: " << ring+2
                         << " layer: " << layer
                         << " roll: " << roll
                         << " type: " << type
                         << endl;
        */

    }

    ev.rpc=0;
    for(size_t nrpc=0; nrpc<cluster_types.size(); nrpc++) {
        ev.rpc_phi[ev.rpc] = cluster_phis[nrpc]/cluster_nhits[nrpc];
        ev.rpc_eta[ev.rpc] = cluster_etas[nrpc]/cluster_nhits[nrpc];
        ev.rpc_region[ev.rpc] = cluster_regions[nrpc];
        ev.rpc_station[ev.rpc] = cluster_stations[nrpc];
        ev.rpc_sector[ev.rpc] = cluster_sectors[nrpc];
        ev.rpc_cscsector[ev.rpc] = sectorRPC2CSC(ev.rpc_phi[ev.rpc]);
        ev.rpc_phibit[ev.rpc] = calcPhiBits(ev.rpc_phi[ev.rpc],ev.rpc_cscsector[ev.rpc]);
        ev.rpc_subsector[ev.rpc] = cluster_subsectors[nrpc];
        ev.rpc_ring[ev.rpc] = cluster_rings[nrpc];
        ev.rpc_layer[ev.rpc] = cluster_layers[nrpc];
        ev.rpc_roll[ev.rpc] = cluster_rolls[nrpc];
        ev.rpc_nhits[ev.rpc] = cluster_nhits[nrpc];
        ev.rpc++;
    }






    ev.csctrk=0;
    for(auto csctftrk = CSCTFTracks->cbegin(); csctftrk < CSCTFTracks->cend(); csctftrk++) {

        // Access the track variables in bit form
        int trkPt_bit         = csctftrk -> pt_packed();
        float trkEta_bit      = csctftrk -> eta_packed();
        float trkPhi_bit      = csctftrk -> phi_packed();
        unsigned long trkMode = csctftrk -> cscMode();
        int trk_charge        = csctftrk -> chargeValue();

        // Access pt/eta in human readable form
        float trkPt  = ptscale[trkPt_bit];
        double trEta = scale -> getRegionalEtaScale(2) -> getCenter(trkEta_bit);
        double trPhi = scale -> getPhiScale() -> getLowEdge(trkPhi_bit);

        //TODO, add above 8 variables in trees


        cout << "\n\n====================== CSCTF Track ==================" << endl;
        cout << "trkPt_bit: " << trkPt_bit
             << " trkEta_bit: " << trkEta_bit
             << " trkPhi_bit: " << trkPhi_bit
             << " trkMode: " << trkMode
             << " trk_charge: " << trk_charge << endl;
        cout << "trk Pt: " << trkPt
             << " trk Eta: " << trEta
             << " trk Phi: " << trPhi
             << endl;



        auto lct_map = csctftrk -> getStubs();
        for( unsigned station = 1; station <= 4; ++station ) {

            const unsigned id = 4*L1TMuon::InternalTrack::kCSC+station-1;
            auto x_LCTs = lct_map[id];

            // Loop over lcts in each station
            ev.csclct=0;
            for ( auto t_lcts = x_LCTs.cbegin(); t_lcts != x_LCTs.cend(); t_lcts++ ) {

                auto lcts = *t_lcts; // dereference the edm:Ref object
                CSCDetId id = lcts->detId<CSCDetId>();


                // choose the objects you'd like to access
                auto trlct_station       = id.station();
                auto trlct_endcap        = id.endcap();
                auto trlct_chamber       = id.chamber();
                auto trlct_ring          = id.ring();
                double trlct_phi         = lcts->getCMSGlobalPhi();
                double trlct_eta         = lcts->getCMSGlobalEta();
                //uint16_t trlct_bx        = lcts->getCSCData().bx;
                int trlct_sector         = CSCTriggerNumbering::triggerSectorFromLabels(id)-1;
                int trlct_subsector      = CSCTriggerNumbering::triggerSubSectorFromLabels(id);
                //uint16_t trlct_bx0       = lcts->getCSCData().bx0;
                uint16_t trlct_cscID     = lcts->getCSCData().cscID;
                uint16_t trlct_strip     = lcts->getCSCData().strip;
                uint16_t trlct_pattern   = lcts->getCSCData().pattern;
                uint16_t trlct_bend      = lcts->getCSCData().bend;
                uint16_t trlct_quality   = lcts->getCSCData().quality;
                uint16_t trlct_keywire   = lcts->getCSCData().keywire;



                // Now access eta/phi in bit form used by CSCTF
                int FPGALct = ( trlct_subsector ? trlct_subsector-1 : trlct_station );
                if (trlct_endcap == 2) trlct_endcap = 0;

                lclphidat lclPhi = srLUTs_[FPGALct][trlct_endcap] -> localPhi(trlct_strip,
                                   trlct_pattern,
                                   trlct_quality,
                                   trlct_bend);

                gblphidat gblPhi = srLUTs_[FPGALct][trlct_endcap] -> globalPhiME(lclPhi.phi_local,
                                   trlct_keywire,
                                   trlct_cscID);
                /*
                        gbletadat gblEta = srLUTs_[FPGALct][trlct_endcap] -> globalEtaME(lclPhi.phi_bend_local,
                                           lclPhi.phi_local,
                                           trlct_keywire,
                                           trlct_cscID);
                */

                auto trlct_phibit = gblPhi.global_phi;
                //auto trlct_etabit = gblEta.global_eta;
                int trlct_etabit = calcEtaBits(trlct_eta);


                ev.csc_lctstation[ev.csctrk][ev.csclct] = trlct_station;
                ev.csc_lctendcap[ev.csctrk][ev.csclct] = trlct_endcap;
                ev.csc_lctchamber[ev.csctrk][ev.csclct] = trlct_chamber;
                ev.csc_lctring[ev.csctrk][ev.csclct] = trlct_ring;
                ev.csc_lctsector[ev.csctrk][ev.csclct] = trlct_sector;
                ev.csc_lctsubsector[ev.csctrk][ev.csclct] = trlct_subsector;
                ev.csc_lctcscID[ev.csctrk][ev.csclct] = trlct_cscID;
                ev.csc_lctphibit[ev.csctrk][ev.csclct] = trlct_phibit;
                ev.csc_lctetabit[ev.csctrk][ev.csclct] = trlct_etabit;

                ev.csc_lctphi[ev.csctrk][ev.csclct] = trlct_phi;
                ev.csc_lcteta[ev.csctrk][ev.csclct] = trlct_eta;


                // These are called global but are really CSCTF bits
                //const double thePhiBinning = CSCTFConstants::SECTOR_RAD/(1<<CSCBitWidths::kGlobalPhiDataBitWidth);
                //const double theEtaBinning = (CSCTFConstants::maxEta - CSCTFConstants::minEta)/(CSCTFConstants::etaBins);

                //float etaG = gblEta.global_eta*theEtaBinning + CSCTFConstants::minEta;
                //if(trlct_endcap==0) etaG *=-1.;
                //float phiG = fmod( gblPhi.global_phi*thePhiBinning+15.0*M_PI/180+(trlct_sector)*60.0*M_PI/180, 2.*M_PI );

                /*
                        cout << "phiG: " << phiG
                             << " etaG: " << etaG
                             << " trlct_phi: " << trlct_phi
                             << " trlct_eta: " << trlct_eta
                             << " trlct_endcap: " << trlct_endcap
                             << " trlct_station: " << trlct_station
                             << " trlct_ring: " << trlct_ring
                             << " trlct_chamber: " << trlct_chamber
                             << " trlct_sector: " << trlct_sector
                             << endl;
                */

                ev.csclct++;
            } //lct END
        } // station END


        ev.csctrk++;
    } // CSC track END









/*
    // rpc
    for(int rpc=0; rpc<ev.rpc; rpc++) {
        // csctrk
        for(int csctrk=0; csctrk<ev.csctrk; csctrk++) {

            // Create a vector to be filled with LCT variables.  Used by PtAddress.h
            std::vector<std::vector<int> > CSChits;
            CSChits.clear();

            for(int csclct=0; csclct<ev.csclct; csclct++) {

                vector<int> cschit;
                // fill hits vector with LCT variables
                cschit.push_back(ev.csc_lctstation[csctrk][csclct]);
                cschit.push_back(ev.csc_lctphibit[csctrk][csclct]);
                cschit.push_back(ev.csc_lcteta[csctrk][csclct]);
                cschit.push_back(ev.csc_lctsector[csctrk][csclct]);
                cschit.push_back(ev.csc_lctsubsector[csctrk][csclct]);
                cschit.push_back(ev.csc_lctcscID[csctrk][csclct]);
                cschit.push_back(ev.csc_lctchamber[csctrk][csclct]);

                // Fill vector to be used in PtAddress.h
                CSChits.push_back(cschit);




                vector<int> rpchits;
                // Fill vector to use in PtAddress.h
                rpchits.push_back(trlct_station);
                rpchits.push_back(ev.rpc_phibit[rpc]);
                rpchits.push_back(trlct_etaBits);
                rpchits.push_back(trlct_sector);
                rpchits.push_back(trlct_subsector);
                rpchits.push_back(trlct_cscID);
                rpchits.push_back(trlct_chamber);


                // Now swap out lct with matched rpc hit
                temp_cschits = csc_hits;
                temp_cschits.at(iLct) = rpchits;

                ptadd address1 = getAddress1(temp_cschits);
                ptadd address0 = getAddress0(temp_cschits);

            }
        }
    }
*/



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


// convert RPC global phi to CSC LUT phi Bits, only for RPC
int
RPCInclusionAnalyzer::calcPhiBits(float GblPhi, int sector)
{
    float phiBit = -999;
    if (sector == 1)
        phiBit = (GblPhi - 0.243) / 1.0835;
    else if (sector == 2)
        phiBit = (GblPhi - 1.2914) / 1.0835;
    else if (sector == 3) {
        if (GblPhi > 0) phiBit = (GblPhi - 2.338) / 1.0835;
        else {
            float sector_distance = abs(GblPhi + 3.1416) + (3.1416 - 2.338);
            phiBit = sector_distance / 1.0835;
        }
    } else if (sector == 4)
        phiBit = (GblPhi + 2.898) / 1.0835;
    else if (sector == 5)
        phiBit = (GblPhi + 1.8507) / 1.0835;
    else if (sector == 6) {
        if (GblPhi < 0) phiBit = (GblPhi + 0.803) / 1.0835;
        else {
            float sector_distance = GblPhi + 0.803;
            phiBit = sector_distance / 1.0835;
        }
    }

    phiBit = phiBit*4096;

    return static_cast<int>(phiBit);
}

// convert global Eta to Eta Bit for pT assignment, for both CSC and RPC
int
RPCInclusionAnalyzer::calcEtaBits(float GblEta)
{
    int etabits = 0;
    for(int iETA = 0; iETA < 15; iETA++) {
        if(fabs(GblEta) >= etabins[iETA] && fabs(GblEta) < etabins[iETA+1] ) {
            etabits = iETA;
            break;
        }
    }
    return etabits;
}

// Cluster sector must be determined manually.  Sometimes the rpc sector and the cluster sector are not the same
int
RPCInclusionAnalyzer::sectorRPC2CSC(float rpcphi)
{
    if (rpcphi >= 0.243 && rpcphi < 1.2914)
        return 1;
    else if (rpcphi >= 1.2914 && rpcphi < 2.338)
        return 2;
    else if (rpcphi >= 2.388 || rpcphi < -2.898)
        return 3;
    else if (rpcphi >= -2.898 && rpcphi < -1.8507)
        return 4;
    else if (rpcphi >= -1.8507 && rpcphi < -0.803)
        return 5;
    else if (rpcphi >= -0.803 && rpcphi < 0.243)
        return 6;
    else
	return -1;
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

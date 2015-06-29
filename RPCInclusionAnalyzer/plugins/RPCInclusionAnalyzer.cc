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

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"

#include "TMath.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <unistd.h>

#include "L1TCSCTFUpgrade/RPCInclusionAnalyzer/interface/SelectionMonitor.h"
#include "L1TCSCTFUpgrade/RPCInclusionAnalyzer/interface/TSelectionMonitor.h"
#include "L1TCSCTFUpgrade/RPCInclusionAnalyzer/interface/DataEvtSummaryHandler.h"
#include "L1TCSCTFUpgrade/RPCInclusionAnalyzer/interface/PtAddress.h"
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
    int scaling(double a);

    bool isMC_;
    edm::InputTag GenParticles_;
    edm::InputTag RPCTPTag_, CSCTFTag_;
    edm::ParameterSet LUTparam_;
    DataEvtSummaryHandler summaryHandler_;
    TSelectionMonitor controlHistos_;

    // Needed for CSCTF LUTs
    CSCSectorReceiverLUT* srLUTs_[5][2][6];
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
    isMC_(              iConfig.getParameter<bool>("isMC")),
    GenParticles_(	iConfig.getUntrackedParameter<edm::InputTag>("GenParticles")),
    controlHistos_(     iConfig.getParameter<std::string>("dtag"))
{
    //now do what ever initialization is needed
    edm::Service<TFileService> fs;
    summaryHandler_.initTree(  fs->make<TTree>("data","Event Summary") );
    TFileDirectory baseDir=fs->mkdir(iConfig.getParameter<std::string>("dtag"));

    Double_t ptbins[14] = {
        2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 16, 20, 50, 140
    };

    controlHistos_.addHistogram("nevents",";nevents; nevents",1,-0.5,0.5);

    controlHistos_.addHistogram("deltaR",";#Delta R; Event",500,0.,0.5);

    //Histograms of changes in pT
    controlHistos_.addHistogram("dpt_front",";#it{p}_{T} (Front) - #it{p}_{T} (Track) [GeV]; Events",300,-150.,150.);
    controlHistos_.addHistogram("dpt_rear",";#it{p}_{T} (Rear) - #it{p}_{T} (Track) [GeV]; Events",300,-150.,150.);
    controlHistos_.addHistogram("dpt",";#it{p}_{T} (Rear/Front) - #it{p}_{T} (Track) [GeV]; Events",100,0.,100.);
    controlHistos_.addHistogram("dpt_rpc",";#it{p}_{T} (CSC) - #it{p}_{T} (RPC Added) [GeV]; Events",300,-150.,150.);
    controlHistos_.addHistogram("dphi_rpc2_csc2_dpt",";#it{p}_{T} (CSC) - #it{p}_{T} (RPC Added) [GeV];#phi(CSC2) - #phi(RPC2) [rad]",300,-150.,150.,300,-.1,.1);

    controlHistos_.addHistogram("rpc_matches",";RPC Station; Events",4,0,4);

    //1D Histograms of dphi distributions
    controlHistos_.addHistogram("dphi_csc1_csc2_all",";#phi(CSC2) - #phi(CSC1) [rad]; Events",300,-1,1);
    controlHistos_.addHistogram("dphi_csc1_csc2",";#phi(CSC2) - #phi(CSC1) [rad]; Events",300,-1,1);

    //example
    controlHistos_.addHistogram("dphi_csc1_rpc2",";#phi(RPC2) - #phi(CSC1) [rad]; Events",300,-.2,.2);
    controlHistos_.addHistogram("dphi_rpc2_csc3",";#phi(CSC3) - #phi(RPC2) [rad]; Events",300,-.2,.2);
    controlHistos_.addHistogram("dphi_rpc2_csc2",";#phi(CSC2) - #phi(RPC2) [rad]; Events",300,-0.2,0.2);

    //2D Histograms of dphi vs pT
    controlHistos_.addHistogram("dphi_csc1_csc2_pt_all",";#it{p}_{T} (3 Hit Address) [GeV];#phi(CSC2) - #phi(CSC1) [rad]",140,0,140,300,-.2,.2);
    controlHistos_.addHistogram("dphi_csc1_csc2_pt",";#it{p}_{T} (3 Hit Address) [GeV];#phi(CSC2) - #phi(CSC1) [rad]",140,0,140,300,-.2,.2);
    controlHistos_.addHistogram("dphi_csc1_rpc2_pt",";#it{p}_{T} (RPC Address) [GeV];#phi(RPC2) - #phi(CSC1) [rad]",140,0,140,300,-.2,.2);
    controlHistos_.addHistogram("dphi_rpc2_csc2_pt",";#it{p}_{T} (RPC Address) [GeV];#phi(CSC2) - #phi(RPC2) [rad]",140,0,140,3000,-.2,.2);
    controlHistos_.addHistogram("dphi_rpc2_csc3_pt",";#it{p}_{T} (RPC Address) [GeV];#phi(CSC3) - #phi(RPC2) [rad]",140,0,140,3000,-.2,.2);

    //Gen pT distribution w/ matching
    controlHistos_.addHistogram("gen_pt_matched",";#it{p}_{T} (Gen) [GeV]; Counts",140,0,140);
    controlHistos_.addHistogram("dr_gen",";dR (Gen - Trk); Counts",300,0,1.3);

    //2D Histograms of dphi vs 1/pT
    controlHistos_.addHistogram("dphi_csc1_csc2_invpt_all",";1 / #it{p}_{T} (3 Hit Address) [1/GeV];CSC2 - CSC1",140,0,.4,300,-.2,.2);
    controlHistos_.addHistogram("dphi_csc1_csc2_invpt",";1 / #it{p}_{T} (3 Hit Address) [1/GeV];#phi(CSC2) - #phi(CSC1) [rad]",140,0,.4,300,-.2,.2);
    controlHistos_.addHistogram("dphi_csc1_rpc2_invpt",";1 / #it{p}_{T} (RPC Address) [1/GeV];#phi(RPC2) - #phi(CSC1) [rad]",140,0,.4,300,-.2,.2);
    controlHistos_.addHistogram("dphi_rpc2_csc2_invpt",";1 / #it{p}_{T} (RPC Address) [1/GeV];#phi(CSC2) - #phi(RPC2) [rad]",140,0,.4,3000,-.2,.2);
    controlHistos_.addHistogram("dphi_rpc2_csc3_invpt",";1 / #it{p}_{T} (RPC Address) [1/GeV];#phi(CSC3) - #phi(RPC2) [rad]",140,0,.4,3000,-.2,.2);

    //2D Histograms with gen pt
    controlHistos_.addHistogram("dphi_csc1_csc2_genpt_all",";#it{p}_{T} (Gen) [GeV];#phi(CSC2) - #phi(CSC1) [rad]",140,0,140,300,-.2,.2);
    controlHistos_.addHistogram("dphi_csc1_csc2_genpt",";#it{p}_{T} (Gen) [GeV];#phi(CSC2) - #phi(CSC1) [rad]",140,0,140,300,-.2,.2);
    controlHistos_.addHistogram("dphi_csc1_rpc2_genpt",";#it{p}_{T} (Gen) [GeV];#phi(RPC2) - #phi(CSC1) [rad]",140,0,140,300,-.2,.2);
    controlHistos_.addHistogram("dphi_rpc2_csc2_genpt",";#it{p}_{T} (Gen) [GeV];#phi(CSC2) - #phi(RPC2) [rad]",140,0,140,3000,-.2,.2);
    controlHistos_.addHistogram("dphi_rpc2_csc3_genpt",";#it{p}_{T} (Gen) [GeV];#phi(CSC3) - #phi(RPC2) [rad]",140,0,140,3000,-.2,.2);
    controlHistos_.addHistogram("dphi_csc1_csc2_invgenpt_all",";1 / #it{p}_{T} (Gen) [1/GeV];CSC2 - CSC1",140,0,.4,300,-.2,.2);
    controlHistos_.addHistogram("dphi_csc1_csc2_invgenpt",";1 / #it{p}_{T} (Gen) [1/GeV];#phi(CSC2) - #phi(CSC1) [rad]",140,0,.4,300,-.2,.2);
    controlHistos_.addHistogram("dphi_csc1_rpc2_invgenpt",";1 / #it{p}_{T} (Gen) [1/GeV];#phi(RPC2) - #phi(CSC1) [rad]",140,0,.4,300,-.2,.2);
    controlHistos_.addHistogram("dphi_rpc2_csc2_invgenpt",";1 / #it{p}_{T} (Gen) [1/GeV];#phi(CSC2) - #phi(RPC2) [rad]",140,0,.4,3000,-.2,.2);
    controlHistos_.addHistogram("dphi_rpc2_csc3_invgenpt",";1 / #it{p}_{T} (Gen) [1/GeV];#phi(CSC3) - #phi(RPC2) [rad]",140,0,.4,3000,-.2,.2);

    //Efficiency histograms
    controlHistos_.addHistogram("pt_turnon_threehit_all",";#it{p}_{T} (Gen) [GeV]; Counts",13,ptbins);
    //controlHistos_.addHistogram("pt_turnon_twohit_all",";#it{p}_{T} (Gen) [GeV]; Counts",140,0,140);
    controlHistos_.addHistogram("pt_turnon_rpc_all",";#it{p}_{T} (Gen) [GeV]; Counts",13,ptbins);
    controlHistos_.addHistogram("pt_turnon_threehit",";#it{p}_{T} (Gen) [GeV]; Counts",13,ptbins);
    //controlHistos_.addHistogram("pt_turnon_twohit",";#it{p}_{T} (Gen) [GeV]; Counts",140,0,140);
    controlHistos_.addHistogram("pt_turnon_rpc",";#it{p}_{T} (Gen) [GeV]; Counts",13,ptbins);

    RPCTPTag_      = iConfig.getParameter<edm::InputTag>("RPCTPTag");
    CSCTFTag_      = iConfig.getParameter<edm::InputTag>("CSCTFTag");
    LUTparam_      = iConfig.getParameter<edm::ParameterSet>("lutParam");


    bzero(srLUTs_ , sizeof(srLUTs_));
    //int sector=1;    // assume SR LUTs are all same for every sector
    bool TMB07=true; // specific TMB firmware
    edm::ParameterSet srLUTset;
    srLUTset.addUntrackedParameter<bool>("ReadLUTs", false);
    srLUTset.addUntrackedParameter<bool>("Binary",   false);
    srLUTset.addUntrackedParameter<std::string>("LUTPath", "./");


    // positive endcap
    int endcap = 1;
    for(int sector=0; sector<6; sector++) {
        for(int station=1,fpga=0; station<=4 && fpga<5; station++) {
            if(station==1)
                for(int subSector=0; subSector<2 && fpga<5; subSector++)
                    srLUTs_[fpga++][1][sector] = new CSCSectorReceiverLUT(endcap,sector+1,subSector+1,
                            station, srLUTset, TMB07);
            else
                srLUTs_[fpga++][1][sector]   = new CSCSectorReceiverLUT(endcap,  sector+1,   0,
                        station, srLUTset, TMB07);
        }
    }

    // negative endcap
    endcap = 2;
    for(int sector=0; sector<6; sector++) {
        for(int station=1,fpga=0; station<=4 && fpga<5; station++) {
            if(station==1)
                for(int subSector=0; subSector<2 && fpga<5; subSector++)
                    srLUTs_[fpga++][0][sector] = new CSCSectorReceiverLUT(endcap,sector+1,subSector+1,
                            station, srLUTset, TMB07);
            else
                srLUTs_[fpga++][0][sector]   = new CSCSectorReceiverLUT(endcap,  sector+1,   0,
                        station, srLUTset, TMB07);
        }
    }

}


RPCInclusionAnalyzer::~RPCInclusionAnalyzer()
{

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    //free the CSCTF array of pointers
    for(unsigned int j=0; j<2; j++)
        for(unsigned int i=0; i<5; i++)
            for(unsigned int s=0; s<6; s++)
                delete srLUTs_[i][j][s];
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
    using namespace reco;

    controlHistos_.fillHisto("nevents","all",0);

    summaryHandler_.resetStruct();
    //event summary to be filled
    DataEvtSummary_t &ev=summaryHandler_.getEvent();


    //event header
    ev.run    = event.id().run();
    ev.lumi   = event.luminosityBlock();
    ev.event  = event.id().event();


    //
    // gen particles
    //
    edm::Handle< std::vector<reco::GenParticle> > genParticles;
    event.getByLabel(GenParticles_, genParticles);
    if(!genParticles.isValid())     cerr << "  WARNING: genParticles is not valid! " << endl;
    else {
        ev.nmcparticles = 0;

        for(size_t i=0; i<genParticles->size(); i++) {
            const Candidate * genParticle = &(*genParticles)[i];
            int status=genParticle->status();
            int pid=genParticle->pdgId();
            if(status!=1 || abs(pid)!=13) continue; // PYTHIA 6 based

            ev.mc_px[ev.nmcparticles] = genParticle->px();
            ev.mc_py[ev.nmcparticles] = genParticle->py();
            ev.mc_pz[ev.nmcparticles] = genParticle->pz();
            ev.mc_en[ev.nmcparticles] = genParticle->energy();
            ev.mc_id[ev.nmcparticles] = genParticle->pdgId();
            ev.mc_mom[ev.nmcparticles] = 0;
            ev.mc_status[ev.nmcparticles] = genParticle->status();

            LorentzVector genCandidate(ev.mc_px[ev.nmcparticles],ev.mc_py[ev.nmcparticles],ev.mc_pz[ev.nmcparticles],ev.mc_en[ev.nmcparticles]);

            ev.mc_phi[ev.nmcparticles] = genCandidate.Phi();
            ev.mc_eta[ev.nmcparticles] = genCandidate.Eta();
            ev.mc_pt[ev.nmcparticles] = genCandidate.Pt();

            ev.nmcparticles++;
        }
    }


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

        if(region==0) continue; //remove barrel hit

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
        ev.rpc_etabit[ev.rpc] = calcEtaBits(ev.rpc_eta[ev.rpc]);
        ev.rpc_subsector[ev.rpc] = cluster_subsectors[nrpc];
        ev.rpc_ring[ev.rpc] = cluster_rings[nrpc];
        ev.rpc_layer[ev.rpc] = cluster_layers[nrpc];
        ev.rpc_roll[ev.rpc] = cluster_rolls[nrpc];
        ev.rpc_nhits[ev.rpc] = cluster_nhits[nrpc];
        ev.rpc++;
    }





    cout << "\n\n====================== CSCTF Track ==================" << endl;

    ev.csctrk=0;
    ev.csclct=0;
    for(auto csctftrk = CSCTFTracks->cbegin(); csctftrk < CSCTFTracks->cend(); csctftrk++) {

        // Access the track variables in bit form
        int trkPt_bit         = csctftrk -> pt_packed();
        float trkEta_bit      = csctftrk -> eta_packed();
        float trkPhi_bit      = csctftrk -> phi_packed();
        unsigned long trkMode = csctftrk -> cscMode();
        int trkCharge        = csctftrk -> chargeValue();

        // Access pt/eta in human readable form
        float trkPt  = ptscale[trkPt_bit];
        double trkEta = scale -> getRegionalEtaScale(2) -> getCenter(trkEta_bit);
        double trkPhi = scale -> getPhiScale() -> getLowEdge(trkPhi_bit);

        // Write variables to tree
        ev.trkPt_bit[ev.csctrk] = trkPt_bit;
        ev.trkEta_bit[ev.csctrk] = trkEta_bit;
        ev.trkPhi_bit[ev.csctrk] = trkPhi_bit;
        ev.trkMode[ev.csctrk] = trkMode;
        ev.trkCharge[ev.csctrk] = trkCharge;
        ev.trkPt[ev.csctrk] = trkPt;
        ev.trkEta[ev.csctrk] = trkEta;
        ev.trkPhi[ev.csctrk] = trkPhi;

        cout << "trkPt_bit: " << trkPt_bit
             << " trkEta_bit: " << trkEta_bit
             << " trkPhi_bit: " << trkPhi_bit
             << " trkMode: " << trkMode
             << " trkCharge: " << trkCharge
             <<  "trk Pt: " << trkPt
             << " trk Eta: " << trkEta
             << " trk Phi: " << trkPhi
             << endl;



        auto lct_map = csctftrk -> getStubs();
        for( unsigned station = 1; station <= 4; ++station ) {

            const unsigned id = 4*L1TMuon::InternalTrack::kCSC+station-1;
            auto x_LCTs = lct_map[id];

            // Loop over lcts in each station
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
                uint16_t trlct_bx        = lcts->getCSCData().bx;
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

                lclphidat lclPhi = srLUTs_[FPGALct][trlct_endcap][trlct_sector] -> localPhi(trlct_strip,
                                   trlct_pattern,
                                   trlct_quality,
                                   trlct_bend);

                gblphidat gblPhi = srLUTs_[FPGALct][trlct_endcap][trlct_sector] -> globalPhiME(lclPhi.phi_local,
                                   trlct_keywire,
                                   trlct_cscID);
                /*
                        gbletadat gblEta = srLUTs_[FPGALct][trlct_endcap][trlct_sector] -> globalEtaME(lclPhi.phi_bend_local,
                                           lclPhi.phi_local,
                                           trlct_keywire,
                                           trlct_cscID);
                */

                auto trlct_phibit = gblPhi.global_phi;
                //auto trlct_etabit = gblEta.global_eta;
                int trlct_etabit = calcEtaBits(trlct_eta);


                ev.csc_lctnthtrk[ev.csclct] = ev.csctrk;
                ev.csc_lctstation[ev.csclct] = trlct_station;
                ev.csc_lctendcap[ev.csclct] = trlct_endcap;
                ev.csc_lctchamber[ev.csclct] = trlct_chamber;
                ev.csc_lctring[ev.csclct] = trlct_ring;
                ev.csc_lctsector[ev.csclct] = trlct_sector;
                ev.csc_lctsubsector[ev.csclct] = trlct_subsector;
                ev.csc_lctcscID[ev.csclct] = trlct_cscID;
                ev.csc_lctphibit[ev.csclct] = trlct_phibit;
                ev.csc_lctetabit[ev.csclct] = trlct_etabit;

                ev.csc_lctphi[ev.csclct] = trlct_phi;
                ev.csc_lcteta[ev.csclct] = trlct_eta;


                // These are called global but are really CSCTF bits
                //const double thePhiBinning = CSCTFConstants::SECTOR_RAD/(1<<CSCBitWidths::kGlobalPhiDataBitWidth);
                //const double theEtaBinning = (CSCTFConstants::maxEta - CSCTFConstants::minEta)/(CSCTFConstants::etaBins);

                //float etaG = gblEta.global_eta*theEtaBinning + CSCTFConstants::minEta;
                //if(trlct_endcap==0) etaG *=-1.;
                //float phiG = fmod( gblPhi.global_phi*thePhiBinning+15.0*M_PI/180+(trlct_sector)*60.0*M_PI/180, 2.*M_PI );

                cout << " LCTs "
                     << " endcap: " << trlct_endcap
                     << " station: " << trlct_station
                     << " ring: " << trlct_ring
                     << " chamber: " << trlct_chamber
                     << " sector: " << trlct_sector
                     << " eta: " << trlct_eta
                     << " phi: " << trlct_phi
                     << " bx: " << trlct_bx
                     << endl;


                ev.csclct++;
            } //lct END
        } // station END


        ev.csctrk++;
    } // CSC track END







    // PT assignment

    // rpc
    for(int rpc=0; rpc<ev.rpc; rpc++) {

        // csctrk
        for(int csctrk=0; csctrk<ev.csctrk; csctrk++) {

            // Create a vector to be filled with LCT variables.  Used by PtAddress.h
            std::vector<std::vector<int> > CSChits;
            CSChits.clear();

            bool hasRPC(false);
            bool hasGen(false);
            bool rear_address(false);
            //hasRPC = true; // always use CSC hits, double-check PT assignment

            float mc_pt = -999.0,dR_Gen_Trk = 999.0;

            //Get gen info and check for matched tracks
            for (int mc = 0; mc < ev.nmcparticles; mc++) {
                float dR_Gen_Trk_temp = deltaR(ev.mc_phi[mc],ev.mc_eta[mc],ev.trkPhi[csctrk],ev.trkEta[csctrk]);
		controlHistos_.fillHisto("dr_gen","all",dR_Gen_Trk_temp);

                if (dR_Gen_Trk_temp < 1.3 && dR_Gen_Trk_temp < dR_Gen_Trk) {
                    dR_Gen_Trk = dR_Gen_Trk_temp;
                    mc_pt  = ev.mc_pt[mc];
                    hasGen = true;
                }
            }

            if (hasGen) {
                controlHistos_.fillHisto("gen_pt_matched","all",mc_pt);
            }


            //Getting the CSC hits w/o RPC info
            for (int csclct = 0; csclct < ev.csclct; csclct++) {
                vector<int> cschit;
                // fill hits vector with LCT variables
                cschit.push_back(ev.csc_lctstation[csclct]);
                cschit.push_back(ev.csc_lctphibit[csclct]);
                cschit.push_back(ev.csc_lctetabit[csclct]);
                cschit.push_back(ev.csc_lctsector[csclct]);
                cschit.push_back(ev.csc_lctsubsector[csclct]);
                cschit.push_back(ev.csc_lctcscID[csclct]);
                cschit.push_back(ev.csc_lctchamber[csclct]);

                // Fill vector to be used in PtAddress.h
                CSChits.push_back(cschit);
            }

            //Ignore tracks that aren't 3 hit with hits in 1-2-3
            if(CSChits.size()!=3) continue;
            if(CSChits[0][0] != 1 || CSChits[1][0] != 2 || CSChits[2][0] != 3) continue;

            cout << "station: " << CSChits[0][0]
                 << " " << CSChits[1][0]
                 << " " << CSChits[2][0]
                 << endl;

            double dphi12_all = ev.csc_lctphi[1] - ev.csc_lctphi[0];
            double dphi12, dphi12_rpc, dphi22_rpc, dphi23_rpc;

            //get the address to compute the original 3 hit AddressPt
            ptadd address1 = getAddress1(CSChits);
            ptadd address0 = getAddress0(CSChits);

            CSCTFPtLUT lut = CSCTFPtLUT(LUTparam_, scale, ptScale);
            float pt_front = scaling(lut.PtReal(address1));
            float pt_rear = scaling(lut.PtReal(address0));
            float pt_best = pt_front;

            cout << "      ----> Calculating pT with PtAddress.h" << endl;
            cout << "       front scaled  = " << pt_front << " dpt: " << pt_front - ev.trkPt[csctrk] << endl;
            cout << "       rear  scaled  = " << pt_rear << " dpt: " << pt_rear - ev.trkPt[csctrk] << endl;
            cout << "       Actual pT = " << ev.trkPt[csctrk] << endl;

            //check whether front or rear address gives better value
            controlHistos_.fillHisto("dpt_front","all", pt_front - ev.trkPt[csctrk]);
            controlHistos_.fillHisto("dpt_rear","all",  pt_rear  - ev.trkPt[csctrk]);
            int minpt = abs(pt_front - ev.trkPt[csctrk]);
            int minpt_rear = abs(pt_rear - ev.trkPt[csctrk]);
            if (minpt_rear < minpt) {
                minpt = minpt_rear;
                rear_address = true;
                pt_best = pt_rear;
            }

            //Fill histograms with better pt value
            controlHistos_.fillHisto("dpt","all",minpt);
            controlHistos_.fillHisto("dphi_csc1_csc2_all","all",dphi12_all);
            controlHistos_.fillHisto("dphi_csc1_csc2_pt_all","all",pt_best,dphi12_all);
            controlHistos_.fillHisto("dphi_csc1_csc2_invpt_all","all",1.0/pt_best,dphi12_all);
            //gen pt histograms
            if (hasGen) {
                controlHistos_.fillHisto("dphi_csc1_csc2_genpt_all","all",mc_pt,dphi12_all);
                controlHistos_.fillHisto("dphi_csc1_csc2_invgenpt_all","all",1.0/mc_pt,dphi12_all);
            }

            CSChits.clear();

            // all lcts
            // now replace one hit w/ rpc hit
            for(int csclct=0; csclct<ev.csclct; csclct++) {
                if(ev.csc_lctnthtrk[csclct] != csctrk) continue;

                double dR_RPC_CSC = deltaR(ev.rpc_eta[rpc], ev.rpc_phi[rpc], ev.csc_lcteta[csclct], ev.csc_lctphi[csclct]);
                controlHistos_.fillHisto("deltaR","all",dR_RPC_CSC);


                if( !hasRPC && dR_RPC_CSC<0.05 && ev.csc_lctstation[csclct] == 2 && ev.rpc_station[rpc] == 2) {

                    bool inSameEndCap(false);
                    if(ev.csc_lctendcap[csclct]==1 && ev.rpc_region[rpc] ==1 ) inSameEndCap = true; // plus endcap
                    else if(ev.csc_lctendcap[csclct]==0 && ev.rpc_region[rpc] == -1 ) inSameEndCap = true;// minus endcap
                    if(!inSameEndCap) continue;

                    vector<int> rpchits;
                    // Fill vector to use in PtAddress.h
                    rpchits.push_back(ev.rpc_station[rpc]);
                    rpchits.push_back(ev.rpc_phibit[rpc]);
                    rpchits.push_back(ev.rpc_etabit[rpc]);
                    rpchits.push_back(ev.csc_lctsector[csclct]);
                    rpchits.push_back(ev.csc_lctsubsector[csclct]);
                    rpchits.push_back(ev.csc_lctcscID[csclct]);
                    rpchits.push_back(ev.csc_lctchamber[csclct]);

                    CSChits.push_back(rpchits);
                    hasRPC = true;
                    cout << "---------------------------------------" << endl;
                    cout << "trk: " << ev.csc_lctnthtrk[csclct]
                         << " lcteta: " << ev.csc_lcteta[csclct]
                         << " lctphi: " << ev.csc_lctphi[csclct]
                         << " lctendcap: " << ev.csc_lctendcap[csclct]
                         << " lctstation: " << ev.csc_lctstation[csclct]
                         << " lctring: " << ev.csc_lctring[csclct]
                         << " lctchamber: " << ev.csc_lctchamber[csclct]
                         << endl;

                    controlHistos_.fillHisto("rpc_matches","all", ev.rpc_station[rpc]-1);

                } else {

                    vector<int> cschit;
                    // fill hits vector with LCT variables
                    cschit.push_back(ev.csc_lctstation[csclct]);
                    cschit.push_back(ev.csc_lctphibit[csclct]);
                    cschit.push_back(ev.csc_lctetabit[csclct]);
                    cschit.push_back(ev.csc_lctsector[csclct]);
                    cschit.push_back(ev.csc_lctsubsector[csclct]);
                    cschit.push_back(ev.csc_lctcscID[csclct]);
                    cschit.push_back(ev.csc_lctchamber[csclct]);

                    // Fill vector to be used in PtAddress.h
                    CSChits.push_back(cschit);
                }

            } // END LCT

            if(!hasRPC) continue;

            // Compute dphi values for various tracks
            dphi12 = ev.csc_lctphi[1] - ev.csc_lctphi[0];
            dphi12_rpc = - ev.rpc_phi[rpc] + ev.csc_lctphi[0];
            dphi22_rpc = ev.csc_lctphi[1] - ev.rpc_phi[rpc];
            dphi23_rpc = ev.csc_lctphi[2] - ev.rpc_phi[rpc];

            // Fill 1D dPhi histograms
            controlHistos_.fillHisto("dphi_csc1_csc2","all",dphi12);
            controlHistos_.fillHisto("dphi_csc1_rpc2","all",dphi12_rpc);
            controlHistos_.fillHisto("dphi_rpc2_csc2","all",dphi22_rpc);
            controlHistos_.fillHisto("dphi_rpc2_csc3","all",dphi23_rpc);

            //get addresses for computing pT values
            ptadd address1_rpc = getAddress1(CSChits);
            ptadd address0_rpc = getAddress0(CSChits);

            float pt_front_rpc = scaling(lut.PtReal(address1_rpc));
            float pt_rear_rpc = scaling(lut.PtReal(address0_rpc));
            float pt_best_rpc = pt_front_rpc;

            /*cout << "      ----> Calculating pT with PtAddress.h" << endl;
            cout << "       front scaled  = " << pt_front << " dpt: " << pt_front - ev.trkPt[csctrk] << endl;
            cout << "       rear  scaled  = " << pt_rear << " dpt: " << pt_rear - ev.trkPt[csctrk] << endl;
            cout << "       Actual pT = " << ev.trkPt[csctrk] << endl;*/

            float dpt = pt_front - pt_front_rpc;

            //check which pT is better
            if (rear_address) {
                dpt = pt_rear - pt_rear_rpc;
                pt_best_rpc = pt_rear_rpc;
            }

            //Fill 2D dPhi/pT histograms
            controlHistos_.fillHisto("dphi_csc1_csc2_pt","all",pt_best,dphi12);
            controlHistos_.fillHisto("dphi_csc1_rpc2_pt","all",pt_best_rpc,dphi12_rpc);
            controlHistos_.fillHisto("dphi_rpc2_csc2_pt","all",pt_best_rpc,dphi22_rpc);
            controlHistos_.fillHisto("dphi_rpc2_csc3_pt","all",pt_best_rpc,dphi23_rpc);
            controlHistos_.fillHisto("dphi_csc1_csc2_invpt","all",1.0/pt_best,dphi12);
            controlHistos_.fillHisto("dphi_csc1_rpc2_invpt","all",1.0/pt_best_rpc,dphi12_rpc);
            controlHistos_.fillHisto("dphi_rpc2_csc2_invpt","all",1.0/pt_best_rpc,dphi22_rpc);
            controlHistos_.fillHisto("dphi_rpc2_csc3_invpt","all",1.0/pt_best_rpc,dphi23_rpc);

            //Fill Gen pT Histograms
            if (hasGen) {
                controlHistos_.fillHisto("dphi_csc1_csc2_genpt","all",mc_pt,dphi12);
                controlHistos_.fillHisto("dphi_csc1_rpc2_genpt","all",mc_pt,dphi12_rpc);
                controlHistos_.fillHisto("dphi_rpc2_csc2_genpt","all",mc_pt,dphi22_rpc);
                controlHistos_.fillHisto("dphi_rpc2_csc3_genpt","all",mc_pt,dphi23_rpc);
                controlHistos_.fillHisto("dphi_csc1_csc2_invgenpt","all",1.0/mc_pt,dphi12);
                controlHistos_.fillHisto("dphi_csc1_rpc2_invgenpt","all",1.0/mc_pt,dphi12_rpc);
                controlHistos_.fillHisto("dphi_rpc2_csc2_invgenpt","all",1.0/mc_pt,dphi22_rpc);
                controlHistos_.fillHisto("dphi_rpc2_csc3_invgenpt","all",1.0/mc_pt,dphi23_rpc);
            }

            //Fill dpt histograms
            controlHistos_.fillHisto("dpt_rpc","all", dpt);
            controlHistos_.fillHisto("dphi_rpc2_csc2_dpt","all",dpt,dphi22_rpc);


            //Fill turn on curves for efficiency histograms
            controlHistos_.fillHisto("pt_turnon_threehit_all","all",mc_pt);
            controlHistos_.fillHisto("pt_turnon_rpc_all","all",mc_pt);

            float thresh = 16.0;

            if (pt_best > thresh) {
                controlHistos_.fillHisto("pt_turnon_threehit","all",mc_pt);
            }
            if (pt_best_rpc > thresh) {
                controlHistos_.fillHisto("pt_turnon_rpc","all",mc_pt);
            }

        } // track

    } //RPC

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

// pT scaling used with Matt's PtAddress.h
int
RPCInclusionAnalyzer::scaling(double a)//original pt scale
{
    if (a >= 0 && a < 1.5) {
        a = 0.0;
    } else if (a >= 1.5 && a < 2.0) {
        a = 1.5;
    } else if (a >= 2.0 && a < 2.5) {
        a = 2.5;
    } else if (a >= 2.5 && a < 3.0) {
        a = 3.0;
    } else if (a >= 3.0 && a < 3.5) {
        a = 3.0;
    } else if (a >= 3.5 && a < 4.0) {
        a = 3.5;
    } else if(a >= 4.0 && a < 4.5) {
        a = 4.0;
    } else if(a >= 4.5 && a < 5.0) {
        a = 4.5;
    } else if (a >= 5.0 && a < 6.0) {
        a = 5.0;
    } else if (a >= 6.0 && a < 7.0) {
        a = 6.0;
    } else if(a >= 7.0 && a < 8.0) {
        a = 7.0;
    } else if(a >= 8.0 && a < 10.0) {
        a = 8.0;
    } else if (a >= 10.0 && a < 12.0) {
        a = 10.0;
    } else if (a >= 12.0 && a < 14.0) {
        a = 12.0;
    } else if(a >= 14.0 && a < 16.0) {
        a = 14.0;
    } else if(a >= 16.0 && a < 18.0) {
        a = 16.0;
    } else if (a >= 18.0 && a < 20.0) {
        a = 18.0;
    } else if (a >= 20.0 && a < 25.0) {
        a = 20.0;
    } else if(a >= 25.0 && a < 30.0) {
        a = 25.0;
    } else if(a >= 30.0 && a < 35.0) {
        a = 30.0;
    } else if (a >= 35.0 && a < 40.0) {
        a = 35.0;
    } else if (a >= 40.0 && a < 45.0) {
        a = 40.0;
    } else if(a >= 45.0 && a < 50.0) {
        a = 45.0;
    } else if(a >= 50.0 && a < 60.0) {
        a = 50.0;
    } else if (a >= 60.0 && a < 70.0) {
        a = 60.0;
    } else if (a >= 70.0 && a < 80.0) {
        a = 70.0;
    } else if(a >= 80.0 && a < 90.0) {
        a = 80.0;
    } else if(a >= 90.0 && a < 100.0) {
        a = 90.0;
    } else if (a >= 100.0 && a < 120.0) {
        a = 100.0;
    } else if (a >= 120.0 && a < 140.0) {
        a = 120.0;
    } else if(a >= 140.0 && a < 1E6) {
        a = 140.0;
    }
    return a;
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

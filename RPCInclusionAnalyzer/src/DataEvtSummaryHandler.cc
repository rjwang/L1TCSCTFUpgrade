#include "L1TCSCTFUpgrade/RPCInclusionAnalyzer/interface/DataEvtSummaryHandler.h"

using namespace std;

//
DataEvtSummaryHandler::DataEvtSummaryHandler()
{
}

//
bool DataEvtSummaryHandler::initTree(TTree *t)
{
    if(t==0) return false;
    t_ = t;


    //event info
    t_->Branch("run",           &evSummary_.run,            "run/I");
    t_->Branch("lumi",          &evSummary_.lumi,           "lumi/I");
    t_->Branch("event",         &evSummary_.event,          "event/I");

    //rpc
    t_->Branch("rpc",           &evSummary_.rpc,             "rpc/I");
    t_->Branch("rpc_region",   evSummary_.rpc_region,                     "rpc_region/I");
    t_->Branch("rpc_station",   evSummary_.rpc_station,                     "rpc_station/I");
    t_->Branch("rpc_sector",    evSummary_.rpc_sector,                     "rpc_sector/I");
    t_->Branch("rpc_cscsector",    evSummary_.rpc_cscsector,                   "rpc_cscsector/I");
    t_->Branch("rpc_phibit",    evSummary_.rpc_phibit,                   "rpc_phibit/I");
    t_->Branch("rpc_subsector",	evSummary_.rpc_subsector,                     "rpc_subsector/I");
    t_->Branch("rpc_ring",      evSummary_.rpc_ring,                     "rpc_ring/I");
    t_->Branch("rpc_layer",     evSummary_.rpc_layer,                     "rpc_layer/I");
    t_->Branch("rpc_roll",   evSummary_.rpc_roll,                     "rpc_roll/I");
    t_->Branch("rpc_nhits",   evSummary_.rpc_nhits,                     "rpc_nhits/I");

    t_->Branch("rpc_phi",       evSummary_.rpc_phi,                     "rpc_phi/F");
    t_->Branch("rpc_eta",       evSummary_.rpc_eta,                     "rpc_eta/F");


    //trk
    t_->Branch("trkPt_bit",        evSummary_.trkPt_bit,           "trkPt_bit/I");
    t_->Branch("trkEta_bit",       evSummary_.trkEta_bit,           "trkEta_bit/F");
    t_->Branch("trkPhi_bit",       evSummary_.trkPhi_bit,           "trkPhi_bit/F");
    t_->Branch("trkMode",       evSummary_.trkMode,           "trkMode/L");
    t_->Branch("trkCharge",       evSummary_.trkCharge,           "trkCharge/I");
    t_->Branch("trkPt",        evSummary_.trkPt,           "trkPt/F");
    t_->Branch("trkEta",       evSummary_.trkEta,           "trkEta/D");
    t_->Branch("trkPhi",       evSummary_.trkPhi,           "trkPhi/D");


    //csc
    t_->Branch("csctrk",           &evSummary_.csctrk,             "csctrk/I");
    t_->Branch("csclct",           &evSummary_.csclct,             "csclct/I");

    t_->Branch("csc_lctstation",       evSummary_.csc_lctstation,                     "csc_lctstation/I");
    t_->Branch("csc_lctendcap",       evSummary_.csc_lctendcap,                     "csc_lctendcap/I");
    t_->Branch("csc_lctchamber",       evSummary_.csc_lctchamber,                     "csc_lctchamber/I");
    t_->Branch("csc_lctring",       evSummary_.csc_lctring,                     "csc_lctring/I");
    t_->Branch("csc_lctsector",       evSummary_.csc_lctsector,                     "csc_lctsector/I");
    t_->Branch("csc_lctsubsector",       evSummary_.csc_lctsubsector,                     "csc_lctsubsector/I");
    t_->Branch("csc_lctcscID",       evSummary_.csc_lctcscID,                     "csc_lctcscID/I");
    t_->Branch("csc_lctphibit",       evSummary_.csc_lctphibit,                     "csc_lctphibit/I");
    t_->Branch("csc_lctetabit",       evSummary_.csc_lctetabit,                     "csc_lctetabit/I");

    t_->Branch("csc_lctphi",       evSummary_.csc_lctphi,                     "csc_lctphi/F");
    t_->Branch("csc_lcteta",       evSummary_.csc_lcteta,                     "csc_lcteta/F");













    return true;
}

//
bool DataEvtSummaryHandler::attachToTree(TTree *t)
{
    if(t==0) return false;
    t_ = t;

    //event info
    t_->SetBranchAddress("run",             &evSummary_.run);
    t_->SetBranchAddress("lumi",            &evSummary_.lumi);
    t_->SetBranchAddress("event",           &evSummary_.event);

    //rpc
    t_->SetBranchAddress("rpc",           	  &evSummary_.rpc);
    t_->SetBranchAddress("rpc_region",           evSummary_.rpc_region);
    t_->SetBranchAddress("rpc_station",           evSummary_.rpc_station);
    t_->SetBranchAddress("rpc_sector",           evSummary_.rpc_sector);
    t_->SetBranchAddress("rpc_cscsector",           evSummary_.rpc_cscsector);
    t_->SetBranchAddress("rpc_phibit",           evSummary_.rpc_phibit);
    t_->SetBranchAddress("rpc_subsector",           evSummary_.rpc_subsector);
    t_->SetBranchAddress("rpc_ring",           evSummary_.rpc_ring);
    t_->SetBranchAddress("rpc_layer",           evSummary_.rpc_layer);
    t_->SetBranchAddress("rpc_roll",           evSummary_.rpc_roll);
    t_->SetBranchAddress("rpc_nhits",           evSummary_.rpc_nhits);

    t_->SetBranchAddress("rpc_phi",           evSummary_.rpc_phi);
    t_->SetBranchAddress("rpc_eta",           evSummary_.rpc_eta);

    //trk
    t_->SetBranchAddress("trkPt_bit",         evSummary_.trkPt_bit);
    t_->SetBranchAddress("trkEta_bit",         evSummary_.trkEta_bit);
    t_->SetBranchAddress("trkPhi_bit",         evSummary_.trkPhi_bit);
    t_->SetBranchAddress("trkMode",         evSummary_.trkMode);
    t_->SetBranchAddress("trkCharge",         evSummary_.trkCharge);
    t_->SetBranchAddress("trkPt",         evSummary_.trkPt);
    t_->SetBranchAddress("trkEta",         evSummary_.trkEta);
    t_->SetBranchAddress("trkPhi",         evSummary_.trkPhi);

    //csc
    t_->SetBranchAddress("csctrk",                   &evSummary_.csctrk);
    t_->SetBranchAddress("csclct",                   &evSummary_.csclct);

    t_->SetBranchAddress("csc_lctstation",           evSummary_.csc_lctstation);
    t_->SetBranchAddress("csc_lctendcap",           evSummary_.csc_lctendcap);
    t_->SetBranchAddress("csc_lctchamber",           evSummary_.csc_lctchamber);
    t_->SetBranchAddress("csc_lctring",           evSummary_.csc_lctring);
    t_->SetBranchAddress("csc_lctsector",           evSummary_.csc_lctsector);
    t_->SetBranchAddress("csc_lctsubsector",           evSummary_.csc_lctsubsector);
    t_->SetBranchAddress("csc_lctcscID",           evSummary_.csc_lctcscID);
    t_->SetBranchAddress("csc_lctphibit",           evSummary_.csc_lctphibit);
    t_->SetBranchAddress("csc_lctetabit",           evSummary_.csc_lctetabit);



    t_->SetBranchAddress("csc_lctphi",           evSummary_.csc_lctphi);
    t_->SetBranchAddress("csc_lcteta",           evSummary_.csc_lcteta);
















    return true;
}


//
void DataEvtSummaryHandler::resetStruct()
{
}

//
void DataEvtSummaryHandler::fillTree()
{
    if(t_) t_->Fill();
}

//
DataEvtSummaryHandler::~DataEvtSummaryHandler()
{
}

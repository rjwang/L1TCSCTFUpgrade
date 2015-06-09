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
    t_->Branch("csc_sector",    evSummary_.csc_sector,                   "csc_sector/I");
    t_->Branch("rpc_subsector",	evSummary_.rpc_subsector,                     "rpc_subsector/I");
    t_->Branch("rpc_ring",      evSummary_.rpc_ring,                     "rpc_ring/I");
    t_->Branch("rpc_layer",     evSummary_.rpc_layer,                     "rpc_layer/I");
    t_->Branch("rpc_roll",   evSummary_.rpc_roll,                     "rpc_roll/I");
    t_->Branch("rpc_nhits",   evSummary_.rpc_nhits,                     "rpc_nhits/I");

    t_->Branch("rpc_phi",       evSummary_.rpc_phi,                     "rpc_phi/F");
    t_->Branch("rpc_eta",       evSummary_.rpc_eta,                     "rpc_eta/F");









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

    t_->SetBranchAddress("rpc",           	  &evSummary_.rpc);
    t_->SetBranchAddress("rpc_region",           evSummary_.rpc_region);
    t_->SetBranchAddress("rpc_station",           evSummary_.rpc_station);
    t_->SetBranchAddress("rpc_sector",           evSummary_.rpc_sector);
    t_->SetBranchAddress("csc_sector",           evSummary_.csc_sector);
    t_->SetBranchAddress("rpc_subsector",           evSummary_.rpc_subsector);
    t_->SetBranchAddress("rpc_ring",           evSummary_.rpc_ring);
    t_->SetBranchAddress("rpc_layer",           evSummary_.rpc_layer);
    t_->SetBranchAddress("rpc_roll",           evSummary_.rpc_roll);
    t_->SetBranchAddress("rpc_nhits",           evSummary_.rpc_nhits);



    t_->SetBranchAddress("rpc_phi",           evSummary_.rpc_phi);
    t_->SetBranchAddress("rpc_eta",           evSummary_.rpc_eta);







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

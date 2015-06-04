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

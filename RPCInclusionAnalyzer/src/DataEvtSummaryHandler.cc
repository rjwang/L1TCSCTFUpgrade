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


    //mc truth
    t_->Branch("nmcparticles",  &evSummary_.nmcparticles,   "nmcparticles/I");
    t_->Branch("mc_px",         evSummary_.mc_px,           "mc_px[nmcparticles]/F");
    t_->Branch("mc_py",         evSummary_.mc_py,           "mc_py[nmcparticles]/F");
    t_->Branch("mc_pz",         evSummary_.mc_pz,           "mc_pz[nmcparticles]/F");
    t_->Branch("mc_en",         evSummary_.mc_en,           "mc_en[nmcparticles]/F");
    t_->Branch("mc_phi",        evSummary_.mc_phi,          "mc_phi[nmcparticles]/F");
    t_->Branch("mc_eta",        evSummary_.mc_eta,          "mc_eta[nmcparticles]/F");
    t_->Branch("mc_pt",         evSummary_.mc_pt,           "mc_pt[nmcparticles]/F");
    t_->Branch("mc_id",         evSummary_.mc_id,           "mc_id[nmcparticles]/I");
    t_->Branch("mc_status",     evSummary_.mc_status,       "mc_status[nmcparticles]/I");
    t_->Branch("mc_mom",        evSummary_.mc_mom,          "mc_mom[nmcparticles]/I");

    //rpc
    t_->Branch("rpc",           &evSummary_.rpc,             "rpc/I");
    t_->Branch("rpc_region",   evSummary_.rpc_region,                     "rpc_region[rpc]/I");
    t_->Branch("rpc_station",   evSummary_.rpc_station,                     "rpc_station[rpc]/I");
    t_->Branch("rpc_sector",    evSummary_.rpc_sector,                     "rpc_sector[rpc]/I");
    t_->Branch("rpc_cscsector",    evSummary_.rpc_cscsector,                   "rpc_cscsector[rpc]/I");
    t_->Branch("rpc_phibit",    evSummary_.rpc_phibit,                   "rpc_phibit[rpc]/I");
    t_->Branch("rpc_etabit",    evSummary_.rpc_etabit,                   "rpc_etabit[rpc]/I");
    t_->Branch("rpc_subsector",	evSummary_.rpc_subsector,                     "rpc_subsector[rpc]/I");
    t_->Branch("rpc_ring",      evSummary_.rpc_ring,                     "rpc_ring[rpc]/I");
    t_->Branch("rpc_layer",     evSummary_.rpc_layer,                     "rpc_layer[rpc]/I");
    t_->Branch("rpc_roll",   evSummary_.rpc_roll,                     "rpc_roll[rpc]/I");
    t_->Branch("rpc_nhits",   evSummary_.rpc_nhits,                     "rpc_nhits[rpc]/I");
    t_->Branch("rpc_phi",       evSummary_.rpc_phi,                     "rpc_phi[rpc]/F");
    t_->Branch("rpc_eta",       evSummary_.rpc_eta,                     "rpc_eta[rpc]/F");


    //csc
    t_->Branch("csctrk",           &evSummary_.csctrk,             "csctrk/I");
    t_->Branch("trkPt_bit",        evSummary_.trkPt_bit,           "trkPt_bit[csctrk]/I");
    t_->Branch("trkEta_bit",       evSummary_.trkEta_bit,           "trkEta_bit[csctrk]/F");
    t_->Branch("trkPhi_bit",       evSummary_.trkPhi_bit,           "trkPhi_bit[csctrk]/F");
    t_->Branch("trkMode",       evSummary_.trkMode,           "trkMode[csctrk]/L");
    t_->Branch("trkCharge",       evSummary_.trkCharge,           "trkCharge[csctrk]/I");
    t_->Branch("trkPt",        evSummary_.trkPt,           "trkPt[csctrk]/F");
    t_->Branch("trkEta",       evSummary_.trkEta,           "trkEta[csctrk]/D");
    t_->Branch("trkPhi",       evSummary_.trkPhi,           "trkPhi[csctrk]/D");


    t_->Branch("csclct",           &evSummary_.csclct,             "csclct/I");
    t_->Branch("csc_lctnthtrk",       evSummary_.csc_lctnthtrk,                     "csc_lctnthtrk[csclct]/I");
    t_->Branch("csc_lctstation",       evSummary_.csc_lctstation,                     "csc_lctstation[csclct]/I");
    t_->Branch("csc_lctendcap",       evSummary_.csc_lctendcap,                     "csc_lctendcap[csclct]/I");
    t_->Branch("csc_lctchamber",       evSummary_.csc_lctchamber,                     "csc_lctchamber[csclct]/I");
    t_->Branch("csc_lctring",       evSummary_.csc_lctring,                     "csc_lctring[csclct]/I");
    t_->Branch("csc_lctsector",       evSummary_.csc_lctsector,                     "csc_lctsector[csclct]/I");
    t_->Branch("csc_lctsubsector",       evSummary_.csc_lctsubsector,                     "csc_lctsubsector[csclct]/I");
    t_->Branch("csc_lctcscID",       evSummary_.csc_lctcscID,                     "csc_lctcscID[csclct]/I");
    t_->Branch("csc_lctphibit",       evSummary_.csc_lctphibit,                     "csc_lctphibit[csclct]/I");
    t_->Branch("csc_lctetabit",       evSummary_.csc_lctetabit,                     "csc_lctetabit[csclct]/I");
    t_->Branch("csc_lctphi",       evSummary_.csc_lctphi,                     "csc_lctphi[csclct]/F");
    t_->Branch("csc_lcteta",       evSummary_.csc_lcteta,                     "csc_lcteta[csclct]/F");













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

    //gen info
    t_->SetBranchAddress("nmcparticles",  &evSummary_.nmcparticles);
    t_->SetBranchAddress("mc_px",         evSummary_.mc_px);
    t_->SetBranchAddress("mc_py",         evSummary_.mc_py);
    t_->SetBranchAddress("mc_pz",         evSummary_.mc_pz);
    t_->SetBranchAddress("mc_en",         evSummary_.mc_en);
    t_->SetBranchAddress("mc_phi",        evSummary_.mc_phi);
    t_->SetBranchAddress("mc_eta",        evSummary_.mc_eta);
    t_->SetBranchAddress("mc_pt",         evSummary_.mc_pt);
    t_->SetBranchAddress("mc_id",         evSummary_.mc_id);
    t_->SetBranchAddress("mc_status",     evSummary_.mc_status);
    t_->SetBranchAddress("mc_mom",        evSummary_.mc_mom);

    //rpc
    t_->SetBranchAddress("rpc",           	  &evSummary_.rpc);
    t_->SetBranchAddress("rpc_region",           evSummary_.rpc_region);
    t_->SetBranchAddress("rpc_station",           evSummary_.rpc_station);
    t_->SetBranchAddress("rpc_sector",           evSummary_.rpc_sector);
    t_->SetBranchAddress("rpc_cscsector",           evSummary_.rpc_cscsector);
    t_->SetBranchAddress("rpc_phibit",           evSummary_.rpc_phibit);
    t_->SetBranchAddress("rpc_etabit",           evSummary_.rpc_etabit);
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

    t_->SetBranchAddress("csc_lctnthtrk",           evSummary_.csc_lctnthtrk);
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

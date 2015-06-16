#ifndef dataevtsummaryhandler_h
#define dataevtsummaryhandler_h

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <string.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <set>
#include <cmath>

#include "Math/LorentzVector.h"
#include "TMath.h"
#include "TVector2.h"
#include "TTree.h"
#include "TLorentzVector.h"

#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;
typedef std::vector<LorentzVector> LorentzVectorCollection;

#define MAXPARTICLES 500

struct DataEvtSummary_t {

    Int_t run,lumi,event;

    //rpc
    Int_t rpc;
    Int_t rpc_region[MAXPARTICLES],rpc_roll[MAXPARTICLES],rpc_nhits[MAXPARTICLES];
    Int_t rpc_station[MAXPARTICLES],rpc_sector[MAXPARTICLES],rpc_cscsector[MAXPARTICLES],rpc_phibit[MAXPARTICLES],rpc_etabit[MAXPARTICLES],rpc_subsector[MAXPARTICLES],rpc_ring[MAXPARTICLES],rpc_layer[MAXPARTICLES];
    Float_t rpc_phi[MAXPARTICLES],rpc_eta[MAXPARTICLES];

    //trk
    Int_t csctrk;
    Int_t trkPt_bit[MAXPARTICLES];
    Float_t trkEta_bit[MAXPARTICLES];
    Float_t trkPhi_bit[MAXPARTICLES];
    Long_t trkMode[MAXPARTICLES];
    Int_t trkCharge[MAXPARTICLES];
    Float_t trkPt[MAXPARTICLES];
    Double_t trkEta[MAXPARTICLES];
    Double_t trkPhi[MAXPARTICLES];

    //csc
    Int_t csclct;
    Int_t csc_lctnthtrk[MAXPARTICLES];
    Int_t csc_lctstation[MAXPARTICLES];
    Int_t csc_lctendcap[MAXPARTICLES];
    Int_t csc_lctchamber[MAXPARTICLES];
    Int_t csc_lctring[MAXPARTICLES];
    Int_t csc_lctsector[MAXPARTICLES];
    Int_t csc_lctsubsector[MAXPARTICLES];
    Int_t csc_lctcscID[MAXPARTICLES];
    Int_t csc_lctphibit[MAXPARTICLES];
    Int_t csc_lctetabit[MAXPARTICLES];
    Float_t csc_lctphi[MAXPARTICLES];
    Float_t csc_lcteta[MAXPARTICLES];









};

class DataEvtSummaryHandler {
public:
    //
    DataEvtSummaryHandler();
    ~DataEvtSummaryHandler();

    //current event
    DataEvtSummary_t evSummary_;
    DataEvtSummary_t &getEvent() {
        return evSummary_;
    }

    //write mode
    bool initTree(TTree *t);
    void fillTree();

    //read mode
    bool attachToTree(TTree *t);
    int getEntries() {
        return (t_ ? t_->GetEntriesFast() : 0);
    }
    void getEntry(int ientry) {
        resetStruct();
        if(t_) t_->GetEntry(ientry);
    }

    void resetStruct();

private:
    //the tree
    TTree *t_;
};

#endif

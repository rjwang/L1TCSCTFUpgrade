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
    Int_t rpc_station[MAXPARTICLES],rpc_sector[MAXPARTICLES],rpc_cscsector[MAXPARTICLES],rpc_phibit[MAXPARTICLES],rpc_subsector[MAXPARTICLES],rpc_ring[MAXPARTICLES],rpc_layer[MAXPARTICLES];
    Float_t rpc_phi[MAXPARTICLES],rpc_eta[MAXPARTICLES];

    //csc
    Int_t csctrk, csclct;

    Int_t csc_lctstation[MAXPARTICLES][MAXPARTICLES];
    Int_t csc_lctendcap[MAXPARTICLES][MAXPARTICLES];
    Int_t csc_lctchamber[MAXPARTICLES][MAXPARTICLES];
    Int_t csc_lctring[MAXPARTICLES][MAXPARTICLES];
    Int_t csc_lctsector[MAXPARTICLES][MAXPARTICLES];
    Int_t csc_lctsubsector[MAXPARTICLES][MAXPARTICLES];
    Int_t csc_lctcscID[MAXPARTICLES][MAXPARTICLES];
    Int_t csc_lctphibit[MAXPARTICLES][MAXPARTICLES];
    Int_t csc_lctetabit[MAXPARTICLES][MAXPARTICLES];



    Float_t csc_lctphi[MAXPARTICLES][MAXPARTICLES];
    Float_t csc_lcteta[MAXPARTICLES][MAXPARTICLES];









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

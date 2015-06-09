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

    //electron
    Int_t rpc;
    Int_t rpc_region[MAXPARTICLES],rpc_roll[MAXPARTICLES],rpc_nhits[MAXPARTICLES];
    Int_t rpc_station[MAXPARTICLES],rpc_sector[MAXPARTICLES],csc_sector[MAXPARTICLES],rpc_subsector[MAXPARTICLES],rpc_ring[MAXPARTICLES],rpc_layer[MAXPARTICLES];
    Float_t rpc_phi[MAXPARTICLES],rpc_eta[MAXPARTICLES];

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

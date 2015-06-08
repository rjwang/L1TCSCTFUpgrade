////Returns a Pt address from integer Phi and eta values to be used in
////pt assignment
////
////
////

#include <CondFormats/DataRecord/interface/L1MuTriggerScalesRcd.h>
#include <CondFormats/DataRecord/interface/L1MuTriggerPtScaleRcd.h>

using namespace std;

ptadd getAddress1(std::vector<std::vector<int>> TP)
{
    int sector_type[3] = {0,0,0};
    unsigned fr_table[16][6]= {{0,1,1,0,1,0},
        {1,0,0,1,0,1},
        {0,1,1,0,1,0},
        {0,0,1,1,1,0},
        {1,1,0,0,0,1},
        {0,0,1,1,1,0},
        {1,1,0,0,0,0},
        {0,0,1,1,0,0},
        {1,1,0,0,0,0},
        {0,0,0,0,1,0},
        {0,0,0,0,0,1},
        {0,0,0,0,1,0},
        {0,0,0,0,0,0},
        {0,0,0,0,0,0},
        {0,0,0,0,0,0},
        {0,0,0,0,0,0}
    };




    int hits[3][6];
    int count = 0,station = 0;
    for(std::vector<std::vector<int>>::iterator i = TP.begin(); i != TP.end(); i++) {

        station += (*i)[0];
        //cout << "station " << count << " = " << (*i)[0] << endl;

        for(int j=0; j<6; j++) {
            hits[count][j] = (*i)[j];
        }
        /*
        hits[count][0] = (*i)[0];//station
        hits[count][1] = (*i)[1];//integer phi
        hits[count][2] = (*i)[2];//integer eta
        hits[count][3] = (*i)[3];//sector
        hits[count][4] = (*i)[4];//subsector
        hits[count][5] = (*i)[5];//cscid
        */
        count++;
    }

    // In this example dphi12, dphi23, and eta are free parameters looped over below. You could fix any of these,
    // and allow other parameters to run free as well.

    // Mode, FR, and Sign are fixed.

    // "dPhi12" in phi-units, generally the range is 0 - 511 units
    int dphi12 = 0; // this is just a placeholder
    dphi12 = dphi12;

    // "dPhi23" in phi-units, generally the range is 0 - 255 units
    int dphi23 = 0; // this is just a placeholder
    dphi23 = dphi23;

    // track-eta, from 0.9 - 2.4
    //float eta = 1.5; // this is just a placeholder

    // Fixed Values //

    ////////// SET THE MODE ////////
    // 3-Station Modes
    // mode = 2 => ME1-2-3
    // mode = 3 => ME1-2-4
    // mode = 4 => ME1-3-4
    // mode = 5 => ME2-3-4
    // Just use 3-station modes for now. Some additional changes have to be made to
    // the code below to work with 2-station modes.
    // Also note, the Pt scaling by +20% has been disabled for the purposes of this excercise.

    unsigned track_mode = 0;
    track_mode = track_mode;
    switch(station) {

        cout << "station sum  = " << station << endl;
    case 6:
        track_mode = 2;
        break;
    case 7:
        track_mode = 3;
        break;
    case 8:
        track_mode = 4;
        break;
    case 9:
        track_mode = 5;
        break;
        //default: std::cout<<"station is out of range"<<std::endl;
    }

    //std::cout<<"Track_Mode = "<<track_mode<<"\n";

    dphi12 = fabs(hits[1][1] - hits[0][1]);
    dphi23 = fabs(hits[2][1] - hits[1][1]);


    //////// SET FR //////////////
    int fr[3] = {0,0,0};
    unsigned track_fr = 1;
    track_fr = track_fr;
    for(int i=0; i<3; i++) {

        int st = hits[i][0];
        switch(st) {

        case 1:
            sector_type[i] = 4 + hits[i][4];
            break;
        case 2:
            sector_type[i] = 1 - (hits[i][3]%2);
            break;
        case 3:
            sector_type[i] = 3 - (hits[i][3]%2);
            break;
        case 4:
            sector_type[i] = 3 - (hits[i][3]%2);
            break;

        }


        fr[i] = fr_table[(hits[i][5])-1][sector_type[i]];
        if(fr[i] == 4) {
            std::cout<<"Station = "<<st<<" and fr = "<<fr[i]<<std::endl;
        }


    }



    ///////// SET dPHI SIGN /////
    // sign = 1 => dphi12 and dphi23 have same sign
    // sign = 0 => dphi12 and dphi23 have opposite sign
    unsigned delta_phi_sign = 1;



    // Array that maps the 5-bit integer dPhi --> dPhi-units. It is assumed that this is used for dPhi23,
    // which has a maximum value of 3.83 degrees (255 units) in the extrapolation units.
    int dPhiNLBMap_5bit[32] =
    {0,1,2,4,5,7,9,11,13,15,18,21,24,28,32,37,41,47,53,60,67,75,84,94,105,117,131,145,162,180,200,222};

    // Array that maps the 7-bit integer dPhi --> dPhi-units. It is assumed that this is used for dPhi12,
    // which has a maximum value of 7.67 degrees (511 units) in the extrapolation units.
    int dPhiNLBMap_7bit[128] =
    {0,1,2,3,4,5,6,8,9,10,11,12,14,15,16,17,19,20,21,23,24,26,27,29,30,32,33,35,37,38,40,42,44,45,47,49,51,53,55,57,59,61,63,65,67,70,72,74,77,79,81,84,86,89,92,94,97,100,103,105,108,111,114,117,121,124,127,130,134,137,141,144,148,151,155,159,163,167,171,175,179,183,188,192,197,201,206,210,215,220,225,230,235,241,246,251,257,263,268,274,280,286,292,299,305,312,318,325,332,339,346,353,361,368,376,383,391,399,408,416,425,433,442,451,460,469,479,489};



    /////////// SET "dPhi12" Word ////////////
    unsigned delta_phi_12 = 0;
    int Nbins = 128;
    for(int iPHI = 0; iPHI < Nbins-1; iPHI++)
        if(fabs(dphi12) >= dPhiNLBMap_7bit[iPHI] && fabs(dphi12) < dPhiNLBMap_7bit[iPHI+1] ) {
            delta_phi_12 = iPHI;
            break;
        }

    /////////// SET "dPhi23" Word ////////////
    unsigned delta_phi_23 = 0;
    Nbins = 32;
    for(int iPHI = 0; iPHI < Nbins-1; iPHI++)
        if(fabs(dphi23) >= dPhiNLBMap_5bit[iPHI] && fabs(dphi23) < dPhiNLBMap_5bit[iPHI+1] ) {
            delta_phi_23 = iPHI;
            break;
        }

    int track_eta = hits[1][2];

    // Now make the full PTLUT Address
    ptadd address;
    // reform "dphi" portion of PTLUT due to using 7 bits + 5 bits rather than the old 8 bits + 4 bits.
    unsigned merged = ((delta_phi_12 & ((1<<7)-1)) | ((delta_phi_23 & ((1<<5)-1)) << 7 ) );
    address.delta_phi_12 = ((1<<8)-1) &  merged;
    address.delta_phi_23 = ((1<<4)-1) & (merged >> 8);
    address.track_eta = track_eta;
    address.track_mode = track_mode;
    address.track_fr = track_fr;
    address.delta_phi_sign = delta_phi_sign;


    return address;

}

ptadd getAddress0(std::vector<std::vector<int>> TP)
{
    int sector_type[3] = {0,0,0};
    unsigned fr_table[16][6]= {{0,1,1,0,1,0},
        {1,0,0,1,0,1},
        {0,1,1,0,1,0},
        {0,0,1,1,1,0},
        {1,1,0,0,0,1},
        {0,0,1,1,1,0},
        {1,1,0,0,0,0},
        {0,0,1,1,0,0},
        {1,1,0,0,0,0},
        {0,0,0,0,1,0},
        {0,0,0,0,0,1},
        {0,0,0,0,1,0},
        {0,0,0,0,0,0},
        {0,0,0,0,0,0},
        {0,0,0,0,0,0},
        {0,0,0,0,0,0}
    };




    int hits[3][6];
    int count = 0,station = 0;
    for(std::vector<std::vector<int>>::iterator i = TP.begin(); i != TP.end(); i++) {

        station += (*i)[0];

        for(int j=0; j<6; j++) {
            hits[count][j] = (*i)[j];
        }
        /*
        hits[count][0] = (*i)[0];//station
        hits[count][1] = (*i)[1];//integer phi
        hits[count][2] = (*i)[2];//integer eta
        hits[count][3] = (*i)[3];//sector
        hits[count][4] = (*i)[4];//subsector
        hits[count][5] = (*i)[5];//cscid
        */
        count++;
    }

    // In this example dphi12, dphi23, and eta are free parameters looped over below. You could fix any of these,
    // and allow other parameters to run free as well.

    // Mode, FR, and Sign are fixed.

    // "dPhi12" in phi-units, generally the range is 0 - 511 units
    int dphi12 = 0; // this is just a placeholder
    dphi12 = dphi12;
    // "dPhi23" in phi-units, generally the range is 0 - 255 units
    int dphi23 = 0; // this is just a placeholder
    dphi23 = dphi23;
    // track-eta, from 0.9 - 2.4
    //float eta = 1.5; // this is just a placeholder

    // Fixed Values //

    ////////// SET THE MODE ////////
    // 3-Station Modes
    // mode = 2 => ME1-2-3
    // mode = 3 => ME1-2-4
    // mode = 4 => ME1-3-4
    // mode = 5 => ME2-3-4
    // Just use 3-station modes for now. Some additional changes have to be made to
    // the code below to work with 2-station modes.
    // Also note, the Pt scaling by +20% has been disabled for the purposes of this excercise.

    unsigned track_mode = 0;
    track_mode = track_mode;
    switch(station) {

    case 6:
        track_mode = 2;
        break;
    case 7:
        track_mode = 3;
        break;
    case 8:
        track_mode = 4;
        break;
    case 9:
        track_mode = 5;
        break;
        //default: std::cout<<"station is out of range"<<std::endl;
    }
    //std::cout<<"Track_Mode = "<<track_mode<<"\n";

    dphi12 = fabs(hits[1][1] - hits[0][1]);
    dphi23 = fabs(hits[2][1] - hits[1][1]);

    //////// SET FR //////////////
    int fr[3] = {0,0,0};
    unsigned track_fr = 0;
    track_fr = track_fr;
    for(int i=0; i<3; i++) {

        int st = hits[i][0];
        switch(st) {

        case 1:
            sector_type[i] = 4 + hits[i][4];
            break;
        case 2:
            sector_type[i] = 1 - (hits[i][3]%2);
            break;
        case 3:
            sector_type[i] = 3 - (hits[i][3]%2);
            break;
        case 4:
            sector_type[i] = 3 - (hits[i][3]%2);
            break;

        }
        fr[i] = fr_table[(hits[i][5])-1][sector_type[i]];
        if(fr[i] == 4) {
            std::cout<<"Station = "<<st<<" and fr = "<<fr[i]<<std::endl;
        }
    }










    ///////// SET dPHI SIGN /////
    // sign = 1 => dphi12 and dphi23 have same sign
    // sign = 0 => dphi12 and dphi23 have opposite sign
    unsigned delta_phi_sign = 1;



    // Array that maps the 5-bit integer dPhi --> dPhi-units. It is assumed that this is used for dPhi23,
    // which has a maximum value of 3.83 degrees (255 units) in the extrapolation units.
    int dPhiNLBMap_5bit[32] =
    {0,1,2,4,5,7,9,11,13,15,18,21,24,28,32,37,41,47,53,60,67,75,84,94,105,117,131,145,162,180,200,222};

    // Array that maps the 7-bit integer dPhi --> dPhi-units. It is assumed that this is used for dPhi12,
    // which has a maximum value of 7.67 degrees (511 units) in the extrapolation units.
    int dPhiNLBMap_7bit[128] =
    {0,1,2,3,4,5,6,8,9,10,11,12,14,15,16,17,19,20,21,23,24,26,27,29,30,32,33,35,37,38,40,42,44,45,47,49,51,53,55,57,59,61,63,65,67,70,72,74,77,79,81,84,86,89,92,94,97,100,103,105,108,111,114,117,121,124,127,130,134,137,141,144,148,151,155,159,163,167,171,175,179,183,188,192,197,201,206,210,215,220,225,230,235,241,246,251,257,263,268,274,280,286,292,299,305,312,318,325,332,339,346,353,361,368,376,383,391,399,408,416,425,433,442,451,460,469,479,489};



    /////////// SET "dPhi12" Word ////////////
    unsigned delta_phi_12 = 0;
    int Nbins = 128;
    for(int iPHI = 0; iPHI < Nbins-1; iPHI++)
        if(fabs(dphi12) >= dPhiNLBMap_7bit[iPHI] && fabs(dphi12) < dPhiNLBMap_7bit[iPHI+1] ) {
            delta_phi_12 = iPHI;
            break;
        }

    /////////// SET "dPhi23" Word ////////////
    unsigned delta_phi_23 = 0;
    Nbins = 32;
    for(int iPHI = 0; iPHI < Nbins-1; iPHI++)
        if(fabs(dphi23) >= dPhiNLBMap_5bit[iPHI] && fabs(dphi23) < dPhiNLBMap_5bit[iPHI+1] ) {
            delta_phi_23 = iPHI;
            break;
        }

    int track_eta = hits[1][2];

    // Now make the full PTLUT Address
    ptadd address;
    // reform "dphi" portion of PTLUT due to using 7 bits + 5 bits rather than the old 8 bits + 4 bits.
    unsigned merged = ((delta_phi_12 & ((1<<7)-1)) | ((delta_phi_23 & ((1<<5)-1)) << 7 ) );
    address.delta_phi_12 = ((1<<8)-1) &  merged;
    address.delta_phi_23 = ((1<<4)-1) & (merged >> 8);
    address.track_eta = track_eta;
    address.track_mode = track_mode;
    address.track_fr = track_fr;
    address.delta_phi_sign = delta_phi_sign;


    return address;

}


#ifndef _SHELL_PAIRS_WRAPPER_H_
#define _SHELL_PAIRS_WRAPPER_H_

#include <array>
#include <vector>
#include <hls_stream.h>
#include "shell_pairs.h"

const int MAX_SHELL_PAIRS = 100;
const int MAX_PRIMITIVES = 50;

// struct shell_pairs
// {
//     int la, lb; // angular momentum
//     int na, nb; // number of shells
//     int ma, mb; // number of primitivesa
//     std::vector<double> Z; // zeta_a + zeta_b
//     std::vector<double> ZA; // zeta_a
//     std::vector<double> K; // kappa constant
//     std::vector<double> S; // Schwarz factor sPrt[(ab|ab)]
//     std::array<std::vector<int>,2> idx; // shell indices (a and b)
//     std::array<std::vector<double>,3> P; // zeta_a*A + zeta_b*B / (zeta_a + zeta_b)
//     std::array<std::vector<double>,3> PA; // P - A
//     std::array<std::vector<double>,3> AB; // A - B
// };

struct shell_pairs_hls
{
    int la, lb; // angular momentum
    int na, nb; // number of shells
    int ma, mb; // number of primitives
    const double *Z; // zeta_a + zeta_b
    const double *ZA; // zeta_a
    const double *K; // kappa constant
    const double *S; // Schwarz factor sPrt[(ab|ab)]
    const int *idx[2]; // shell indices (a and b)
    const double *P[3]; // zeta_a*A + zeta_b*B / (zeta_a + zeta_b)
    const double *PA[3]; // P - A
    const double *AB[3]; // A - B

    shell_pairs_hls *ptr;

#pragma HLS INTERFACE m_axi port=I_ offset=slave bundle=gmem1 max_read_burst_length=256 max_write_burst_length=256 depth= 1
// scaler inputs
#pragma HLS INTERFACE s_axilite port=la 
#pragma HLS INTERFACE s_axilite port=lb 
#pragma HLS INTERFACE s_axilite port=na 
#pragma HLS INTERFACE s_axilite port=nb 
#pragma HLS INTERFACE s_axilite port=ma 
#pragma HLS INTERFACE s_axilite port=mb
// array inputs
#pragma HLS INTERFACE m_axi port=Z  offset=slave bundle=gmem2 max_read_burst_length=256 max_write_burst_length=256 depth= 6  
#pragma HLS INTERFACE s_axilite port=Z 
#pragma HLS INTERFACE m_axi port=ZA offset=slave bundle=gmem3 max_read_burst_length=256 max_write_burst_length=256 depth= 6
#pragma HLS INTERFACE s_axilite port=ZA 
#pragma HLS INTERFACE m_axi port=K  offset=slave bundle=gmem4 max_read_burst_length=256 max_write_burst_length=256 depth= 6
#pragma HLS INTERFACE s_axilite port=K 
#pragma HLS INTERFACE m_axi port=S  offset=slave bundle=gmem5 max_read_burst_length=256 max_write_burst_length=256 depth= 6
#pragma HLS INTERFACE s_axilite port=S 
#pragma HLS INTERFACE m_axi port=idx  offset=slave bundle=gmem6 max_read_burst_length=256 max_write_burst_length=256 depth= 4
#pragma HLS INTERFACE s_axilite port=idx 
// 2d array inputs
#pragma HLS INTERFACE m_axi port=P  offset=slave bundle=gmem7 max_read_burst_length=256 max_write_burst_length=256 depth= 18
#pragma HLS INTERFACE s_axilite port=P
#pragma HLS INTERFACE m_axi port=PA offset=slave bundle=gmem8 max_read_burst_length=256 max_write_burst_length=256 depth= 18
#pragma HLS INTERFACE s_axilite port=PA 
#pragma HLS INTERFACE m_axi port=AB offset=slave bundle=gmem9 max_read_burst_length=256 max_write_burst_length=256 depth= 18
#pragma HLS INTERFACE s_axilite port=AB 

#pragma HLS INTERFACE s_axilite port=return

    // default constructor
    shell_pairs_hls() {
        la = 0;
        lb = 0;
        na = 0;
        nb = 0;
        ma = 0;
        mb = 0;
        Z = nullptr;
        ZA = nullptr;
        K = nullptr;
        S = nullptr;
        idx[0] = nullptr;
        idx[1] = nullptr;
        P[0] = nullptr;
        P[1] = nullptr;
        P[2] = nullptr;
        PA[0] = nullptr;
        PA[1] = nullptr;
        PA[2] = nullptr;
        AB[0] = nullptr;
        AB[1] = nullptr;
        AB[2] = nullptr;
    }

    // copy Constructor
    shell_pairs_hls(const shell_pairs &og):la(og.la), lb(og.lb), na(og.na), nb(og.nb), ma(og.ma), mb(og.mb), Z(og.Z.data()), ZA(og.ZA.data()), K(og.K.data()), S(og.S.data()) {
        idx[0] = og.idx[0].data();
        idx[1] = og.idx[1].data();
        P[0] = og.P[0].data();
        P[1] = og.P[1].data();
        P[2] = og.P[2].data();
        PA[0] = og.PA[0].data();
        PA[1] = og.PA[1].data();
        PA[2] = og.PA[2].data();
        AB[0] = og.AB[0].data();
        AB[1] = og.AB[1].data();
        AB[2] = og.AB[2].data();
    }
    
    // destructor
    ~shell_pairs_hls() {
        if (ptr != nullptr) {
            delete ptr;
        }
    }

    operator shell_pairs_hls*() {
        if (ptr == nullptr) {
            ptr = new shell_pairs_hls();
        }
        return ptr;
    }

};


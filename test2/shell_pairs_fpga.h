#ifndef _SHELL_PAIRS_WRAPPER_H_
#define _SHELL_PAIRS_WRAPPER_H_

#include <array>
#include <vector>
#include <hls_stream.h>

const int MAX_SHELL_PAIRS = 100;
const int MAX_PRIMITIVES = 50;

// struct shell_pairs_original
// {
//     int la, lb; // angular momentum
//     int na, nb; // number of shells
//     int ma, mb; // number of primitives
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
    double Z[MAX_SHELL_PAIRS]; // zeta_a + zeta_b
    double ZA[MAX_SHELL_PAIRS]; // zeta_a
    double K[MAX_SHELL_PAIRS]; // kappa constant
    double S[MAX_SHELL_PAIRS]; // Schwarz factor sPrt[(ab|ab)]
    int idx[2][MAX_SHELL_PAIRS]; // shell indices (a and b)
    double P[3][MAX_SHELL_PAIRS]; // zeta_a*A + zeta_b*B / (zeta_a + zeta_b)
    double PA[3][MAX_SHELL_PAIRS]; // P - A
    double AB[3][MAX_SHELL_PAIRS]; // A - B

    // Constructor
    shell_pairs_hls() {
#pragma HLS ARRAY_PARTITION variable=Z complete
#pragma HLS ARRAY_PARTITION variable=ZA complete
#pragma HLS ARRAY_PARTITION variable=K complete
#pragma HLS ARRAY_PARTITION variable=S complete
#pragma HLS ARRAY_PARTITION variable=idx complete
#pragma HLS ARRAY_PARTITION variable=P complete
#pragma HLS ARRAY_PARTITION variable=PA complete
#pragma HLS ARRAY_PARTITION variable=AB complete
    }
};

class ShellPairsWrapper {
public:
    shell_pairs_hls sp_hls;

    ShellPairsWrapper() {
        init_shell_pairs_hls();
    }

    void init_shell_pairs_hls() {
        sp_hls.la = 0;
        sp_hls.lb = 0;
        sp_hls.na = 0;
        sp_hls.nb = 0;
        sp_hls.ma = 0;
        sp_hls.mb = 0;
        
        for (int i = 0; i < MAX_SHELL_PAIRS; ++i) {
#pragma HLS UNROLL
            sp_hls.Z[i] = 0.0;
            sp_hls.ZA[i] = 0.0;
            sp_hls.K[i] = 0.0;
            sp_hls.S[i] = 0.0;
            sp_hls.idx[0][i] = 0;
            sp_hls.idx[1][i] = 0;
            sp_hls.P[0][i] = 0.0;
            sp_hls.P[1][i] = 0.0;
            sp_hls.P[2][i] = 0.0;
            sp_hls.PA[0][i] = 0.0;
            sp_hls.PA[1][i] = 0.0;
            sp_hls.PA[2][i] = 0.0;
            sp_hls.AB[0][i] = 0.0;
            sp_hls.AB[1][i] = 0.0;
            sp_hls.AB[2][i] = 0.0;
        }
    }

    void populate_from_original(const shell_pairs_original& sp_orig) {
        sp_hls.la = sp_orig.la;
        sp_hls.lb = sp_orig.lb;
        sp_hls.na = sp_orig.na;
        sp_hls.nb = sp_orig.nb;
        sp_hls.ma = sp_orig.ma;
        sp_hls.mb = sp_orig.mb;

        int size = std::min(static_cast<int>(sp_orig.Z.size()), MAX_SHELL_PAIRS);
        
        for (int i = 0; i < size; ++i) {
#pragma HLS UNROLL
            sp_hls.Z[i] = sp_orig.Z[i];
            sp_hls.ZA[i] = sp_orig.ZA[i];
            sp_hls.K[i] = sp_orig.K[i];
            sp_hls.S[i] = sp_orig.S[i];
            sp_hls.idx[0][i] = sp_orig.idx[0][i];
            sp_hls.idx[1][i] = sp_orig.idx[1][i];
            sp_hls.P[0][i] = sp_orig.P[0][i];
            sp_hls.P[1][i] = sp_orig.P[1][i];
            sp_hls.P[2][i] = sp_orig.P[2][i];
            sp_hls.PA[0][i] = sp_orig.PA[0][i];
            sp_hls.PA[1][i] = sp_orig.PA[1][i];
            sp_hls.PA[2][i] = sp_orig.PA[2][i];
            sp_hls.AB[0][i] = sp_orig.AB[0][i];
            sp_hls.AB[1][i] = sp_orig.AB[1][i];
            sp_hls.AB[2][i] = sp_orig.AB[2][i];
        }
    }

    // AXI Stream function to write data
    void write_axi_stream(hls::stream<shell_pairs_hls>& stream) {
#pragma HLS INLINE
        stream.write(sp_hls);
    }

    // AXI Stream function to read data
    void read_axi_stream(hls::stream<shell_pairs_hls>& stream) {
#pragma HLS INLINE
        sp_hls = stream.read();
    }

};

#endif // _SHELL_PAIRS_WRAPPER_H_

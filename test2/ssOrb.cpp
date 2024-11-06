#include "ssOrb.h"
#include <cmath>
#include <array>

void hrr_ssss(  int la, int lb, int na, int nb, int ma, int mb, std::array<double, 6> abZ, std::array<double, 6> abZA, std::array<double, 6> abK, std::array<double, 6> abS, std::array<std::array<int,2>,2> abidx, std::array<std::array<double,6>,3> abP, std::array<std::array<double,6>,3> abPA, std::array<std::array<double,6>,3> abAB,
                int lc, int ld, int nc, int nd, int mc, int md, std::array<double, 6> cdZ, std::array<double, 6> cdZA, std::array<double, 6> cdK, std::array<double, 6> cdS, std::array<std::array<int,2>,2> cdidx, std::array<std::array<double,6>,3> cdP, std::array<std::array<double,6>,3> cdPA, std::array<std::array<double,6>,3> cdAB,
                double* I_ ) {

#pragma HLS INTERFACE m_axi port=I_ offset=slave bundle=gmem1 max_read_burst_length=256 max_write_burst_length=256 depth= 1

// scaler inputs
// ab
#pragma HLS INTERFACE s_axilite port=la 
#pragma HLS INTERFACE s_axilite port=lb 
#pragma HLS INTERFACE s_axilite port=na 
#pragma HLS INTERFACE s_axilite port=nb 
#pragma HLS INTERFACE s_axilite port=ma 
#pragma HLS INTERFACE s_axilite port=mb
// cd
#pragma HLS INTERFACE s_axilite port=lc 
#pragma HLS INTERFACE s_axilite port=ld 
#pragma HLS INTERFACE s_axilite port=nc 
#pragma HLS INTERFACE s_axilite port=nd 
#pragma HLS INTERFACE s_axilite port=mc 
#pragma HLS INTERFACE s_axilite port=md

// array inputs
// ab
#pragma HLS INTERFACE m_axi port=abZ  offset=slave bundle=gmem2 max_read_burst_length=256 max_write_burst_length=256 depth= 6  
#pragma HLS INTERFACE s_axilite port=abZ 
#pragma HLS INTERFACE m_axi port=abZA offset=slave bundle=gmem3 max_read_burst_length=256 max_write_burst_length=256 depth= 6
#pragma HLS INTERFACE s_axilite port=abZA 
#pragma HLS INTERFACE m_axi port=abK  offset=slave bundle=gmem4 max_read_burst_length=256 max_write_burst_length=256 depth= 6
#pragma HLS INTERFACE s_axilite port=abK 
#pragma HLS INTERFACE m_axi port=abS  offset=slave bundle=gmem5 max_read_burst_length=256 max_write_burst_length=256 depth= 6
#pragma HLS INTERFACE s_axilite port=abS 
#pragma HLS INTERFACE m_axi port=abidx  offset=slave bundle=gmem6 max_read_burst_length=256 max_write_burst_length=256 depth= 4
#pragma HLS INTERFACE s_axilite port=abidx 
// cd
#pragma HLS INTERFACE m_axi port=cdZ  offset=slave bundle=gmem7 max_read_burst_length=256 max_write_burst_length=256 depth= 6  
#pragma HLS INTERFACE s_axilite port=cdZ 
#pragma HLS INTERFACE m_axi port=cdZA offset=slave bundle=gmem8 max_read_burst_length=256 max_write_burst_length=256 depth= 6
#pragma HLS INTERFACE s_axilite port=cdZA 
#pragma HLS INTERFACE m_axi port=cdK  offset=slave bundle=gmem9 max_read_burst_length=256 max_write_burst_length=256 depth= 6
#pragma HLS INTERFACE s_axilite port=cdK 
#pragma HLS INTERFACE m_axi port=cdS  offset=slave bundle=gmem10 max_read_burst_length=256 max_write_burst_length=256 depth= 6
#pragma HLS INTERFACE s_axilite port=cdS 
#pragma HLS INTERFACE m_axi port=cdidx  offset=slave bundle=gmem11 max_read_burst_length=256 max_write_burst_length=256 depth= 4
#pragma HLS INTERFACE s_axilite port=cdidx 

// 2d array inputs
// ab
#pragma HLS INTERFACE m_axi port=adP  offset=slave bundle=gmem12 max_read_burst_length=256 max_write_burst_length=256 depth= 18
#pragma HLS INTERFACE s_axilite port=adP
#pragma HLS INTERFACE m_axi port=adPA offset=slave bundle=gmem13 max_read_burst_length=256 max_write_burst_length=256 depth= 18
#pragma HLS INTERFACE s_axilite port=adPA 
#pragma HLS INTERFACE m_axi port=adAB offset=slave bundle=gmem14 max_read_burst_length=256 max_write_burst_length=256 depth= 18
#pragma HLS INTERFACE s_axilite port=adAB 
// cd
#pragma HLS INTERFACE m_axi port=cdP  offset=slave bundle=gmem15 max_read_burst_length=256 max_write_burst_length=256 depth= 18
#pragma HLS INTERFACE s_axilite port=cdP
#pragma HLS INTERFACE m_axi port=cdPA offset=slave bundle=gmem16 max_read_burst_length=256 max_write_burst_length=256 depth= 18
#pragma HLS INTERFACE s_axilite port=cdPA 
#pragma HLS INTERFACE m_axi port=cdAB offset=slave bundle=gmem17 max_read_burst_length=256 max_write_burst_length=256 depth= 18
#pragma HLS INTERFACE s_axilite port=cdAB

#pragma HLS INTERFACE s_axilite port=return

    int AB_la;
    int AB_lb;
    int AB_na;
    int AB_nb;
    int AB_ma;
    int AB_mb; 
    double AB_Z[6]; // zeta_a + zeta_b
    double AB_ZA[6]; // zeta_a
    double AB_K[6]; // kappa constant
    double AB_S[6]; // Schwarz factor sPrt[(ab|ab)]
    int AB_idx[2][2]; // shell indices (a and b)
    double AB_P[3][6]; // zeta_a*A + zeta_b*B / (zeta_a + zeta_b)
    double AB_PA[3][6]; // P - A
    double AB_AB[3][6]; // A - B

    int CD_la;
    int CD_lb;
    int CD_na;
    int CD_nb;
    int CD_ma;
    int CD_mb; 
    double CD_Z[6]; // zeta_a + zeta_b
    double CD_ZA[6]; // zeta_a
    double CD_K[6]; // kappa constant
    double CD_S[6]; // Schwarz factor sPrt[(ab|ab)]
    int CD_idx[2][2]; // shell indices (a and b)
    double CD_P[3][6]; // zeta_a*A + zeta_b*B / (zeta_a + zeta_b)
    double CD_PA[3][6]; // P - A
    double CD_AB[3][6]; // A - B
    
    AB_la = la;
    AB_lb = lb;
    AB_na = na;
    AB_nb = nb;
    AB_ma = ma;
    AB_mb = mb;
    
    CD_la = lc;
    CD_lb = ld;
    CD_na = nc;
    CD_nb = nd;
    CD_ma = mc;
    CD_mb = md;

    for (int i = 0; i < 6; i++) {
        AB_Z[i] = abZ[i];
        AB_ZA[i] = abZA[i];
        AB_K[i] = abK[i];
        AB_S[i] = abS[i];

        CD_Z[i] = cdZ[i];
        CD_ZA[i] = cdZA[i];
        CD_K[i] = cdK[i];
        CD_S[i] = cdS[i];
    }

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            AB_idx[i][j] = abidx[i][j];
            CD_idx[i][j] = cdidx[i][j];
        }
    }

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 6; j++) {
            AB_P[i][j] = abP[i][j];
            AB_PA[i][j] = abPA[i][j];
            AB_AB[i][j] = abAB[i][j];

            CD_P[i][j] = cdP[i][j];
            CD_PA[i][j] = cdPA[i][j];
            CD_AB[i][j] = cdAB[i][j];
        }
    }


    
    int nab = AB_na * AB_nb;
    int ncd = CD_na * CD_nb; 
    for (int idx = 0; idx < nab; idx++) {
        for (int idy = 0; idy < ncd; idy++) {
            #pragma HLS PIPELINE II=1

            double AB[3] = {
                AB_AB[0][idx],
                AB_AB[1][idx],
                AB_AB[2][idx]
            };
            double CD[3] = {
                CD_AB[0][idy],
                CD_AB[1][idy],
                CD_AB[2][idy]
            };

            double s000_s000_s000_s000_0_con = 0.0;

            int mab = AB_ma * AB_mb;
            int mcd = CD_ma * CD_mb;

            for (int pab = 0; pab < mab; pab++) {
                for (int pcd = 0; pcd < mcd; pcd++) {
                    double iab = idx + pab * nab;
                    double icd = idy + pcd * ncd;

                    double zab = AB_Z[(int)iab];
                    double zcd = CD_Z[(int)icd];

                    double zab_inv = 1 / zab;
                    double zcd_inv = 1 / zcd;
                    double zabcd_inv_sqrt = 1.0 / (sqrt(zab + zcd));
                    double zabcd_inv = zabcd_inv_sqrt * zabcd_inv_sqrt;
                    double rho = zab * zcd * zabcd_inv;

                    double P[3] = {
                        AB_P[0][(int)iab],
                        AB_P[1][(int)iab],
                        AB_P[2][(int)iab],
                    };

                    double Q[3] = {
                        CD_P[0][(int)icd],
                        CD_P[1][(int)icd],
                        CD_P[2][(int)icd],
                    };

                    double T = rho * ((P[0] - Q[0]) * (P[0] - Q[0]) + (P[1] - Q[1]) * (P[1] - Q[1]) + (P[2] - Q[2]) * (P[2] - Q[2]));
                    double K = zabcd_inv_sqrt * (AB_K[(int)iab]) * CD_K[(int)icd];

                    double fm[1];
                    // vgamma<double>(0, T, fm);

                    for (int i = 0; i < 1; i++) {
                        fm[i] *= K;
                    }

                    double s000_s000_s000_s000_0 = fm[0];
                    s000_s000_s000_s000_0_con += s000_s000_s000_s000_0;
                }
            }

            I_[idx + idy * (int)nab] = s000_s000_s000_s000_0_con;
        }
    }
}




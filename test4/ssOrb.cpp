#include "ssOrb.h"
#include <cmath>
#include <array>

void hrr_ssss(  int la, int lb, int na, int nb, int ma, int mb, std::array<double, 6> Z, std::array<double, 6> ZA, std::array<double, 6> K, std::array<double, 6> S, std::array<std::array<int,2>,2> idx, std::array<std::array<double,6>,3> P, std::array<std::array<double,6>,3> PA, std::array<std::array<double,6>,3> AB,
                double* I_ ) {

#pragma HLS INTERFACE m_axi port=I_ offset=slave bundle=gmem1 max_read_burst_length=256 max_write_burst_length=256 depth= 1
// ab
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
    CD_la = la;
    CD_lb = lb;
    CD_na = na;
    CD_nb = nb;
    CD_ma = ma;
    CD_mb = mb;

    for (int i = 0; i < 6; i++) {
        AB_Z[i] = Z[i];
        AB_ZA[i] = ZA[i];
        AB_K[i] = K[i];
        AB_S[i] = S[i];

        CD_Z[i] = Z[i];
        CD_ZA[i] = ZA[i];
        CD_K[i] = K[i];
        CD_S[i] = S[i];
    }

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            AB_idx[i][j] = idx[i][j];
            CD_idx[i][j] = idx[i][j];
        }
    }

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 6; j++) {
            AB_P[i][j] = P[i][j];
            AB_PA[i][j] = PA[i][j];
            AB_AB[i][j] = AB[i][j];

            CD_P[i][j] = P[i][j];
            CD_PA[i][j] = PA[i][j];
            CD_AB[i][j] = AB[i][j];
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


void hrr_psss(  int la, int lb, int na, int nb, int ma, int mb, std::array<double, 6> Z, std::array<double, 6> ZA, std::array<double, 6> K, std::array<double, 6> S, std::array<std::array<int,2>,2> idx, std::array<std::array<double,6>,3> P, std::array<std::array<double,6>,3> PA, std::array<std::array<double,6>,3> AB,
                double* I_ ) {

#pragma HLS INTERFACE m_axi port=I_ offset=slave bundle=gmem1 max_read_burst_length=256 max_write_burst_length=256 depth= 1
// ab
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
    CD_la = la;
    CD_lb = lb;
    CD_na = na;
    CD_nb = nb;
    CD_ma = ma;
    CD_mb = mb;

    for (int i = 0; i < 6; i++) {
        AB_Z[i] = Z[i];
        AB_ZA[i] = ZA[i];
        AB_K[i] = K[i];
        AB_S[i] = S[i];

        CD_Z[i] = Z[i];
        CD_ZA[i] = ZA[i];
        CD_K[i] = K[i];
        CD_S[i] = S[i];
    }

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            AB_idx[i][j] = idx[i][j];
            CD_idx[i][j] = idx[i][j];
        }
    }

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 6; j++) {
            AB_P[i][j] = P[i][j];
            AB_PA[i][j] = PA[i][j];
            AB_AB[i][j] = AB[i][j];

            CD_P[i][j] = P[i][j];
            CD_PA[i][j] = PA[i][j];
            CD_AB[i][j] = AB[i][j];
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

            double p100_s000_s000_s000_0_con = 0.0;
            double p010_s000_s000_s000_0_con = 0.0;
            double p001_s000_s000_s000_0_con = 0.0;

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

                    double PA[3] = {
                        AB_PA[0][(int)iab],
                        AB_PA[1][(int)iab],
                        AB_PA[2][(int)iab],
                    };

                    double PB[3] = {
                        AB_PA[0][(int)iab] + AB_AB[0][(int)iab],
                        AB_PA[1][(int)iab] + AB_AB[1][(int)iab],
                        AB_PA[2][(int)iab] + AB_AB[2][(int)iab],
                    };

                    double QC[3] = {
                        CD_PA[0][(int)icd],
                        CD_PA[1][(int)icd],
                        CD_PA[2][(int)icd],
                    };

                    double QD[3] = {
                        CD_PA[0][(int)icd] + CD_AB[0][(int)icd],
                        CD_PA[1][(int)icd] + CD_AB[1][(int)icd],
                        CD_PA[2][(int)icd] + CD_AB[2][(int)icd],
                    };

                    double W[3] = {
                        (zab * P[0] + zcd * Q[0]) * zabcd_inv,
                        (zab * P[1] + zcd * Q[1]) * zabcd_inv,
                        (zab * P[2] + zcd * Q[2]) * zabcd_inv,
                    };

                    double WP[3] = {
                        W[0] - PA[0],
                        W[1] - PA[1],
                        W[2] - PA[2],
                    };

                    double WQ[3] = {
                        W[0] - Q[0],
                        W[1] - Q[1],
                        W[2] - Q[2],
                    };

                    double T = rho * ((P[0] - Q[0]) * (P[0] - Q[0]) + (P[1] - Q[1]) * (P[1] - Q[1]) + (P[2] - Q[2]) * (P[2] - Q[2]));
                    double K = zabcd_inv_sqrt * (AB_K[(int)iab]) * CD_K[(int)icd];

                    double fm[2];
                    // vgamma<double>(1, T, fm);

                    for (int i = 0; i < 2; i++) {
                        fm[i] *= K;
                    }

                    double s000_s000_s000_s000_0 = fm[0];
                    double s000_s000_s000_s000_1 = fm[1] ;
                    double p100_s000_s000_s000_0 = PA[0] * s000_s000_s000_s000_0 + WP[0] * s000_s000_s000_s000_1;
                    double p010_s000_s000_s000_0 = PA[1] * s000_s000_s000_s000_0 + WP[1] * s000_s000_s000_s000_1;
                    double p001_s000_s000_s000_0 = PA[2] * s000_s000_s000_s000_0 + WP[2] * s000_s000_s000_s000_1;
                    p100_s000_s000_s000_0_con += p100_s000_s000_s000_0;
                    p010_s000_s000_s000_0_con += p010_s000_s000_s000_0;
                    p001_s000_s000_s000_0_con += p001_s000_s000_s000_0;
                }
            }

            I_[idx + (0 + 0 * 3 + 0 * 3 * 1 + 0 * 3 * 1 * 1) * nab + idy * nab * 3] =                 1 * p100_s000_s000_s000_0_con ;
            I_[idx + (1 + 0 * 3 + 0 * 3 * 1 + 0 * 3 * 1 * 1) * nab + idy * nab * 3] =                 1 * p010_s000_s000_s000_0_con ;
            I_[idx + (2 + 0 * 3 + 0 * 3 * 1 + 0 * 3 * 1 * 1) * nab + idy * nab * 3] =                 1 * p001_s000_s000_s000_0_con ;
        }
    }   
}


void hrr_psps(  int la, int lb, int na, int nb, int ma, int mb, std::array<double, 6> Z, std::array<double, 6> ZA, std::array<double, 6> K, std::array<double, 6> S, std::array<std::array<int,2>,2> idx, std::array<std::array<double,6>,3> P, std::array<std::array<double,6>,3> PA, std::array<std::array<double,6>,3> AB,
                double* I_ ) {

#pragma HLS INTERFACE m_axi port=I_ offset=slave bundle=gmem1 max_read_burst_length=256 max_write_burst_length=256 depth= 1
// ab
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
    CD_la = la;
    CD_lb = lb;
    CD_na = na;
    CD_nb = nb;
    CD_ma = ma;
    CD_mb = mb;

    for (int i = 0; i < 6; i++) {
        AB_Z[i] = Z[i];
        AB_ZA[i] = ZA[i];
        AB_K[i] = K[i];
        AB_S[i] = S[i];

        CD_Z[i] = Z[i];
        CD_ZA[i] = ZA[i];
        CD_K[i] = K[i];
        CD_S[i] = S[i];
    }

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            AB_idx[i][j] = idx[i][j];
            CD_idx[i][j] = idx[i][j];
        }
    }

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 6; j++) {
            AB_P[i][j] = P[i][j];
            AB_PA[i][j] = PA[i][j];
            AB_AB[i][j] = AB[i][j];

            CD_P[i][j] = P[i][j];
            CD_PA[i][j] = PA[i][j];
            CD_AB[i][j] = AB[i][j];
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
            
            double p100_s000_p100_s000_0_con = 0.0;
            double p100_s000_p010_s000_0_con = 0.0;
            double p100_s000_p001_s000_0_con = 0.0;
            double p010_s000_p100_s000_0_con = 0.0;
            double p010_s000_p010_s000_0_con = 0.0;
            double p010_s000_p001_s000_0_con = 0.0;
            double p001_s000_p100_s000_0_con = 0.0;
            double p001_s000_p010_s000_0_con = 0.0;
            double p001_s000_p001_s000_0_con = 0.0;


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

                    double PA[3] = {
                        AB_PA[0][(int)iab],
                        AB_PA[1][(int)iab],
                        AB_PA[2][(int)iab],
                    };

                    double PB[3] = {
                        AB_PA[0][(int)iab] + AB_AB[0][(int)iab],
                        AB_PA[1][(int)iab] + AB_AB[1][(int)iab],
                        AB_PA[2][(int)iab] + AB_AB[2][(int)iab],
                    };

                    double QC[3] = {
                        CD_PA[0][(int)icd],
                        CD_PA[1][(int)icd],
                        CD_PA[2][(int)icd],
                    };

                    double QD[3] = {
                        CD_PA[0][(int)icd] + CD_AB[0][(int)icd],
                        CD_PA[1][(int)icd] + CD_AB[1][(int)icd],
                        CD_PA[2][(int)icd] + CD_AB[2][(int)icd],
                    };

                    double W[3] = {
                        (zab * P[0] + zcd * Q[0]) * zabcd_inv,
                        (zab * P[1] + zcd * Q[1]) * zabcd_inv,
                        (zab * P[2] + zcd * Q[2]) * zabcd_inv,
                    };

                    double WP[3] = {
                        W[0] - PA[0],
                        W[1] - PA[1],
                        W[2] - PA[2],
                    };

                    double WQ[3] = {
                        W[0] - Q[0],
                        W[1] - Q[1],
                        W[2] - Q[2],
                    };

                    double T = rho * ((P[0] - Q[0]) * (P[0] - Q[0]) + (P[1] - Q[1]) * (P[1] - Q[1]) + (P[2] - Q[2]) * (P[2] - Q[2]));
                    double K = zabcd_inv_sqrt * (AB_K[(int)iab]) * CD_K[(int)icd];

                    double fm[3];
                    // vgamma<double>(2, T, fm);

                    for (int i = 0; i < 3; i++) {
                        fm[i] *= K;
                    }

                    double s000_s000_s000_s000_0 = fm[0] ;
                    double s000_s000_s000_s000_1 = fm[1] ;
                    double s000_s000_s000_s000_2 = fm[2] ;
                    double s000_s000_p100_s000_0 = QC[0] * s000_s000_s000_s000_0 + WQ[0] * s000_s000_s000_s000_1;
                    double s000_s000_p100_s000_1 = QC[0] * s000_s000_s000_s000_1 + WQ[0] * s000_s000_s000_s000_2;
                    double s000_s000_p010_s000_0 = QC[1] * s000_s000_s000_s000_0 + WQ[1] * s000_s000_s000_s000_1;
                    double s000_s000_p010_s000_1 = QC[1] * s000_s000_s000_s000_1 + WQ[1] * s000_s000_s000_s000_2;
                    double s000_s000_p001_s000_0 = QC[2] * s000_s000_s000_s000_0 + WQ[2] * s000_s000_s000_s000_1;
                    double s000_s000_p001_s000_1 = QC[2] * s000_s000_s000_s000_1 + WQ[2] * s000_s000_s000_s000_2;
                    double p100_s000_p100_s000_0 = PA[0] * s000_s000_p100_s000_0 + WP[0] * s000_s000_p100_s000_1 + 0.5 * zabcd_inv * 1 * s000_s000_s000_s000_1;
                    double p100_s000_p010_s000_0 = PA[0] * s000_s000_p010_s000_0 + WP[0] * s000_s000_p010_s000_1;
                    double p100_s000_p001_s000_0 = PA[0] * s000_s000_p001_s000_0 + WP[0] * s000_s000_p001_s000_1;
                    double p010_s000_p100_s000_0 = PA[1] * s000_s000_p100_s000_0 + WP[1] * s000_s000_p100_s000_1;
                    double p010_s000_p010_s000_0 = PA[1] * s000_s000_p010_s000_0 + WP[1] * s000_s000_p010_s000_1 + 0.5 * zabcd_inv * 1 * s000_s000_s000_s000_1;
                    double p010_s000_p001_s000_0 = PA[1] * s000_s000_p001_s000_0 + WP[1] * s000_s000_p001_s000_1;
                    double p001_s000_p100_s000_0 = PA[2] * s000_s000_p100_s000_0 + WP[2] * s000_s000_p100_s000_1;
                    double p001_s000_p010_s000_0 = PA[2] * s000_s000_p010_s000_0 + WP[2] * s000_s000_p010_s000_1;
                    double p001_s000_p001_s000_0 = PA[2] * s000_s000_p001_s000_0 + WP[2] * s000_s000_p001_s000_1 + 0.5 * zabcd_inv * 1 * s000_s000_s000_s000_1;
                    p100_s000_p100_s000_0_con += p100_s000_p100_s000_0;
                    p100_s000_p010_s000_0_con += p100_s000_p010_s000_0;
                    p100_s000_p001_s000_0_con += p100_s000_p001_s000_0;
                    p010_s000_p100_s000_0_con += p010_s000_p100_s000_0;
                    p010_s000_p010_s000_0_con += p010_s000_p010_s000_0;
                    p010_s000_p001_s000_0_con += p010_s000_p001_s000_0;
                    p001_s000_p100_s000_0_con += p001_s000_p100_s000_0;
                    p001_s000_p010_s000_0_con += p001_s000_p010_s000_0;
                    p001_s000_p001_s000_0_con += p001_s000_p001_s000_0;
                    
                    
                }
            }
            I_[idx + (0 + 0 * 3 + 0 * 3 * 1 + 0 * 3 * 1 * 3) * nab + idy * nab * 9] =                 1 * p100_s000_p100_s000_0_con ;
            I_[idx + (1 + 0 * 3 + 0 * 3 * 1 + 0 * 3 * 1 * 3) * nab + idy * nab * 9] =                 1 * p010_s000_p100_s000_0_con ;
            I_[idx + (2 + 0 * 3 + 0 * 3 * 1 + 0 * 3 * 1 * 3) * nab + idy * nab * 9] =                 1 * p001_s000_p100_s000_0_con ;
            I_[idx + (0 + 0 * 3 + 1 * 3 * 1 + 0 * 3 * 1 * 3) * nab + idy * nab * 9] =                 1 * p100_s000_p010_s000_0_con ;
            I_[idx + (1 + 0 * 3 + 1 * 3 * 1 + 0 * 3 * 1 * 3) * nab + idy * nab * 9] =                 1 * p010_s000_p010_s000_0_con ;
            I_[idx + (2 + 0 * 3 + 1 * 3 * 1 + 0 * 3 * 1 * 3) * nab + idy * nab * 9] =                 1 * p001_s000_p010_s000_0_con ;
            I_[idx + (0 + 0 * 3 + 2 * 3 * 1 + 0 * 3 * 1 * 3) * nab + idy * nab * 9] =                 1 * p100_s000_p001_s000_0_con ;
            I_[idx + (1 + 0 * 3 + 2 * 3 * 1 + 0 * 3 * 1 * 3) * nab + idy * nab * 9] =                 1 * p010_s000_p001_s000_0_con ;
            I_[idx + (2 + 0 * 3 + 2 * 3 * 1 + 0 * 3 * 1 * 3) * nab + idy * nab * 9] =                 1 * p001_s000_p001_s000_0_con ;
        }
    }   
}


void hrr_ppss(  int la, int lb, int na, int nb, int ma, int mb, std::array<double, 6> Z, std::array<double, 6> ZA, std::array<double, 6> K, std::array<double, 6> S, std::array<std::array<int,2>,2> idx, std::array<std::array<double,6>,3> P, std::array<std::array<double,6>,3> PA, std::array<std::array<double,6>,3> AB,
                double* I_ ) {

#pragma HLS INTERFACE m_axi port=I_ offset=slave bundle=gmem1 max_read_burst_length=256 max_write_burst_length=256 depth= 1
// ab
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
    CD_la = la;
    CD_lb = lb;
    CD_na = na;
    CD_nb = nb;
    CD_ma = ma;
    CD_mb = mb;

    for (int i = 0; i < 6; i++) {
        AB_Z[i] = Z[i];
        AB_ZA[i] = ZA[i];
        AB_K[i] = K[i];
        AB_S[i] = S[i];

        CD_Z[i] = Z[i];
        CD_ZA[i] = ZA[i];
        CD_K[i] = K[i];
        CD_S[i] = S[i];
    }

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            AB_idx[i][j] = idx[i][j];
            CD_idx[i][j] = idx[i][j];
        }
    }

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 6; j++) {
            AB_P[i][j] = P[i][j];
            AB_PA[i][j] = PA[i][j];
            AB_AB[i][j] = AB[i][j];

            CD_P[i][j] = P[i][j];
            CD_PA[i][j] = PA[i][j];
            CD_AB[i][j] = AB[i][j];
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
            
            double d200_s000_s000_s000_0_con = 0.0;
            double d110_s000_s000_s000_0_con = 0.0;
            double d101_s000_s000_s000_0_con = 0.0;
            double d020_s000_s000_s000_0_con = 0.0;
            double d011_s000_s000_s000_0_con = 0.0;
            double d002_s000_s000_s000_0_con = 0.0;
            double p100_s000_s000_s000_0_con = 0.0;
            double p010_s000_s000_s000_0_con = 0.0;
            double p001_s000_s000_s000_0_con = 0.0;



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

                    double PA[3] = {
                        AB_PA[0][(int)iab],
                        AB_PA[1][(int)iab],
                        AB_PA[2][(int)iab],
                    };

                    double PB[3] = {
                        AB_PA[0][(int)iab] + AB_AB[0][(int)iab],
                        AB_PA[1][(int)iab] + AB_AB[1][(int)iab],
                        AB_PA[2][(int)iab] + AB_AB[2][(int)iab],
                    };

                    double QC[3] = {
                        CD_PA[0][(int)icd],
                        CD_PA[1][(int)icd],
                        CD_PA[2][(int)icd],
                    };

                    double QD[3] = {
                        CD_PA[0][(int)icd] + CD_AB[0][(int)icd],
                        CD_PA[1][(int)icd] + CD_AB[1][(int)icd],
                        CD_PA[2][(int)icd] + CD_AB[2][(int)icd],
                    };

                    double W[3] = {
                        (zab * P[0] + zcd * Q[0]) * zabcd_inv,
                        (zab * P[1] + zcd * Q[1]) * zabcd_inv,
                        (zab * P[2] + zcd * Q[2]) * zabcd_inv,
                    };

                    double WP[3] = {
                        W[0] - PA[0],
                        W[1] - PA[1],
                        W[2] - PA[2],
                    };

                    double WQ[3] = {
                        W[0] - Q[0],
                        W[1] - Q[1],
                        W[2] - Q[2],
                    };

                    double T = rho * ((P[0] - Q[0]) * (P[0] - Q[0]) + (P[1] - Q[1]) * (P[1] - Q[1]) + (P[2] - Q[2]) * (P[2] - Q[2]));
                    double K = zabcd_inv_sqrt * (AB_K[(int)iab]) * CD_K[(int)icd];

                    double fm[3];
                    // vgamma<double>(2, T, fm);

                    for (int i = 0; i < 3; i++) {
                        fm[i] *= K;
                    }

                    double s000_s000_s000_s000_0 = fm[0] ;
                    double s000_s000_s000_s000_1 = fm[1] ;
                    double s000_s000_s000_s000_2 = fm[2] ;
                    double p100_s000_s000_s000_0 = PA[0] * s000_s000_s000_s000_0 + WP[0] * s000_s000_s000_s000_1;
                    double p100_s000_s000_s000_1 = PA[0] * s000_s000_s000_s000_1 + WP[0] * s000_s000_s000_s000_2;
                    double p010_s000_s000_s000_0 = PA[1] * s000_s000_s000_s000_0 + WP[1] * s000_s000_s000_s000_1;
                    double p010_s000_s000_s000_1 = PA[1] * s000_s000_s000_s000_1 + WP[1] * s000_s000_s000_s000_2;
                    double p001_s000_s000_s000_0 = PA[2] * s000_s000_s000_s000_0 + WP[2] * s000_s000_s000_s000_1;
                    double p001_s000_s000_s000_1 = PA[2] * s000_s000_s000_s000_1 + WP[2] * s000_s000_s000_s000_2;
                    double d200_s000_s000_s000_0 = PA[0] * p100_s000_s000_s000_0 + WP[0] * p100_s000_s000_s000_1 + 0.5 * zab_inv * 1 * s000_s000_s000_s000_0 + (- rho) * zab_inv * 0.5 * zab_inv * 1 * s000_s000_s000_s000_1;
                    double d110_s000_s000_s000_0 = PA[0] * p010_s000_s000_s000_0 + WP[0] * p010_s000_s000_s000_1;
                    double d101_s000_s000_s000_0 = PA[0] * p001_s000_s000_s000_0 + WP[0] * p001_s000_s000_s000_1;
                    double d020_s000_s000_s000_0 = PA[1] * p010_s000_s000_s000_0 + WP[1] * p010_s000_s000_s000_1 + 0.5 * zab_inv * 1 * s000_s000_s000_s000_0 + (- rho) * zab_inv * 0.5 * zab_inv * 1 * s000_s000_s000_s000_1;
                    double d011_s000_s000_s000_0 = PA[1] * p001_s000_s000_s000_0 + WP[1] * p001_s000_s000_s000_1;
                    double d002_s000_s000_s000_0 = PA[2] * p001_s000_s000_s000_0 + WP[2] * p001_s000_s000_s000_1 + 0.5 * zab_inv * 1 * s000_s000_s000_s000_0 + (- rho) * zab_inv * 0.5 * zab_inv * 1 * s000_s000_s000_s000_1;
                    d200_s000_s000_s000_0_con += d200_s000_s000_s000_0;
                    d110_s000_s000_s000_0_con += d110_s000_s000_s000_0;
                    d101_s000_s000_s000_0_con += d101_s000_s000_s000_0;
                    d020_s000_s000_s000_0_con += d020_s000_s000_s000_0;
                    d011_s000_s000_s000_0_con += d011_s000_s000_s000_0;
                    d002_s000_s000_s000_0_con += d002_s000_s000_s000_0;
                    p100_s000_s000_s000_0_con += p100_s000_s000_s000_0;
                    p010_s000_s000_s000_0_con += p010_s000_s000_s000_0;
                    p001_s000_s000_s000_0_con += p001_s000_s000_s000_0;
                    
                }
            }

            double p100_p100_s000_s000_0_con = d200_s000_s000_s000_0_con + AB[0] * p100_s000_s000_s000_0_con ; 
            double p100_p010_s000_s000_0_con = d110_s000_s000_s000_0_con + AB[1] * p100_s000_s000_s000_0_con ; 
            double p100_p001_s000_s000_0_con = d101_s000_s000_s000_0_con + AB[2] * p100_s000_s000_s000_0_con ; 
            double p010_p100_s000_s000_0_con = d110_s000_s000_s000_0_con + AB[0] * p010_s000_s000_s000_0_con ; 
            double p010_p010_s000_s000_0_con = d020_s000_s000_s000_0_con + AB[1] * p010_s000_s000_s000_0_con ; 
            double p010_p001_s000_s000_0_con = d011_s000_s000_s000_0_con + AB[2] * p010_s000_s000_s000_0_con ; 
            double p001_p100_s000_s000_0_con = d101_s000_s000_s000_0_con + AB[0] * p001_s000_s000_s000_0_con ; 
            double p001_p010_s000_s000_0_con = d011_s000_s000_s000_0_con + AB[1] * p001_s000_s000_s000_0_con ; 
            double p001_p001_s000_s000_0_con = d002_s000_s000_s000_0_con + AB[2] * p001_s000_s000_s000_0_con ; 
            I_[idx + (0 + 0 * 3 + 0 * 3 * 3 + 0 * 3 * 3 * 1) * nab + idy * nab * 9] =                 1 * p100_p100_s000_s000_0_con ;
            I_[idx + (1 + 0 * 3 + 0 * 3 * 3 + 0 * 3 * 3 * 1) * nab + idy * nab * 9] =                 1 * p010_p100_s000_s000_0_con ;
            I_[idx + (2 + 0 * 3 + 0 * 3 * 3 + 0 * 3 * 3 * 1) * nab + idy * nab * 9] =                 1 * p001_p100_s000_s000_0_con ;
            I_[idx + (0 + 1 * 3 + 0 * 3 * 3 + 0 * 3 * 3 * 1) * nab + idy * nab * 9] =                 1 * p100_p010_s000_s000_0_con ;
            I_[idx + (1 + 1 * 3 + 0 * 3 * 3 + 0 * 3 * 3 * 1) * nab + idy * nab * 9] =                 1 * p010_p010_s000_s000_0_con ;
            I_[idx + (2 + 1 * 3 + 0 * 3 * 3 + 0 * 3 * 3 * 1) * nab + idy * nab * 9] =                 1 * p001_p010_s000_s000_0_con ;
            I_[idx + (0 + 2 * 3 + 0 * 3 * 3 + 0 * 3 * 3 * 1) * nab + idy * nab * 9] =                 1 * p100_p001_s000_s000_0_con ;
            I_[idx + (1 + 2 * 3 + 0 * 3 * 3 + 0 * 3 * 3 * 1) * nab + idy * nab * 9] =                 1 * p010_p001_s000_s000_0_con ;
            I_[idx + (2 + 2 * 3 + 0 * 3 * 3 + 0 * 3 * 3 * 1) * nab + idy * nab * 9] =                 1 * p001_p001_s000_s000_0_con ;
        }
    }   
}

void hrr_ppps(  int la, int lb, int na, int nb, int ma, int mb, std::array<double, 6> Z, std::array<double, 6> ZA, std::array<double, 6> K, std::array<double, 6> S, std::array<std::array<int,2>,2> idx, std::array<std::array<double,6>,3> P, std::array<std::array<double,6>,3> PA, std::array<std::array<double,6>,3> AB,
                double* I_ ) {

#pragma HLS INTERFACE m_axi port=I_ offset=slave bundle=gmem1 max_read_burst_length=256 max_write_burst_length=256 depth= 1
// ab
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
    CD_la = la;
    CD_lb = lb;
    CD_na = na;
    CD_nb = nb;
    CD_ma = ma;
    CD_mb = mb;

    for (int i = 0; i < 6; i++) {
        AB_Z[i] = Z[i];
        AB_ZA[i] = ZA[i];
        AB_K[i] = K[i];
        AB_S[i] = S[i];

        CD_Z[i] = Z[i];
        CD_ZA[i] = ZA[i];
        CD_K[i] = K[i];
        CD_S[i] = S[i];
    }

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            AB_idx[i][j] = idx[i][j];
            CD_idx[i][j] = idx[i][j];
        }
    }

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 6; j++) {
            AB_P[i][j] = P[i][j];
            AB_PA[i][j] = PA[i][j];
            AB_AB[i][j] = AB[i][j];

            CD_P[i][j] = P[i][j];
            CD_PA[i][j] = PA[i][j];
            CD_AB[i][j] = AB[i][j];
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
            
            double d200_s000_p100_s000_0_con = 0.0;
            double d200_s000_p010_s000_0_con = 0.0;
            double d200_s000_p001_s000_0_con = 0.0;
            double d110_s000_p100_s000_0_con = 0.0;
            double d110_s000_p010_s000_0_con = 0.0;
            double d110_s000_p001_s000_0_con = 0.0;
            double d101_s000_p100_s000_0_con = 0.0;
            double d101_s000_p010_s000_0_con = 0.0;
            double d101_s000_p001_s000_0_con = 0.0;
            double d020_s000_p100_s000_0_con = 0.0;
            double d020_s000_p010_s000_0_con = 0.0;
            double d020_s000_p001_s000_0_con = 0.0;
            double d011_s000_p100_s000_0_con = 0.0;
            double d011_s000_p010_s000_0_con = 0.0;
            double d011_s000_p001_s000_0_con = 0.0;
            double d002_s000_p100_s000_0_con = 0.0;
            double d002_s000_p010_s000_0_con = 0.0;
            double d002_s000_p001_s000_0_con = 0.0;
            double p100_s000_p100_s000_0_con = 0.0;
            double p100_s000_p010_s000_0_con = 0.0;
            double p100_s000_p001_s000_0_con = 0.0;
            double p010_s000_p100_s000_0_con = 0.0;
            double p010_s000_p010_s000_0_con = 0.0;
            double p010_s000_p001_s000_0_con = 0.0;
            double p001_s000_p100_s000_0_con = 0.0;
            double p001_s000_p010_s000_0_con = 0.0;
            double p001_s000_p001_s000_0_con = 0.0;

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

                    double PA[3] = {
                        AB_PA[0][(int)iab],
                        AB_PA[1][(int)iab],
                        AB_PA[2][(int)iab],
                    };

                    double PB[3] = {
                        AB_PA[0][(int)iab] + AB_AB[0][(int)iab],
                        AB_PA[1][(int)iab] + AB_AB[1][(int)iab],
                        AB_PA[2][(int)iab] + AB_AB[2][(int)iab],
                    };

                    double QC[3] = {
                        CD_PA[0][(int)icd],
                        CD_PA[1][(int)icd],
                        CD_PA[2][(int)icd],
                    };

                    double QD[3] = {
                        CD_PA[0][(int)icd] + CD_AB[0][(int)icd],
                        CD_PA[1][(int)icd] + CD_AB[1][(int)icd],
                        CD_PA[2][(int)icd] + CD_AB[2][(int)icd],
                    };

                    double W[3] = {
                        (zab * P[0] + zcd * Q[0]) * zabcd_inv,
                        (zab * P[1] + zcd * Q[1]) * zabcd_inv,
                        (zab * P[2] + zcd * Q[2]) * zabcd_inv,
                    };

                    double WP[3] = {
                        W[0] - PA[0],
                        W[1] - PA[1],
                        W[2] - PA[2],
                    };

                    double WQ[3] = {
                        W[0] - Q[0],
                        W[1] - Q[1],
                        W[2] - Q[2],
                    };

                    double T = rho * ((P[0] - Q[0]) * (P[0] - Q[0]) + (P[1] - Q[1]) * (P[1] - Q[1]) + (P[2] - Q[2]) * (P[2] - Q[2]));
                    double K = zabcd_inv_sqrt * (AB_K[(int)iab]) * CD_K[(int)icd];

                    double fm[4];
                    // vgamma<double>(3, T, fm);

                    for (int i = 0; i < 4; i++) {
                        fm[i] *= K;
                    }

                    double s000_s000_s000_s000_0 = fm[0] ;
                    double s000_s000_s000_s000_1 = fm[1] ;
                    double s000_s000_s000_s000_2 = fm[2] ;
                    double s000_s000_s000_s000_3 = fm[3] ;
                    double s000_s000_p100_s000_0 = QC[0] * s000_s000_s000_s000_0 + WQ[0] * s000_s000_s000_s000_1;
                    double s000_s000_p100_s000_1 = QC[0] * s000_s000_s000_s000_1 + WQ[0] * s000_s000_s000_s000_2;
                    double s000_s000_p010_s000_0 = QC[1] * s000_s000_s000_s000_0 + WQ[1] * s000_s000_s000_s000_1;
                    double s000_s000_p010_s000_1 = QC[1] * s000_s000_s000_s000_1 + WQ[1] * s000_s000_s000_s000_2;
                    double s000_s000_p001_s000_0 = QC[2] * s000_s000_s000_s000_0 + WQ[2] * s000_s000_s000_s000_1;
                    double s000_s000_p001_s000_1 = QC[2] * s000_s000_s000_s000_1 + WQ[2] * s000_s000_s000_s000_2;
                    double p100_s000_s000_s000_0 = PA[0] * s000_s000_s000_s000_0 + WP[0] * s000_s000_s000_s000_1;
                    double p100_s000_s000_s000_1 = PA[0] * s000_s000_s000_s000_1 + WP[0] * s000_s000_s000_s000_2;
                    double p100_s000_s000_s000_2 = PA[0] * s000_s000_s000_s000_2 + WP[0] * s000_s000_s000_s000_3;
                    double p010_s000_s000_s000_0 = PA[1] * s000_s000_s000_s000_0 + WP[1] * s000_s000_s000_s000_1;
                    double p010_s000_s000_s000_1 = PA[1] * s000_s000_s000_s000_1 + WP[1] * s000_s000_s000_s000_2;
                    double p010_s000_s000_s000_2 = PA[1] * s000_s000_s000_s000_2 + WP[1] * s000_s000_s000_s000_3;
                    double p001_s000_s000_s000_0 = PA[2] * s000_s000_s000_s000_0 + WP[2] * s000_s000_s000_s000_1;
                    double p001_s000_s000_s000_1 = PA[2] * s000_s000_s000_s000_1 + WP[2] * s000_s000_s000_s000_2;
                    double p001_s000_s000_s000_2 = PA[2] * s000_s000_s000_s000_2 + WP[2] * s000_s000_s000_s000_3;
                    double p100_s000_p100_s000_0 = PA[0] * s000_s000_p100_s000_0 + WP[0] * s000_s000_p100_s000_1 + 0.5 * zabcd_inv * 1 * s000_s000_s000_s000_1;
                    double p100_s000_p010_s000_0 = PA[0] * s000_s000_p010_s000_0 + WP[0] * s000_s000_p010_s000_1;
                    double p100_s000_p001_s000_0 = PA[0] * s000_s000_p001_s000_0 + WP[0] * s000_s000_p001_s000_1;
                    double p010_s000_p100_s000_0 = PA[1] * s000_s000_p100_s000_0 + WP[1] * s000_s000_p100_s000_1;
                    double p010_s000_p010_s000_0 = PA[1] * s000_s000_p010_s000_0 + WP[1] * s000_s000_p010_s000_1 + 0.5 * zabcd_inv * 1 * s000_s000_s000_s000_1;
                    double p010_s000_p001_s000_0 = PA[1] * s000_s000_p001_s000_0 + WP[1] * s000_s000_p001_s000_1;
                    double p001_s000_p100_s000_0 = PA[2] * s000_s000_p100_s000_0 + WP[2] * s000_s000_p100_s000_1;
                    double p001_s000_p010_s000_0 = PA[2] * s000_s000_p010_s000_0 + WP[2] * s000_s000_p010_s000_1;
                    double p001_s000_p001_s000_0 = PA[2] * s000_s000_p001_s000_0 + WP[2] * s000_s000_p001_s000_1 + 0.5 * zabcd_inv * 1 * s000_s000_s000_s000_1;
                    double d200_s000_s000_s000_0 = PA[0] * p100_s000_s000_s000_0 + WP[0] * p100_s000_s000_s000_1 + 0.5 * zab_inv * 1 * s000_s000_s000_s000_0 + (- rho) * zab_inv * 0.5 * zab_inv * 1 * s000_s000_s000_s000_1;
                    double d200_s000_s000_s000_1 = PA[0] * p100_s000_s000_s000_1 + WP[0] * p100_s000_s000_s000_2 + 0.5 * zab_inv * 1 * s000_s000_s000_s000_1 + (- rho) * zab_inv * 0.5 * zab_inv * 1 * s000_s000_s000_s000_2;
                    double d110_s000_s000_s000_0 = PA[0] * p010_s000_s000_s000_0 + WP[0] * p010_s000_s000_s000_1;
                    double d110_s000_s000_s000_1 = PA[0] * p010_s000_s000_s000_1 + WP[0] * p010_s000_s000_s000_2;
                    double d101_s000_s000_s000_0 = PA[0] * p001_s000_s000_s000_0 + WP[0] * p001_s000_s000_s000_1;
                    double d101_s000_s000_s000_1 = PA[0] * p001_s000_s000_s000_1 + WP[0] * p001_s000_s000_s000_2;
                    double d020_s000_s000_s000_0 = PA[1] * p010_s000_s000_s000_0 + WP[1] * p010_s000_s000_s000_1 + 0.5 * zab_inv * 1 * s000_s000_s000_s000_0 + (- rho) * zab_inv * 0.5 * zab_inv * 1 * s000_s000_s000_s000_1;
                    double d020_s000_s000_s000_1 = PA[1] * p010_s000_s000_s000_1 + WP[1] * p010_s000_s000_s000_2 + 0.5 * zab_inv * 1 * s000_s000_s000_s000_1 + (- rho) * zab_inv * 0.5 * zab_inv * 1 * s000_s000_s000_s000_2;
                    double d011_s000_s000_s000_0 = PA[1] * p001_s000_s000_s000_0 + WP[1] * p001_s000_s000_s000_1;
                    double d011_s000_s000_s000_1 = PA[1] * p001_s000_s000_s000_1 + WP[1] * p001_s000_s000_s000_2;
                    double d002_s000_s000_s000_0 = PA[2] * p001_s000_s000_s000_0 + WP[2] * p001_s000_s000_s000_1 + 0.5 * zab_inv * 1 * s000_s000_s000_s000_0 + (- rho) * zab_inv * 0.5 * zab_inv * 1 * s000_s000_s000_s000_1;
                    double d002_s000_s000_s000_1 = PA[2] * p001_s000_s000_s000_1 + WP[2] * p001_s000_s000_s000_2 + 0.5 * zab_inv * 1 * s000_s000_s000_s000_1 + (- rho) * zab_inv * 0.5 * zab_inv * 1 * s000_s000_s000_s000_2;
                    double d200_s000_p100_s000_0 = QC[0] * d200_s000_s000_s000_0 + WQ[0] * d200_s000_s000_s000_1 +  0.0 * p100_s000_s000_s000_0 + 0.5 * zabcd_inv * 2 * p100_s000_s000_s000_1;
                    double d200_s000_p010_s000_0 = QC[1] * d200_s000_s000_s000_0 + WQ[1] * d200_s000_s000_s000_1;
                    double d200_s000_p001_s000_0 = QC[2] * d200_s000_s000_s000_0 + WQ[2] * d200_s000_s000_s000_1;
                    double d110_s000_p100_s000_0 = QC[0] * d110_s000_s000_s000_0 + WQ[0] * d110_s000_s000_s000_1 +  0.0 * p010_s000_s000_s000_0 + 0.5 * zabcd_inv * 1 * p010_s000_s000_s000_1;
                    double d110_s000_p010_s000_0 = QC[1] * d110_s000_s000_s000_0 + WQ[1] * d110_s000_s000_s000_1 +  0.0 * p100_s000_s000_s000_0 + 0.5 * zabcd_inv * 1 * p100_s000_s000_s000_1;
                    double d110_s000_p001_s000_0 = QC[2] * d110_s000_s000_s000_0 + WQ[2] * d110_s000_s000_s000_1;
                    double d101_s000_p100_s000_0 = QC[0] * d101_s000_s000_s000_0 + WQ[0] * d101_s000_s000_s000_1 +  0.0 * p001_s000_s000_s000_0 + 0.5 * zabcd_inv * 1 * p001_s000_s000_s000_1;
                    double d101_s000_p010_s000_0 = QC[1] * d101_s000_s000_s000_0 + WQ[1] * d101_s000_s000_s000_1;
                    double d101_s000_p001_s000_0 = QC[2] * d101_s000_s000_s000_0 + WQ[2] * d101_s000_s000_s000_1 +  0.0 * p100_s000_s000_s000_0 + 0.5 * zabcd_inv * 1 * p100_s000_s000_s000_1;
                    double d020_s000_p100_s000_0 = QC[0] * d020_s000_s000_s000_0 + WQ[0] * d020_s000_s000_s000_1;
                    double d020_s000_p010_s000_0 = QC[1] * d020_s000_s000_s000_0 + WQ[1] * d020_s000_s000_s000_1 +  0.0 * p010_s000_s000_s000_0 + 0.5 * zabcd_inv * 2 * p010_s000_s000_s000_1;
                    double d020_s000_p001_s000_0 = QC[2] * d020_s000_s000_s000_0 + WQ[2] * d020_s000_s000_s000_1;
                    double d011_s000_p100_s000_0 = QC[0] * d011_s000_s000_s000_0 + WQ[0] * d011_s000_s000_s000_1;
                    double d011_s000_p010_s000_0 = QC[1] * d011_s000_s000_s000_0 + WQ[1] * d011_s000_s000_s000_1 +  0.0 * p001_s000_s000_s000_0 + 0.5 * zabcd_inv * 1 * p001_s000_s000_s000_1;
                    double d011_s000_p001_s000_0 = QC[2] * d011_s000_s000_s000_0 + WQ[2] * d011_s000_s000_s000_1 +  0.0 * p010_s000_s000_s000_0 + 0.5 * zabcd_inv * 1 * p010_s000_s000_s000_1;
                    double d002_s000_p100_s000_0 = QC[0] * d002_s000_s000_s000_0 + WQ[0] * d002_s000_s000_s000_1;
                    double d002_s000_p010_s000_0 = QC[1] * d002_s000_s000_s000_0 + WQ[1] * d002_s000_s000_s000_1;
                    double d002_s000_p001_s000_0 = QC[2] * d002_s000_s000_s000_0 + WQ[2] * d002_s000_s000_s000_1 +  0.0 * p001_s000_s000_s000_0 + 0.5 * zabcd_inv * 2 * p001_s000_s000_s000_1;
                    d200_s000_p100_s000_0_con += d200_s000_p100_s000_0;
                    d200_s000_p010_s000_0_con += d200_s000_p010_s000_0;
                    d200_s000_p001_s000_0_con += d200_s000_p001_s000_0;
                    d110_s000_p100_s000_0_con += d110_s000_p100_s000_0;
                    d110_s000_p010_s000_0_con += d110_s000_p010_s000_0;
                    d110_s000_p001_s000_0_con += d110_s000_p001_s000_0;
                    d101_s000_p100_s000_0_con += d101_s000_p100_s000_0;
                    d101_s000_p010_s000_0_con += d101_s000_p010_s000_0;
                    d101_s000_p001_s000_0_con += d101_s000_p001_s000_0;
                    d020_s000_p100_s000_0_con += d020_s000_p100_s000_0;
                    d020_s000_p010_s000_0_con += d020_s000_p010_s000_0;
                    d020_s000_p001_s000_0_con += d020_s000_p001_s000_0;
                    d011_s000_p100_s000_0_con += d011_s000_p100_s000_0;
                    d011_s000_p010_s000_0_con += d011_s000_p010_s000_0;
                    d011_s000_p001_s000_0_con += d011_s000_p001_s000_0;
                    d002_s000_p100_s000_0_con += d002_s000_p100_s000_0;
                    d002_s000_p010_s000_0_con += d002_s000_p010_s000_0;
                    d002_s000_p001_s000_0_con += d002_s000_p001_s000_0;
                    p100_s000_p100_s000_0_con += p100_s000_p100_s000_0;
                    p100_s000_p010_s000_0_con += p100_s000_p010_s000_0;
                    p100_s000_p001_s000_0_con += p100_s000_p001_s000_0;
                    p010_s000_p100_s000_0_con += p010_s000_p100_s000_0;
                    p010_s000_p010_s000_0_con += p010_s000_p010_s000_0;
                    p010_s000_p001_s000_0_con += p010_s000_p001_s000_0;
                    p001_s000_p100_s000_0_con += p001_s000_p100_s000_0;
                    p001_s000_p010_s000_0_con += p001_s000_p010_s000_0;
                    p001_s000_p001_s000_0_con += p001_s000_p001_s000_0;  
                }
            }

            double p100_p100_p100_s000_0_con = d200_s000_p100_s000_0_con + AB[0] * p100_s000_p100_s000_0_con ; 
            double p100_p100_p010_s000_0_con = d200_s000_p010_s000_0_con + AB[0] * p100_s000_p010_s000_0_con ; 
            double p100_p100_p001_s000_0_con = d200_s000_p001_s000_0_con + AB[0] * p100_s000_p001_s000_0_con ; 
            double p100_p010_p100_s000_0_con = d110_s000_p100_s000_0_con + AB[1] * p100_s000_p100_s000_0_con ; 
            double p100_p010_p010_s000_0_con = d110_s000_p010_s000_0_con + AB[1] * p100_s000_p010_s000_0_con ; 
            double p100_p010_p001_s000_0_con = d110_s000_p001_s000_0_con + AB[1] * p100_s000_p001_s000_0_con ; 
            double p100_p001_p100_s000_0_con = d101_s000_p100_s000_0_con + AB[2] * p100_s000_p100_s000_0_con ; 
            double p100_p001_p010_s000_0_con = d101_s000_p010_s000_0_con + AB[2] * p100_s000_p010_s000_0_con ; 
            double p100_p001_p001_s000_0_con = d101_s000_p001_s000_0_con + AB[2] * p100_s000_p001_s000_0_con ; 
            double p010_p100_p100_s000_0_con = d110_s000_p100_s000_0_con + AB[0] * p010_s000_p100_s000_0_con ; 
            double p010_p100_p010_s000_0_con = d110_s000_p010_s000_0_con + AB[0] * p010_s000_p010_s000_0_con ; 
            double p010_p100_p001_s000_0_con = d110_s000_p001_s000_0_con + AB[0] * p010_s000_p001_s000_0_con ; 
            double p010_p010_p100_s000_0_con = d020_s000_p100_s000_0_con + AB[1] * p010_s000_p100_s000_0_con ; 
            double p010_p010_p010_s000_0_con = d020_s000_p010_s000_0_con + AB[1] * p010_s000_p010_s000_0_con ; 
            double p010_p010_p001_s000_0_con = d020_s000_p001_s000_0_con + AB[1] * p010_s000_p001_s000_0_con ; 
            double p010_p001_p100_s000_0_con = d011_s000_p100_s000_0_con + AB[2] * p010_s000_p100_s000_0_con ; 
            double p010_p001_p010_s000_0_con = d011_s000_p010_s000_0_con + AB[2] * p010_s000_p010_s000_0_con ; 
            double p010_p001_p001_s000_0_con = d011_s000_p001_s000_0_con + AB[2] * p010_s000_p001_s000_0_con ; 
            double p001_p100_p100_s000_0_con = d101_s000_p100_s000_0_con + AB[0] * p001_s000_p100_s000_0_con ; 
            double p001_p100_p010_s000_0_con = d101_s000_p010_s000_0_con + AB[0] * p001_s000_p010_s000_0_con ; 
            double p001_p100_p001_s000_0_con = d101_s000_p001_s000_0_con + AB[0] * p001_s000_p001_s000_0_con ; 
            double p001_p010_p100_s000_0_con = d011_s000_p100_s000_0_con + AB[1] * p001_s000_p100_s000_0_con ; 
            double p001_p010_p010_s000_0_con = d011_s000_p010_s000_0_con + AB[1] * p001_s000_p010_s000_0_con ; 
            double p001_p010_p001_s000_0_con = d011_s000_p001_s000_0_con + AB[1] * p001_s000_p001_s000_0_con ; 
            double p001_p001_p100_s000_0_con = d002_s000_p100_s000_0_con + AB[2] * p001_s000_p100_s000_0_con ; 
            double p001_p001_p010_s000_0_con = d002_s000_p010_s000_0_con + AB[2] * p001_s000_p010_s000_0_con ; 
            double p001_p001_p001_s000_0_con = d002_s000_p001_s000_0_con + AB[2] * p001_s000_p001_s000_0_con ; 
            I_[idx + (0 + 0 * 3 + 0 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 27] =                 1 * p100_p100_p100_s000_0_con ;
            I_[idx + (1 + 0 * 3 + 0 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 27] =                 1 * p010_p100_p100_s000_0_con ;
            I_[idx + (2 + 0 * 3 + 0 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 27] =                 1 * p001_p100_p100_s000_0_con ;
            I_[idx + (0 + 1 * 3 + 0 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 27] =                 1 * p100_p010_p100_s000_0_con ;
            I_[idx + (1 + 1 * 3 + 0 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 27] =                 1 * p010_p010_p100_s000_0_con ;
            I_[idx + (2 + 1 * 3 + 0 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 27] =                 1 * p001_p010_p100_s000_0_con ;
            I_[idx + (0 + 2 * 3 + 0 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 27] =                 1 * p100_p001_p100_s000_0_con ;
            I_[idx + (1 + 2 * 3 + 0 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 27] =                 1 * p010_p001_p100_s000_0_con ;
            I_[idx + (2 + 2 * 3 + 0 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 27] =                 1 * p001_p001_p100_s000_0_con ;
            I_[idx + (0 + 0 * 3 + 1 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 27] =                 1 * p100_p100_p010_s000_0_con ;
            I_[idx + (1 + 0 * 3 + 1 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 27] =                 1 * p010_p100_p010_s000_0_con ;
            I_[idx + (2 + 0 * 3 + 1 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 27] =                 1 * p001_p100_p010_s000_0_con ;
            I_[idx + (0 + 1 * 3 + 1 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 27] =                 1 * p100_p010_p010_s000_0_con ;
            I_[idx + (1 + 1 * 3 + 1 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 27] =                 1 * p010_p010_p010_s000_0_con ;
            I_[idx + (2 + 1 * 3 + 1 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 27] =                 1 * p001_p010_p010_s000_0_con ;
            I_[idx + (0 + 2 * 3 + 1 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 27] =                 1 * p100_p001_p010_s000_0_con ;
            I_[idx + (1 + 2 * 3 + 1 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 27] =                 1 * p010_p001_p010_s000_0_con ;
            I_[idx + (2 + 2 * 3 + 1 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 27] =                 1 * p001_p001_p010_s000_0_con ;
            I_[idx + (0 + 0 * 3 + 2 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 27] =                 1 * p100_p100_p001_s000_0_con ;
            I_[idx + (1 + 0 * 3 + 2 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 27] =                 1 * p010_p100_p001_s000_0_con ;
            I_[idx + (2 + 0 * 3 + 2 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 27] =                 1 * p001_p100_p001_s000_0_con ;
            I_[idx + (0 + 1 * 3 + 2 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 27] =                 1 * p100_p010_p001_s000_0_con ;
            I_[idx + (1 + 1 * 3 + 2 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 27] =                 1 * p010_p010_p001_s000_0_con ;
            I_[idx + (2 + 1 * 3 + 2 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 27] =                 1 * p001_p010_p001_s000_0_con ;
            I_[idx + (0 + 2 * 3 + 2 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 27] =                 1 * p100_p001_p001_s000_0_con ;
            I_[idx + (1 + 2 * 3 + 2 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 27] =                 1 * p010_p001_p001_s000_0_con ;
            I_[idx + (2 + 2 * 3 + 2 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 27] =                 1 * p001_p001_p001_s000_0_con ;
        }
    } 
}


void hrr_pppp(  int la, int lb, int na, int nb, int ma, int mb, std::array<double, 6> Z, std::array<double, 6> ZA, std::array<double, 6> K, std::array<double, 6> S, std::array<std::array<int,2>,2> idx, std::array<std::array<double,6>,3> P, std::array<std::array<double,6>,3> PA, std::array<std::array<double,6>,3> AB,
                double* I_ ) {

#pragma HLS INTERFACE m_axi port=I_ offset=slave bundle=gmem1 max_read_burst_length=256 max_write_burst_length=256 depth= 1
// ab
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
    CD_la = la;
    CD_lb = lb;
    CD_na = na;
    CD_nb = nb;
    CD_ma = ma;
    CD_mb = mb;

    for (int i = 0; i < 6; i++) {
        AB_Z[i] = Z[i];
        AB_ZA[i] = ZA[i];
        AB_K[i] = K[i];
        AB_S[i] = S[i];

        CD_Z[i] = Z[i];
        CD_ZA[i] = ZA[i];
        CD_K[i] = K[i];
        CD_S[i] = S[i];
    }

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            AB_idx[i][j] = idx[i][j];
            CD_idx[i][j] = idx[i][j];
        }
    }

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 6; j++) {
            AB_P[i][j] = P[i][j];
            AB_PA[i][j] = PA[i][j];
            AB_AB[i][j] = AB[i][j];

            CD_P[i][j] = P[i][j];
            CD_PA[i][j] = PA[i][j];
            CD_AB[i][j] = AB[i][j];
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
            
            double d200_s000_d200_s000_0_con = 0.0;
            double d200_s000_d110_s000_0_con = 0.0;
            double d200_s000_d101_s000_0_con = 0.0;
            double d200_s000_d020_s000_0_con = 0.0;
            double d200_s000_d011_s000_0_con = 0.0;
            double d200_s000_d002_s000_0_con = 0.0;
            double d110_s000_d200_s000_0_con = 0.0;
            double d110_s000_d110_s000_0_con = 0.0;
            double d110_s000_d101_s000_0_con = 0.0;
            double d110_s000_d020_s000_0_con = 0.0;
            double d110_s000_d011_s000_0_con = 0.0;
            double d110_s000_d002_s000_0_con = 0.0;
            double d101_s000_d200_s000_0_con = 0.0;
            double d101_s000_d110_s000_0_con = 0.0;
            double d101_s000_d101_s000_0_con = 0.0;
            double d101_s000_d020_s000_0_con = 0.0;
            double d101_s000_d011_s000_0_con = 0.0;
            double d101_s000_d002_s000_0_con = 0.0;
            double d020_s000_d200_s000_0_con = 0.0;
            double d020_s000_d110_s000_0_con = 0.0;
            double d020_s000_d101_s000_0_con = 0.0;
            double d020_s000_d020_s000_0_con = 0.0;
            double d020_s000_d011_s000_0_con = 0.0;
            double d020_s000_d002_s000_0_con = 0.0;
            double d011_s000_d200_s000_0_con = 0.0;
            double d011_s000_d110_s000_0_con = 0.0;
            double d011_s000_d101_s000_0_con = 0.0;
            double d011_s000_d020_s000_0_con = 0.0;
            double d011_s000_d011_s000_0_con = 0.0;
            double d011_s000_d002_s000_0_con = 0.0;
            double d002_s000_d200_s000_0_con = 0.0;
            double d002_s000_d110_s000_0_con = 0.0;
            double d002_s000_d101_s000_0_con = 0.0;
            double d002_s000_d020_s000_0_con = 0.0;
            double d002_s000_d011_s000_0_con = 0.0;
            double d002_s000_d002_s000_0_con = 0.0;
            double d200_s000_p100_s000_0_con = 0.0;
            double d200_s000_p010_s000_0_con = 0.0;
            double d200_s000_p001_s000_0_con = 0.0;
            double d110_s000_p100_s000_0_con = 0.0;
            double d110_s000_p010_s000_0_con = 0.0;
            double d110_s000_p001_s000_0_con = 0.0;
            double d101_s000_p100_s000_0_con = 0.0;
            double d101_s000_p010_s000_0_con = 0.0;
            double d101_s000_p001_s000_0_con = 0.0;
            double d020_s000_p100_s000_0_con = 0.0;
            double d020_s000_p010_s000_0_con = 0.0;
            double d020_s000_p001_s000_0_con = 0.0;
            double d011_s000_p100_s000_0_con = 0.0;
            double d011_s000_p010_s000_0_con = 0.0;
            double d011_s000_p001_s000_0_con = 0.0;
            double d002_s000_p100_s000_0_con = 0.0;
            double d002_s000_p010_s000_0_con = 0.0;
            double d002_s000_p001_s000_0_con = 0.0;
            double p100_s000_d200_s000_0_con = 0.0;
            double p100_s000_d110_s000_0_con = 0.0;
            double p100_s000_d101_s000_0_con = 0.0;
            double p100_s000_d020_s000_0_con = 0.0;
            double p100_s000_d011_s000_0_con = 0.0;
            double p100_s000_d002_s000_0_con = 0.0;
            double p010_s000_d200_s000_0_con = 0.0;
            double p010_s000_d110_s000_0_con = 0.0;
            double p010_s000_d101_s000_0_con = 0.0;
            double p010_s000_d020_s000_0_con = 0.0;
            double p010_s000_d011_s000_0_con = 0.0;
            double p010_s000_d002_s000_0_con = 0.0;
            double p001_s000_d200_s000_0_con = 0.0;
            double p001_s000_d110_s000_0_con = 0.0;
            double p001_s000_d101_s000_0_con = 0.0;
            double p001_s000_d020_s000_0_con = 0.0;
            double p001_s000_d011_s000_0_con = 0.0;
            double p001_s000_d002_s000_0_con = 0.0;
            double p100_s000_p100_s000_0_con = 0.0;
            double p100_s000_p010_s000_0_con = 0.0;
            double p100_s000_p001_s000_0_con = 0.0;
            double p010_s000_p100_s000_0_con = 0.0;
            double p010_s000_p010_s000_0_con = 0.0;
            double p010_s000_p001_s000_0_con = 0.0;
            double p001_s000_p100_s000_0_con = 0.0;
            double p001_s000_p010_s000_0_con = 0.0;
            double p001_s000_p001_s000_0_con = 0.0;

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

                    double PA[3] = {
                        AB_PA[0][(int)iab],
                        AB_PA[1][(int)iab],
                        AB_PA[2][(int)iab],
                    };

                    double PB[3] = {
                        AB_PA[0][(int)iab] + AB_AB[0][(int)iab],
                        AB_PA[1][(int)iab] + AB_AB[1][(int)iab],
                        AB_PA[2][(int)iab] + AB_AB[2][(int)iab],
                    };

                    double QC[3] = {
                        CD_PA[0][(int)icd],
                        CD_PA[1][(int)icd],
                        CD_PA[2][(int)icd],
                    };

                    double QD[3] = {
                        CD_PA[0][(int)icd] + CD_AB[0][(int)icd],
                        CD_PA[1][(int)icd] + CD_AB[1][(int)icd],
                        CD_PA[2][(int)icd] + CD_AB[2][(int)icd],
                    };

                    double W[3] = {
                        (zab * P[0] + zcd * Q[0]) * zabcd_inv,
                        (zab * P[1] + zcd * Q[1]) * zabcd_inv,
                        (zab * P[2] + zcd * Q[2]) * zabcd_inv,
                    };

                    double WP[3] = {
                        W[0] - PA[0],
                        W[1] - PA[1],
                        W[2] - PA[2],
                    };

                    double WQ[3] = {
                        W[0] - Q[0],
                        W[1] - Q[1],
                        W[2] - Q[2],
                    };

                    double T = rho * ((P[0] - Q[0]) * (P[0] - Q[0]) + (P[1] - Q[1]) * (P[1] - Q[1]) + (P[2] - Q[2]) * (P[2] - Q[2]));
                    double K = zabcd_inv_sqrt * (AB_K[(int)iab]) * CD_K[(int)icd];

                    double fm[5];
                    // vgamma<double>(4, T, fm);

                    for (int i = 0; i < 5; i++) {
                        fm[i] *= K;
                    }

                    double s000_s000_s000_s000_0 = fm[0] ;
                    double s000_s000_s000_s000_1 = fm[1] ;
                    double s000_s000_s000_s000_2 = fm[2] ;
                    double s000_s000_s000_s000_3 = fm[3] ;
                    double s000_s000_s000_s000_4 = fm[4] ;
                    double s000_s000_p100_s000_0 = QC[0] * s000_s000_s000_s000_0 + WQ[0] * s000_s000_s000_s000_1;
                    double s000_s000_p100_s000_1 = QC[0] * s000_s000_s000_s000_1 + WQ[0] * s000_s000_s000_s000_2;
                    double s000_s000_p100_s000_2 = QC[0] * s000_s000_s000_s000_2 + WQ[0] * s000_s000_s000_s000_3;
                    double s000_s000_p100_s000_3 = QC[0] * s000_s000_s000_s000_3 + WQ[0] * s000_s000_s000_s000_4;
                    double s000_s000_p010_s000_0 = QC[1] * s000_s000_s000_s000_0 + WQ[1] * s000_s000_s000_s000_1;
                    double s000_s000_p010_s000_1 = QC[1] * s000_s000_s000_s000_1 + WQ[1] * s000_s000_s000_s000_2;
                    double s000_s000_p010_s000_2 = QC[1] * s000_s000_s000_s000_2 + WQ[1] * s000_s000_s000_s000_3;
                    double s000_s000_p010_s000_3 = QC[1] * s000_s000_s000_s000_3 + WQ[1] * s000_s000_s000_s000_4;
                    double s000_s000_p001_s000_0 = QC[2] * s000_s000_s000_s000_0 + WQ[2] * s000_s000_s000_s000_1;
                    double s000_s000_p001_s000_1 = QC[2] * s000_s000_s000_s000_1 + WQ[2] * s000_s000_s000_s000_2;
                    double s000_s000_p001_s000_2 = QC[2] * s000_s000_s000_s000_2 + WQ[2] * s000_s000_s000_s000_3;
                    double s000_s000_p001_s000_3 = QC[2] * s000_s000_s000_s000_3 + WQ[2] * s000_s000_s000_s000_4;
                    double s000_s000_d200_s000_0 = QC[0] * s000_s000_p100_s000_0 + WQ[0] * s000_s000_p100_s000_1 + 0.5 * zcd_inv * 1 * s000_s000_s000_s000_0 + (- rho) * zcd_inv * 0.5 * zcd_inv * 1 * s000_s000_s000_s000_1;
                    double s000_s000_d200_s000_1 = QC[0] * s000_s000_p100_s000_1 + WQ[0] * s000_s000_p100_s000_2 + 0.5 * zcd_inv * 1 * s000_s000_s000_s000_1 + (- rho) * zcd_inv * 0.5 * zcd_inv * 1 * s000_s000_s000_s000_2;
                    double s000_s000_d200_s000_2 = QC[0] * s000_s000_p100_s000_2 + WQ[0] * s000_s000_p100_s000_3 + 0.5 * zcd_inv * 1 * s000_s000_s000_s000_2 + (- rho) * zcd_inv * 0.5 * zcd_inv * 1 * s000_s000_s000_s000_3;
                    double s000_s000_d110_s000_0 = QC[0] * s000_s000_p010_s000_0 + WQ[0] * s000_s000_p010_s000_1;
                    double s000_s000_d110_s000_1 = QC[0] * s000_s000_p010_s000_1 + WQ[0] * s000_s000_p010_s000_2;
                    double s000_s000_d110_s000_2 = QC[0] * s000_s000_p010_s000_2 + WQ[0] * s000_s000_p010_s000_3;
                    double s000_s000_d101_s000_0 = QC[0] * s000_s000_p001_s000_0 + WQ[0] * s000_s000_p001_s000_1;
                    double s000_s000_d101_s000_1 = QC[0] * s000_s000_p001_s000_1 + WQ[0] * s000_s000_p001_s000_2;
                    double s000_s000_d101_s000_2 = QC[0] * s000_s000_p001_s000_2 + WQ[0] * s000_s000_p001_s000_3;
                    double s000_s000_d020_s000_0 = QC[1] * s000_s000_p010_s000_0 + WQ[1] * s000_s000_p010_s000_1 + 0.5 * zcd_inv * 1 * s000_s000_s000_s000_0 + (- rho) * zcd_inv * 0.5 * zcd_inv * 1 * s000_s000_s000_s000_1;
                    double s000_s000_d020_s000_1 = QC[1] * s000_s000_p010_s000_1 + WQ[1] * s000_s000_p010_s000_2 + 0.5 * zcd_inv * 1 * s000_s000_s000_s000_1 + (- rho) * zcd_inv * 0.5 * zcd_inv * 1 * s000_s000_s000_s000_2;
                    double s000_s000_d020_s000_2 = QC[1] * s000_s000_p010_s000_2 + WQ[1] * s000_s000_p010_s000_3 + 0.5 * zcd_inv * 1 * s000_s000_s000_s000_2 + (- rho) * zcd_inv * 0.5 * zcd_inv * 1 * s000_s000_s000_s000_3;
                    double s000_s000_d011_s000_0 = QC[1] * s000_s000_p001_s000_0 + WQ[1] * s000_s000_p001_s000_1;
                    double s000_s000_d011_s000_1 = QC[1] * s000_s000_p001_s000_1 + WQ[1] * s000_s000_p001_s000_2;
                    double s000_s000_d011_s000_2 = QC[1] * s000_s000_p001_s000_2 + WQ[1] * s000_s000_p001_s000_3;
                    double s000_s000_d002_s000_0 = QC[2] * s000_s000_p001_s000_0 + WQ[2] * s000_s000_p001_s000_1 + 0.5 * zcd_inv * 1 * s000_s000_s000_s000_0 + (- rho) * zcd_inv * 0.5 * zcd_inv * 1 * s000_s000_s000_s000_1;
                    double s000_s000_d002_s000_1 = QC[2] * s000_s000_p001_s000_1 + WQ[2] * s000_s000_p001_s000_2 + 0.5 * zcd_inv * 1 * s000_s000_s000_s000_1 + (- rho) * zcd_inv * 0.5 * zcd_inv * 1 * s000_s000_s000_s000_2;
                    double s000_s000_d002_s000_2 = QC[2] * s000_s000_p001_s000_2 + WQ[2] * s000_s000_p001_s000_3 + 0.5 * zcd_inv * 1 * s000_s000_s000_s000_2 + (- rho) * zcd_inv * 0.5 * zcd_inv * 1 * s000_s000_s000_s000_3;
                    double p100_s000_s000_s000_0 = PA[0] * s000_s000_s000_s000_0 + WP[0] * s000_s000_s000_s000_1;
                    double p100_s000_s000_s000_1 = PA[0] * s000_s000_s000_s000_1 + WP[0] * s000_s000_s000_s000_2;
                    double p100_s000_s000_s000_2 = PA[0] * s000_s000_s000_s000_2 + WP[0] * s000_s000_s000_s000_3;
                    double p010_s000_s000_s000_0 = PA[1] * s000_s000_s000_s000_0 + WP[1] * s000_s000_s000_s000_1;
                    double p010_s000_s000_s000_1 = PA[1] * s000_s000_s000_s000_1 + WP[1] * s000_s000_s000_s000_2;
                    double p010_s000_s000_s000_2 = PA[1] * s000_s000_s000_s000_2 + WP[1] * s000_s000_s000_s000_3;
                    double p001_s000_s000_s000_0 = PA[2] * s000_s000_s000_s000_0 + WP[2] * s000_s000_s000_s000_1;
                    double p001_s000_s000_s000_1 = PA[2] * s000_s000_s000_s000_1 + WP[2] * s000_s000_s000_s000_2;
                    double p001_s000_s000_s000_2 = PA[2] * s000_s000_s000_s000_2 + WP[2] * s000_s000_s000_s000_3;
                    double p100_s000_p100_s000_0 = PA[0] * s000_s000_p100_s000_0 + WP[0] * s000_s000_p100_s000_1 + 0.5 * zabcd_inv * 1 * s000_s000_s000_s000_1;
                    double p100_s000_p100_s000_1 = PA[0] * s000_s000_p100_s000_1 + WP[0] * s000_s000_p100_s000_2 + 0.5 * zabcd_inv * 1 * s000_s000_s000_s000_2;
                    double p100_s000_p010_s000_0 = PA[0] * s000_s000_p010_s000_0 + WP[0] * s000_s000_p010_s000_1;
                    double p100_s000_p010_s000_1 = PA[0] * s000_s000_p010_s000_1 + WP[0] * s000_s000_p010_s000_2;
                    double p100_s000_p001_s000_0 = PA[0] * s000_s000_p001_s000_0 + WP[0] * s000_s000_p001_s000_1;
                    double p100_s000_p001_s000_1 = PA[0] * s000_s000_p001_s000_1 + WP[0] * s000_s000_p001_s000_2;
                    double p010_s000_p100_s000_0 = PA[1] * s000_s000_p100_s000_0 + WP[1] * s000_s000_p100_s000_1;
                    double p010_s000_p100_s000_1 = PA[1] * s000_s000_p100_s000_1 + WP[1] * s000_s000_p100_s000_2;
                    double p010_s000_p010_s000_0 = PA[1] * s000_s000_p010_s000_0 + WP[1] * s000_s000_p010_s000_1 + 0.5 * zabcd_inv * 1 * s000_s000_s000_s000_1;
                    double p010_s000_p010_s000_1 = PA[1] * s000_s000_p010_s000_1 + WP[1] * s000_s000_p010_s000_2 + 0.5 * zabcd_inv * 1 * s000_s000_s000_s000_2;
                    double p010_s000_p001_s000_0 = PA[1] * s000_s000_p001_s000_0 + WP[1] * s000_s000_p001_s000_1;
                    double p010_s000_p001_s000_1 = PA[1] * s000_s000_p001_s000_1 + WP[1] * s000_s000_p001_s000_2;
                    double p001_s000_p100_s000_0 = PA[2] * s000_s000_p100_s000_0 + WP[2] * s000_s000_p100_s000_1;
                    double p001_s000_p100_s000_1 = PA[2] * s000_s000_p100_s000_1 + WP[2] * s000_s000_p100_s000_2;
                    double p001_s000_p010_s000_0 = PA[2] * s000_s000_p010_s000_0 + WP[2] * s000_s000_p010_s000_1;
                    double p001_s000_p010_s000_1 = PA[2] * s000_s000_p010_s000_1 + WP[2] * s000_s000_p010_s000_2;
                    double p001_s000_p001_s000_0 = PA[2] * s000_s000_p001_s000_0 + WP[2] * s000_s000_p001_s000_1 + 0.5 * zabcd_inv * 1 * s000_s000_s000_s000_1;
                    double p001_s000_p001_s000_1 = PA[2] * s000_s000_p001_s000_1 + WP[2] * s000_s000_p001_s000_2 + 0.5 * zabcd_inv * 1 * s000_s000_s000_s000_2;
                    double p100_s000_d200_s000_0 = PA[0] * s000_s000_d200_s000_0 + WP[0] * s000_s000_d200_s000_1 + 0.5 * zabcd_inv * 2 * s000_s000_p100_s000_1;
                    double p100_s000_d200_s000_1 = PA[0] * s000_s000_d200_s000_1 + WP[0] * s000_s000_d200_s000_2 + 0.5 * zabcd_inv * 2 * s000_s000_p100_s000_2;
                    double p100_s000_d110_s000_0 = PA[0] * s000_s000_d110_s000_0 + WP[0] * s000_s000_d110_s000_1 + 0.5 * zabcd_inv * 1 * s000_s000_p010_s000_1;
                    double p100_s000_d110_s000_1 = PA[0] * s000_s000_d110_s000_1 + WP[0] * s000_s000_d110_s000_2 + 0.5 * zabcd_inv * 1 * s000_s000_p010_s000_2;
                    double p100_s000_d101_s000_0 = PA[0] * s000_s000_d101_s000_0 + WP[0] * s000_s000_d101_s000_1 + 0.5 * zabcd_inv * 1 * s000_s000_p001_s000_1;
                    double p100_s000_d101_s000_1 = PA[0] * s000_s000_d101_s000_1 + WP[0] * s000_s000_d101_s000_2 + 0.5 * zabcd_inv * 1 * s000_s000_p001_s000_2;
                    double p100_s000_d020_s000_0 = PA[0] * s000_s000_d020_s000_0 + WP[0] * s000_s000_d020_s000_1;
                    double p100_s000_d020_s000_1 = PA[0] * s000_s000_d020_s000_1 + WP[0] * s000_s000_d020_s000_2;
                    double p100_s000_d011_s000_0 = PA[0] * s000_s000_d011_s000_0 + WP[0] * s000_s000_d011_s000_1;
                    double p100_s000_d011_s000_1 = PA[0] * s000_s000_d011_s000_1 + WP[0] * s000_s000_d011_s000_2;
                    double p100_s000_d002_s000_0 = PA[0] * s000_s000_d002_s000_0 + WP[0] * s000_s000_d002_s000_1;
                    double p100_s000_d002_s000_1 = PA[0] * s000_s000_d002_s000_1 + WP[0] * s000_s000_d002_s000_2;
                    double p010_s000_d200_s000_0 = PA[1] * s000_s000_d200_s000_0 + WP[1] * s000_s000_d200_s000_1;
                    double p010_s000_d200_s000_1 = PA[1] * s000_s000_d200_s000_1 + WP[1] * s000_s000_d200_s000_2;
                    double p010_s000_d110_s000_0 = PA[1] * s000_s000_d110_s000_0 + WP[1] * s000_s000_d110_s000_1 + 0.5 * zabcd_inv * 1 * s000_s000_p100_s000_1;
                    double p010_s000_d110_s000_1 = PA[1] * s000_s000_d110_s000_1 + WP[1] * s000_s000_d110_s000_2 + 0.5 * zabcd_inv * 1 * s000_s000_p100_s000_2;
                    double p010_s000_d101_s000_0 = PA[1] * s000_s000_d101_s000_0 + WP[1] * s000_s000_d101_s000_1;
                    double p010_s000_d101_s000_1 = PA[1] * s000_s000_d101_s000_1 + WP[1] * s000_s000_d101_s000_2;
                    double p010_s000_d020_s000_0 = PA[1] * s000_s000_d020_s000_0 + WP[1] * s000_s000_d020_s000_1 + 0.5 * zabcd_inv * 2 * s000_s000_p010_s000_1;
                    double p010_s000_d020_s000_1 = PA[1] * s000_s000_d020_s000_1 + WP[1] * s000_s000_d020_s000_2 + 0.5 * zabcd_inv * 2 * s000_s000_p010_s000_2;
                    double p010_s000_d011_s000_0 = PA[1] * s000_s000_d011_s000_0 + WP[1] * s000_s000_d011_s000_1 + 0.5 * zabcd_inv * 1 * s000_s000_p001_s000_1;
                    double p010_s000_d011_s000_1 = PA[1] * s000_s000_d011_s000_1 + WP[1] * s000_s000_d011_s000_2 + 0.5 * zabcd_inv * 1 * s000_s000_p001_s000_2;
                    double p010_s000_d002_s000_0 = PA[1] * s000_s000_d002_s000_0 + WP[1] * s000_s000_d002_s000_1;
                    double p010_s000_d002_s000_1 = PA[1] * s000_s000_d002_s000_1 + WP[1] * s000_s000_d002_s000_2;
                    double p001_s000_d200_s000_0 = PA[2] * s000_s000_d200_s000_0 + WP[2] * s000_s000_d200_s000_1;
                    double p001_s000_d200_s000_1 = PA[2] * s000_s000_d200_s000_1 + WP[2] * s000_s000_d200_s000_2;
                    double p001_s000_d110_s000_0 = PA[2] * s000_s000_d110_s000_0 + WP[2] * s000_s000_d110_s000_1;
                    double p001_s000_d110_s000_1 = PA[2] * s000_s000_d110_s000_1 + WP[2] * s000_s000_d110_s000_2;
                    double p001_s000_d101_s000_0 = PA[2] * s000_s000_d101_s000_0 + WP[2] * s000_s000_d101_s000_1 + 0.5 * zabcd_inv * 1 * s000_s000_p100_s000_1;
                    double p001_s000_d101_s000_1 = PA[2] * s000_s000_d101_s000_1 + WP[2] * s000_s000_d101_s000_2 + 0.5 * zabcd_inv * 1 * s000_s000_p100_s000_2;
                    double p001_s000_d020_s000_0 = PA[2] * s000_s000_d020_s000_0 + WP[2] * s000_s000_d020_s000_1;
                    double p001_s000_d020_s000_1 = PA[2] * s000_s000_d020_s000_1 + WP[2] * s000_s000_d020_s000_2;
                    double p001_s000_d011_s000_0 = PA[2] * s000_s000_d011_s000_0 + WP[2] * s000_s000_d011_s000_1 + 0.5 * zabcd_inv * 1 * s000_s000_p010_s000_1;
                    double p001_s000_d011_s000_1 = PA[2] * s000_s000_d011_s000_1 + WP[2] * s000_s000_d011_s000_2 + 0.5 * zabcd_inv * 1 * s000_s000_p010_s000_2;
                    double p001_s000_d002_s000_0 = PA[2] * s000_s000_d002_s000_0 + WP[2] * s000_s000_d002_s000_1 + 0.5 * zabcd_inv * 2 * s000_s000_p001_s000_1;
                    double p001_s000_d002_s000_1 = PA[2] * s000_s000_d002_s000_1 + WP[2] * s000_s000_d002_s000_2 + 0.5 * zabcd_inv * 2 * s000_s000_p001_s000_2;
                    double d200_s000_s000_s000_0 = PA[0] * p100_s000_s000_s000_0 + WP[0] * p100_s000_s000_s000_1 + 0.5 * zab_inv * 1 * s000_s000_s000_s000_0 + (- rho) * zab_inv * 0.5 * zab_inv * 1 * s000_s000_s000_s000_1;
                    double d200_s000_s000_s000_1 = PA[0] * p100_s000_s000_s000_1 + WP[0] * p100_s000_s000_s000_2 + 0.5 * zab_inv * 1 * s000_s000_s000_s000_1 + (- rho) * zab_inv * 0.5 * zab_inv * 1 * s000_s000_s000_s000_2;
                    double d110_s000_s000_s000_0 = PA[0] * p010_s000_s000_s000_0 + WP[0] * p010_s000_s000_s000_1;
                    double d110_s000_s000_s000_1 = PA[0] * p010_s000_s000_s000_1 + WP[0] * p010_s000_s000_s000_2;
                    double d101_s000_s000_s000_0 = PA[0] * p001_s000_s000_s000_0 + WP[0] * p001_s000_s000_s000_1;
                    double d101_s000_s000_s000_1 = PA[0] * p001_s000_s000_s000_1 + WP[0] * p001_s000_s000_s000_2;
                    double d020_s000_s000_s000_0 = PA[1] * p010_s000_s000_s000_0 + WP[1] * p010_s000_s000_s000_1 + 0.5 * zab_inv * 1 * s000_s000_s000_s000_0 + (- rho) * zab_inv * 0.5 * zab_inv * 1 * s000_s000_s000_s000_1;
                    double d020_s000_s000_s000_1 = PA[1] * p010_s000_s000_s000_1 + WP[1] * p010_s000_s000_s000_2 + 0.5 * zab_inv * 1 * s000_s000_s000_s000_1 + (- rho) * zab_inv * 0.5 * zab_inv * 1 * s000_s000_s000_s000_2;
                    double d011_s000_s000_s000_0 = PA[1] * p001_s000_s000_s000_0 + WP[1] * p001_s000_s000_s000_1;
                    double d011_s000_s000_s000_1 = PA[1] * p001_s000_s000_s000_1 + WP[1] * p001_s000_s000_s000_2;
                    double d002_s000_s000_s000_0 = PA[2] * p001_s000_s000_s000_0 + WP[2] * p001_s000_s000_s000_1 + 0.5 * zab_inv * 1 * s000_s000_s000_s000_0 + (- rho) * zab_inv * 0.5 * zab_inv * 1 * s000_s000_s000_s000_1;
                    double d002_s000_s000_s000_1 = PA[2] * p001_s000_s000_s000_1 + WP[2] * p001_s000_s000_s000_2 + 0.5 * zab_inv * 1 * s000_s000_s000_s000_1 + (- rho) * zab_inv * 0.5 * zab_inv * 1 * s000_s000_s000_s000_2;
                    double d200_s000_p100_s000_0 = QC[0] * d200_s000_s000_s000_0 + WQ[0] * d200_s000_s000_s000_1 +  0.0 * p100_s000_s000_s000_0 + 0.5 * zabcd_inv * 2 * p100_s000_s000_s000_1;
                    double d200_s000_p010_s000_0 = QC[1] * d200_s000_s000_s000_0 + WQ[1] * d200_s000_s000_s000_1;
                    double d200_s000_p001_s000_0 = QC[2] * d200_s000_s000_s000_0 + WQ[2] * d200_s000_s000_s000_1;
                    double d110_s000_p100_s000_0 = QC[0] * d110_s000_s000_s000_0 + WQ[0] * d110_s000_s000_s000_1 +  0.0 * p010_s000_s000_s000_0 + 0.5 * zabcd_inv * 1 * p010_s000_s000_s000_1;
                    double d110_s000_p010_s000_0 = QC[1] * d110_s000_s000_s000_0 + WQ[1] * d110_s000_s000_s000_1 +  0.0 * p100_s000_s000_s000_0 + 0.5 * zabcd_inv * 1 * p100_s000_s000_s000_1;
                    double d110_s000_p001_s000_0 = QC[2] * d110_s000_s000_s000_0 + WQ[2] * d110_s000_s000_s000_1;
                    double d101_s000_p100_s000_0 = QC[0] * d101_s000_s000_s000_0 + WQ[0] * d101_s000_s000_s000_1 +  0.0 * p001_s000_s000_s000_0 + 0.5 * zabcd_inv * 1 * p001_s000_s000_s000_1;
                    double d101_s000_p010_s000_0 = QC[1] * d101_s000_s000_s000_0 + WQ[1] * d101_s000_s000_s000_1;
                    double d101_s000_p001_s000_0 = QC[2] * d101_s000_s000_s000_0 + WQ[2] * d101_s000_s000_s000_1 +  0.0 * p100_s000_s000_s000_0 + 0.5 * zabcd_inv * 1 * p100_s000_s000_s000_1;
                    double d020_s000_p100_s000_0 = QC[0] * d020_s000_s000_s000_0 + WQ[0] * d020_s000_s000_s000_1;
                    double d020_s000_p010_s000_0 = QC[1] * d020_s000_s000_s000_0 + WQ[1] * d020_s000_s000_s000_1 +  0.0 * p010_s000_s000_s000_0 + 0.5 * zabcd_inv * 2 * p010_s000_s000_s000_1;
                    double d020_s000_p001_s000_0 = QC[2] * d020_s000_s000_s000_0 + WQ[2] * d020_s000_s000_s000_1;
                    double d011_s000_p100_s000_0 = QC[0] * d011_s000_s000_s000_0 + WQ[0] * d011_s000_s000_s000_1;
                    double d011_s000_p010_s000_0 = QC[1] * d011_s000_s000_s000_0 + WQ[1] * d011_s000_s000_s000_1 +  0.0 * p001_s000_s000_s000_0 + 0.5 * zabcd_inv * 1 * p001_s000_s000_s000_1;
                    double d011_s000_p001_s000_0 = QC[2] * d011_s000_s000_s000_0 + WQ[2] * d011_s000_s000_s000_1 +  0.0 * p010_s000_s000_s000_0 + 0.5 * zabcd_inv * 1 * p010_s000_s000_s000_1;
                    double d002_s000_p100_s000_0 = QC[0] * d002_s000_s000_s000_0 + WQ[0] * d002_s000_s000_s000_1;
                    double d002_s000_p010_s000_0 = QC[1] * d002_s000_s000_s000_0 + WQ[1] * d002_s000_s000_s000_1;
                    double d002_s000_p001_s000_0 = QC[2] * d002_s000_s000_s000_0 + WQ[2] * d002_s000_s000_s000_1 +  0.0 * p001_s000_s000_s000_0 + 0.5 * zabcd_inv * 2 * p001_s000_s000_s000_1;
                    double d200_s000_d200_s000_0 = PA[0] * p100_s000_d200_s000_0 + WP[0] * p100_s000_d200_s000_1 + 0.5 * zab_inv * 1 * s000_s000_d200_s000_0 + (- rho) * zab_inv * 0.5 * zab_inv * 1 * s000_s000_d200_s000_1 + 0.5 * zabcd_inv * 2 * p100_s000_p100_s000_1;
                    double d200_s000_d110_s000_0 = PA[0] * p100_s000_d110_s000_0 + WP[0] * p100_s000_d110_s000_1 + 0.5 * zab_inv * 1 * s000_s000_d110_s000_0 + (- rho) * zab_inv * 0.5 * zab_inv * 1 * s000_s000_d110_s000_1 + 0.5 * zabcd_inv * 1 * p100_s000_p010_s000_1;
                    double d200_s000_d101_s000_0 = PA[0] * p100_s000_d101_s000_0 + WP[0] * p100_s000_d101_s000_1 + 0.5 * zab_inv * 1 * s000_s000_d101_s000_0 + (- rho) * zab_inv * 0.5 * zab_inv * 1 * s000_s000_d101_s000_1 + 0.5 * zabcd_inv * 1 * p100_s000_p001_s000_1;
                    double d200_s000_d020_s000_0 = PA[0] * p100_s000_d020_s000_0 + WP[0] * p100_s000_d020_s000_1 + 0.5 * zab_inv * 1 * s000_s000_d020_s000_0 + (- rho) * zab_inv * 0.5 * zab_inv * 1 * s000_s000_d020_s000_1;
                    double d200_s000_d011_s000_0 = PA[0] * p100_s000_d011_s000_0 + WP[0] * p100_s000_d011_s000_1 + 0.5 * zab_inv * 1 * s000_s000_d011_s000_0 + (- rho) * zab_inv * 0.5 * zab_inv * 1 * s000_s000_d011_s000_1;
                    double d200_s000_d002_s000_0 = PA[0] * p100_s000_d002_s000_0 + WP[0] * p100_s000_d002_s000_1 + 0.5 * zab_inv * 1 * s000_s000_d002_s000_0 + (- rho) * zab_inv * 0.5 * zab_inv * 1 * s000_s000_d002_s000_1;
                    double d110_s000_d200_s000_0 = PA[0] * p010_s000_d200_s000_0 + WP[0] * p010_s000_d200_s000_1 + 0.5 * zabcd_inv * 2 * p010_s000_p100_s000_1;
                    double d110_s000_d110_s000_0 = PA[0] * p010_s000_d110_s000_0 + WP[0] * p010_s000_d110_s000_1 + 0.5 * zabcd_inv * 1 * p010_s000_p010_s000_1;
                    double d110_s000_d101_s000_0 = PA[0] * p010_s000_d101_s000_0 + WP[0] * p010_s000_d101_s000_1 + 0.5 * zabcd_inv * 1 * p010_s000_p001_s000_1;
                    double d110_s000_d020_s000_0 = PA[0] * p010_s000_d020_s000_0 + WP[0] * p010_s000_d020_s000_1;
                    double d110_s000_d011_s000_0 = PA[0] * p010_s000_d011_s000_0 + WP[0] * p010_s000_d011_s000_1;
                    double d110_s000_d002_s000_0 = PA[0] * p010_s000_d002_s000_0 + WP[0] * p010_s000_d002_s000_1;
                    double d101_s000_d200_s000_0 = PA[0] * p001_s000_d200_s000_0 + WP[0] * p001_s000_d200_s000_1 + 0.5 * zabcd_inv * 2 * p001_s000_p100_s000_1;
                    double d101_s000_d110_s000_0 = PA[0] * p001_s000_d110_s000_0 + WP[0] * p001_s000_d110_s000_1 + 0.5 * zabcd_inv * 1 * p001_s000_p010_s000_1;
                    double d101_s000_d101_s000_0 = PA[0] * p001_s000_d101_s000_0 + WP[0] * p001_s000_d101_s000_1 + 0.5 * zabcd_inv * 1 * p001_s000_p001_s000_1;
                    double d101_s000_d020_s000_0 = PA[0] * p001_s000_d020_s000_0 + WP[0] * p001_s000_d020_s000_1;
                    double d101_s000_d011_s000_0 = PA[0] * p001_s000_d011_s000_0 + WP[0] * p001_s000_d011_s000_1;
                    double d101_s000_d002_s000_0 = PA[0] * p001_s000_d002_s000_0 + WP[0] * p001_s000_d002_s000_1;
                    double d020_s000_d200_s000_0 = PA[1] * p010_s000_d200_s000_0 + WP[1] * p010_s000_d200_s000_1 + 0.5 * zab_inv * 1 * s000_s000_d200_s000_0 + (- rho) * zab_inv * 0.5 * zab_inv * 1 * s000_s000_d200_s000_1;
                    double d020_s000_d110_s000_0 = PA[1] * p010_s000_d110_s000_0 + WP[1] * p010_s000_d110_s000_1 + 0.5 * zab_inv * 1 * s000_s000_d110_s000_0 + (- rho) * zab_inv * 0.5 * zab_inv * 1 * s000_s000_d110_s000_1 + 0.5 * zabcd_inv * 1 * p010_s000_p100_s000_1;
                    double d020_s000_d101_s000_0 = PA[1] * p010_s000_d101_s000_0 + WP[1] * p010_s000_d101_s000_1 + 0.5 * zab_inv * 1 * s000_s000_d101_s000_0 + (- rho) * zab_inv * 0.5 * zab_inv * 1 * s000_s000_d101_s000_1;
                    double d020_s000_d020_s000_0 = PA[1] * p010_s000_d020_s000_0 + WP[1] * p010_s000_d020_s000_1 + 0.5 * zab_inv * 1 * s000_s000_d020_s000_0 + (- rho) * zab_inv * 0.5 * zab_inv * 1 * s000_s000_d020_s000_1 + 0.5 * zabcd_inv * 2 * p010_s000_p010_s000_1;
                    double d020_s000_d011_s000_0 = PA[1] * p010_s000_d011_s000_0 + WP[1] * p010_s000_d011_s000_1 + 0.5 * zab_inv * 1 * s000_s000_d011_s000_0 + (- rho) * zab_inv * 0.5 * zab_inv * 1 * s000_s000_d011_s000_1 + 0.5 * zabcd_inv * 1 * p010_s000_p001_s000_1;
                    double d020_s000_d002_s000_0 = PA[1] * p010_s000_d002_s000_0 + WP[1] * p010_s000_d002_s000_1 + 0.5 * zab_inv * 1 * s000_s000_d002_s000_0 + (- rho) * zab_inv * 0.5 * zab_inv * 1 * s000_s000_d002_s000_1;
                    double d011_s000_d200_s000_0 = PA[1] * p001_s000_d200_s000_0 + WP[1] * p001_s000_d200_s000_1;
                    double d011_s000_d110_s000_0 = PA[1] * p001_s000_d110_s000_0 + WP[1] * p001_s000_d110_s000_1 + 0.5 * zabcd_inv * 1 * p001_s000_p100_s000_1;
                    double d011_s000_d101_s000_0 = PA[1] * p001_s000_d101_s000_0 + WP[1] * p001_s000_d101_s000_1;
                    double d011_s000_d020_s000_0 = PA[1] * p001_s000_d020_s000_0 + WP[1] * p001_s000_d020_s000_1 + 0.5 * zabcd_inv * 2 * p001_s000_p010_s000_1;
                    double d011_s000_d011_s000_0 = PA[1] * p001_s000_d011_s000_0 + WP[1] * p001_s000_d011_s000_1 + 0.5 * zabcd_inv * 1 * p001_s000_p001_s000_1;
                    double d011_s000_d002_s000_0 = PA[1] * p001_s000_d002_s000_0 + WP[1] * p001_s000_d002_s000_1;
                    double d002_s000_d200_s000_0 = PA[2] * p001_s000_d200_s000_0 + WP[2] * p001_s000_d200_s000_1 + 0.5 * zab_inv * 1 * s000_s000_d200_s000_0 + (- rho) * zab_inv * 0.5 * zab_inv * 1 * s000_s000_d200_s000_1;
                    double d002_s000_d110_s000_0 = PA[2] * p001_s000_d110_s000_0 + WP[2] * p001_s000_d110_s000_1 + 0.5 * zab_inv * 1 * s000_s000_d110_s000_0 + (- rho) * zab_inv * 0.5 * zab_inv * 1 * s000_s000_d110_s000_1;
                    double d002_s000_d101_s000_0 = PA[2] * p001_s000_d101_s000_0 + WP[2] * p001_s000_d101_s000_1 + 0.5 * zab_inv * 1 * s000_s000_d101_s000_0 + (- rho) * zab_inv * 0.5 * zab_inv * 1 * s000_s000_d101_s000_1 + 0.5 * zabcd_inv * 1 * p001_s000_p100_s000_1;
                    double d002_s000_d020_s000_0 = PA[2] * p001_s000_d020_s000_0 + WP[2] * p001_s000_d020_s000_1 + 0.5 * zab_inv * 1 * s000_s000_d020_s000_0 + (- rho) * zab_inv * 0.5 * zab_inv * 1 * s000_s000_d020_s000_1;
                    double d002_s000_d011_s000_0 = PA[2] * p001_s000_d011_s000_0 + WP[2] * p001_s000_d011_s000_1 + 0.5 * zab_inv * 1 * s000_s000_d011_s000_0 + (- rho) * zab_inv * 0.5 * zab_inv * 1 * s000_s000_d011_s000_1 + 0.5 * zabcd_inv * 1 * p001_s000_p010_s000_1;
                    double d002_s000_d002_s000_0 = PA[2] * p001_s000_d002_s000_0 + WP[2] * p001_s000_d002_s000_1 + 0.5 * zab_inv * 1 * s000_s000_d002_s000_0 + (- rho) * zab_inv * 0.5 * zab_inv * 1 * s000_s000_d002_s000_1 + 0.5 * zabcd_inv * 2 * p001_s000_p001_s000_1;
                    d200_s000_d200_s000_0_con += d200_s000_d200_s000_0;
                    d200_s000_d110_s000_0_con += d200_s000_d110_s000_0;
                    d200_s000_d101_s000_0_con += d200_s000_d101_s000_0;
                    d200_s000_d020_s000_0_con += d200_s000_d020_s000_0;
                    d200_s000_d011_s000_0_con += d200_s000_d011_s000_0;
                    d200_s000_d002_s000_0_con += d200_s000_d002_s000_0;
                    d110_s000_d200_s000_0_con += d110_s000_d200_s000_0;
                    d110_s000_d110_s000_0_con += d110_s000_d110_s000_0;
                    d110_s000_d101_s000_0_con += d110_s000_d101_s000_0;
                    d110_s000_d020_s000_0_con += d110_s000_d020_s000_0;
                    d110_s000_d011_s000_0_con += d110_s000_d011_s000_0;
                    d110_s000_d002_s000_0_con += d110_s000_d002_s000_0;
                    d101_s000_d200_s000_0_con += d101_s000_d200_s000_0;
                    d101_s000_d110_s000_0_con += d101_s000_d110_s000_0;
                    d101_s000_d101_s000_0_con += d101_s000_d101_s000_0;
                    d101_s000_d020_s000_0_con += d101_s000_d020_s000_0;
                    d101_s000_d011_s000_0_con += d101_s000_d011_s000_0;
                    d101_s000_d002_s000_0_con += d101_s000_d002_s000_0;
                    d020_s000_d200_s000_0_con += d020_s000_d200_s000_0;
                    d020_s000_d110_s000_0_con += d020_s000_d110_s000_0;
                    d020_s000_d101_s000_0_con += d020_s000_d101_s000_0;
                    d020_s000_d020_s000_0_con += d020_s000_d020_s000_0;
                    d020_s000_d011_s000_0_con += d020_s000_d011_s000_0;
                    d020_s000_d002_s000_0_con += d020_s000_d002_s000_0;
                    d011_s000_d200_s000_0_con += d011_s000_d200_s000_0;
                    d011_s000_d110_s000_0_con += d011_s000_d110_s000_0;
                    d011_s000_d101_s000_0_con += d011_s000_d101_s000_0;
                    d011_s000_d020_s000_0_con += d011_s000_d020_s000_0;
                    d011_s000_d011_s000_0_con += d011_s000_d011_s000_0;
                    d011_s000_d002_s000_0_con += d011_s000_d002_s000_0;
                    d002_s000_d200_s000_0_con += d002_s000_d200_s000_0;
                    d002_s000_d110_s000_0_con += d002_s000_d110_s000_0;
                    d002_s000_d101_s000_0_con += d002_s000_d101_s000_0;
                    d002_s000_d020_s000_0_con += d002_s000_d020_s000_0;
                    d002_s000_d011_s000_0_con += d002_s000_d011_s000_0;
                    d002_s000_d002_s000_0_con += d002_s000_d002_s000_0;
                    d200_s000_p100_s000_0_con += d200_s000_p100_s000_0;
                    d200_s000_p010_s000_0_con += d200_s000_p010_s000_0;
                    d200_s000_p001_s000_0_con += d200_s000_p001_s000_0;
                    d110_s000_p100_s000_0_con += d110_s000_p100_s000_0;
                    d110_s000_p010_s000_0_con += d110_s000_p010_s000_0;
                    d110_s000_p001_s000_0_con += d110_s000_p001_s000_0;
                    d101_s000_p100_s000_0_con += d101_s000_p100_s000_0;
                    d101_s000_p010_s000_0_con += d101_s000_p010_s000_0;
                    d101_s000_p001_s000_0_con += d101_s000_p001_s000_0;
                    d020_s000_p100_s000_0_con += d020_s000_p100_s000_0;
                    d020_s000_p010_s000_0_con += d020_s000_p010_s000_0;
                    d020_s000_p001_s000_0_con += d020_s000_p001_s000_0;
                    d011_s000_p100_s000_0_con += d011_s000_p100_s000_0;
                    d011_s000_p010_s000_0_con += d011_s000_p010_s000_0;
                    d011_s000_p001_s000_0_con += d011_s000_p001_s000_0;
                    d002_s000_p100_s000_0_con += d002_s000_p100_s000_0;
                    d002_s000_p010_s000_0_con += d002_s000_p010_s000_0;
                    d002_s000_p001_s000_0_con += d002_s000_p001_s000_0;
                    p100_s000_d200_s000_0_con += p100_s000_d200_s000_0;
                    p100_s000_d110_s000_0_con += p100_s000_d110_s000_0;
                    p100_s000_d101_s000_0_con += p100_s000_d101_s000_0;
                    p100_s000_d020_s000_0_con += p100_s000_d020_s000_0;
                    p100_s000_d011_s000_0_con += p100_s000_d011_s000_0;
                    p100_s000_d002_s000_0_con += p100_s000_d002_s000_0;
                    p010_s000_d200_s000_0_con += p010_s000_d200_s000_0;
                    p010_s000_d110_s000_0_con += p010_s000_d110_s000_0;
                    p010_s000_d101_s000_0_con += p010_s000_d101_s000_0;
                    p010_s000_d020_s000_0_con += p010_s000_d020_s000_0;
                    p010_s000_d011_s000_0_con += p010_s000_d011_s000_0;
                    p010_s000_d002_s000_0_con += p010_s000_d002_s000_0;
                    p001_s000_d200_s000_0_con += p001_s000_d200_s000_0;
                    p001_s000_d110_s000_0_con += p001_s000_d110_s000_0;
                    p001_s000_d101_s000_0_con += p001_s000_d101_s000_0;
                    p001_s000_d020_s000_0_con += p001_s000_d020_s000_0;
                    p001_s000_d011_s000_0_con += p001_s000_d011_s000_0;
                    p001_s000_d002_s000_0_con += p001_s000_d002_s000_0;
                    p100_s000_p100_s000_0_con += p100_s000_p100_s000_0;
                    p100_s000_p010_s000_0_con += p100_s000_p010_s000_0;
                    p100_s000_p001_s000_0_con += p100_s000_p001_s000_0;
                    p010_s000_p100_s000_0_con += p010_s000_p100_s000_0;
                    p010_s000_p010_s000_0_con += p010_s000_p010_s000_0;
                    p010_s000_p001_s000_0_con += p010_s000_p001_s000_0;
                    p001_s000_p100_s000_0_con += p001_s000_p100_s000_0;
                    p001_s000_p010_s000_0_con += p001_s000_p010_s000_0;
                    p001_s000_p001_s000_0_con += p001_s000_p001_s000_0;  
                }
            }

            double p100_p100_d200_s000_0_con = d200_s000_d200_s000_0_con + AB[0] * p100_s000_d200_s000_0_con ; 
            double p100_p100_d110_s000_0_con = d200_s000_d110_s000_0_con + AB[0] * p100_s000_d110_s000_0_con ; 
            double p100_p100_d101_s000_0_con = d200_s000_d101_s000_0_con + AB[0] * p100_s000_d101_s000_0_con ; 
            double p100_p100_d020_s000_0_con = d200_s000_d020_s000_0_con + AB[0] * p100_s000_d020_s000_0_con ; 
            double p100_p100_d011_s000_0_con = d200_s000_d011_s000_0_con + AB[0] * p100_s000_d011_s000_0_con ; 
            double p100_p100_d002_s000_0_con = d200_s000_d002_s000_0_con + AB[0] * p100_s000_d002_s000_0_con ; 
            double p100_p010_d200_s000_0_con = d110_s000_d200_s000_0_con + AB[1] * p100_s000_d200_s000_0_con ; 
            double p100_p010_d110_s000_0_con = d110_s000_d110_s000_0_con + AB[1] * p100_s000_d110_s000_0_con ; 
            double p100_p010_d101_s000_0_con = d110_s000_d101_s000_0_con + AB[1] * p100_s000_d101_s000_0_con ; 
            double p100_p010_d020_s000_0_con = d110_s000_d020_s000_0_con + AB[1] * p100_s000_d020_s000_0_con ; 
            double p100_p010_d011_s000_0_con = d110_s000_d011_s000_0_con + AB[1] * p100_s000_d011_s000_0_con ; 
            double p100_p010_d002_s000_0_con = d110_s000_d002_s000_0_con + AB[1] * p100_s000_d002_s000_0_con ; 
            double p100_p001_d200_s000_0_con = d101_s000_d200_s000_0_con + AB[2] * p100_s000_d200_s000_0_con ; 
            double p100_p001_d110_s000_0_con = d101_s000_d110_s000_0_con + AB[2] * p100_s000_d110_s000_0_con ; 
            double p100_p001_d101_s000_0_con = d101_s000_d101_s000_0_con + AB[2] * p100_s000_d101_s000_0_con ; 
            double p100_p001_d020_s000_0_con = d101_s000_d020_s000_0_con + AB[2] * p100_s000_d020_s000_0_con ; 
            double p100_p001_d011_s000_0_con = d101_s000_d011_s000_0_con + AB[2] * p100_s000_d011_s000_0_con ; 
            double p100_p001_d002_s000_0_con = d101_s000_d002_s000_0_con + AB[2] * p100_s000_d002_s000_0_con ; 
            double p010_p100_d200_s000_0_con = d110_s000_d200_s000_0_con + AB[0] * p010_s000_d200_s000_0_con ; 
            double p010_p100_d110_s000_0_con = d110_s000_d110_s000_0_con + AB[0] * p010_s000_d110_s000_0_con ; 
            double p010_p100_d101_s000_0_con = d110_s000_d101_s000_0_con + AB[0] * p010_s000_d101_s000_0_con ; 
            double p010_p100_d020_s000_0_con = d110_s000_d020_s000_0_con + AB[0] * p010_s000_d020_s000_0_con ; 
            double p010_p100_d011_s000_0_con = d110_s000_d011_s000_0_con + AB[0] * p010_s000_d011_s000_0_con ; 
            double p010_p100_d002_s000_0_con = d110_s000_d002_s000_0_con + AB[0] * p010_s000_d002_s000_0_con ; 
            double p010_p010_d200_s000_0_con = d020_s000_d200_s000_0_con + AB[1] * p010_s000_d200_s000_0_con ; 
            double p010_p010_d110_s000_0_con = d020_s000_d110_s000_0_con + AB[1] * p010_s000_d110_s000_0_con ; 
            double p010_p010_d101_s000_0_con = d020_s000_d101_s000_0_con + AB[1] * p010_s000_d101_s000_0_con ; 
            double p010_p010_d020_s000_0_con = d020_s000_d020_s000_0_con + AB[1] * p010_s000_d020_s000_0_con ; 
            double p010_p010_d011_s000_0_con = d020_s000_d011_s000_0_con + AB[1] * p010_s000_d011_s000_0_con ; 
            double p010_p010_d002_s000_0_con = d020_s000_d002_s000_0_con + AB[1] * p010_s000_d002_s000_0_con ; 
            double p010_p001_d200_s000_0_con = d011_s000_d200_s000_0_con + AB[2] * p010_s000_d200_s000_0_con ; 
            double p010_p001_d110_s000_0_con = d011_s000_d110_s000_0_con + AB[2] * p010_s000_d110_s000_0_con ; 
            double p010_p001_d101_s000_0_con = d011_s000_d101_s000_0_con + AB[2] * p010_s000_d101_s000_0_con ; 
            double p010_p001_d020_s000_0_con = d011_s000_d020_s000_0_con + AB[2] * p010_s000_d020_s000_0_con ; 
            double p010_p001_d011_s000_0_con = d011_s000_d011_s000_0_con + AB[2] * p010_s000_d011_s000_0_con ; 
            double p010_p001_d002_s000_0_con = d011_s000_d002_s000_0_con + AB[2] * p010_s000_d002_s000_0_con ; 
            double p001_p100_d200_s000_0_con = d101_s000_d200_s000_0_con + AB[0] * p001_s000_d200_s000_0_con ; 
            double p001_p100_d110_s000_0_con = d101_s000_d110_s000_0_con + AB[0] * p001_s000_d110_s000_0_con ; 
            double p001_p100_d101_s000_0_con = d101_s000_d101_s000_0_con + AB[0] * p001_s000_d101_s000_0_con ; 
            double p001_p100_d020_s000_0_con = d101_s000_d020_s000_0_con + AB[0] * p001_s000_d020_s000_0_con ; 
            double p001_p100_d011_s000_0_con = d101_s000_d011_s000_0_con + AB[0] * p001_s000_d011_s000_0_con ; 
            double p001_p100_d002_s000_0_con = d101_s000_d002_s000_0_con + AB[0] * p001_s000_d002_s000_0_con ; 
            double p001_p010_d200_s000_0_con = d011_s000_d200_s000_0_con + AB[1] * p001_s000_d200_s000_0_con ; 
            double p001_p010_d110_s000_0_con = d011_s000_d110_s000_0_con + AB[1] * p001_s000_d110_s000_0_con ; 
            double p001_p010_d101_s000_0_con = d011_s000_d101_s000_0_con + AB[1] * p001_s000_d101_s000_0_con ; 
            double p001_p010_d020_s000_0_con = d011_s000_d020_s000_0_con + AB[1] * p001_s000_d020_s000_0_con ; 
            double p001_p010_d011_s000_0_con = d011_s000_d011_s000_0_con + AB[1] * p001_s000_d011_s000_0_con ; 
            double p001_p010_d002_s000_0_con = d011_s000_d002_s000_0_con + AB[1] * p001_s000_d002_s000_0_con ; 
            double p001_p001_d200_s000_0_con = d002_s000_d200_s000_0_con + AB[2] * p001_s000_d200_s000_0_con ; 
            double p001_p001_d110_s000_0_con = d002_s000_d110_s000_0_con + AB[2] * p001_s000_d110_s000_0_con ; 
            double p001_p001_d101_s000_0_con = d002_s000_d101_s000_0_con + AB[2] * p001_s000_d101_s000_0_con ; 
            double p001_p001_d020_s000_0_con = d002_s000_d020_s000_0_con + AB[2] * p001_s000_d020_s000_0_con ; 
            double p001_p001_d011_s000_0_con = d002_s000_d011_s000_0_con + AB[2] * p001_s000_d011_s000_0_con ; 
            double p001_p001_d002_s000_0_con = d002_s000_d002_s000_0_con + AB[2] * p001_s000_d002_s000_0_con ; 
            double p100_p100_p100_s000_0_con = d200_s000_p100_s000_0_con + AB[0] * p100_s000_p100_s000_0_con ; 
            double p100_p100_p010_s000_0_con = d200_s000_p010_s000_0_con + AB[0] * p100_s000_p010_s000_0_con ; 
            double p100_p100_p001_s000_0_con = d200_s000_p001_s000_0_con + AB[0] * p100_s000_p001_s000_0_con ; 
            double p100_p010_p100_s000_0_con = d110_s000_p100_s000_0_con + AB[1] * p100_s000_p100_s000_0_con ; 
            double p100_p010_p010_s000_0_con = d110_s000_p010_s000_0_con + AB[1] * p100_s000_p010_s000_0_con ; 
            double p100_p010_p001_s000_0_con = d110_s000_p001_s000_0_con + AB[1] * p100_s000_p001_s000_0_con ; 
            double p100_p001_p100_s000_0_con = d101_s000_p100_s000_0_con + AB[2] * p100_s000_p100_s000_0_con ; 
            double p100_p001_p010_s000_0_con = d101_s000_p010_s000_0_con + AB[2] * p100_s000_p010_s000_0_con ; 
            double p100_p001_p001_s000_0_con = d101_s000_p001_s000_0_con + AB[2] * p100_s000_p001_s000_0_con ; 
            double p010_p100_p100_s000_0_con = d110_s000_p100_s000_0_con + AB[0] * p010_s000_p100_s000_0_con ; 
            double p010_p100_p010_s000_0_con = d110_s000_p010_s000_0_con + AB[0] * p010_s000_p010_s000_0_con ; 
            double p010_p100_p001_s000_0_con = d110_s000_p001_s000_0_con + AB[0] * p010_s000_p001_s000_0_con ; 
            double p010_p010_p100_s000_0_con = d020_s000_p100_s000_0_con + AB[1] * p010_s000_p100_s000_0_con ; 
            double p010_p010_p010_s000_0_con = d020_s000_p010_s000_0_con + AB[1] * p010_s000_p010_s000_0_con ; 
            double p010_p010_p001_s000_0_con = d020_s000_p001_s000_0_con + AB[1] * p010_s000_p001_s000_0_con ; 
            double p010_p001_p100_s000_0_con = d011_s000_p100_s000_0_con + AB[2] * p010_s000_p100_s000_0_con ; 
            double p010_p001_p010_s000_0_con = d011_s000_p010_s000_0_con + AB[2] * p010_s000_p010_s000_0_con ; 
            double p010_p001_p001_s000_0_con = d011_s000_p001_s000_0_con + AB[2] * p010_s000_p001_s000_0_con ; 
            double p001_p100_p100_s000_0_con = d101_s000_p100_s000_0_con + AB[0] * p001_s000_p100_s000_0_con ; 
            double p001_p100_p010_s000_0_con = d101_s000_p010_s000_0_con + AB[0] * p001_s000_p010_s000_0_con ; 
            double p001_p100_p001_s000_0_con = d101_s000_p001_s000_0_con + AB[0] * p001_s000_p001_s000_0_con ; 
            double p001_p010_p100_s000_0_con = d011_s000_p100_s000_0_con + AB[1] * p001_s000_p100_s000_0_con ; 
            double p001_p010_p010_s000_0_con = d011_s000_p010_s000_0_con + AB[1] * p001_s000_p010_s000_0_con ; 
            double p001_p010_p001_s000_0_con = d011_s000_p001_s000_0_con + AB[1] * p001_s000_p001_s000_0_con ; 
            double p001_p001_p100_s000_0_con = d002_s000_p100_s000_0_con + AB[2] * p001_s000_p100_s000_0_con ; 
            double p001_p001_p010_s000_0_con = d002_s000_p010_s000_0_con + AB[2] * p001_s000_p010_s000_0_con ; 
            double p001_p001_p001_s000_0_con = d002_s000_p001_s000_0_con + AB[2] * p001_s000_p001_s000_0_con ; 
            double p100_p100_p100_p100_0_con = p100_p100_d200_s000_0_con + CD[0] * p100_p100_p100_s000_0_con ; 
            double p100_p100_p100_p010_0_con = p100_p100_d110_s000_0_con + CD[1] * p100_p100_p100_s000_0_con ; 
            double p100_p100_p100_p001_0_con = p100_p100_d101_s000_0_con + CD[2] * p100_p100_p100_s000_0_con ; 
            double p100_p100_p010_p100_0_con = p100_p100_d110_s000_0_con + CD[0] * p100_p100_p010_s000_0_con ; 
            double p100_p100_p010_p010_0_con = p100_p100_d020_s000_0_con + CD[1] * p100_p100_p010_s000_0_con ; 
            double p100_p100_p010_p001_0_con = p100_p100_d011_s000_0_con + CD[2] * p100_p100_p010_s000_0_con ; 
            double p100_p100_p001_p100_0_con = p100_p100_d101_s000_0_con + CD[0] * p100_p100_p001_s000_0_con ; 
            double p100_p100_p001_p010_0_con = p100_p100_d011_s000_0_con + CD[1] * p100_p100_p001_s000_0_con ; 
            double p100_p100_p001_p001_0_con = p100_p100_d002_s000_0_con + CD[2] * p100_p100_p001_s000_0_con ; 
            double p100_p010_p100_p100_0_con = p100_p010_d200_s000_0_con + CD[0] * p100_p010_p100_s000_0_con ; 
            double p100_p010_p100_p010_0_con = p100_p010_d110_s000_0_con + CD[1] * p100_p010_p100_s000_0_con ; 
            double p100_p010_p100_p001_0_con = p100_p010_d101_s000_0_con + CD[2] * p100_p010_p100_s000_0_con ; 
            double p100_p010_p010_p100_0_con = p100_p010_d110_s000_0_con + CD[0] * p100_p010_p010_s000_0_con ; 
            double p100_p010_p010_p010_0_con = p100_p010_d020_s000_0_con + CD[1] * p100_p010_p010_s000_0_con ; 
            double p100_p010_p010_p001_0_con = p100_p010_d011_s000_0_con + CD[2] * p100_p010_p010_s000_0_con ; 
            double p100_p010_p001_p100_0_con = p100_p010_d101_s000_0_con + CD[0] * p100_p010_p001_s000_0_con ; 
            double p100_p010_p001_p010_0_con = p100_p010_d011_s000_0_con + CD[1] * p100_p010_p001_s000_0_con ; 
            double p100_p010_p001_p001_0_con = p100_p010_d002_s000_0_con + CD[2] * p100_p010_p001_s000_0_con ; 
            double p100_p001_p100_p100_0_con = p100_p001_d200_s000_0_con + CD[0] * p100_p001_p100_s000_0_con ; 
            double p100_p001_p100_p010_0_con = p100_p001_d110_s000_0_con + CD[1] * p100_p001_p100_s000_0_con ; 
            double p100_p001_p100_p001_0_con = p100_p001_d101_s000_0_con + CD[2] * p100_p001_p100_s000_0_con ; 
            double p100_p001_p010_p100_0_con = p100_p001_d110_s000_0_con + CD[0] * p100_p001_p010_s000_0_con ; 
            double p100_p001_p010_p010_0_con = p100_p001_d020_s000_0_con + CD[1] * p100_p001_p010_s000_0_con ; 
            double p100_p001_p010_p001_0_con = p100_p001_d011_s000_0_con + CD[2] * p100_p001_p010_s000_0_con ; 
            double p100_p001_p001_p100_0_con = p100_p001_d101_s000_0_con + CD[0] * p100_p001_p001_s000_0_con ; 
            double p100_p001_p001_p010_0_con = p100_p001_d011_s000_0_con + CD[1] * p100_p001_p001_s000_0_con ; 
            double p100_p001_p001_p001_0_con = p100_p001_d002_s000_0_con + CD[2] * p100_p001_p001_s000_0_con ; 
            double p010_p100_p100_p100_0_con = p010_p100_d200_s000_0_con + CD[0] * p010_p100_p100_s000_0_con ; 
            double p010_p100_p100_p010_0_con = p010_p100_d110_s000_0_con + CD[1] * p010_p100_p100_s000_0_con ; 
            double p010_p100_p100_p001_0_con = p010_p100_d101_s000_0_con + CD[2] * p010_p100_p100_s000_0_con ; 
            double p010_p100_p010_p100_0_con = p010_p100_d110_s000_0_con + CD[0] * p010_p100_p010_s000_0_con ; 
            double p010_p100_p010_p010_0_con = p010_p100_d020_s000_0_con + CD[1] * p010_p100_p010_s000_0_con ; 
            double p010_p100_p010_p001_0_con = p010_p100_d011_s000_0_con + CD[2] * p010_p100_p010_s000_0_con ; 
            double p010_p100_p001_p100_0_con = p010_p100_d101_s000_0_con + CD[0] * p010_p100_p001_s000_0_con ; 
            double p010_p100_p001_p010_0_con = p010_p100_d011_s000_0_con + CD[1] * p010_p100_p001_s000_0_con ; 
            double p010_p100_p001_p001_0_con = p010_p100_d002_s000_0_con + CD[2] * p010_p100_p001_s000_0_con ; 
            double p010_p010_p100_p100_0_con = p010_p010_d200_s000_0_con + CD[0] * p010_p010_p100_s000_0_con ; 
            double p010_p010_p100_p010_0_con = p010_p010_d110_s000_0_con + CD[1] * p010_p010_p100_s000_0_con ; 
            double p010_p010_p100_p001_0_con = p010_p010_d101_s000_0_con + CD[2] * p010_p010_p100_s000_0_con ; 
            double p010_p010_p010_p100_0_con = p010_p010_d110_s000_0_con + CD[0] * p010_p010_p010_s000_0_con ; 
            double p010_p010_p010_p010_0_con = p010_p010_d020_s000_0_con + CD[1] * p010_p010_p010_s000_0_con ; 
            double p010_p010_p010_p001_0_con = p010_p010_d011_s000_0_con + CD[2] * p010_p010_p010_s000_0_con ; 
            double p010_p010_p001_p100_0_con = p010_p010_d101_s000_0_con + CD[0] * p010_p010_p001_s000_0_con ; 
            double p010_p010_p001_p010_0_con = p010_p010_d011_s000_0_con + CD[1] * p010_p010_p001_s000_0_con ; 
            double p010_p010_p001_p001_0_con = p010_p010_d002_s000_0_con + CD[2] * p010_p010_p001_s000_0_con ; 
            double p010_p001_p100_p100_0_con = p010_p001_d200_s000_0_con + CD[0] * p010_p001_p100_s000_0_con ; 
            double p010_p001_p100_p010_0_con = p010_p001_d110_s000_0_con + CD[1] * p010_p001_p100_s000_0_con ; 
            double p010_p001_p100_p001_0_con = p010_p001_d101_s000_0_con + CD[2] * p010_p001_p100_s000_0_con ; 
            double p010_p001_p010_p100_0_con = p010_p001_d110_s000_0_con + CD[0] * p010_p001_p010_s000_0_con ; 
            double p010_p001_p010_p010_0_con = p010_p001_d020_s000_0_con + CD[1] * p010_p001_p010_s000_0_con ; 
            double p010_p001_p010_p001_0_con = p010_p001_d011_s000_0_con + CD[2] * p010_p001_p010_s000_0_con ; 
            double p010_p001_p001_p100_0_con = p010_p001_d101_s000_0_con + CD[0] * p010_p001_p001_s000_0_con ; 
            double p010_p001_p001_p010_0_con = p010_p001_d011_s000_0_con + CD[1] * p010_p001_p001_s000_0_con ; 
            double p010_p001_p001_p001_0_con = p010_p001_d002_s000_0_con + CD[2] * p010_p001_p001_s000_0_con ; 
            double p001_p100_p100_p100_0_con = p001_p100_d200_s000_0_con + CD[0] * p001_p100_p100_s000_0_con ; 
            double p001_p100_p100_p010_0_con = p001_p100_d110_s000_0_con + CD[1] * p001_p100_p100_s000_0_con ; 
            double p001_p100_p100_p001_0_con = p001_p100_d101_s000_0_con + CD[2] * p001_p100_p100_s000_0_con ; 
            double p001_p100_p010_p100_0_con = p001_p100_d110_s000_0_con + CD[0] * p001_p100_p010_s000_0_con ; 
            double p001_p100_p010_p010_0_con = p001_p100_d020_s000_0_con + CD[1] * p001_p100_p010_s000_0_con ; 
            double p001_p100_p010_p001_0_con = p001_p100_d011_s000_0_con + CD[2] * p001_p100_p010_s000_0_con ; 
            double p001_p100_p001_p100_0_con = p001_p100_d101_s000_0_con + CD[0] * p001_p100_p001_s000_0_con ; 
            double p001_p100_p001_p010_0_con = p001_p100_d011_s000_0_con + CD[1] * p001_p100_p001_s000_0_con ; 
            double p001_p100_p001_p001_0_con = p001_p100_d002_s000_0_con + CD[2] * p001_p100_p001_s000_0_con ; 
            double p001_p010_p100_p100_0_con = p001_p010_d200_s000_0_con + CD[0] * p001_p010_p100_s000_0_con ; 
            double p001_p010_p100_p010_0_con = p001_p010_d110_s000_0_con + CD[1] * p001_p010_p100_s000_0_con ; 
            double p001_p010_p100_p001_0_con = p001_p010_d101_s000_0_con + CD[2] * p001_p010_p100_s000_0_con ; 
            double p001_p010_p010_p100_0_con = p001_p010_d110_s000_0_con + CD[0] * p001_p010_p010_s000_0_con ; 
            double p001_p010_p010_p010_0_con = p001_p010_d020_s000_0_con + CD[1] * p001_p010_p010_s000_0_con ; 
            double p001_p010_p010_p001_0_con = p001_p010_d011_s000_0_con + CD[2] * p001_p010_p010_s000_0_con ; 
            double p001_p010_p001_p100_0_con = p001_p010_d101_s000_0_con + CD[0] * p001_p010_p001_s000_0_con ; 
            double p001_p010_p001_p010_0_con = p001_p010_d011_s000_0_con + CD[1] * p001_p010_p001_s000_0_con ; 
            double p001_p010_p001_p001_0_con = p001_p010_d002_s000_0_con + CD[2] * p001_p010_p001_s000_0_con ; 
            double p001_p001_p100_p100_0_con = p001_p001_d200_s000_0_con + CD[0] * p001_p001_p100_s000_0_con ; 
            double p001_p001_p100_p010_0_con = p001_p001_d110_s000_0_con + CD[1] * p001_p001_p100_s000_0_con ; 
            double p001_p001_p100_p001_0_con = p001_p001_d101_s000_0_con + CD[2] * p001_p001_p100_s000_0_con ; 
            double p001_p001_p010_p100_0_con = p001_p001_d110_s000_0_con + CD[0] * p001_p001_p010_s000_0_con ; 
            double p001_p001_p010_p010_0_con = p001_p001_d020_s000_0_con + CD[1] * p001_p001_p010_s000_0_con ; 
            double p001_p001_p010_p001_0_con = p001_p001_d011_s000_0_con + CD[2] * p001_p001_p010_s000_0_con ; 
            double p001_p001_p001_p100_0_con = p001_p001_d101_s000_0_con + CD[0] * p001_p001_p001_s000_0_con ; 
            double p001_p001_p001_p010_0_con = p001_p001_d011_s000_0_con + CD[1] * p001_p001_p001_s000_0_con ; 
            double p001_p001_p001_p001_0_con = p001_p001_d002_s000_0_con + CD[2] * p001_p001_p001_s000_0_con ; 
            I_[idx + (0 + 0 * 3 + 0 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p100_p100_p100_p100_0_con ;
            I_[idx + (1 + 0 * 3 + 0 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p010_p100_p100_p100_0_con ;
            I_[idx + (2 + 0 * 3 + 0 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p001_p100_p100_p100_0_con ;
            I_[idx + (0 + 1 * 3 + 0 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p100_p010_p100_p100_0_con ;
            I_[idx + (1 + 1 * 3 + 0 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p010_p010_p100_p100_0_con ;
            I_[idx + (2 + 1 * 3 + 0 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p001_p010_p100_p100_0_con ;
            I_[idx + (0 + 2 * 3 + 0 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p100_p001_p100_p100_0_con ;
            I_[idx + (1 + 2 * 3 + 0 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p010_p001_p100_p100_0_con ;
            I_[idx + (2 + 2 * 3 + 0 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p001_p001_p100_p100_0_con ;
            I_[idx + (0 + 0 * 3 + 1 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p100_p100_p010_p100_0_con ;
            I_[idx + (1 + 0 * 3 + 1 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p010_p100_p010_p100_0_con ;
            I_[idx + (2 + 0 * 3 + 1 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p001_p100_p010_p100_0_con ;
            I_[idx + (0 + 1 * 3 + 1 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p100_p010_p010_p100_0_con ;
            I_[idx + (1 + 1 * 3 + 1 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p010_p010_p010_p100_0_con ;
            I_[idx + (2 + 1 * 3 + 1 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p001_p010_p010_p100_0_con ;
            I_[idx + (0 + 2 * 3 + 1 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p100_p001_p010_p100_0_con ;
            I_[idx + (1 + 2 * 3 + 1 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p010_p001_p010_p100_0_con ;
            I_[idx + (2 + 2 * 3 + 1 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p001_p001_p010_p100_0_con ;
            I_[idx + (0 + 0 * 3 + 2 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p100_p100_p001_p100_0_con ;
            I_[idx + (1 + 0 * 3 + 2 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p010_p100_p001_p100_0_con ;
            I_[idx + (2 + 0 * 3 + 2 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p001_p100_p001_p100_0_con ;
            I_[idx + (0 + 1 * 3 + 2 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p100_p010_p001_p100_0_con ;
            I_[idx + (1 + 1 * 3 + 2 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p010_p010_p001_p100_0_con ;
            I_[idx + (2 + 1 * 3 + 2 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p001_p010_p001_p100_0_con ;
            I_[idx + (0 + 2 * 3 + 2 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p100_p001_p001_p100_0_con ;
            I_[idx + (1 + 2 * 3 + 2 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p010_p001_p001_p100_0_con ;
            I_[idx + (2 + 2 * 3 + 2 * 3 * 3 + 0 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p001_p001_p001_p100_0_con ;
            I_[idx + (0 + 0 * 3 + 0 * 3 * 3 + 1 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p100_p100_p100_p010_0_con ;
            I_[idx + (1 + 0 * 3 + 0 * 3 * 3 + 1 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p010_p100_p100_p010_0_con ;
            I_[idx + (2 + 0 * 3 + 0 * 3 * 3 + 1 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p001_p100_p100_p010_0_con ;
            I_[idx + (0 + 1 * 3 + 0 * 3 * 3 + 1 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p100_p010_p100_p010_0_con ;
            I_[idx + (1 + 1 * 3 + 0 * 3 * 3 + 1 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p010_p010_p100_p010_0_con ;
            I_[idx + (2 + 1 * 3 + 0 * 3 * 3 + 1 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p001_p010_p100_p010_0_con ;
            I_[idx + (0 + 2 * 3 + 0 * 3 * 3 + 1 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p100_p001_p100_p010_0_con ;
            I_[idx + (1 + 2 * 3 + 0 * 3 * 3 + 1 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p010_p001_p100_p010_0_con ;
            I_[idx + (2 + 2 * 3 + 0 * 3 * 3 + 1 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p001_p001_p100_p010_0_con ;
            I_[idx + (0 + 0 * 3 + 1 * 3 * 3 + 1 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p100_p100_p010_p010_0_con ;
            I_[idx + (1 + 0 * 3 + 1 * 3 * 3 + 1 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p010_p100_p010_p010_0_con ;
            I_[idx + (2 + 0 * 3 + 1 * 3 * 3 + 1 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p001_p100_p010_p010_0_con ;
            I_[idx + (0 + 1 * 3 + 1 * 3 * 3 + 1 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p100_p010_p010_p010_0_con ;
            I_[idx + (1 + 1 * 3 + 1 * 3 * 3 + 1 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p010_p010_p010_p010_0_con ;
            I_[idx + (2 + 1 * 3 + 1 * 3 * 3 + 1 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p001_p010_p010_p010_0_con ;
            I_[idx + (0 + 2 * 3 + 1 * 3 * 3 + 1 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p100_p001_p010_p010_0_con ;
            I_[idx + (1 + 2 * 3 + 1 * 3 * 3 + 1 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p010_p001_p010_p010_0_con ;
            I_[idx + (2 + 2 * 3 + 1 * 3 * 3 + 1 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p001_p001_p010_p010_0_con ;
            I_[idx + (0 + 0 * 3 + 2 * 3 * 3 + 1 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p100_p100_p001_p010_0_con ;
            I_[idx + (1 + 0 * 3 + 2 * 3 * 3 + 1 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p010_p100_p001_p010_0_con ;
            I_[idx + (2 + 0 * 3 + 2 * 3 * 3 + 1 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p001_p100_p001_p010_0_con ;
            I_[idx + (0 + 1 * 3 + 2 * 3 * 3 + 1 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p100_p010_p001_p010_0_con ;
            I_[idx + (1 + 1 * 3 + 2 * 3 * 3 + 1 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p010_p010_p001_p010_0_con ;
            I_[idx + (2 + 1 * 3 + 2 * 3 * 3 + 1 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p001_p010_p001_p010_0_con ;
            I_[idx + (0 + 2 * 3 + 2 * 3 * 3 + 1 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p100_p001_p001_p010_0_con ;
            I_[idx + (1 + 2 * 3 + 2 * 3 * 3 + 1 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p010_p001_p001_p010_0_con ;
            I_[idx + (2 + 2 * 3 + 2 * 3 * 3 + 1 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p001_p001_p001_p010_0_con ;
            I_[idx + (0 + 0 * 3 + 0 * 3 * 3 + 2 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p100_p100_p100_p001_0_con ;
            I_[idx + (1 + 0 * 3 + 0 * 3 * 3 + 2 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p010_p100_p100_p001_0_con ;
            I_[idx + (2 + 0 * 3 + 0 * 3 * 3 + 2 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p001_p100_p100_p001_0_con ;
            I_[idx + (0 + 1 * 3 + 0 * 3 * 3 + 2 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p100_p010_p100_p001_0_con ;
            I_[idx + (1 + 1 * 3 + 0 * 3 * 3 + 2 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p010_p010_p100_p001_0_con ;
            I_[idx + (2 + 1 * 3 + 0 * 3 * 3 + 2 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p001_p010_p100_p001_0_con ;
            I_[idx + (0 + 2 * 3 + 0 * 3 * 3 + 2 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p100_p001_p100_p001_0_con ;
            I_[idx + (1 + 2 * 3 + 0 * 3 * 3 + 2 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p010_p001_p100_p001_0_con ;
            I_[idx + (2 + 2 * 3 + 0 * 3 * 3 + 2 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p001_p001_p100_p001_0_con ;
            I_[idx + (0 + 0 * 3 + 1 * 3 * 3 + 2 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p100_p100_p010_p001_0_con ;
            I_[idx + (1 + 0 * 3 + 1 * 3 * 3 + 2 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p010_p100_p010_p001_0_con ;
            I_[idx + (2 + 0 * 3 + 1 * 3 * 3 + 2 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p001_p100_p010_p001_0_con ;
            I_[idx + (0 + 1 * 3 + 1 * 3 * 3 + 2 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p100_p010_p010_p001_0_con ;
            I_[idx + (1 + 1 * 3 + 1 * 3 * 3 + 2 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p010_p010_p010_p001_0_con ;
            I_[idx + (2 + 1 * 3 + 1 * 3 * 3 + 2 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p001_p010_p010_p001_0_con ;
            I_[idx + (0 + 2 * 3 + 1 * 3 * 3 + 2 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p100_p001_p010_p001_0_con ;
            I_[idx + (1 + 2 * 3 + 1 * 3 * 3 + 2 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p010_p001_p010_p001_0_con ;
            I_[idx + (2 + 2 * 3 + 1 * 3 * 3 + 2 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p001_p001_p010_p001_0_con ;
            I_[idx + (0 + 0 * 3 + 2 * 3 * 3 + 2 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p100_p100_p001_p001_0_con ;
            I_[idx + (1 + 0 * 3 + 2 * 3 * 3 + 2 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p010_p100_p001_p001_0_con ;
            I_[idx + (2 + 0 * 3 + 2 * 3 * 3 + 2 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p001_p100_p001_p001_0_con ;
            I_[idx + (0 + 1 * 3 + 2 * 3 * 3 + 2 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p100_p010_p001_p001_0_con ;
            I_[idx + (1 + 1 * 3 + 2 * 3 * 3 + 2 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p010_p010_p001_p001_0_con ;
            I_[idx + (2 + 1 * 3 + 2 * 3 * 3 + 2 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p001_p010_p001_p001_0_con ;
            I_[idx + (0 + 2 * 3 + 2 * 3 * 3 + 2 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p100_p001_p001_p001_0_con ;
            I_[idx + (1 + 2 * 3 + 2 * 3 * 3 + 2 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p010_p001_p001_p001_0_con ;
            I_[idx + (2 + 2 * 3 + 2 * 3 * 3 + 2 * 3 * 3 * 3) * nab + idy * nab * 81] =                 1 * p001_p001_p001_p001_0_con ;
                
        }
    } 
}






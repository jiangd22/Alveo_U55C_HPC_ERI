#include "ssOrb.h"
#include <cmath>
#include <array>

void hrr_ssss(  int AB_la, int AB_lb, int AB_na, int AB_nb, int AB_ma, int AB_mb, std::array<double, 6> AB_Z, std::array<double, 6> AB_ZA, std::array<double, 6> AB_K, std::array<double, 6> AB_S, std::array<std::array<int,2>,2> AB_idx, std::array<std::array<double,6>,3> AB_P, std::array<std::array<double,6>,3> AB_PA, std::array<std::array<double,6>,3> AB_AB,
                int CD_la, int CD_lb, int CD_na, int CD_nb, int CD_ma, int CD_mb, std::array<double, 6> CD_Z, std::array<double, 6> CD_ZA, std::array<double, 6> CD_K, std::array<double, 6> CD_S, std::array<std::array<int,2>,2> CD_idx, std::array<std::array<double,6>,3> CD_P, std::array<std::array<double,6>,3> CD_PA, std::array<std::array<double,6>,3> CD_AB,
                double* I_ ) {

    #pragma HLS INTERFACE m_axi port=I_ offset=slave bundle=gmem
    #pragma HLS INTERFACE s_axilite port=AB_la bundle=control
    #pragma HLS INTERFACE s_axilite port=AB_lb bundle=control
    #pragma HLS INTERFACE s_axilite port=AB_na bundle=control
    #pragma HLS INTERFACE s_axilite port=AB_nb bundle=control
    #pragma HLS INTERFACE s_axilite port=AB_ma bundle=control
    #pragma HLS INTERFACE s_axilite port=AB_mb bundle=control
    #pragma HLS INTERFACE s_axilite port=AB_Z bundle=control
    #pragma HLS INTERFACE s_axilite port=AB_ZA bundle=control
    #pragma HLS INTERFACE s_axilite port=AB_K bundle=control
    #pragma HLS INTERFACE s_axilite port=AB_S bundle=control
    #pragma HLS INTERFACE s_axilite port=AB_idx bundle=control
    #pragma HLS INTERFACE s_axilite port=AB_P bundle=control
    #pragma HLS INTERFACE s_axilite port=AB_PA bundle=control
    #pragma HLS INTERFACE s_axilite port=AB_AB bundle=control

    #pragma HLS INTERFACE s_axilite port=CD_la bundle=control
    #pragma HLS INTERFACE s_axilite port=CD_lb bundle=control
    #pragma HLS INTERFACE s_axilite port=CD_na bundle=control
    #pragma HLS INTERFACE s_axilite port=CD_nb bundle=control
    #pragma HLS INTERFACE s_axilite port=CD_ma bundle=control
    #pragma HLS INTERFACE s_axilite port=CD_mb bundle=control
    #pragma HLS INTERFACE s_axilite port=CD_Z bundle=control
    #pragma HLS INTERFACE s_axilite port=CD_ZA bundle=control
    #pragma HLS INTERFACE s_axilite port=CD_K bundle=control
    #pragma HLS INTERFACE s_axilite port=CD_S bundle=control
    #pragma HLS INTERFACE s_axilite port=CD_idx bundle=control
    #pragma HLS INTERFACE s_axilite port=CD_P bundle=control
    #pragma HLS INTERFACE s_axilite port=CD_PA bundle=control
    #pragma HLS INTERFACE s_axilite port=CD_AB bundle=control

    #pragma HLS INTERFACE s_axilite port=return bundle=control

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



// void hrr_ssss(shell_pairs ab, shell_pairs cd,  double* I_ )
// {
//     double idx = threadIdx.x + blockDim.x * blockIdx.x;
//     double idy = threadIdx.y + blockDim.y * blockIdx.y;


//     double nab = ab.nab;
//     double ncd = cd.nab;

//     double mab = ab.ma * ab.mb; 
//     double mcd = cd.ma * cd.mb; 

//     if (idx < nab  and idy < ncd)
//     {
//     double AB[3] = 
//     {
//         ab.AB[0][idx],
//         ab.AB[1][idx],
//         ab.AB[2][idx]
//     };
//     double CD[3] = 
//     {
//         cd.AB[0][idy],
//         cd.AB[1][idy],
//         cd.AB[2][idy]
//     };

//     double s000_s000_s000_s000_0_con = 0.0;

//     for (int pab = 0; pab < mab; pab++)
//     for (int pcd = 0; pcd < mcd; pcd++)
//     {
//         double iab = idx + pab * nab; 
//         double icd = idy + pcd * ncd; 

//         double zab = ab.Z[iab];
//         double zcd = cd.Z[icd];
                              
//         double zab_inv = 1 / zab;
//         double zcd_inv = 1 / zcd;
//         double zabcd_inv_sqrt = 1.0/(sqrt(zab + zcd));
//         double zabcd_inv = zabcd_inv_sqrt * zabcd_inv_sqrt;
//         double rho = zab * zcd * zabcd_inv;
                                         
//         double P[3] = 
//         {
//             ab.P[0][iab],
//             ab.P[1][iab],
//             ab.P[2][iab],
//         };
        
//         double Q[3] = 
//         {
//             cd.P[0][icd],
//             cd.P[1][icd],
//             cd.P[2][icd],
//         };
        
//         double T = rho * ((P[0]-Q[0])*(P[0]-Q[0]) + (P[1]-Q[1])*(P[1]-Q[1]) + (P[2]-Q[2])*(P[2]-Q[2]));
//         double K = zabcd_inv_sqrt * (ab.K[iab]) * cd.K[icd]; 
        
//         double fm[1]; 
//         // vgamma<double>(0, T, fm);
        
//         for (auto i = 0; i <1 ; i++) 
//             fm[i] *= K;

//         auto s000_s000_s000_s000_0 = fm[0] ;
//         s000_s000_s000_s000_0_con += s000_s000_s000_s000_0;
//         }
//          //****************//
//     I_[idx + (0 + 0 * 1 + 0 * 1 * 1 + 0 * 1 * 1 * 1) * nab + idy * nab * 1] =                 1 * s000_s000_s000_s000_0_con ;
//     }
// }

// void hrr_psss(shell_pairs ab, shell_pairs cd,  double* I_ )
// {
//     double idx = threadIdx.x + blockDim.x * blockIdx.x; // idx for shell_pairs ab 
//     double idy = threadIdx.y + blockDim.y * blockIdx.y; // idy for shell_pairs cd 

//     double nab = ab.nab;
//     double ncd = cd.nab;

//     double mab = ab.ma * ab.mb; 
//     double mcd = cd.ma * cd.mb; 

//     if (idx < nab  and idy < ncd)
//     {
//     double AB[3] = 
//     {
//         ab.AB[0][idx],
//         ab.AB[1][idx],
//         ab.AB[2][idx]
//     };
//     double CD[3] = 
//     {
//         cd.AB[0][idy],
//         cd.AB[1][idy],
//         cd.AB[2][idy]
//     };

//     double p100_s000_s000_s000_0_con = 0.0;
//     double p010_s000_s000_s000_0_con = 0.0;
//     double p001_s000_s000_s000_0_con = 0.0;

//     for (int pab = 0; pab < mab; pab++)
//     for (int pcd = 0; pcd < mcd; pcd++)
//     {
//         double iab = idx + pab * nab; 
//         double icd = idy + pcd * ncd; 
                                    
//         double zab = ab.Z[iab];
//         double zcd = cd.Z[icd];
                              
//         double zab_inv = 1 / zab;
//         double zcd_inv = 1 / zcd;
//         double zabcd_inv_sqrt = 1.0/sqrt(zab + zcd);
//         double zabcd_inv = zabcd_inv_sqrt * zabcd_inv_sqrt;
//         double rho = zab * zcd * zabcd_inv;
                                         
//         double P[3] = 
//         {
//             ab.P[0][iab],
//             ab.P[1][iab],
//             ab.P[2][iab],
//         };
        
//         double Q[3] = 
//         {
//             cd.P[0][icd],
//             cd.P[1][icd],
//             cd.P[2][icd],
//         };
        
//         double PA[3] =
//         {
//             ab.PA[0][iab],
//             ab.PA[1][iab],
//             ab.PA[2][iab],
//         };
        
//         double PB[3] =
//         {
//             ab.PA[0][iab] + ab.AB[0][iab],
//             ab.PA[1][iab] + ab.AB[1][iab],
//             ab.PA[2][iab] + ab.AB[2][iab],
//         };
        
//         double QC[3] = 
//         {
//             cd.PA[0][icd],
//             cd.PA[1][icd],
//             cd.PA[2][icd],
//         };
        
//         double QD[3] =
//         {
//             cd.PA[0][icd] + cd.AB[0][icd],
//             cd.PA[1][icd] + cd.AB[1][icd],
//             cd.PA[2][icd] + cd.AB[2][icd],
//         };
        
//         double W[3] =
//         {
//             (zab * P[0] + zcd * Q[0]) * zabcd_inv,
//             (zab * P[1] + zcd * Q[1]) * zabcd_inv,
//             (zab * P[2] + zcd * Q[2]) * zabcd_inv,
//         };
//         double WP[3] =
//         {
//             W[0] - P[0],
//             W[1] - P[1],
//             W[2] - P[2]
//         };
        
//         double WQ[3] =
//         {
//             W[0] - Q[0],
//             W[1] - Q[1],
//             W[2] - Q[2]
//         };
        
//         double T = rho * ((P[0]-Q[0])*(P[0]-Q[0]) + (P[1]-Q[1])*(P[1]-Q[1]) + (P[2]-Q[2])*(P[2]-Q[2]));
//         double K = zabcd_inv_sqrt * (ab.K[iab]) * cd.K[icd]; 
        
//         double fm[2]; 
//         // vgamma<double>(1, T, fm);
        
//         for (auto i = 0; i <2 ; i++) 
//             fm[i] *= K;

//         double s000_s000_s000_s000_0 = fm[0] ;
//         double s000_s000_s000_s000_1 = fm[1] ;
//         double p100_s000_s000_s000_0 = PA[0] * s000_s000_s000_s000_0 + WP[0] * s000_s000_s000_s000_1;
//         double p010_s000_s000_s000_0 = PA[1] * s000_s000_s000_s000_0 + WP[1] * s000_s000_s000_s000_1;
//         double p001_s000_s000_s000_0 = PA[2] * s000_s000_s000_s000_0 + WP[2] * s000_s000_s000_s000_1;
//         p100_s000_s000_s000_0_con += p100_s000_s000_s000_0;
//         p010_s000_s000_s000_0_con += p010_s000_s000_s000_0;
//         p001_s000_s000_s000_0_con += p001_s000_s000_s000_0;
//         }
//          //****************//
//     I_[idx + (0 + 0 * 3 + 0 * 3 * 1 + 0 * 3 * 1 * 1) * nab + idy * nab * 3] =                 1 * p100_s000_s000_s000_0_con ;
//     I_[idx + (1 + 0 * 3 + 0 * 3 * 1 + 0 * 3 * 1 * 1) * nab + idy * nab * 3] =                 1 * p010_s000_s000_s000_0_con ;
//     I_[idx + (2 + 0 * 3 + 0 * 3 * 1 + 0 * 3 * 1 * 1) * nab + idy * nab * 3] =                 1 * p001_s000_s000_s000_0_con ;
//     }
// }


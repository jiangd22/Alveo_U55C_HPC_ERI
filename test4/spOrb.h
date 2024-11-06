#include <cmath>
// #include "hls_math.h"
#include <vector>
#include <array>
// #include "shell_pairs.h"


void driver (int order, 
            int la, int lb, int na, int nb, int ma, int mb, std::array<double, 6> abZ, std::array<double, 6> abZA, std::array<double, 6> abK, std::array<double, 6> abS, std::array<std::array<int,2>,2> abidx, std::array<std::array<double,6>,3> abP, std::array<std::array<double,6>,3> abPA, std::array<std::array<double,6>,3> abAB, 
            int lc, int ld, int nc, int nd, int mc, int md, std::array<double, 6> cdZ, std::array<double, 6> cdZA, std::array<double, 6> cdK, std::array<double, 6> cdS, std::array<std::array<int,2>,2> cdidx, std::array<std::array<double,6>,3> cdP, std::array<std::array<double,6>,3> cdPA, std::array<std::array<double,6>,3> cdAB,
            double* I_ );

void hrr_ssss(  int AB_la, int AB_lb, int AB_na, int AB_nb, int AB_ma, int AB_mb, double AB_Z[6], double AB_ZA[6], double AB_K[6], double AB_S[6], int AB_idx[2][2], double AB_P[3][6], double AB_PA[3][6], double AB_AB[3][6],
                int CD_la, int CD_lb, int CD_na, int CD_nb, int CD_ma, int CD_mb, double CD_Z[6], double CD_ZA[6], double CD_K[6], double CD_S[6], int CD_idx[2][2], double CD_P[3][6], double CD_PA[3][6], double CD_AB[3][6],
                double* I_ );

void hrr_psss(  int AB_la, int AB_lb, int AB_na, int AB_nb, int AB_ma, int AB_mb, double AB_Z[6], double AB_ZA[6], double AB_K[6], double AB_S[6], int AB_idx[2][2], double AB_P[3][6], double AB_PA[3][6], double AB_AB[3][6],
                int CD_la, int CD_lb, int CD_na, int CD_nb, int CD_ma, int CD_mb, double CD_Z[6], double CD_ZA[6], double CD_K[6], double CD_S[6], int CD_idx[2][2], double CD_P[3][6], double CD_PA[3][6], double CD_AB[3][6],
                double* I_ );

void hrr_psps(  int AB_la, int AB_lb, int AB_na, int AB_nb, int AB_ma, int AB_mb, double AB_Z[6], double AB_ZA[6], double AB_K[6], double AB_S[6], int AB_idx[2][2], double AB_P[3][6], double AB_PA[3][6], double AB_AB[3][6],
                int CD_la, int CD_lb, int CD_na, int CD_nb, int CD_ma, int CD_mb, double CD_Z[6], double CD_ZA[6], double CD_K[6], double CD_S[6], int CD_idx[2][2], double CD_P[3][6], double CD_PA[3][6], double CD_AB[3][6],
                double* I_ );


void hrr_ppss(  int AB_la, int AB_lb, int AB_na, int AB_nb, int AB_ma, int AB_mb, double AB_Z[6], double AB_ZA[6], double AB_K[6], double AB_S[6], int AB_idx[2][2], double AB_P[3][6], double AB_PA[3][6], double AB_AB[3][6],
                int CD_la, int CD_lb, int CD_na, int CD_nb, int CD_ma, int CD_mb, double CD_Z[6], double CD_ZA[6], double CD_K[6], double CD_S[6], int CD_idx[2][2], double CD_P[3][6], double CD_PA[3][6], double CD_AB[3][6],
                double* I_ );

void hrr_ppps(  int AB_la, int AB_lb, int AB_na, int AB_nb, int AB_ma, int AB_mb, double AB_Z[6], double AB_ZA[6], double AB_K[6], double AB_S[6], int AB_idx[2][2], double AB_P[3][6], double AB_PA[3][6], double AB_AB[3][6],
                int CD_la, int CD_lb, int CD_na, int CD_nb, int CD_ma, int CD_mb, double CD_Z[6], double CD_ZA[6], double CD_K[6], double CD_S[6], int CD_idx[2][2], double CD_P[3][6], double CD_PA[3][6], double CD_AB[3][6],
                double* I_ );

void hrr_pppp(  int AB_la, int AB_lb, int AB_na, int AB_nb, int AB_ma, int AB_mb, double AB_Z[6], double AB_ZA[6], double AB_K[6], double AB_S[6], int AB_idx[2][2], double AB_P[3][6], double AB_PA[3][6], double AB_AB[3][6],
                int CD_la, int CD_lb, int CD_na, int CD_nb, int CD_ma, int CD_mb, double CD_Z[6], double CD_ZA[6], double CD_K[6], double CD_S[6], int CD_idx[2][2], double CD_P[3][6], double CD_PA[3][6], double CD_AB[3][6],
                double* I_ );  

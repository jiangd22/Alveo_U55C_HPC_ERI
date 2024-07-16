#include <cmath>
// #include "hls_math.h"
#include <vector>
#include <array>
// #include "shell_pairs.h"


// void hrr_ssss(shell_pairs ab, shell_pairs cd,  double* I_ );
void hrr_ssss(  int AB_la,
                int AB_lb,
                int AB_na,
                int AB_nb, 
                int AB_ma, 
                int AB_mb, 
                std::array<double, 6> AB_Z, 
                std::array<double, 6> AB_ZA, 
                std::array<double, 6> AB_K, 
                std::array<double, 6> AB_S, 
                std::array<std::array<int,2>,2> AB_idx, 
                std::array<std::array<double,6>,3> AB_P, 
                std::array<std::array<double,6>,3> AB_PA, 
                std::array<std::array<double,6>,3> AB_AB,

                int CD_la, 
                int CD_lb, 
                int CD_na, 
                int CD_nb, 
                int CD_ma, 
                int CD_mb, 
                std::array<double, 6> CD_Z, 
                std::array<double, 6> CD_ZA, 
                std::array<double, 6> CD_K, 
                std::array<double, 6> CD_S, 
                std::array<std::array<int,2>,2> CD_idx, 
                std::array<std::array<double,6>,3> CD_P, 
                std::array<std::array<double,6>,3> CD_PA, 
                std::array<std::array<double,6>,3> CD_AB,

                double* I_ );

// void hrr_psss(shell_pairs ab, shell_pairs cd,  double* I_ );
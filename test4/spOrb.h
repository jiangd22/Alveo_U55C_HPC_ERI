#include <cmath>
// #include "hls_math.h"
#include <vector>
#include <array>
// #include "shell_pairs.h"


// void hrr_ssss(shell_pairs ab, shell_pairs cd,  double* I_ );
void hrr_ssss(  int la,
                int lb,
                int na,
                int nb,
                int ma,
                int mb,
                std::array<double, 6> Z,
                std::array<double, 6> ZA,
                std::array<double, 6> K,
                std::array<double, 6> S,
                std::array<std::array<int,2>,2> idx,
                std::array<std::array<double,6>,3> P,
                std::array<std::array<double,6>,3> PA,
                std::array<std::array<double,6>,3> AB,
                double* I_ 
                );

void hrr_psss(  int la,
                int lb,
                int na,
                int nb,
                int ma,
                int mb,
                std::array<double, 6> Z,
                std::array<double, 6> ZA,
                std::array<double, 6> K,
                std::array<double, 6> S,
                std::array<std::array<int,2>,2> idx,
                std::array<std::array<double,6>,3> P,
                std::array<std::array<double,6>,3> PA,
                std::array<std::array<double,6>,3> AB,
                double* I_ 
                );

void hrr_psps(  int la,
                int lb,
                int na,
                int nb,
                int ma,
                int mb,
                std::array<double, 6> Z,
                std::array<double, 6> ZA,
                std::array<double, 6> K,
                std::array<double, 6> S,
                std::array<std::array<int,2>,2> idx,
                std::array<std::array<double,6>,3> P,
                std::array<std::array<double,6>,3> PA,
                std::array<std::array<double,6>,3> AB,
                double* I_ 
                );


void hrr_ppss(  int la,
                int lb,
                int na,
                int nb,
                int ma,
                int mb,
                std::array<double, 6> Z,
                std::array<double, 6> ZA,
                std::array<double, 6> K,
                std::array<double, 6> S,
                std::array<std::array<int,2>,2> idx,
                std::array<std::array<double,6>,3> P,
                std::array<std::array<double,6>,3> PA,
                std::array<std::array<double,6>,3> AB,
                double* I_ 
                );

void hrr_ppps(  int la,
                int lb,
                int na,
                int nb,
                int ma,
                int mb,
                std::array<double, 6> Z,
                std::array<double, 6> ZA,
                std::array<double, 6> K,
                std::array<double, 6> S,
                std::array<std::array<int,2>,2> idx,
                std::array<std::array<double,6>,3> P,
                std::array<std::array<double,6>,3> PA,
                std::array<std::array<double,6>,3> AB,
                double* I_ 
                );

void hrr_pppp(  int la,
                int lb,
                int na,
                int nb,
                int ma,
                int mb,
                std::array<double, 6> Z,
                std::array<double, 6> ZA,
                std::array<double, 6> K,
                std::array<double, 6> S,
                std::array<std::array<int,2>,2> idx,
                std::array<std::array<double,6>,3> P,
                std::array<std::array<double,6>,3> PA,
                std::array<std::array<double,6>,3> AB,
                double* I_ 
                );  

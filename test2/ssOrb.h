#include <cmath>
// #include "hls_math.h"
#include <vector>
#include <array>
// #include "shell_pairs.h"


// void hrr_ssss(shell_pairs ab, shell_pairs cd,  double* I_ );
void hrr_ssss(  int la, int lb, 
                int na, int nb, 
                int ma, int mb, 
                std::array<double, 6> abZ, 
                std::array<double, 6> abZA, 
                std::array<double, 6> abK, 
                std::array<double, 6> abS, 
                std::array<std::array<int,2>,2> abidx, 
                std::array<std::array<double,6>,3> abP, 
                std::array<std::array<double,6>,3> abPA, 
                std::array<std::array<double,6>,3> abAB,
                int lc, int ld, 
                int nc, int nd, 
                int mc, int md, 
                std::array<double, 6> cdZ, 
                std::array<double, 6> cdZA, 
                std::array<double, 6> cdK, 
                std::array<double, 6> cdS, 
                std::array<std::array<int,2>,2> cdidx, 
                std::array<std::array<double,6>,3> cdP, 
                std::array<std::array<double,6>,3> cdPA, 
                std::array<std::array<double,6>,3> cdAB,
                double* I_ ) ;

// void hrr_psss(shell_pairs ab, shell_pairs cd,  double* I_ );
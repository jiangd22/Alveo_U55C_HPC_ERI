#include "ssOrb.h"
// #include "shell_pairs.h"
#include <iostream>
#include <array>
#include <chrono>  // Include for timing

int main() {

    int la = 0;
    int lb = 0; // angular momentum
    int na = 2;
    int nb = 3; // number of shells
    int ma = 1;
    int mb = 1; // number of primitives
    std::array<double, 6> Z = {6.342, 6.29243, 13.0011, 6.33731, 6.94737, 6.4369}; // zeta_a + zeta_b
    std::array<double, 6> ZA = {6.342, 6.29243, 13.0011, 6.33731, 6.94737, 6.4369}; // zeta_a
    std::array<double, 6> K = {0.110144, 0.00938326, 0.00124135, 0.203043, 0.0265948, 0.0490335}; // kappa constant
    std::array<double,6> S = {1,}; // Schwarz factor sPrt[(ab|ab)]
    std::array<std::array<int,2>,2> idx = {{{ 0, }, { 0, }}}; // shell indices (a and b)
    std::array<std::array<double,6>,3> P = {{
        { 0.00829901, 0.00581651, 0.169694, 0.00806565, 0.0357558, 0.0129443, },
        { -0.0288298, -0.0202059, -0.589495, -0.0280191, -0.124211, -0.0449669, },
        { 0.00455825, 0.00319474, 0.0932046, 0.00443008, 0.019639, 0.00710968, },
    }} ; // zeta_a*A + zeta_b*B / (zeta_a + zeta_b)
    std::array<std::array<double, 6>,3> PA = {{
        { -0.315103, -0.317585, -0.153708, -0.315336, -0.287646, -0.310458, },
        { 1.09463, 1.10325, 0.533965, 1.09544, 0.999248, 1.07849, },
        { -0.173071, -0.174435, -0.0844248, -0.173199, -0.15799, -0.17052, },
    }};// P - A
    std::array<std::array<double, 6>,3> AB = {{
        { 0.323402, 0.323402, 0.323402, 0.323402, 0.323402, 0.323402, },
        { -1.12346, -1.12346, -1.12346, -1.12346, -1.12346, -1.12346, },
        { 0.177629, 0.177629, 0.177629, 0.177629, 0.177629, 0.177629, },
    }}; // A - B

    double I_ssss[1];
    
    // Start the timer
    auto start = std::chrono::high_resolution_clock::now();

    hrr_ssss(   la, lb, na, nb, ma, mb, Z, ZA, K, S, idx, P, PA, AB, 
                la, lb, na, nb, ma, mb, Z, ZA, K, S, idx, P, PA, AB,
                I_ssss);

    // Stop the timer
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate the elapsed time in microseconds (or other units if needed)
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    // Print the result and the timing
    std::cout << "Result: " << I_ssss[0] << std::endl;
    std::cout << "Elapsed Time: " << duration << " microseconds" << std::endl;

    return 0;
}

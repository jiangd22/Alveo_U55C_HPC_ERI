#ifndef _SHELL_PAIRS_H_
#define _SHELL_PAIRS_H_

// #include "molecule.h"
class molecule;

#include <vector>
#include <array>

struct shell_pairs
{
    int la, lb; // angular momentum
    int na, nb; // number of shells
    int ma, mb; // number of primitives
    std::vector<double> Z; // zeta_a + zeta_b
    std::vector<double> ZA; // zeta_a
    std::vector<double> K; // kappa constant
    std::vector<double> S; // Schwarz factor sPrt[(ab|ab)]
    std::array<std::vector<int>,2> idx; // shell indices (a and b)
    std::array<std::vector<double>,3> P; // zeta_a*A + zeta_b*B / (zeta_a + zeta_b)
    std::array<std::vector<double>,3> PA; // P - A
    std::array<std::vector<double>,3> AB; // A - B

};

std::vector<shell_pairs> shell_pairs_for_molecule(const molecule& mol, double tol = 0.0);

#endif

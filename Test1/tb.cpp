#include "expTest.h"
#include <iostream>

int main() {
    data_t input = 0;
    data_t cMath_output[N];
    data_t hlsMath_output[N];
    data_t time[2*N];
    expTest(input, cMath_output, hlsMath_output);

    for (int i = 0; i < N; i++) {
        std::cout << "cMath_output[" << i << "] = " << cMath_output[i] << std::endl;
        std::cout << "hlsMath_output[" << i << "] = " << hlsMath_output[i] << std::endl;
    }

    return 0;
}
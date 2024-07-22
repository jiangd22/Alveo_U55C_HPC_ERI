#include "expTest.h"

typedef int data_t;

data_t cpp_math(data_t imput) {
    data_t temp = exp(imput);
    return temp;
}

data_t cpp_math_hls(data_t imput) {
    data_t temp = hls::exp(imput);
    return temp;
}

void expTest(data_t input, data_t cMath_output[N], data_t hlsMath_output[N]) {
    for (int i = 0; i < N; i++) {
        cMath_output[i] = cpp_math(input);
        hlsMath_output[i] = cpp_math_hls(input);
        input++;
    }
}


    

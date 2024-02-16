#include <CL/sycl.hpp>
#include <iostream>
#include <vector>

namespace sycl = cl::sycl;

constexpr size_t N = 1024;

int main() {
    // Create two vectors
    std::vector<int> vecA(N, 1);
    std::vector<int> vecB(N, 2);
    std::vector<int> result(N);

    // Create a SYCL queue targeting FPGA device
    // sycl::queue q(sycl::default_selector{});

    // Submit a command group to the queue
    q.submit([&](sycl::handler& cgh) {
        // Accessors for input and output buffers
        auto accessorA = sycl::accessor(vecA, cgh, sycl::read_only);
        auto accessorB = sycl::accessor(vecB, cgh, sycl::read_only);
        auto accessorResult = sycl::accessor(result, cgh, sycl::write_only);

        // Kernel function to calculate sum of two vectors
        cgh.parallel_for<class vector_addition>(sycl::range<1>(N), [=](sycl::id<1> idx) {
            accessorResult[idx] = accessorA[idx] + accessorB[idx];
        });
    });

    // Wait for command group to finish and then print the result
    q.wait();
    
    std::cout << "Result:" << std::endl;
    for (int i = 0; i < N; ++i) {
        std::cout << result[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}

#include <chrono>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <memory>

#include "hilbert.hpp"

int main() {
  using Int = std::uint16_t;
  using UInt = std::size_t;

  constexpr std::size_t N = 2;
  constexpr std::size_t K = 14;

  constexpr std::size_t bytes = sizeof(Int) * N << N * K;

  auto time = [&](const char* desc, auto&& f) {
    auto start = std::chrono::high_resolution_clock::now();
    f();
    auto duration = std::chrono::high_resolution_clock::now() - start;
    double count = std::chrono::duration<double>(duration).count();
    double throughput = bytes / count / 1024 / 1024 / 1024;
    std::cout << desc << ":\t" << count << " s (" << throughput << " GiB/s)\n";
  };

  std::unique_ptr<Int[]> curve;
  time("new[]", [&]() { curve = std::make_unique<Int[]>(N << N * K); });
  // memset() is not necessary, but just gives a frame of reference
  // for how fast we can sequentially write to main memory.
  time("memset", [&]() { std::memset(curve.get(), 0, bytes); });
  time("Curve", [&]() { Hilbert<Int, UInt>::Curve<N, K>(curve.get()); });

  return 0;
}

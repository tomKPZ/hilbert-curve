#include <chrono>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <memory>

#include "hilbert.hpp"

int main() {
  using ITy = std::uint32_t;
  using ViTy = std::uint16_t;

  constexpr uint N = 2;
  constexpr uint K = 14;

  constexpr uint bytes = sizeof(ViTy) * N << N * K;
  auto time = [&](const char* desc, auto&& f) {
    auto start = std::chrono::high_resolution_clock::now();
    f();
    auto duration = std::chrono::high_resolution_clock::now() - start;
    double count = std::chrono::duration<double>(duration).count();
    double throughput = bytes / count / 1024 / 1024 / 1024;
    std::cout << desc << ":\t" << count << " s (" << throughput << " GiB/s)\n";
  };

  std::unique_ptr<ViTy[]> vs;
  time("new[]", [&]() { vs = std::make_unique<ViTy[]>(N << N * K); });
  // memset() is not necessary, but just gives a frame of reference
  // for how fast we can sequentially write to main memory.
  time("memset", [&]() { std::memset(vs.get(), 0, bytes); });
  time("IsToVs", [&]() { Hilbert<N, K, ViTy, ITy>::IsToVs(vs.get()); });
  std::unique_ptr<ITy[]> is;
  time("new[]", [&]() { is = std::make_unique<ITy[]>(1 << N * K); });
  time("memset", [&]() { std::memset(vs.get(), 0, bytes); });
  time("VsToIs", [&]() { Hilbert<N, K, ViTy, ITy>::VsToIs(is.get()); });

  return 0;
}

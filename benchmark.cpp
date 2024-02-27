#include <chrono>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <memory>

#include "hilbert.hpp"

int main() {
  using ITy = std::uint32_t;
  using ViTy = std::uint16_t;
  using STy = unsigned int;

  constexpr STy N = 2;
  constexpr STy K = 14;

  auto time_impl = [&](STy bytes, const char* desc, auto&& f) {
    auto start = std::chrono::high_resolution_clock::now();
    f();
    auto duration = std::chrono::high_resolution_clock::now() - start;
    double count = std::chrono::duration<double>(duration).count();
    double throughput = bytes / count / 1024 / 1024 / 1024;
    std::cout << desc << ":\t" << count << " s (" << throughput << " GiB/s)\n";
  };

  constexpr STy bytes = sizeof(ViTy) * N << N * K;
  auto time = [&](const char* desc, auto&& f) { time_impl(bytes, desc, f); };
  std::unique_ptr<ViTy[]> vs;
  time("new[]", [&]() { vs = std::make_unique<ViTy[]>(N << N * K); });
  // memset() is not necessary, but just gives a frame of reference
  // for how fast we can sequentially write to main memory.
  time("memset", [&]() { std::memset(vs.get(), 0, bytes); });
  time("IsToVs", [&]() {
    for (ITy i = 0; i < 1U << (N * K); i++) {
      Hilbert<N, K, ViTy, ITy>::IToV(i, &vs[N * i]);
    }
  });
  time("VsToIs", [&]() {
    for (ITy i = 0; i < 1U << (N * K); i++) {
      Hilbert<N, K, ViTy, ITy>::VToI(&vs[N * i]);
    }
  });

  return 0;
}

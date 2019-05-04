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

  {
    constexpr STy bytes = sizeof(ViTy) * N << N * K;
    auto time = [&](const char* desc, auto&& f) { time_impl(bytes, desc, f); };
    std::unique_ptr<ViTy[]> vs;
    time("new[]", [&]() { vs = std::make_unique<ViTy[]>(N << N * K); });
    // memset() is not necessary, but just gives a frame of reference
    // for how fast we can sequentially write to main memory.
    time("memset", [&]() { std::memset(vs.get(), 0, bytes); });
    time("IsToVs", [&]() { Hilbert<ViTy, ITy>::IsToVs<N, K>(vs.get()); });
  }

  {
    constexpr STy bytes = sizeof(ITy) << N * K;
    auto time = [&](const char* desc, auto&& f) { time_impl(bytes, desc, f); };
    std::unique_ptr<ITy[]> is;
    time("new[]", [&]() { is = std::make_unique<ITy[]>(1U << N * K); });
    time("memset", [&]() { std::memset(is.get(), 0, bytes); });
    time("IsToVs", [&]() { Hilbert<ViTy, ITy>::VsToIs<N, K>(is.get()); });
  }

  return 0;
}

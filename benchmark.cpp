#include <chrono>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <memory>

#include "hilbert.hpp"

int main() {
  using Int = std::uint16_t;
  using UInt = std::uint32_t;

  constexpr std::size_t N = 2;
  constexpr std::size_t K = 14;

  auto time_impl = [&](std::size_t bytes, const char* desc, auto&& f) {
    auto start = std::chrono::high_resolution_clock::now();
    f();
    auto duration = std::chrono::high_resolution_clock::now() - start;
    double count = std::chrono::duration<double>(duration).count();
    double throughput = bytes / count / 1024 / 1024 / 1024;
    std::cout << desc << ":\t" << count << " s (" << throughput << " GiB/s)\n";
  };

  {
    constexpr std::size_t bytes = sizeof(Int) * N << N * K;
    auto time = [&](const char* desc, auto&& f) { time_impl(bytes, desc, f); };
    std::unique_ptr<Int[]> vs;
    time("new[]", [&]() { vs = std::make_unique<Int[]>(N << N * K); });
    // memset() is not necessary, but just gives a frame of reference
    // for how fast we can sequentially write to main memory.
    time("memset", [&]() { std::memset(vs.get(), 0, bytes); });
    time("IsToVs", [&]() { Hilbert<Int, UInt>::IsToVs<N, K>(vs.get()); });
  }

  {
    constexpr std::size_t bytes = sizeof(UInt) << N * K;
    auto time = [&](const char* desc, auto&& f) { time_impl(bytes, desc, f); };
    std::unique_ptr<UInt[]> is;
    time("new[]", [&]() { is = std::make_unique<UInt[]>(1 << N * K); });
    time("memset", [&]() { std::memset(is.get(), 0, bytes); });
    time("IsToVs", [&]() { Hilbert<Int, UInt>::VsToIs<N, K>(is.get()); });
  }

  return 0;
}

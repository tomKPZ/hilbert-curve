#include "hilbert.hpp"

#include <memory>

int main(void) {
  constexpr unsigned int N = 2;
  constexpr unsigned int K = 14;

  auto curve = std::make_unique<std::array<short, N>[]>(1 << N * K);
  Hilbert<short, unsigned int>::Curve<N, K>(curve[0].data());
  // Hilbert<short, unsigned int>::Curve(N, K, curve[0].data());

  // short v[N];
  // for (unsigned int i = 0; i < 1 << N * K; i++) {
  //   Hilbert<short, unsigned int>::IToV(N, K, i, v);
  // }
  return 0;
}

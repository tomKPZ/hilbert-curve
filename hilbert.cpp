#include "hilbert.hpp"

#include <memory>

int main(void) {
  using Int = short;
  constexpr unsigned int N = 2;
  constexpr unsigned int K = 14;

  auto curve = std::make_unique<std::array<Int, N>[]>(1 << N * K);
  Hilbert<Int, unsigned int>::Curve(N, K, curve[0].data());

  // Int v[N];
  // for (unsigned int i = 0; i < 1 << N * K; i++) {
  //   Hilbert<Int, unsigned int>::IToV(N, K, i, v);
  // }
  return 0;
}

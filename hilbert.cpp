#include "hilbert.hpp"

#include <memory>

int main(void) {
  using Int = short;
  using UInt = std::size_t;

  constexpr UInt N = 2;
  constexpr UInt K = 14;

  auto curve = std::make_unique<Int[]>(N << N * K);
  Hilbert<Int, UInt>::Curve<N>(K, curve.get());

  // Int v[N];
  // for (UInt i = 0; i < 1 << N * K; i++) {
  //   Hilbert<Int, UInt>::IToV(N, K, i, v);
  // }
  return 0;
}

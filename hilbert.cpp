#include "hilbert.hpp"

#include <memory>

int main(void) {
  using Int = short;
  using UInt = unsigned int;

  constexpr UInt N = 12;
  constexpr UInt K = 2;

  auto curve = std::make_unique<Int[]>(N << N * K);
  Hilbert<Int, UInt>::Curve<N, K>(curve.get());

  // Int v[N];
  // for (UInt i = 0; i < 1 << N * K; i++) {
  //   Hilbert<Int, UInt>::IToV(N, K, i, v);
  // }
  return 0;
}

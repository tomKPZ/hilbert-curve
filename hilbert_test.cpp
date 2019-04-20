#include "hilbert.hpp"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iterator>
#include <memory>

template <bool Write = false> void OpTestData(std::size_t N, std::size_t K) {
  std::string fname =
      "test_data/" + std::to_string(N) + '_' + std::to_string(K);
  std::fstream f;
  f.open(fname, std::ios::binary | (Write ? std::ios::out : std::ios::in));
  auto curve = std::make_unique<int[]>(N << N * K);
  Hilbert<>::Curve(N, K, curve.get());
  for (std::size_t i = 0; i < 1 << (N * K); i++) {
    int center[N];
    Hilbert<>::IToV(N, K, i, center);
    assert(std::equal(center, center + N, curve.get() + N * i));

    int copy[N];
    std::copy(curve.get() + N * i, curve.get() + N * (i + 1), copy);
    auto i2 = Hilbert<>::VToI(N, K, copy);
    assert(i == i2);
    for (int x : copy) {
      assert(x == 0);
    }

    int offset[N];
    Hilbert<>::OffsetV(N, K, curve.get() + N * i, offset);
    int recenter[N];
    Hilbert<>::CenterV(N, K, offset, recenter);
    assert(std::equal(recenter, recenter + N, curve.get() + N * i));

    for (int x : offset) {
      uint8_t bytes[2];
      if constexpr (Write) {
        bytes[0] = (x & 0xff00) >> 8;
        bytes[1] = x & 0xff;
        f.write(reinterpret_cast<const char*>(bytes), sizeof(bytes));
      } else {
        f.read(reinterpret_cast<char*>(bytes), sizeof(bytes));
        assert(f);
        assert(x == (bytes[0] << 8) + bytes[1]);
      }
    }
  }
  if constexpr (!Write) {
    assert(f.peek() == EOF);
  }
}

int main() {
  static constexpr struct {
    std::size_t N;
    std::size_t K;
  } tests[] = {
      {0, 0}, {0, 1}, {1, 0}, {1, 1},  {1, 2}, {1, 3}, {1, 4}, {1, 5},  {1, 6},
      {1, 7}, {1, 8}, {1, 9}, {1, 10}, {2, 1}, {2, 2}, {2, 3}, {2, 4},  {2, 5},
      {2, 6}, {2, 7}, {3, 1}, {3, 2},  {3, 3}, {3, 4}, {4, 1}, {4, 2},  {4, 3},
      {5, 1}, {5, 2}, {6, 1}, {6, 2},  {7, 1}, {8, 1}, {9, 1}, {10, 1},
  };

  for (const auto& test : tests) {
    OpTestData(test.N, test.K);
  }

  return 0;
}

#include "hilbert.hpp"

#include <cassert>
#include <fstream>

template <std::size_t N, std::size_t K, bool Write = false>
void OpTestData() {
  std::string fname =
      "test_data/" + std::to_string(N) + '_' + std::to_string(K);
  std::fstream f;
  f.open(fname, std::ios::binary | (Write ? std::ios::out : std::ios::in));
  auto buf = Hilbert<>::Curve<N>(K);
  for (std::size_t i = 0; i < 1 << (N * K); i++) {
    assert(Hilbert<>::IToV<N>(K, i) == buf[i]);
    assert(Hilbert<>::VToI<N>(K, buf[i]) == i);
    const auto v = Hilbert<>::OffsetV<N>(K, buf[i]);
    assert(Hilbert<>::CenterV<N>(K, v) == buf[i]);
    for (int x : v) {
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
  OpTestData<0, 0>();
  OpTestData<0, 1>();
  OpTestData<1, 0>();

  OpTestData<1, 1>();
  OpTestData<1, 2>();
  OpTestData<1, 3>();
  OpTestData<1, 4>();
  OpTestData<1, 5>();
  OpTestData<1, 6>();
  OpTestData<1, 7>();
  OpTestData<1, 8>();
  OpTestData<1, 9>();
  OpTestData<1, 10>();

  OpTestData<2, 1>();
  OpTestData<2, 2>();
  OpTestData<2, 3>();
  OpTestData<2, 4>();
  OpTestData<2, 5>();
  OpTestData<2, 6>();
  OpTestData<2, 7>();

  OpTestData<3, 1>();
  OpTestData<3, 2>();
  OpTestData<3, 3>();
  OpTestData<3, 4>();

  OpTestData<4, 1>();
  OpTestData<4, 2>();
  OpTestData<4, 3>();

  OpTestData<5, 1>();
  OpTestData<5, 2>();

  OpTestData<6, 1>();
  OpTestData<6, 2>();

  OpTestData<7, 1>();

  OpTestData<8, 1>();

  OpTestData<9, 1>();

  OpTestData<10, 1>();

  return 0;
}

// hilbert-curve: Optimized N-dimensional Hilbert curve generation in C++
// Copyright (C) 2019 <tomKPZ@gmail.com>
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA

#include "hilbert.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <vector>

template <typename T>
T& check_aux(T&& t, const char* file, std::size_t line, const char* message) {
  if (!t) {
    std::cerr << "CHECK failed: " << file << ':' << line << ": " << message
              << std::endl;
    std::abort();
  }
  return t;
}

#define CHECK(x) check_aux(x, __FILE__, __LINE__, #x)

using Int = unsigned int;
using UInt = std::size_t;
std::vector<UInt> VsToIs(std::size_t N, std::size_t K) {
  std::vector<UInt> prev(1);
  prev[0] = 0;
  if (N == 0 || K == 0) {
    return prev;
  }
  for (std::size_t k = 0; k < K; ++k) {
    std::vector<UInt> is(1U << N * (k + 1));
    for (std::size_t i = 0; i < (1U << N); ++i) {
      UInt orthant = 0;
      bool parity = 0;
      for (std::size_t j = N; j-- > 0;) {
        parity ^= i & (1 << j);
        orthant |= parity << j;
      }

      std::size_t rotate = 0;
      if (i != 0 && i != (1U << N) - 1) {
        UInt j = (i - 1) >> 1;
        for (UInt bits = ~j & (j + 1); bits != 0; bits >>= 1) {
          ++rotate;
        }
      }

      UInt gray = ((i - 1) >> 1) ^ (i - 1);
      for (std::size_t j = 0; j < (1 << (N * k)); ++j) {
        UInt src = prev[j] + orthant * (1U << (N * k));
        UInt dest = 0;
        for (std::size_t vi = 0; vi < N; ++vi) {
          std::size_t nvi = (vi + rotate) % N;
          std::size_t mask = ((1 << k) - 1);
          std::size_t mask_shifted = mask << (vi * k);
          std::size_t value_shifted = mask_shifted & j;
          std::size_t value = value_shifted >> (vi * k);
          bool reflect = nvi == 0 ? (i == 0 || (i + 1) & 2) : gray & (1U << nvi);
          if (reflect) {
            value = (1 << k) - value - 1;
          }
          if (i & (1U << vi)) {
            value += (1 << k);
          }
          dest += value << (nvi * (k + 1));
        }
        std::cout << dest << " <- " << src << std::endl;
        is[dest] = src;
      }
    }
    prev = is;
    std::cout << std::endl;
  }
  return prev;
}

void RunTest(std::size_t N, std::size_t K) {
  std::string fname =
      "test_data/" + std::to_string(N) + '_' + std::to_string(K);
  std::fstream f;
  f.open(fname, std::ios::binary | std::ios::in);
  auto curve = std::make_unique<unsigned int[]>(N << N * K);
  auto curve_inverse = VsToIs(N, K);
  Hilbert<>::Curve(N, K, curve.get());
  for (std::size_t i = 0; i < 1U << (N * K); i++) {
    for (std::size_t j = 0; j < N; j++) {
      int x = curve[N * i + j];
      uint8_t bytes[2];
      f.read(reinterpret_cast<char*>(bytes), sizeof(bytes));
      CHECK(f);
      CHECK(x == (bytes[0] << 8) + bytes[1]);
    }

    unsigned int v[N];
    Hilbert<>::IToV(N, K, i, v);
    CHECK(std::equal(v, v + N, curve.get() + N * i));

    std::size_t index = 0;
    for (std::size_t j = 0; j < N; j++) {
      index |= (v[j]) << (j * K);
    }
    std::cout << "checking " << curve_inverse[index] << " vs " << i
              << std::endl;
    CHECK(curve_inverse[index] == i);

    CHECK(i == Hilbert<>::VToI(N, K, v));
    for (int x : v) {
      CHECK(x == 0);
    }
  }
  CHECK(f.peek() == EOF);
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
    RunTest(test.N, test.K);
  }

  return 0;
}

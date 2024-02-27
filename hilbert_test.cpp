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

#include <iostream>
#include <unordered_set>
#include <vector>

#define CHECK(x) check_aux(x, __FILE__, __LINE__, #x)

using ITy = unsigned int;
using ViTy = unsigned int;
using uint = unsigned int;

template <typename T>
T& check_aux(T&& t, const char* file, std::size_t line, const char* message) {
  if (!t) {
    std::cerr << file << ':' << line << ": CHECK failed: " << message
              << std::endl;
    std::abort();
  }
  return t;
}

std::string VToString(ViTy v[], uint N) {
  std::string s;
  for (uint i = 0; i < N; ++i) {
    s += std::to_string(v[i]) + ",";
  }
  return s;
}

template <uint N, uint K> void TestIToV() {
  std::unordered_set<std::string> seen;
  std::vector<ViTy> v0(N), v1(N);

  for (ITy i = 0; i < 1U << N * K; ++i) {
    Hilbert<N, K>::IToV(i, &v0[0]);
    // Check 1: v[N] should be unique
    auto vStr = VToString(&v0[0], N);
    CHECK(seen.find(vStr) == seen.end());
    seen.insert(vStr);

    // Check 2: All v[j] should be in the range [0, K-1]
    for (uint j = 0; j < N; ++j) {
      CHECK(v0[j] >= 0 && v0[j] < 1 << K);
    }

    if (i < (1U << N * K) - 1) {
      Hilbert<N, K>::IToV(i + 1, &v1[0]);
      int diffCount = 0;
      for (uint k = 0; k < N; ++k) {
        if (v0[k] != v1[k]) {
          // Check 3: v1 should be the same as v0 except for exactly one k
          CHECK(v1[k] == v0[k] + 1 || v1[k] == v0[k] - 1);
          diffCount++;
        }
      }
      CHECK(diffCount == 1);
    }
  }
}

template <uint N, uint K> void TestVToI() {
  std::vector<ViTy> v(N);
  for (ITy i = 0; i < 1U << N * K; ++i) {
    Hilbert<N, K>::IToV(i, &v[0]);
    CHECK((i == Hilbert<N, K>::VToI(&v[0])));
  }
}

template <uint N, uint K> void RunTest() {
  TestIToV<N, K>();
  TestVToI<N, K>();

  if constexpr (N == 1) {
    for (ITy i = 0; i < 1 << K; ++i) {
      ViTy v;
      Hilbert<1, K>::IToV(i, &v);
      CHECK(v == i);
    }
  }
}

int main() {
  RunTest<0, 0>();
  RunTest<0, 1>();
  RunTest<1, 0>();
  RunTest<1, 1>();
  RunTest<1, 2>();
  RunTest<1, 3>();
  RunTest<1, 4>();
  RunTest<1, 5>();
  RunTest<1, 6>();
  RunTest<1, 7>();
  RunTest<1, 8>();
  RunTest<1, 9>();
  RunTest<1, 10>();
  RunTest<2, 1>();
  RunTest<2, 2>();
  RunTest<2, 3>();
  RunTest<2, 4>();
  RunTest<2, 5>();
  RunTest<2, 6>();
  RunTest<2, 7>();
  RunTest<3, 1>();
  RunTest<3, 2>();
  RunTest<3, 3>();
  RunTest<3, 4>();
  RunTest<4, 1>();
  RunTest<4, 2>();
  RunTest<4, 3>();
  RunTest<5, 1>();
  RunTest<5, 2>();
  RunTest<6, 1>();
  RunTest<6, 2>();
  RunTest<7, 1>();
  RunTest<8, 1>();
  RunTest<9, 1>();
  RunTest<10, 1>();

  return 0;
}

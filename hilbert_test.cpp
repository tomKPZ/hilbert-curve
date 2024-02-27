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

#include <fstream>
#include <iostream>
#include <memory>
#include <unordered_set>
#include <vector>

#define CHECK(x) check_aux(x, __FILE__, __LINE__, #x)

using ITy = unsigned int;
using ViTy = unsigned int;
using STy = unsigned int;

template <typename T>
T& check_aux(T&& t, const char* file, std::size_t line, const char* message) {
  if (!t) {
    std::cerr << "CHECK failed: " << file << ':' << line << ": " << message
              << std::endl;
    std::abort();
  }
  return t;
}

std::string VToString(ViTy v[], STy N) {
  std::string s;
  for (STy i = 0; i < N; ++i) {
    s += std::to_string(v[i]) + ",";
  }
  return s;
}

void TestIToV(STy N, STy K) {
  ITy range = 1ULL << (N * K);
  std::unordered_set<std::string> seen;
  std::vector<ViTy> v0(N), v1(N);

  for (ITy i = 0; i < range; ++i) {
    Hilbert<>::IToV(N, K, i, &v0[0]);
    // Check 1: v[N] should be unique
    auto vStr = VToString(&v0[0], N);
    CHECK(seen.find(vStr) == seen.end());
    seen.insert(vStr);

    // Check 2: All v[j] should be in the range [0, K-1]
    for (STy j = 0; j < N; ++j) {
      CHECK(v0[j] >= 0 && v0[j] < 1 << K);
    }

    if (i < range - 1) {
      Hilbert<>::IToV(N, K, i + 1, &v1[0]);
      int diffCount = 0;
      for (STy k = 0; k < N; ++k) {
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

void RunTest(STy N, STy K) {
  TestIToV(N, K);
  std::string fname =
      "test_data/" + std::to_string(N) + '_' + std::to_string(K);
  std::fstream f;
  f.open(fname, std::ios::binary | std::ios::in);
  auto vs = std::make_unique<ViTy[]>(N << N * K);
  Hilbert<>::IsToVs(N, K, vs.get());
  auto is = std::make_unique<ITy[]>(1U << N * K);
  Hilbert<>::VsToIs(N, K, is.get());
  for (STy i = 0; i < 1U << (N * K); ++i) {
    for (STy j = 0; j < N; ++j) {
      ViTy x = vs[N * i + j];
      uint8_t bytes[2];
      f.read(reinterpret_cast<char*>(bytes), sizeof(bytes));
      CHECK(f);
      CHECK(x == static_cast<ViTy>((bytes[0] << 8) + bytes[1]));
    }

    ViTy v[N];
    for (STy k = 0; k < N; k++) {
      v[k] = vs[N * i + k];
    }

    STy index = 0;
    for (STy j = 0; j < N; ++j) {
      index |= v[j] << j * K;
    }
    CHECK(is[index] == i);

    CHECK(i == Hilbert<>::VToI(N, K, v));
    for (ViTy x : v) {
      CHECK(x == 0);
    }
  }
  CHECK(f.peek() == EOF);
}

int main() {
  static constexpr struct {
    STy N;
    STy K;
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

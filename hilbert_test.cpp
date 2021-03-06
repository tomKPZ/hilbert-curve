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

void RunTest(STy N, STy K) {
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
    Hilbert<>::IToV(N, K, i, v);
    CHECK(std::equal(v, v + N, vs.get() + N * i));

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

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

#include <iostream>
#include <memory>

#include "hilbert.hpp"

using ITy = unsigned int;
using ViTy = unsigned int;
using STy = unsigned int;

// Dimension of curve used in examples.
constexpr STy N = 2;

// Iteration of curve used in examples.
constexpr STy K = 2;

// Compute and print the points of the curve.
void BasicIsToVsExample() {
  ViTy vs[1U << N * K][N];
  Hilbert<>::IsToVs<N, K>(vs[0]);
  for (STy i = 0; i < 1U << N * K; ++i) {
    for (STy j = 0; j < N; ++j) {
      std::cout << vs[i][j] << '\t';
    }
    std::cout << std::endl;
  }
}

// Compute and point the points of the curve one-at-a-time and in
// reverse.
void IToVExample() {
  for (ITy i = (1U << N * K) - 1; i > 0; --i) {
    ViTy v[N];
    Hilbert<>::IToV<N, K>(i, v);
    for (STy j = 0; j < N; ++j) {
      std::cout << v[j] << '\t';
    }
    std::cout << std::endl;
  }
}

// Helper function for VToIExample().
void VToIExampleImpl(ViTy v[], STy j) {
  if (j == N) {
    std::cout << '(';
    for (STy k = 0; k < N; ++k) {
      std::cout << v[k];
      if (k != N - 1) {
        std::cout << ",\t";
      }
    }
    ViTy copy[N];
    for (STy k = 0; k < N; ++k) {
      copy[k] = v[k];
    }
    auto i = Hilbert<>::VToI<N, K>(copy);
    std::cout << "):\t" << i << std::endl;
  } else {
    for (STy k = 0; k < 1U << K; ++k) {
      v[j] = k;
      VToIExampleImpl(v, j + 1);
    }
  }
}

// Compute the indices that points would have given their coordinates.
void VToIExample() {
  ViTy v[N];
  VToIExampleImpl(v, 0);
}

// Helper function for ConstexprIsToVsExample().
constexpr std::array<ViTy, N << N * K> ConstexprIsToVsImpl() {
  std::array<ViTy, N << N * K> vs{};
  Hilbert<>::IsToVs<N, K>(vs.data());
  return vs;
}

// Compute the points of the curve at compile-time and print them.
void ConstexprIsToVsExample() {
  constexpr auto vs = ConstexprIsToVsImpl();
  for (STy i = 0; i < 1U << N * K; ++i) {
    for (STy j = 0; j < N; ++j) {
      std::cout << vs[i * N + j] << '\t';
    }
    std::cout << std::endl;
  }
}

int main(void) {
  std::cout << "BasicIsToVsExample" << std::endl;
  BasicIsToVsExample();

  std::cout << "IToVExample" << std::endl;
  IToVExample();

  std::cout << "VToIExample" << std::endl;
  VToIExample();

  std::cout << "ConstexprIsToVsExample" << std::endl;
  ConstexprIsToVsExample();
  return 0;
}

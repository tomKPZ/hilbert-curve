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

#include "hilbert.hpp"

using ITy = unsigned int;
using ViTy = unsigned int;
using STy = unsigned int;

// Dimension of curve used in examples.
constexpr STy N = 2;

// Iteration of curve used in examples.
constexpr STy K = 6;

// Compute and print the points of the curve.
void BasicIsToVsExample() {
  ViTy vs[1U << N * K][N];
  for (STy i = 0; i < 1U << N * K; ++i) {
    Hilbert<N, K>::IToV(i, vs[i]);
    std::cout << "[" << vs[i][0] << ',' << vs[i][1] << "]," << std::endl;
  }
}
int main(void) {
  BasicIsToVsExample();
  return 0;
}

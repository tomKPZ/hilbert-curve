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

////////////////////////////////////////////////////////////////////////////////
// Temporary scratch file for testing.
////////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include "hilbert.hpp"

constexpr int N = 2;
constexpr int K = 3;

void VToIExample2(const unsigned int orthant[], unsigned int v[], int j) {
  if (j == N) {
    std::cout << '(';
    for (int k = 0; k < N; k++) {
      std::cout << v[k];
      if (k != N - 1) {
        std::cout << ",\t";
      }
    }
    unsigned int copy[N];
    for (int k = 0; k < N; k++) {
      copy[k] = v[k] + orthant[k] * (1 << (K - 1));
    }
    auto i = Hilbert<>::VToI<N, K>(copy);
    std::cout << "):\t" << i % (1 << (N * (K - 1))) << std::endl;
  } else {
    for (unsigned int k = 0; k < 1 << (K - 1); ++k) {
      v[j] = k;
      VToIExample2(orthant, v, j + 1);
    }
  }
}

void VToIExample1(unsigned int orthant[], int j) {
  if (j == N) {
    unsigned int v[N];
    VToIExample2(orthant, v, 0);
    std::cout << std::endl;
  } else {
    for (unsigned int k = 0; k < 2; ++k) {
      orthant[j] = k;
      VToIExample1(orthant, j + 1);
    }
  }
}

void VToIExample() {
  unsigned int v[N];
  VToIExample1(v, 0);
}

int main(void) {
  VToIExample();
  return 0;
}

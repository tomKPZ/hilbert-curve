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
#include <vector>

#include "hilbert.hpp"

constexpr int N = 2;
constexpr int K = 2;

using Int = unsigned int;
using UInt = std::size_t;

void PrintPrevVs(unsigned int v[], int j) {
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
      copy[k] = v[k];
    }
    auto i = Hilbert<>::VToI<N, (K - 1)>(copy);
    std::cout << "):\t" << i << std::endl;
  } else {
    for (unsigned int k = 0; k < 1 << (K - 1); ++k) {
      v[j] = k;
      PrintPrevVs(v, j + 1);
    }
  }
}

void PrintPrevVs() {
  unsigned int v[N];
  PrintPrevVs(v, 0);
}

void PrintVs(const unsigned int orthant[], unsigned int v[], int j) {
  if (j == N) {
    unsigned int copy[N];
    for (int k = 0; k < N; k++) {
      copy[k] = v[k] + orthant[k] * (1 << (K - 1));
    }
    std::cout << '(';
    for (int k = 0; k < N; k++) {
      // std::cout << v[k];
      std::cout << v[k] + (orthant[k] << (K - 1));
      if (k != N - 1) {
        std::cout << ",\t";
      }
    }
    auto i = Hilbert<>::VToI<N, K>(copy);
    // std::cout << "):\t" << i % (1 << (N * (K - 1))) << std::endl;
    std::cout << "):\t" << i << std::endl;
  } else {
    for (unsigned int k = 0; k < 1 << (K - 1); ++k) {
      v[j] = k;
      PrintVs(orthant, v, j + 1);
    }
  }
}

void PrintVs(unsigned int orthant[], int j) {
  if (j == N) {
    unsigned int v[N];
    PrintVs(orthant, v, 0);
    std::cout << std::endl;
  } else {
    for (unsigned int k = 0; k < 2; ++k) {
      orthant[j] = k;
      PrintVs(orthant, j + 1);
    }
  }
}

void PrintVs() {
  unsigned int v[N];
  PrintVs(v, 0);
}

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

      std::size_t rotate = N - 1;
      if (orthant != 0 && orthant != (1U << N) - 1) {
        UInt j = (orthant - 1) >> 1;
        for (UInt bits = ~j & (j + 1); bits != 0; bits >>= 1) {
          --rotate;
        }
      }

      UInt gray = ((orthant - 1) >> 1) ^ (orthant - 1);
      UInt orthant_is = i * (1 << (N * k));
      for (std::size_t j = 0; j < (1 << (N * k)); ++j) {
        UInt src = prev[j] + orthant * (1 << (N * k));
        UInt dest = orthant_is;
        bool reflect = !(orthant == 0 || (orthant + 1) & 2);
        for (std::size_t vi = 0; vi < N; ++vi) {
          std::size_t mask = ((1 << k) - 1);
          std::size_t mask_shifted = mask << (vi * k);
          std::size_t value_shifted = mask_shifted & j;
          std::size_t value = value_shifted >> (vi * k);
          std::size_t nvi = (vi + rotate) % N;
          if (reflect) {
            value = (1 << k) - value - 1;
          }
          dest += value << (nvi * k);
          reflect = gray & (1U << (vi + 1));
        }
        is[dest] = src;
      }
    }
    prev = is;
  }
  return prev;
}

void PrintIs(const std::vector<UInt>& is) {
  Int v[N];
  for (std::size_t i = 0; i < 1 << (N * K); ++i) {
    for (std::size_t j = 0; j < N; ++j) {
      std::size_t mask = ((1 << K) - 1);
      std::size_t mask_shifted = mask << (j * K);
      std::size_t vj_shifted = mask_shifted & i;
      std::size_t vj = vj_shifted >> (j * K);
      v[j] = vj;
    }
    std::cout << '(';
    for (int j = 0; j < N; ++j) {
      std::cout << v[j];
      if (j != N - 1) {
        std::cout << ",\t";
      }
    }
    std::cout << "):\t" << is[i] << std::endl;
  }
}

int main(void) {
  {
    PrintPrevVs();
    std::cout << std::endl;
    PrintVs();
  }
  {
    std::vector<UInt> is = VsToIs(N, K);
    PrintIs(is);
  }
  return 0;
}

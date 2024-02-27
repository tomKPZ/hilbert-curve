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

#ifndef HILBERT_HPP
#define HILBERT_HPP

template <unsigned int N,
          unsigned int K,
          typename ViTy = unsigned int,
          typename ITy = unsigned int>
struct Hilbert {
  using uint = unsigned int;

  Hilbert() = delete;

  static constexpr void Scatter(ITy i, ViTy v[]) {
    for (uint n = 0; n < N; ++n) {
      v[n] = 0;
    }
    for (uint k = 0; k < K; ++k) {
      for (uint n = 0; n < N; ++n) {
        ViTy bit = (i >> (k * N + n)) & 1;
        v[N - n - 1] |= bit << k;
      }
    }
  }

  static constexpr ITy Gather(ViTy v[]) {
    ITy i = 0;
    for (uint k = 0; k < K; ++k) {
      for (uint n = 0; n < N; ++n) {
        ViTy bit = (v[N - n - 1] >> k) & 1;
        i |= bit << (k * N + n);
      }
    }
    return i;
  }

  static constexpr void IToV(ITy i, ViTy v[]) {
    Scatter(i ^ (i >> 1), v);
    for (unsigned int q = 1; q < K; ++q) {
      ViTy Q = 1U << q;
      ViTy P = Q - 1;
      for (int i = N; i-- > 0;) {
        if (v[i] & Q) {
          // invert
          v[0] ^= P;
        } else {
          // exchange
          ViTy t = (v[0] ^ v[i]) & P;
          v[0] ^= t;
          v[i] ^= t;
        }
      }
    }
  }

  static constexpr ITy VToI(ViTy v[]) {
    for (int q = K; q-- > 1;) {
      ViTy Q = 1U << q;
      ViTy P = Q - 1;
      for (uint i = 0; i < N; ++i) {
        if (v[i] & Q) {
          // invert
          v[0] ^= P;
        } else {
          // exchange
          ViTy t = (v[0] ^ v[i]) & P;
          v[0] ^= t;
          v[i] ^= t;
        }
      }
    }

    // Gray encode
    for (uint i = 1; i < N; ++i) {
      v[i] ^= v[i - 1];
    }
    ViTy t = 0;
    for (int q = K; q-- > 1;) {
      ViTy Q = 1U << q;
      if (v[N - 1] & Q) {
        t ^= Q - 1;
      }
    }
    for (uint i = 0; i < N; ++i) {
      v[i] ^= t;
    }

    return Gather(v);
  }
};

#endif  // HILBERT_HPP

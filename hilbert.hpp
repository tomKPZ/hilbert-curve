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

#include <cstdint>

#if defined(__GNUC__)
#define HB_ALWAYS_INLINE inline __attribute__((__always_inline__))
#elif defined(_MSC_VER)
#define HB_ALWAYS_INLINE __forceinline
#else
#define HB_ALWAYS_INLINE inline
#endif

template <typename Int = unsigned int, typename UInt = std::size_t>
class Hilbert {
 public:
  // Computes the K'th step of an N dimensional Hilbert curve.  The
  // result will be stored in curve, which must have space for 2^(N*K)
  // N-dimensional vectors, for a total of N*2^(N*K) integers.
  static constexpr void Curve(std::size_t N, std::size_t K, Int curve[]) {
    CurveImpl(N, K, curve);
  }
  template <std::size_t N>
  static constexpr void CurveN(std::size_t K, Int curve[]) {
    CurveImpl(N, K, curve);
  }
  template <std::size_t K>
  static constexpr void CurveK(std::size_t N, Int curve[]) {
    CurveImpl(N, K, curve);
  }
  template <std::size_t N, std::size_t K>
  static constexpr void Curve(Int curve[]) {
    CurveImpl(N, K, curve);
  }

  // Computes the i'th vector in Curve(N, K).  The result will be
  // stored in v.  Requires i < 2^(N*K).
  static constexpr void IToV(std::size_t N, std::size_t K, UInt i, Int v[]) {
    IToVImpl(N, K, i, v);
  }
  template <std::size_t N>
  static constexpr void IToVN(std::size_t K, UInt i, Int v[]) {
    IToVImpl(N, K, i, v);
  }
  template <std::size_t K>
  static constexpr void IToVK(std::size_t N, UInt i, Int v[]) {
    IToVImpl(N, K, i, v);
  }
  template <std::size_t N, std::size_t K>
  static constexpr void IToV(UInt i, Int v[]) {
    IToVImpl(N, K, i, v);
  }

  // Computes the index that v would have in the K'th step of an N
  // dimensional Hilbert curve.  Behavior is undefined if v is not a
  // point on the curve.  v will be zeroed when this function returns.
  static constexpr UInt VToI(std::size_t N, std::size_t K, Int v[]) {
    return VToIImpl(N, K, v);
  }
  template <std::size_t N> static constexpr UInt VToIN(std::size_t K, Int v[]) {
    return VToIImpl(N, K, v);
  }
  template <std::size_t K> static constexpr UInt VToIK(std::size_t N, Int v[]) {
    return VToIImpl(N, K, v);
  }
  template <std::size_t N, std::size_t K> static constexpr UInt VToI(Int v[]) {
    return VToIImpl(N, K, v);
  }

 private:
  Hilbert() = delete;

  HB_ALWAYS_INLINE static constexpr void CurveImpl(std::size_t N,
                                                   std::size_t K,
                                                   Int curve[]) {
    for (std::size_t i = 0; i < N; ++i) {
      curve[i] = 0;
    }
    for (std::size_t k = 1; k <= K; ++k) {
      for (std::size_t i = 1; i < (1U << N); ++i) {
        std::size_t rotate = N - 1;
        if (i != (1U << N) - 1) {
          std::size_t j = (i - 1) >> 1;
          for (std::size_t bits = ~j & (j + 1); bits != 0; bits >>= 1) {
            --rotate;
          }
        }

        std::size_t gray = ((i - 1) >> 1) ^ (i - 1);
        for (std::size_t p = 0; p < 1U << N * (k - 1); ++p) {
          std::size_t write_base = N * ((i << N * (k - 1)) + p);
          std::size_t read_base = N * p;
          bool reflect = !((i + 1) & 2);
          for (std::size_t j = 0; j < N; ++j) {
            std::size_t order = rotate + j >= N ? rotate + j - N : rotate + j;
            std::size_t coord = (i + (1U << j)) & (1U << (j + 1));
            Int offset = coord ? (1U << (k - 1)) : 0;
            Int temp = curve[read_base + order];
            temp = reflect ? ~temp & ((1U << (k - 1)) - 1) : temp;
            curve[write_base + j] = temp + offset;
            reflect = gray & (1U << (j + 1));
          }
        }
      }

      if (N > 0) {
        for (std::size_t p = 0; p < 1U << N * (k - 1); ++p) {
          Int temp = curve[N * (p + 1) - 1];
          for (std::size_t i = N - 1; i > 0; --i) {
            curve[N * p + i] = curve[N * p + i - 1];
          }
          curve[N * p] = temp;
        }
      }
    }
  }

  HB_ALWAYS_INLINE static constexpr void IToVImpl(std::size_t N,
                                                  std::size_t K,
                                                  UInt i,
                                                  Int v[]) {
    for (std::size_t j = 0; j < N; ++j) {
      v[j] = 0;
    }
    for (std::size_t k = 1; k <= K; ++k) {
      UInt orthant_i = i & ((1U << N * k) - 1);
      UInt orthant = orthant_i >> (N * (k - 1));

      std::size_t rotate = N - 1;
      if (orthant != 0 && orthant != (1U << N) - 1) {
        UInt j = (orthant - 1) >> 1;
        for (UInt bits = ~j & (j + 1); bits != 0; bits >>= 1) {
          --rotate;
        }
      }

      UInt gray = ((orthant - 1) >> 1) ^ (orthant - 1);
      bool reflect = !(orthant == 0 || (orthant + 1) & 2);
      for (std::size_t write = 0; write < N;) {
        for (std::size_t read = rotate; read < N; ++read) {
          if (rotate == write) {
            rotate = read;
          }

          bool coord = (orthant + (1U << write)) & (1U << (write + 1));
          Int offset = coord ? 1U << (k - 1) : 0;

          Int temp = v[read];
          temp = reflect ? ~temp & ((1U << (k - 1)) - 1) : temp;
          v[read] = v[write];
          v[write] = temp + offset;

          ++write;
          reflect = gray & (1U << write);
        }
      }
    }
  }

  HB_ALWAYS_INLINE static constexpr UInt VToIImpl(std::size_t N,
                                                  std::size_t K,
                                                  Int v[]) {
    UInt i = 0;
    for (std::size_t k = K; k > 0; --k) {
      UInt orthant = 0;
      bool parity = 0;
      for (std::size_t j = 0; j < N; ++j) {
        std::size_t jr = N - j - 1;
        parity ^= v[jr] >= (1U << (k - 1));
        orthant |= parity << jr;
      }

      std::size_t rotate = N - 1;
      if (orthant != 0 && orthant != (1U << N) - 1) {
        UInt j = (orthant - 1) >> 1;
        for (UInt bits = ~j & (j + 1); bits != 0; bits >>= 1) {
          --rotate;
        }
      }

      UInt gray = ((orthant - 1) >> 1) ^ (orthant - 1);
      bool reflect = !(orthant == 0 || (orthant + 1) & 2);
      std::size_t write = rotate;
      std::size_t rotate_outer = rotate;
      for (std::size_t order = 0; order < N;) {
        std::size_t count = rotate == 0 ? N : rotate;
        std::size_t read = rotate_outer - rotate;
        std::size_t rotate_inner = rotate;
        for (std::size_t r = 0; r < count; ++r) {
          if (rotate == N - order) {
            rotate = rotate_inner - r;
          }

          Int temp = v[read];
          temp = temp >= (1U << (k - 1)) ? temp - (1U << (k - 1)) : temp;
          v[read] = v[write];
          v[write] = reflect ? ~temp & ((1U << (k - 1)) - 1) : temp;

          ++order;
          reflect = gray & (1U << order);
          write = write + 1 == N ? 0 : write + 1;
          read = read + 1 == N ? 0 : read + 1;
        }
      }

      i += orthant << N * (k - 1);
    }
    return i;
  }
};

#endif  // HILBERT_HPP

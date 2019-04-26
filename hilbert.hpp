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

template <typename Int = int, typename UInt = std::size_t> class Hilbert {
 public:
  // Computes the K'th step of an N dimensional Hilbert curve.  The
  // result will be stored in curve, which must have space for 2^(N*K)
  // N-dimensional vectors, for a total of N*2^(N*K) integers.
  static constexpr void Curve(UInt N, UInt K, Int curve[]) {
    CurveImpl(N, K, curve);
  }
  template <UInt N> static constexpr void CurveN(UInt K, Int curve[]) {
    CurveImpl(N, K, curve);
  }
  template <UInt K> static constexpr void CurveK(UInt N, Int curve[]) {
    CurveImpl(N, K, curve);
  }
  template <UInt N, UInt K> static constexpr void Curve(Int curve[]) {
    CurveImpl(N, K, curve);
  }

  // Computes the i'th vector in Curve(N, K).  The result will be
  // stored in v.  Requires i < 2^(N*K).
  static constexpr void IToV(UInt N, UInt K, UInt i, Int v[]) {
    IToVImpl(N, K, i, v);
  }
  template <UInt N> static constexpr void IToVN(UInt K, UInt i, Int v[]) {
    IToVImpl(N, K, i, v);
  }
  template <UInt K> static constexpr void IToVK(UInt N, UInt i, Int v[]) {
    IToVImpl(N, K, i, v);
  }
  template <UInt N, UInt K> static constexpr void IToV(UInt i, Int v[]) {
    IToVImpl(N, K, i, v);
  }

  // Computes the index that v would have in the K'th step of an N
  // dimensional Hilbert curve.  Behavior is undefined if v is not a
  // point on the curve.  v will be zeroed when this function returns.
  static constexpr UInt VToI(UInt N, UInt K, Int v[]) {
    return VToIImpl(N, K, v);
  }
  template <UInt N> static constexpr UInt VToIN(UInt K, Int v[]) {
    return VToIImpl(N, K, v);
  }
  template <UInt K> static constexpr UInt VToIK(UInt N, Int v[]) {
    return VToIImpl(N, K, v);
  }
  template <UInt N, UInt K> static constexpr UInt VToI(Int v[]) {
    return VToIImpl(N, K, v);
  }

 private:
  Hilbert() = delete;

  HB_ALWAYS_INLINE static constexpr void CurveImpl(UInt N,
                                                   UInt K,
                                                   Int curve[]) {
    for (UInt i = 0; i < N; ++i) {
      curve[i] = 0;
    }
    for (UInt k = 1; k <= K; ++k) {
      for (UInt i = 1; i < (1U << N); ++i) {
        UInt rotate = N - 1;
        if (i != (1U << N) - 1) {
          UInt j = (i - 1) >> 1;
          for (UInt bits = ~j & (j + 1); bits != 0; bits >>= 1) {
            --rotate;
          }
        }

        UInt gray = ~(((i - 1) >> 1) ^ (i - 1));
        for (UInt p = 0; p < 1U << N * (k - 1); ++p) {
          UInt write_base = N * ((i << N * (k - 1)) + p);
          UInt read_base = N * p;
          UInt not_reflect = (i + 1) & 2;
          for (UInt j = 0; j < N; ++j) {
            UInt order = rotate + j >= N ? rotate + j - N : rotate + j;
            UInt coord = (i + (1 << j)) & (1 << (j + 1));
            Int offset = coord ? (1U << (k - 1)) : 0;
            Int temp = curve[read_base + order];
            temp = not_reflect ? temp : (1U << (k - 1)) - temp - 1;
            curve[write_base + j] = temp + offset;
            not_reflect = gray & (1 << (j + 1));
          }
        }
      }

      for (UInt p = 0; p < 1U << N * (k - 1); ++p) {
        for (UInt write = 0; write < N; ++write) {
          Int temp = curve[N * (p + 1) - 1];
          curve[N * (p + 1) - 1] = curve[N * p + write];
          curve[N * p + write] = temp;
        }
      }
    }
  }

  HB_ALWAYS_INLINE static constexpr void IToVImpl(UInt N,
                                                  UInt K,
                                                  UInt i,
                                                  Int v[]) {
    for (UInt j = 0; j < N; ++j) {
      v[j] = 0;
    }
    for (UInt k = 1; k <= K; ++k) {
      UInt orthant_i = i & ((1 << N * k) - 1);
      UInt orthant = orthant_i >> (N * (k - 1));

      UInt rotate = N - 1;
      if (orthant != 0 && orthant != (1U << N) - 1) {
        UInt j = (orthant - 1) >> 1;
        for (UInt bits = ~j & (j + 1); bits != 0; bits >>= 1) {
          --rotate;
        }
      }

      UInt gray = ~(((orthant - 1) >> 1) ^ (orthant - 1));
      UInt not_reflect = orthant == 0 || (orthant + 1) & 2;
      for (UInt write = 0; write < N;) {
        for (UInt read = rotate; read < N; ++read) {
          if (rotate == write) {
            rotate = read;
          }

          UInt coord = (orthant + (1 << write)) & (1 << (write + 1));
          Int offset = coord ? 1 << (k - 1) : 0;

          Int temp = v[read];
          temp = not_reflect ? temp : (1U << (k - 1)) - temp - 1;
          v[read] = v[write];
          v[write] = temp + offset;

          ++write;
          not_reflect = gray & (1 << write);
        }
      }
    }
  }

  HB_ALWAYS_INLINE static constexpr UInt VToIImpl(UInt N, UInt K, Int v[]) {
    UInt i = 0;
    for (UInt k = K; k > 0; --k) {
      UInt orthant = 0;
      UInt parity = 0;
      for (UInt j = 0; j < N; ++j) {
        UInt jr = N - j - 1;
        parity ^= v[jr] >= (1 << (k - 1));
        orthant |= parity << jr;
      }

      UInt rotate = N - 1;
      if (orthant != 0 && orthant != (1U << N) - 1) {
        UInt j = (orthant - 1) >> 1;
        for (UInt bits = ~j & (j + 1); bits != 0; bits >>= 1) {
          --rotate;
        }
      }

      UInt gray = ~(((orthant - 1) >> 1) ^ (orthant - 1));
      UInt not_reflect = orthant == 0 || (orthant + 1) & 2;
      UInt write = rotate;
      UInt rotate_outer = rotate;
      for (UInt order = 0; order < N;) {
        UInt count = rotate == 0 ? N : rotate;
        UInt read = rotate_outer - rotate;
        UInt rotate_inner = rotate;
        for (UInt r = 0; r < count; ++r) {
          if (rotate == N - order) {
            rotate = rotate_inner - r;
          }

          Int temp = v[read];
          temp = temp >= (1 << (k - 1)) ? temp - (1 << (k - 1)) : temp;
          v[read] = v[write];
          v[write] = not_reflect ? temp : (1U << (k - 1)) - temp - 1;

          ++order;
          not_reflect = gray & (1 << order);
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

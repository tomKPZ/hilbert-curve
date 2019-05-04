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

template <typename Int = unsigned int, typename UInt = unsigned int>
class Hilbert {
 public:
  static constexpr void IsToVs(std::size_t N, std::size_t K, Int vs[]) {
    IsToVsImpl(N, K, vs);
  }
  template <std::size_t N>
  static constexpr void IsToVsN(std::size_t K, Int vs[]) {
    IsToVsImpl(N, K, vs);
  }
  template <std::size_t K>
  static constexpr void IsToVsK(std::size_t N, Int vs[]) {
    IsToVsImpl(N, K, vs);
  }
  template <std::size_t N, std::size_t K>
  static constexpr void IsToVs(Int vs[]) {
    IsToVsImpl(N, K, vs);
  }

  static constexpr void VsToIs(std::size_t N, std::size_t K, UInt is[]) {
    VsToIsImpl(N, K, is);
  }
  template <std::size_t N>
  static constexpr void VsToIsN(std::size_t K, UInt is[]) {
    VsToIsImpl(N, K, is);
  }
  template <std::size_t K>
  static constexpr void VsToIsK(std::size_t N, UInt is[]) {
    VsToIsImpl(N, K, is);
  }
  template <std::size_t N, std::size_t K>
  static constexpr void VsToIs(UInt is[]) {
    VsToIsImpl(N, K, is);
  }

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

  HB_ALWAYS_INLINE static constexpr void IsToVsImpl(std::size_t N,
                                                    std::size_t K,
                                                    Int vs[]) {
    for (std::size_t i = 0; i < N; ++i) {
      vs[i] = 0;
    }
    if (N == 0 || K == 0) {
      return;
    }
    for (std::size_t k = 0; k < K; ++k) {
      for (std::size_t i = 1; i < (1U << N); ++i) {
        std::size_t rotate = N - 1;
        if (i != (1U << N) - 1) {
          std::size_t j = (i - 1) >> 1;
          for (std::size_t bits = ~j & (j + 1); bits != 0; bits >>= 1) {
            --rotate;
          }
        }

        std::size_t gray = ((i - 1) >> 1) ^ (i - 1);
        for (std::size_t p = 0; p < 1U << N * k; ++p) {
          std::size_t write_base = N * ((i << N * k) + p);
          std::size_t read_base = N * p;
          bool reflect = !((i + 1) & 2);
          for (std::size_t j = 0; j < N; ++j) {
            std::size_t order = rotate + j >= N ? rotate + j - N : rotate + j;
            std::size_t coord = (i + (1U << j)) & (1U << (j + 1));
            Int offset = coord ? (1U << k) : 0;
            Int temp = vs[read_base + order];
            temp = reflect ? ~temp & ((1U << k) - 1) : temp;
            vs[write_base + j] = temp + offset;
            reflect = gray & (1U << (j + 1));
          }
        }
      }

      for (std::size_t p = 0; p < 1U << N * k; ++p) {
        Int temp = vs[N * (p + 1) - 1];
        for (std::size_t i = N - 1; i > 0; --i) {
          vs[N * p + i] = vs[N * p + i - 1];
        }
        vs[N * p] = temp;
      }
    }
  }

  HB_ALWAYS_INLINE static constexpr void VsToIsImpl(std::size_t N,
                                                    std::size_t K,
                                                    UInt is[]) {
    is[0] = 0;
    if (N == 0 || K == 0) {
      return;
    }
    for (std::size_t k = 0; k < K; ++k) {
      for (std::size_t j = 0; j < (1U << (N * k)); ++j) {
        is[j | (1U << (N * (k + 1) - 1))] = is[j];
      }
      for (std::size_t j = 0; j < (1U << (N * k)); ++j) {
        UInt dest = 0;
        for (std::size_t nvi = 0; nvi < N; ++nvi) {
          std::size_t dsh = (nvi == 0 ? N - 1 : nvi - 1) * k;
          dest |= (((((1U << k) - 1) << dsh) & j) >> dsh) << (nvi * (k + 1));
        }
        is[dest] = is[j | (1U << (N * (k + 1) - 1))];
      }

      for (std::size_t i = 1; i < (1U << N); ++i) {
        UInt orthant = 0;
        bool parity = 0;
        for (std::size_t j = 0; j < N; ++j) {
          parity ^= (i & (1U << j)) >> j;
          orthant |= parity << (N - j - 1);
        }

        std::size_t rotate = N - 1;
        if (orthant != 0 && orthant != (1U << N) - 1) {
          UInt j = (orthant - 1) >> 1;
          for (UInt bits = ~j & (j + 1); bits != 0; bits >>= 1) {
            --rotate;
          }
        }

        UInt gray = ((orthant - 1) >> 1) ^ (orthant - 1);
        for (std::size_t j = 0; j < (1U << (N * k)); ++j) {
          UInt src = 0;
          UInt dest = 0;
          bool reflect = !(orthant == 0 || (orthant + 1) & 2);
          for (std::size_t nvi = 0; nvi < N; ++nvi) {
            std::size_t dsh = nvi + rotate;
            dsh = (dsh >= N ? dsh - N : dsh) * k;
            std::size_t di = ((((1U << k) - 1) << dsh) & j) >> dsh;
            di = reflect ? ~di & ((1U << k) - 1) : di;
            di |= (i & (1U << (N - nvi - 1))) ? 1U << k : 0;
            dest |= di << (nvi * (k + 1));

            std::size_t ssh = (nvi == 0 ? N - 1 : nvi - 1) * k;
            src |= (((((1U << k) - 1) << ssh) & j) >> ssh) << (nvi * (k + 1));

            reflect = gray & (1U << (nvi + 1));
          }
          is[dest] = is[src] + orthant * (1U << (N * k));
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
    if (N == 0 || K == 0) {
      return;
    }
    for (std::size_t k = 0; k < K; ++k) {
      UInt orthant_i = i & ((1U << N * (k + 1)) - 1);
      UInt orthant = orthant_i >> (N * k);

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
          Int offset = coord ? 1U << k : 0;

          Int temp = v[read];
          temp = reflect ? ~temp & ((1U << k) - 1) : temp;
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
    if (N == 0 || K == 0) {
      return i;
    }
    for (std::size_t k = K; k-- > 0;) {
      UInt orthant = 0;
      bool parity = 0;
      for (std::size_t j = N; j-- > 0;) {
        parity ^= v[j] >= (1U << k);
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
          temp = temp >= (1U << k) ? temp - (1U << k) : temp;
          v[read] = v[write];
          v[write] = reflect ? ~temp & ((1U << k) - 1) : temp;

          ++order;
          reflect = gray & (1U << order);
          write = write + 1 == N ? 0 : write + 1;
          read = read + 1 == N ? 0 : read + 1;
        }
      }

      i += orthant << N * k;
    }
    return i;
  }
};

#endif  // HILBERT_HPP

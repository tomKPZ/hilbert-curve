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

template <typename ViTy = unsigned int,
          typename ITy = unsigned int,
          typename STy = unsigned int>
struct Hilbert {
  Hilbert() = delete;

  static constexpr void IsToVs(STy N, STy K, ViTy vs[]) {
    for (STy i = 0; i < N; ++i) {
      vs[i] = 0;
    }
    if (N == 0 || K == 0) {
      return;
    }
    for (STy k = 0; k < K; ++k) {
      for (STy i = 1; i < (1U << N); ++i) {
        STy rotate = N - 1;
        if (i != (1U << N) - 1) {
          STy j = (i - 1) >> 1;
          for (STy bits = ~j & (j + 1); bits != 0; bits >>= 1) {
            --rotate;
          }
        }

        STy gray = ((i - 1) >> 1) ^ (i - 1);
        for (STy p = 0; p < 1U << N * k; ++p) {
          STy write_base = N * ((i << N * k) + p);
          STy read_base = N * p;
          bool reflect = !((i + 1) & 2);
          for (STy j = 0; j < N; ++j) {
            STy order = rotate + j >= N ? rotate + j - N : rotate + j;
            STy coord = (i + (1U << j)) & (1U << (j + 1));
            ViTy offset = coord ? (1U << k) : 0;
            ViTy temp = vs[read_base + order];
            temp = reflect ? ~temp & ((1U << k) - 1) : temp;
            vs[write_base + j] = temp + offset;
            reflect = gray & (1U << (j + 1));
          }
        }
      }

      for (STy p = 0; p < 1U << N * k; ++p) {
        ViTy temp = vs[N * (p + 1) - 1];
        for (STy i = N - 1; i > 0; --i) {
          vs[N * p + i] = vs[N * p + i - 1];
        }
        vs[N * p] = temp;
      }
    }
  }

  static constexpr void VsToIs(STy N, STy K, ITy is[]) {
    is[0] = 0;
    if (N == 0 || K == 0) {
      return;
    }
    for (STy k = 0; k < K; ++k) {
      for (STy j = 0; j < (1U << (N * k)); ++j) {
        is[j | (1U << (N * (k + 1) - 1))] = is[j];
      }
      for (STy j = 0; j < (1U << (N * k)); ++j) {
        ITy dest = 0;
        for (STy nvi = 0; nvi < N; ++nvi) {
          STy dsh = (nvi == 0 ? N - 1 : nvi - 1) * k;
          dest |= (((((1U << k) - 1) << dsh) & j) >> dsh) << (nvi * (k + 1));
        }
        is[dest] = is[j | (1U << (N * (k + 1) - 1))];
      }

      for (STy i = 1; i < (1U << N); ++i) {
        ITy orthant = 0;
        bool parity = 0;
        for (STy j = 0; j < N; ++j) {
          parity ^= (i & (1U << j)) >> j;
          orthant |= parity << (N - j - 1);
        }

        STy rotate = N - 1;
        if (orthant != 0 && orthant != (1U << N) - 1) {
          ITy j = (orthant - 1) >> 1;
          for (ITy bits = ~j & (j + 1); bits != 0; bits >>= 1) {
            --rotate;
          }
        }

        ITy gray = ((orthant - 1) >> 1) ^ (orthant - 1);
        for (STy j = 0; j < (1U << (N * k)); ++j) {
          ITy src = 0;
          ITy dest = 0;
          bool reflect = !(orthant == 0 || (orthant + 1) & 2);
          for (STy nvi = 0; nvi < N; ++nvi) {
            STy dsh = nvi + rotate;
            dsh = (dsh >= N ? dsh - N : dsh) * k;
            STy di = ((((1U << k) - 1) << dsh) & j) >> dsh;
            di = reflect ? ~di & ((1U << k) - 1) : di;
            di |= (i & (1U << (N - nvi - 1))) ? 1U << k : 0;
            dest |= di << (nvi * (k + 1));

            STy ssh = (nvi == 0 ? N - 1 : nvi - 1) * k;
            src |= (((((1U << k) - 1) << ssh) & j) >> ssh) << (nvi * (k + 1));

            reflect = gray & (1U << (nvi + 1));
          }
          is[dest] = is[src] + orthant * (1U << (N * k));
        }
      }
    }
  }

  static constexpr void Scatter(STy N, STy K, ITy i, ViTy v[]) {
    for (STy n = 0; n < N; ++n) {
      v[n] = 0;
    }
    for (STy k = 0; k < K; ++k) {
      for (STy n = 0; n < N; ++n) {
        ViTy bit = (i >> (k * N + n)) & 1;
        // STy vi = n;
        STy vi = N - 1 - n;
        STy ki = k;
        // STy ki = K - 1 - k;
        v[vi] |= bit << ki;
      }
    }
  }

  static constexpr ITy Gather(STy N, STy K, ViTy v[]) {
    ITy i = 0;
    for (STy k = 0; k < K; ++k) {
      for (STy n = 0; n < N; ++n) {
        // STy vi = n;
        STy vi = N - 1 - n;
        STy ki = k;
        // STy ki = K - 1 - k;
        ViTy bit = (v[vi] >> ki) & 1;
        i |= bit << (k * N + n);
      }
    }
    return i;
  }

  static constexpr void IToV(STy N, STy K, ITy i, ViTy v[]) {
    Scatter(N, K, i ^ (i >> 1), v);
    for (ViTy Q = 2; Q < 1U << K; Q <<= 1) {
      ViTy P = Q - 1;
      for (int i = N - 1; i >= 0; i--) {
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

  static constexpr ITy VToI(STy N, STy K, ViTy v[]) {
    if (N == 0 || K == 0) {
      return 0;
    }
    ViTy M = 1 << (K - 1);
    for (ViTy Q = M; Q > 1; Q >>= 1) {
      ViTy P = Q - 1;
      for (STy i = 0; i < N; i++) {
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
    for (STy i = 1; i < N; i++) {
      v[i] ^= v[i - 1];
    }
    ViTy t = 0;
    for (ViTy Q = M; Q > 1; Q >>= 1) {
      if (v[N - 1] & Q) {
        t ^= Q - 1;
      }
    }
    for (STy i = 0; i < N; i++) {
      v[i] ^= t;
    }

    return Gather(N, K, v);
  }
};

#endif  // HILBERT_HPP

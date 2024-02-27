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

  static constexpr void IToV(STy N, STy K, ITy i, ViTy v[]) {
    for (STy j = 0; j < N; ++j) {
      v[j] = 0;
    }
    if (N == 0 || K == 0) {
      return;
    }
    for (STy k = 0; k < K; ++k) {
      ITy orthant_i = i & ((1U << N * (k + 1)) - 1);
      ITy orthant = orthant_i >> (N * k);

      STy rotate = N - 1;
      if (orthant != 0 && orthant != (1U << N) - 1) {
        ITy j = (orthant - 1) >> 1;
        for (ITy bits = ~j & (j + 1); bits != 0; bits >>= 1) {
          --rotate;
        }
      }

      ITy gray = ((orthant - 1) >> 1) ^ (orthant - 1);
      bool reflect = !(orthant == 0 || (orthant + 1) & 2);
      for (STy write = 0; write < N;) {
        for (STy read = rotate; read < N; ++read) {
          if (rotate == write) {
            rotate = read;
          }

          bool coord = (orthant + (1U << write)) & (1U << (write + 1));
          ViTy offset = coord ? 1U << k : 0;

          ViTy temp = v[read];
          temp = reflect ? ~temp & ((1U << k) - 1) : temp;
          v[read] = v[write];
          v[write] = temp + offset;

          ++write;
          reflect = gray & (1U << write);
        }
      }
    }
  }

  static constexpr ITy VToI(STy N, STy K, ViTy v[]) {
    ITy i = 0;
    if (N == 0 || K == 0) {
      return i;
    }
    for (STy k = K; k-- > 0;) {
      ITy orthant = 0;
      bool parity = 0;
      for (STy j = N; j-- > 0;) {
        parity ^= v[j] >= (1U << k);
        orthant |= parity << j;
      }

      STy rotate = N - 1;
      if (orthant != 0 && orthant != (1U << N) - 1) {
        ITy j = (orthant - 1) >> 1;
        for (ITy bits = ~j & (j + 1); bits != 0; bits >>= 1) {
          --rotate;
        }
      }

      ITy gray = ((orthant - 1) >> 1) ^ (orthant - 1);
      bool reflect = !(orthant == 0 || (orthant + 1) & 2);
      STy write = rotate;
      STy rotate_outer = rotate;
      for (STy order = 0; order < N;) {
        STy count = rotate == 0 ? N : rotate;
        STy read = rotate_outer - rotate;
        STy rotate_inner = rotate;
        for (STy r = 0; r < count; ++r) {
          if (rotate == N - order) {
            rotate = rotate_inner - r;
          }

          ViTy temp = v[read];
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

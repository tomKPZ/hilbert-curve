#ifndef HILBERT_HPP
#define HILBERT_HPP

#include <type_traits>

template <typename Int = int, typename UInt = unsigned int> class Hilbert {
 public:
  // Curve: Computes the K'th step of an N dimensional Hilbert curve.
  template <UInt N, UInt K> static constexpr void Curve(Int curve[]);
  template <UInt N> static constexpr void Curve(UInt K, Int curve[]);
  static constexpr void Curve(UInt N, UInt K, Int curve[]);

  // IToV: Computes the i'th vector in Curve(N, K).
  template <UInt N, UInt K> static constexpr void IToV(UInt i, Int v[]);
  template <UInt N> static constexpr void IToV(UInt K, UInt i, Int v[]);
  static constexpr void IToV(UInt N, UInt K, UInt i, Int v[]);

  // IToV: Computes the index that v would have in Curve(N, K).
  // TODO: Make v const
  template <UInt N, UInt K> static constexpr UInt VToI(Int v[]);
  template <UInt N> static constexpr UInt VToI(UInt K, Int v[]);
  static constexpr UInt VToI(UInt N, UInt K, Int v[]);

  // Curve(), IToV(), and VToI() all operate on hilbert curves centered
  // at the origin with points separated a distance of 2.  For example,
  // the 2nd iteration of a 1D hilbert curve would have points at [{-3},
  // {-1}, {1}, {3}].  Sometimes this data is more useful based at 0
  // with a distance 1 between points.
  template <UInt N, UInt K>
  static constexpr void OffsetV(const Int center_v[], Int offset_v[]);
  template <UInt N>
  static constexpr void OffsetV(UInt K, const Int center_v[], Int offset_v[]);
  static constexpr void OffsetV(UInt N,
                                UInt K,
                                const Int center_v[],
                                Int offset_v[]);

  // The inverse operation of Offset() described above.
  template <UInt N, UInt K>
  static constexpr void CenterV(const Int offset_v[], Int center_v[]);
  template <UInt N>
  static constexpr void CenterV(UInt K, const Int offset_v[], Int center_v[]);
  static constexpr void CenterV(UInt N,
                                UInt K,
                                const Int offset_v[],
                                Int center_v[]);

 private:
  static_assert(std::is_signed_v<Int>);
  static_assert(std::is_unsigned_v<UInt>);

  Hilbert() = delete;

  // TODO: Make sure this gets inlined.
  static constexpr void CurveImpl(UInt N, UInt K, Int* vs) {
    if (K == 0) {
      for (UInt j = 0; j < N; j++) {
        vs[j] = 0;
      }
      return;
    }

    for (UInt i = 1; i < (1U << N); i++) {
      UInt d = N - 1;
      if (i != (1U << N) - 1) {
        UInt j = (i - 1) >> 1;
        for (UInt bits = ~j & (j + 1); bits != 0; bits >>= 1) {
          d--;
        }
      }
      for (UInt j = 0; j < N; j++) {
        UInt order = d + j >= N ? d + j - N : d + j;
        UInt c = i + (j == 0 ? 1 : (1 << (j + 2)) - (1 << j) - 1);
        bool sign = i == 0 && j == 0 ? 1 : c & (1 << (j + 1));
        UInt coord = (i + (1 << j)) & (1 << (j + 1));
        UInt offset = (coord ? 1 : -1) * (1 << (K - 1));
        for (UInt k = 0; k < 1U << N * (K - 1); k++) {
          vs[N * ((i << N * (K - 1)) + k) + j] =
              vs[N * k + order] * (sign ? 1 : -1) + offset;
        }
      }
    }

    UInt order = N - 1;
    for (UInt write = 0; write < N;) {
      for (UInt read = order; read < N; read++) {
        if (order == write) {
          order = read;
        }

        UInt c = (read == 0 ? 1 : (1 << (read + 2)) - (1 << read) - 1);
        bool sign = read == 0 ? 1 : c & (1 << (read + 1));
        UInt coord = (1 << read) & (1 << (read + 1));
        UInt offset = (coord ? 1 : -1) * (1 << (K - 1));
        for (UInt k = 0; k < 1U << N * (K - 1); k++) {
          Int temp = vs[N * k + order] * (sign ? 1 : -1) + offset;
          vs[N * k + order] = vs[N * k + write];
          vs[N * k + write] = temp;
        }

        write++;
      }
    }
  }

  static constexpr void IToVImpl(UInt N,
                                 UInt K,
                                 UInt i,
                                 Int* v,
                                 Int* orthant_v) {
    if (K == 0) {
      for (UInt j = 0; j < N; j++) {
        v[j] = 0;
      }
      return;
    }

    UInt orthant = i >> (N * (K - 1));
    UInt orthant_i = i & ((1 << N * (K - 1)) - 1);
    IToVImpl(N, K - 1, orthant_i, orthant_v, v);

    UInt d = N - 1;
    if (orthant != 0 && orthant != (1U << N) - 1) {
      UInt j = (orthant - 1) >> 1;
      for (UInt bits = ~j & (j + 1); bits != 0; bits >>= 1) {
        d--;
      }
    }

    for (UInt j = 0; j < N; j++) {
      UInt order = d + j >= N ? d + j - N : d + j;
      UInt c = orthant + (j == 0 ? 1 : (1 << (j + 2)) - (1 << j) - 1);
      bool sign = orthant == 0 && j == 0 ? 1 : c & (1 << (j + 1));
      UInt coord = (orthant + (1 << j)) & (1 << (j + 1));
      UInt offset = (coord ? 1 : -1) * (1 << (K - 1));
      v[j] = orthant_v[order] * (sign ? 1 : -1) + offset;
    }
  }

  static constexpr UInt VToIImpl(UInt N,
                                 UInt K,
                                 Int* v,
                                 Int* transformed,
                                 Int* orthant_v) {
    if (K == 0) {
      return 0;
    }

    UInt orthant = 0;
    UInt parity = 0;
    for (UInt j = 0; j < N; j++) {
      parity ^= v[N - j - 1] > 0;
      orthant |= parity << (N - j - 1);
      transformed[j] = v[j] + ((v[j] > 0 ? -1 : 1) << (K - 1));
    }

    UInt d = 1;
    if (orthant != 0 && orthant != (1U << N) - 1) {
      UInt j = (orthant - 1) >> 1;
      for (UInt bits = ~j & (j + 1); bits != 0; bits >>= 1) {
        d++;
      }
    }
    d = d == N ? 0 : d;

    for (UInt j = 0; j < N; j++) {
      UInt order = d + j >= N ? d + j - N : d + j;
      UInt c =
          orthant + (order == 0 ? 1 : (1 << (order + 2)) - (1 << order) - 1);
      bool sign = orthant == 0 && j == N - 1 ? 1 : c & (1 << (order + 1));
      orthant_v[j] = transformed[order] * (sign ? 1 : -1);
    }

    UInt offset = orthant * (1 << N * (K - 1));
    return offset + VToIImpl(N, K - 1, orthant_v, transformed, v);
  }

  static constexpr void OffsetVImpl(UInt N,
                                    UInt K,
                                    const Int center_v[],
                                    Int offset_v[]) {
    for (UInt i = 0; i < N; i++) {
      offset_v[i] = (center_v[i] + (1 << K)) >> 1;
    }
  }

  static constexpr void CenterVImpl(UInt N,
                                    UInt K,
                                    const Int offset_v[],
                                    Int center_v[]) {
    for (UInt i = 0; i < N; i++) {
      center_v[i] = offset_v[i] * 2 - (1 << K) + 1;
    }
  }
};

template <typename Int, typename UInt> template <UInt N, UInt K>
constexpr void Hilbert<Int, UInt>::Curve(Int curve[]) {
  if constexpr (K == 0) {
    CurveImpl(N, K, curve);
  } else {
    Curve<N, K - 1>(curve);
    CurveImpl(N, K, curve);
  }
}

template <typename Int, typename UInt> template <UInt N>
constexpr void Hilbert<Int, UInt>::Curve(UInt K, Int curve[]) {
  if (K == 0) {
    CurveImpl(N, K, curve);
  } else {
    Curve<N>(K - 1, curve);
    CurveImpl(N, K, curve);
  }
}

template <typename Int, typename UInt>
constexpr void Hilbert<Int, UInt>::Curve(UInt N, UInt K, Int curve[]) {
  if (K == 0) {
    CurveImpl(N, K, curve);
  } else {
    Curve(N, K - 1, curve);
    CurveImpl(N, K, curve);
  }
}

template <typename Int, typename UInt> template <UInt N, UInt K>
constexpr void Hilbert<Int, UInt>::IToV(UInt i, Int v[]) {
  IToV(N, K, i, v);
}

template <typename Int, typename UInt> template <UInt N>
constexpr void Hilbert<Int, UInt>::IToV(UInt K, UInt i, Int v[]) {
  IToV(N, K, i, v);
}

template <typename Int, typename UInt>
constexpr void Hilbert<Int, UInt>::IToV(UInt N, UInt K, UInt i, Int v[]) {
  Int orthant_v[N];
  IToVImpl(N, K, i, v, orthant_v);
}

template <typename Int, typename UInt> template <UInt N, UInt K>
constexpr UInt Hilbert<Int, UInt>::VToI(Int v[]) {
  return VToI(N, K, v);
}

template <typename Int, typename UInt> template <UInt N>
constexpr UInt Hilbert<Int, UInt>::VToI(UInt K, Int v[]) {
  return VToI(N, K, v);
}

template <typename Int, typename UInt>
constexpr UInt Hilbert<Int, UInt>::VToI(UInt N, UInt K, Int v[]) {
  Int transformed[N];
  Int orthant_v[N];
  return VToIImpl(N, K, v, transformed, orthant_v);
}

template <typename Int, typename UInt> template <UInt N, UInt K>
constexpr void Hilbert<Int, UInt>::OffsetV(const Int center_v[],
                                           Int offset_v[]) {
  OffsetV(N, K, center_v, offset_v);
}

template <typename Int, typename UInt> template <UInt N>
constexpr void Hilbert<Int, UInt>::OffsetV(UInt K,
                                           const Int center_v[],
                                           Int offset_v[]) {
  OffsetV(N, K, center_v, offset_v);
}

template <typename Int, typename UInt>
constexpr void Hilbert<Int, UInt>::OffsetV(UInt N,
                                           UInt K,
                                           const Int center_v[],
                                           Int offset_v[]) {
  OffsetVImpl(N, K, center_v, offset_v);
}

template <typename Int, typename UInt> template <UInt N, UInt K>
constexpr void Hilbert<Int, UInt>::CenterV(const Int offset_v[],
                                           Int center_v[]) {
  CenterV(N, K, offset_v, center_v);
}

template <typename Int, typename UInt> template <UInt N>
constexpr void Hilbert<Int, UInt>::CenterV(UInt K,
                                           const Int offset_v[],
                                           Int center_v[]) {
  CenterV(N, K, offset_v, center_v);
}

template <typename Int, typename UInt>
constexpr void Hilbert<Int, UInt>::CenterV(UInt N,
                                           UInt K,
                                           const Int offset_v[],
                                           Int center_v[]) {
  CenterVImpl(N, K, offset_v, center_v);
}

#endif  // HILBERT_HPP

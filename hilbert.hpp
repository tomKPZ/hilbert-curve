#ifndef HILBERT_HPP
#define HILBERT_HPP

#include <type_traits>

template <typename Int = int, typename UInt = unsigned int> class Hilbert {
 public:
  // Computes the K'th step of an N dimensional Hilbert curve.  The
  // result will be stored in curve, which must have space for 2^(N*K)
  // N-dimensional vectors, for a total of N*2^(N*K) integers.
  static constexpr void Curve(UInt N, UInt K, Int curve[]);

  // Computes the i'th vector in Curve(N, K).  The result will be
  // stored in v.  Requires i < 2^(N*K).
  static constexpr void IToV(UInt N, UInt K, UInt i, Int v[]);

  // Computes the index that v would have in the K'th step of an N
  // dimensional Hilbert curve.  Behavior is undefined if v is not a
  // point on the curve.  v will be zeroed when this function returns.
  static constexpr UInt VToI(UInt N, UInt K, Int v[]);

  // Curve(), IToV(), and VToI() all operate on Hilbert curves
  // centered at the origin with points separated a distance of 2.
  // For example, the 2nd iteration of a 1D Hilbert curve would have
  // points at [{-3}, {-1}, {1}, {3}].  Sometimes this data is more
  // useful based at 0 with a distance 1 between points.  In this
  // format, the curve would have points at [{0}, {1}, {2}, {3}].
  // OffsetV() and CenterV() converts between these formats.

  // Scales center_v down by a factor of 2 and shifts it to lie in the
  // first orthant.  Stores the result in offset_v.  center_v may
  // point to the same vector as offset_v.
  static constexpr void OffsetV(UInt N,
                                UInt K,
                                const Int center_v[],
                                Int offset_v[]);

  // Shifts offset_v to be centered at the origin and scales it up by
  // a factor of 2.  Stores the result in center_v.  offset_v may
  // point to the same vector as center_v.
  static constexpr void CenterV(UInt N,
                                UInt K,
                                const Int offset_v[],
                                Int center_v[]);

 private:
  static_assert(std::is_signed_v<Int>);
  static_assert(std::is_unsigned_v<UInt>);

  Hilbert() = delete;
};

template <typename Int, typename UInt>
constexpr void Hilbert<Int, UInt>::Curve(UInt N, UInt K, Int curve[]) {
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
      for (UInt ip = 0; ip < 1U << N * (k - 1); ++ip) {
        for (UInt j = 0; j < N; ++j) {
          UInt order = rotate + j >= N ? rotate + j - N : rotate + j;
          UInt c = i + (j == 0 ? 1 : (1 << (j + 2)) - (1 << j) - 1);
          UInt sign = i == 0 && j == 0 ? 1 : c & (1 << (j + 1));
          UInt coord = (i + (1 << j)) & (1 << (j + 1));
          UInt offset = (coord ? 1 : -1) * (1 << (k - 1));
          curve[N * ((i << N * (k - 1)) + ip) + j] =
              curve[N * ip + order] * (sign ? 1 : -1) + offset;
        }
      }
    }

    UInt sign = ((1 << (N + 1)) - (1 << (N - 1)) - 1) & (1 << N);
    UInt coord = (1 << (N - 1)) & (1 << N);
    UInt offset = (coord ? 1 : -1) * (1 << (k - 1));
    for (UInt ip = 0; ip < 1U << N * (k - 1); ++ip) {
      for (UInt write = 0; write < N; ++write) {
        Int temp = curve[N * (ip + 1) - 1] * (sign ? 1 : -1) + offset;
        curve[N * (ip + 1) - 1] = curve[N * ip + write];
        curve[N * ip + write] = temp;
      }
    }
  }
}

template <typename Int, typename UInt>
constexpr void Hilbert<Int, UInt>::IToV(UInt N, UInt K, UInt i, Int v[]) {
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

    for (UInt write = 0; write < N;) {
      UInt c =
          orthant + (write == 0 ? 1 : (1 << (write + 2)) - (1 << write) - 1);
      UInt sign = orthant == 0 && write == 0 ? 1 : c & (1 << (write + 1));
      UInt coord = (orthant + (1 << write)) & (1 << (write + 1));
      UInt offset = (coord ? 1 : -1) * (1 << (k - 1));
      for (UInt read = rotate; read < N; ++read) {
        if (rotate == write) {
          rotate = read;
        }

        Int temp = v[read] * (sign ? 1 : -1) + offset;
        v[read] = v[write];
        v[write] = temp;

        ++write;
      }
    }
  }
}

template <typename Int, typename UInt>
constexpr UInt Hilbert<Int, UInt>::VToI(UInt N, UInt K, Int v[]) {
  UInt i = 0;
  for (UInt k = K; k > 0; --k) {
    UInt orthant = 0;
    UInt parity = 0;
    for (UInt j = 0; j < N; ++j) {
      parity ^= v[N - j - 1] > 0;
      orthant |= parity << (N - j - 1);
      v[N - j - 1] = v[N - j - 1] + ((v[N - j - 1] > 0 ? -1 : 1) << (k - 1));
    }

    UInt rotate = 1;
    if (orthant != 0 && orthant != (1U << N) - 1) {
      UInt j = (orthant - 1) >> 1;
      for (UInt bits = ~j & (j + 1); bits != 0; bits >>= 1) {
        ++rotate;
      }
    }
    rotate = rotate == N ? 0 : rotate;

    UInt original_rotate = rotate;
    for (UInt write = 0; write < N;) {
      UInt order = original_rotate + write >= N ? original_rotate + write - N
                                                : original_rotate + write;
      UInt c =
          orthant + (order == 0 ? 1 : (1 << (order + 2)) - (1 << order) - 1);
      UInt sign = orthant == 0 && write == N - 1 ? 1 : c & (1 << (order + 1));
      for (UInt read = rotate; read < N; ++read) {
        if (rotate == write) {
          rotate = read;
        }

        Int temp = v[read] * (sign ? 1 : -1);
        v[read] = v[write];
        v[write] = temp;

        ++write;
      }
    }

    i += (orthant << N * (k - 1));
  }
  return i;
}

template <typename Int, typename UInt>
constexpr void Hilbert<Int, UInt>::OffsetV(UInt N,
                                           UInt K,
                                           const Int center_v[],
                                           Int offset_v[]) {
  for (UInt i = 0; i < N; ++i) {
    offset_v[i] = (center_v[i] + (1 << K)) >> 1;
  }
}

template <typename Int, typename UInt>
constexpr void Hilbert<Int, UInt>::CenterV(UInt N,
                                           UInt K,
                                           const Int offset_v[],
                                           Int center_v[]) {
  for (UInt i = 0; i < N; ++i) {
    center_v[i] = offset_v[i] * 2 - (1 << K) + 1;
  }
}

#endif  // HILBERT_HPP

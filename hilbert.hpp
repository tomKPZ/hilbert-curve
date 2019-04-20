#ifndef HILBERT_HPP
#define HILBERT_HPP

#include <type_traits>

template <typename Int = int, typename UInt = unsigned int> class Hilbert {
 public:
  // Curve: Computes the K'th step of an N dimensional Hilbert curve.
  static constexpr void Curve(UInt N, UInt K, Int curve[]);

  // IToV: Computes the i'th vector in Curve(N, K).
  static constexpr void IToV(UInt N, UInt K, UInt i, Int v[]);

  // IToV: Computes the index that v would have in Curve(N, K).
  // TODO: Make v const
  static constexpr UInt VToI(UInt N, UInt K, Int v[]);

  // Curve(), IToV(), and VToI() all operate on hilbert curves centered
  // at the origin with points separated a distance of 2.  For example,
  // the 2nd iteration of a 1D hilbert curve would have points at [{-3},
  // {-1}, {1}, {3}].  Sometimes this data is more useful based at 0
  // with a distance 1 between points.
  static constexpr void OffsetV(UInt N,
                                UInt K,
                                const Int center_v[],
                                Int offset_v[]);

  // The inverse operation of Offset() described above.
  static constexpr void CenterV(UInt N,
                                UInt K,
                                const Int offset_v[],
                                Int center_v[]);

 private:
  static_assert(std::is_signed_v<Int>);
  static_assert(std::is_unsigned_v<UInt>);

  Hilbert() = delete;

  static constexpr UInt VToIImpl(UInt N, UInt K, Int* v) {
    UInt orthant = 0;
    UInt parity = 0;
    for (UInt j = 0; j < N; j++) {
      parity ^= v[N - j - 1] > 0;
      orthant |= parity << (N - j - 1);
      v[N - j - 1] = v[N - j - 1] + ((v[N - j - 1] > 0 ? -1 : 1) << (K - 1));
    }

    UInt d = 1;
    if (orthant != 0 && orthant != (1U << N) - 1) {
      UInt j = (orthant - 1) >> 1;
      for (UInt bits = ~j & (j + 1); bits != 0; bits >>= 1) {
        d++;
      }
    }
    d = d == N ? 0 : d;

    UInt di = d;
    for (UInt write = 0; write < N;) {
      for (UInt read = d; read < N; read++) {
        if (d == write) {
          d = read;
        }

        UInt order = di + write >= N ? di + write - N : di + write;
        UInt c =
            orthant + (order == 0 ? 1 : (1 << (order + 2)) - (1 << order) - 1);
        bool sign = orthant == 0 && write == N - 1 ? 1 : c & (1 << (order + 1));

        Int temp = v[read] * (sign ? 1 : -1);
        v[read] = v[write];
        v[write] = temp;

        write++;
      }
    }

    return orthant * (1 << N * (K - 1));
  }
};

template <typename Int, typename UInt>
constexpr void Hilbert<Int, UInt>::Curve(UInt N, UInt K, Int curve[]) {
  for (UInt i = 0; i < N; i++) {
    curve[i] = 0;
  }
  for (UInt k = 1; k <= K; k++) {
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
        UInt offset = (coord ? 1 : -1) * (1 << (k - 1));
        for (UInt ip = 0; ip < 1U << N * (k - 1); ip++) {
          curve[N * ((i << N * (k - 1)) + ip) + j] =
              curve[N * ip + order] * (sign ? 1 : -1) + offset;
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
        UInt offset = (coord ? 1 : -1) * (1 << (k - 1));
        for (UInt ip = 0; ip < 1U << N * (k - 1); ip++) {
          Int temp = curve[N * ip + order] * (sign ? 1 : -1) + offset;
          curve[N * ip + order] = curve[N * ip + write];
          curve[N * ip + write] = temp;
        }

        write++;
      }
    }
  }
}

template <typename Int, typename UInt>
constexpr void Hilbert<Int, UInt>::IToV(UInt N, UInt K, UInt i, Int v[]) {
  for (UInt j = 0; j < N; j++) {
    v[j] = 0;
  }
  for (UInt k = 1; k <= K; k++) {
    UInt orthant_i = i & ((1 << N * k) - 1);
    UInt orthant = orthant_i >> (N * (k - 1));

    UInt d = N - 1;
    if (orthant != 0 && orthant != (1U << N) - 1) {
      UInt j = (orthant - 1) >> 1;
      for (UInt bits = ~j & (j + 1); bits != 0; bits >>= 1) {
        d--;
      }
    }

    for (UInt write = 0; write < N;) {
      for (UInt read = d; read < N; read++) {
        if (d == write) {
          d = read;
        }

        UInt c =
            orthant + (write == 0 ? 1 : (1 << (write + 2)) - (1 << write) - 1);
        bool sign = orthant == 0 && write == 0 ? 1 : c & (1 << (write + 1));
        UInt coord = (orthant + (1 << write)) & (1 << (write + 1));
        UInt offset = (coord ? 1 : -1) * (1 << (k - 1));

        Int temp = v[read] * (sign ? 1 : -1) + offset;
        v[read] = v[write];
        v[write] = temp;

        write++;
      }
    }
  }
}

template <typename Int, typename UInt>
constexpr UInt Hilbert<Int, UInt>::VToI(UInt N, UInt K, Int v[]) {
  UInt i = 0;
  for (UInt j = K; j > 0; j--) {
    i += VToIImpl(N, j, v);
  }
  return i;
}

template <typename Int, typename UInt>
constexpr void Hilbert<Int, UInt>::OffsetV(UInt N,
                                           UInt K,
                                           const Int center_v[],
                                           Int offset_v[]) {
  for (UInt i = 0; i < N; i++) {
    offset_v[i] = (center_v[i] + (1 << K)) >> 1;
  }
}

template <typename Int, typename UInt>
constexpr void Hilbert<Int, UInt>::CenterV(UInt N,
                                           UInt K,
                                           const Int offset_v[],
                                           Int center_v[]) {
  for (UInt i = 0; i < N; i++) {
    center_v[i] = offset_v[i] * 2 - (1 << K) + 1;
  }
}

#endif  // HILBERT_HPP

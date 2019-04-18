#ifndef HILBERT_HPP
#define HILBERT_HPP

#include <array>
#include <memory>
#include <type_traits>

template <typename Int = int, typename UInt = unsigned int>
class Hilbert {
 public:
  template <UInt N>
  using Vec = std::array<Int, N>;

  // Curve: Computes the K'th step of an N dimensional Hilbert curve.

  template <UInt N, UInt K>
  static constexpr std::array<Vec<N>, 1 << N * K> Curve();
  template <UInt N, UInt K>
  static std::unique_ptr<std::array<Vec<N>, 1 << N * K>> Curve(std::nullptr_t);
  template <UInt N, UInt K>
  static void Curve(std::array<Vec<N>, 1 << N * K>* curve);

  template <UInt N>
  static std::unique_ptr<Vec<N>[]> Curve(UInt K);
  template <UInt N>
  static constexpr void Curve(UInt K, Vec<N> curve[]);

  static std::unique_ptr<Int[]> Curve(UInt K);
  static constexpr void Curve(UInt N, UInt K, Int curve[]);

  // IToV: Computes the i'th vector in Curve(N, K).

  template <UInt N, UInt K>
  static constexpr Vec<N> IToV(UInt i);
  template <UInt N, UInt K>
  static constexpr void IToV(UInt i, Vec<N>*);

  template <UInt N>
  static constexpr Vec<N> IToV(UInt K, UInt i);
  template <UInt N>
  static constexpr void IToV(UInt K, UInt i, Vec<N>*);

  static constexpr void IToV(UInt N, UInt K, UInt i, Int v[]);

  // IToV: Computes the index that v would have in Curve(N, K).

  template <UInt N, UInt K>
  static constexpr UInt VToI(const Vec<N>& v);

  template <UInt N>
  static constexpr UInt VToI(UInt K, const Vec<N>& v);

  static constexpr UInt VToI(UInt N, UInt K, const Int v[]);

  // Curve(), IToV(), and VToI() all operate on hilbert curves centered
  // at the origin with points separated a distance of 2.  For example,
  // the 2nd iteration of a 1D hilbert curve would have points at [{-3},
  // {-1}, {1}, {3}].  Sometimes this data is more useful based at 0
  // with a distance 1 between points.

  template <UInt N, UInt K>
  static constexpr Vec<N> OffsetV(const Vec<N>& center_v);
  template <UInt N, UInt K>
  static constexpr void OffsetV(const Vec<N>& center_v, Vec<N>* offset_v);

  template <UInt N>
  static constexpr Vec<N> OffsetV(UInt K, const Vec<N>& center_v);
  template <UInt N>
  static constexpr void OffsetV(UInt K,
                                const Vec<N>& center_v,
                                Vec<N>* offset_v);

  static constexpr void OffsetV(UInt K, const Int center_v[], Int offset_v[]);

  // The inverse operation of Offset() described above.

  template <UInt N, UInt K>
  static constexpr Vec<N> CenterV(const Vec<N>& offset_v);
  template <UInt N, UInt K>
  static constexpr void CenterV(const Vec<N>& offset_v, Vec<N>* center_v);

  template <UInt N>
  static constexpr Vec<N> CenterV(UInt K, const Vec<N>& offset_v);
  template <UInt N>
  static constexpr void CenterV(UInt K,
                                const Vec<N>& offset_v,
                                Vec<N>* center_v);

  static constexpr void CenterV(UInt K, const Int offset_v[], Int center_v[]);

 private:
  static_assert(std::is_signed_v<Int>);
  static_assert(std::is_unsigned_v<UInt>);

  Hilbert() = delete;

  // TODO: Make sure this gets inlined.
  static constexpr void Curve(UInt N, UInt K, Int* vs, const Int* pvs) {
    if (K == 0) {
      for (UInt j = 0; j < N; j++) {
        vs[j] = 0;
      }
      return;
    }

    for (UInt i = 0; i < (1U << N); i++) {
      UInt d = N - 1;
      if (i != 0 && i != (1U << N) - 1) {
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
              pvs[N * k + order] * (sign ? 1 : -1) + offset;
        }
      }
    }
  }

  static constexpr void IToV(UInt N, UInt K, UInt i, Int* v, Int* orthant_v) {
    if (K == 0) {
      for (UInt j = 0; j < N; j++) {
        v[i] = 0;
      }
      return;
    }

    UInt orthant = i >> (N * (K - 1));
    UInt orthant_i = i & ((1 << N * (K - 1)) - 1);
    IToV(N, K - 1, orthant_i, orthant_v, v);

    UInt d = N - 1;
    if (orthant != 0 && orthant != (1 << N) - 1) {
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

  static constexpr UInt VToI(UInt N,
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
    if (orthant != 0 && orthant != (1 << N) - 1) {
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
    return offset + VToI(N, K - 1, orthant_v, transformed, v);
  }
};

template <typename Int, typename UInt>
template <UInt N, UInt K>
constexpr std::array<std::array<Int, N>, 1 << N * K>
Hilbert<Int, UInt>::Curve() {
  std::array<Vec<N>, 1 << N * K> curve{};
  Curve(K, curve.data());
  return curve;
}

template <typename Int, typename UInt>
template <UInt N>
std::unique_ptr<std::array<Int, N>[]> Hilbert<Int, UInt>::Curve(UInt K) {
  std::unique_ptr<Vec<N>[]> curve(new Vec<N>[1 << N * K]);
  if (K == 0) {
    Curve(N, K, curve.get()[0].data(), nullptr);
    return curve;
  }
  auto prev = Curve<N>(K - 1);
  Curve(N, K, curve.get()[0].data(), prev.get()[0].data());
  return curve;
}

template <typename Int, typename UInt>
template <UInt N>
constexpr void Hilbert<Int, UInt>::Curve(UInt K, Vec<N>* vs) {
  Vec<N> v{};
  Curve(N, K, vs[0].data(), v.data());
}

template <typename Int, typename UInt>
template <UInt N>
constexpr std::array<Int, N> Hilbert<Int, UInt>::IToV(UInt K, UInt i) {
  Vec<N> v{};
  Vec<N> orthant_v{};
  IToV(N, K, i, v.data(), orthant_v.data());
  return v;
}

template <typename Int, typename UInt>
template <UInt N>
constexpr UInt Hilbert<Int, UInt>::VToI(UInt K, const Vec<N>& v) {
  Vec<N> vec = v;
  Vec<N> transformed{};
  Vec<N> orthant_v{};
  return VToI(N, K, vec.data(), transformed.data(), orthant_v.data());
}

template <typename Int, typename UInt>
template <UInt N>
constexpr std::array<Int, N> Hilbert<Int, UInt>::OffsetV(
    UInt K,
    const Vec<N>& center_v) {
  Vec<N> offset_v{};
  for (UInt i = 0; i < N; i++) {
    offset_v[i] = (center_v[i] + (1 << K)) >> 1;
  }
  return offset_v;
}

template <typename Int, typename UInt>
template <UInt N>
constexpr std::array<Int, N> Hilbert<Int, UInt>::CenterV(
    UInt K,
    const Vec<N>& offset_v) {
  Vec<N> center_v{};
  for (UInt i = 0; i < N; i++) {
    center_v[i] = offset_v[i] * 2 - (1 << K) + 1;
  }
  return center_v;
}

#endif  // HILBERT_HPP

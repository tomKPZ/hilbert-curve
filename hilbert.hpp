#ifndef HILBERT_HPP
#define HILBERT_HPP

#include <array>
#include <cstdint>
#include <memory>
#include <type_traits>

template <std::size_t N, typename Int = int>
class Hilbert {
 public:
  using Vec = std::array<Int, N>;

  // Computes the K'th iteration of the Hilbert curve.  Returns an
  // array of 2^(N*K) Vec's.  Use this version of Curve() when K is a
  // small constant known at compile time.  Example:
  //     // Compute the 3rd iteration of a 2D Hilbert curve.
  //     constexpr auto curve = Hilbert<2>::Curve<3>();
  //     for (const auto& v : curve) { ... }
  template <std::size_t K>
  static constexpr std::array<Vec, 1 << N * K> Curve();

  // Computes the K'th iteration of the Hilbert curve.  Returns a
  // heap-allocated array of 2^(N*K) Vec's.  Use this version of
  // Curve() when K is large or not known at compile time; and you
  // haven't already allocated memory for the result.  Example:
  //     // Compute the 3rd iteration of a 2D Hilbert curve.
  //     auto curve = Hilbert<2>::Curve(3);
  //     for (std::size_t i = 0; i < 1 << (2*3); i++) {
  //       const auto& v = curve[i];
  //       ...
  //     }
  static std::unique_ptr<Vec[]> Curve(std::size_t K);

  // Computes the K'th iteration of the Hilbert curve.  Fills vs
  // with 2^(N*K) Vec's.  Use this version of Curve() when K is large
  // or not known at compile time; and you have already allocated
  // memory for the result.  Example:
  //     // Compute the 3rd iteration of a 2D Hilbert curve.
  //     Hilbert<2>::Vec curve[1 << (2*3)];
  //     Hilbert<2>::Curve(curve, 3);
  //     for (const auto& v : curve) { ... }
  static constexpr void Curve(std::size_t K, Vec* vs);

  // Returns the i'th vector of Curve(K).  Example:
  //     // Compute the 5th vector in the 3rd iteration of a 2D
  //     // Hilbert curve.
  //     auto v = Hilbert<2>::IToV(3, 5);  // v is {-7, -1}.
  static constexpr Vec IToV(std::size_t K, std::size_t i);

  // Returns the index that v would have in Curve(K).  Example:
  //     // Compute the index of the vector {-7, -1} in the 3rd
  //     // iteration of a 2D Hilbert curve.
  //     auto i = Hilbert<2>::VToI(3, {-7, -1});  // i is 5.
  static constexpr std::size_t VToI(std::size_t K, const Vec& v);

  // Curve(), IToV(), and VToI() all operate on hilbert curves
  // centered at the origin with points separated a distance of 2.
  // For example, the 2nd iteration of a 1D hilbert curve would have
  // points at [{-3}, {-1}, {1}, {3}].  Sometimes this data is more
  // useful based at 0 with a distance 1 between points.  Example:
  //     // Offset the points of the 2nd iteration of a 1D Hilbert
  //     // curve.
  //     auto v0 = Hilbert<1>::Offset(2, {-3});  // v0 is {0}.
  //     auto v1 = Hilbert<1>::Offset(2, {-1});  // v1 is {1}.
  //     auto v2 = Hilbert<1>::Offset(2, {+1});  // v1 is {2}.
  //     auto v3 = Hilbert<1>::Offset(2, {+3});  // v1 is {3}.

  static constexpr Vec OffsetV(std::size_t K, const Vec& center_v);

  // The inverse operation of Offset() described above.  Example:
  //     // Center the points of the 2nd iteration of a 1D Hilbert
  //     // curve.
  //     auto v0 = Hilbert<1>::Center(2, {0});  // v0 is {-3}.
  //     auto v1 = Hilbert<1>::Center(2, {1});  // v1 is {-1}.
  //     auto v2 = Hilbert<1>::Center(2, {2});  // v1 is {+1}.
  //     auto v3 = Hilbert<1>::Center(2, {3});  // v1 is {+3}.
  static constexpr Vec CenterV(std::size_t K, const Vec& offset_v);

 private:
  static_assert(std::is_signed_v<Int>);

  Hilbert() = delete;
};

// static
template <std::size_t N, typename Int>
template <std::size_t K>
constexpr std::array<typename Hilbert<N, Int>::Vec, 1 << N * K>
Hilbert<N, Int>::Curve() {
  std::array<Vec, 1 << N * K> ret{};
  Curve(K, &ret[0]);
  return ret;
}

// static
template <std::size_t N, typename Int>
std::unique_ptr<typename Hilbert<N, Int>::Vec[]> Hilbert<N, Int>::Curve(
    std::size_t K) {
  std::unique_ptr<Hilbert<N, Int>::Vec[]> ret(
      new Hilbert<N, Int>::Vec[1 << (N * K)]);
  Curve(K, ret.get());
  return ret;
}

template <typename Int>
constexpr void Curve(std::size_t N, std::size_t K, Int* vs, Int* v) {
  if (K == 0) {
    for (std::size_t j = 0; j < N; j++) {
      vs[j] = 0;
    }
    return;
  }

  Int* prev_end = vs + N * (1 << N * K);
  Int* prev_start = prev_end - N * (1 << N * (K - 1));
  Curve(N, K - 1, prev_start, v);
  size_t current = 0;
  for (std::size_t i = 0; i < (1U << N); i++) {
    for (const Int* p = prev_start; p != prev_end; p += N) {
      // TODO: Avoid copy for orthants 0 to 2^N - 2.
      for (std::size_t j = 0; j < N; j++) {
        v[j] = p[j];
      }

      std::size_t d = N - 1;
      if (i != 0 && i != (1U << N) - 1) {
        std::size_t j = (i - 1) >> 1;
        j = ~j & (j + 1);
        while (j != 0) {
          j >>= 1;
          d--;
        }
      }

      Int* v2 = vs + N * current++;
      for (std::size_t j = 0; j < N; j++) {
        std::size_t order = d + j >= N ? d + j - N : d + j;
        std::size_t c = i + (j == 0 ? 1 : (1 << (j + 2)) - (1 << j) - 1);
        bool sign = i == 0 && j == 0 ? 1 : c & (1 << (j + 1));
        std::size_t coord = (i + (1 << j)) & (1 << (j + 1));
        std::size_t offset = (coord ? 1 : -1) * (1 << (K - 1));
        v2[j] = v[order] * (sign ? 1 : -1) + offset;
      }
    }
  }
}

// static
template <std::size_t N, typename Int>
constexpr void Hilbert<N, Int>::Curve(std::size_t K, Vec* vs) {
  Vec v{};
  ::Curve<Int>(N, K, vs[0].data(), v.data());
}

template <typename Int>
constexpr void IToV(std::size_t N,
                    std::size_t K,
                    std::size_t i,
                    Int* v,
                    Int* orthant_v) {
  if (K == 0) {
    for (std::size_t j = 0; j < N; j++) {
      v[i] = 0;
    }
    return;
  }

  std::size_t orthant = i >> (N * (K - 1));
  std::size_t orthant_i = i & ((1 << N * (K - 1)) - 1);
  IToV<Int>(N, K - 1, orthant_i, orthant_v, v);

  std::size_t d = N - 1;
  if (orthant != 0 && orthant != (1 << N) - 1) {
    std::size_t j = (orthant - 1) >> 1;
    for (std::size_t bits = ~j & (j + 1); bits != 0; bits >>= 1) {
      d--;
    }
  }

  for (std::size_t j = 0; j < N; j++) {
    std::size_t order = d + j >= N ? d + j - N : d + j;
    std::size_t c = orthant + (j == 0 ? 1 : (1 << (j + 2)) - (1 << j) - 1);
    bool sign = orthant == 0 && j == 0 ? 1 : c & (1 << (j + 1));
    std::size_t coord = (orthant + (1 << j)) & (1 << (j + 1));
    std::size_t offset = (coord ? 1 : -1) * (1 << (K - 1));
    v[j] = orthant_v[order] * (sign ? 1 : -1) + offset;
  }
}

// static
template <std::size_t N, typename Int>
constexpr typename Hilbert<N, Int>::Vec Hilbert<N, Int>::IToV(std::size_t K,
                                                              std::size_t i) {
  Vec v{};
  Vec orthant_v{};
  ::IToV<Int>(N, K, i, v.data(), orthant_v.data());
  return v;
}

template <typename Int>
constexpr std::size_t VToI(std::size_t N,
                           std::size_t K,
                           Int* v,
                           Int* transformed,
                           Int* orthant_v) {
  if (K == 0) {
    return 0;
  }

  std::size_t orthant = 0;
  std::size_t parity = 0;
  for (std::size_t j = 0; j < N; j++) {
    parity ^= v[N - j - 1] > 0;
    orthant |= parity << (N - j - 1);
    transformed[j] = v[j] + ((v[j] > 0 ? -1 : 1) << (K - 1));
  }

  std::size_t d = 0;
  if (orthant != 0 && orthant != (1 << N) - 1) {
    std::size_t j = (orthant - 1) >> 1;
    for (std::size_t bits = ~j & (j + 1); bits != 0; bits >>= 1) {
      d++;
    }
  }
  d = d == N - 1 ? 0 : d + 1;

  for (std::size_t j = 0; j < N; j++) {
    std::size_t order = d + j >= N ? d + j - N : d + j;
    std::size_t c =
        orthant + (order == 0 ? 1 : (1 << (order + 2)) - (1 << order) - 1);
    bool sign = orthant == 0 && j == N - 1 ? 1 : c & (1 << (order + 1));
    orthant_v[j] = transformed[order] * (sign ? 1 : -1);
  }

  std::size_t offset = orthant * (1 << N * (K - 1));
  return offset + VToI<Int>(N, K - 1, orthant_v, transformed, v);
}

// static
template <std::size_t N, typename Int>
constexpr std::size_t Hilbert<N, Int>::VToI(std::size_t K, const Vec& v) {
  Vec vec = v;
  Vec transformed{};
  Vec orthant_v{};
  return ::VToI<Int>(N, K, vec.data(), transformed.data(), orthant_v.data());
}

// static
template <std::size_t N, typename Int>
constexpr typename Hilbert<N, Int>::Vec Hilbert<N, Int>::OffsetV(
    std::size_t K,
    const Vec& center_v) {
  Vec offset_v{};
  for (std::size_t i = 0; i < N; i++) {
    offset_v[i] = (center_v[i] + (1 << K)) >> 1;
  }
  return offset_v;
}

// static
template <std::size_t N, typename Int>
constexpr typename Hilbert<N, Int>::Vec Hilbert<N, Int>::CenterV(
    std::size_t K,
    const Vec& offset_v) {
  Vec center_v{};
  for (std::size_t i = 0; i < N; i++) {
    center_v[i] = offset_v[i] * 2 - (1 << K) + 1;
  }
  return center_v;
}

#endif  // HILBERT_HPP

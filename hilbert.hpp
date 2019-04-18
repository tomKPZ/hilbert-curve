#ifndef HILBERT_HPP
#define HILBERT_HPP

#include <array>
#include <cstdint>
#include <memory>

namespace hilbert {

// Computes the K'th iteration of the Hilbert curve.  Returns an array
// of 2^(N*K) Vec's.  Use this version of Curve() when K is a small
// constant known at compile time.  Example:
//     // Compute the 3rd iteration of a 2D Hilbert curve.
//     constexpr auto curve = Hilbert<2>::Curve<3>();
//     for (const auto& v : curve) { ... }
template <std::size_t N, std::size_t K, typename Int = int>
constexpr std::array<std::array<Int, N>, 1 << N * K> Curve();

// Computes the K'th iteration of the Hilbert curve.  Returns a
// heap-allocated array of 2^(N*K) Vec's.  Use this version of Curve()
// when K is large or not known at compile time; and you haven't
// already allocated memory for the result.  Example:
//     // Compute the 3rd iteration of a 2D Hilbert curve.
//     auto curve = Hilbert<2>::Curve(3);
//     for (std::size_t i = 0; i < 1 << (2*3); i++) {
//       const auto& v = curve[i];
//       ...
//     }
template <std::size_t N, typename Int = int>
std::unique_ptr<std::array<Int, N>[]> Curve(std::size_t K);

// Computes the K'th iteration of the Hilbert curve.  Fills vs with
// 2^(N*K) Vec's.  Use this version of Curve() when K is large or not
// known at compile time; and you have already allocated memory for
// the result.  Example:
//     // Compute the 3rd iteration of a 2D Hilbert curve.
//     Hilbert<2>::Vec curve[1 << (2*3)];
//     Hilbert<2>::Curve(curve, 3);
//     for (const auto& v : curve) { ... }
template <std::size_t N, typename Int = int>
constexpr void Curve(std::size_t K, std::array<Int, N>* vs);

// Returns the i'th vector of Curve(K).  Example:
//     // Compute the 5th vector in the 3rd iteration of a 2D
//     // Hilbert curve.
//     auto v = Hilbert<2>::IToV(3, 5);  // v is {-7, -1}.
template <std::size_t N, typename Int = int>
constexpr std::array<Int, N> IToV(std::size_t K, std::size_t i);

// Returns the index that v would have in Curve(K).  Example:
//     // Compute the index of the vector {-7, -1} in the 3rd
//     // iteration of a 2D Hilbert curve.
//     auto i = Hilbert<2>::VToI(3, {-7, -1});  // i is 5.
template <std::size_t N, typename Int = int>
constexpr std::size_t VToI(std::size_t K, const std::array<Int, N>& v);

// Curve(), IToV(), and VToI() all operate on hilbert curves centered
// at the origin with points separated a distance of 2.  For example,
// the 2nd iteration of a 1D hilbert curve would have points at [{-3},
// {-1}, {1}, {3}].  Sometimes this data is more useful based at 0
// with a distance 1 between points.  Example:
//     // Offset the points of the 2nd iteration of a 1D Hilbert
//     // curve.
//     auto v0 = Hilbert<1>::Offset(2, {-3});  // v0 is {0}.
//     auto v1 = Hilbert<1>::Offset(2, {-1});  // v1 is {1}.
//     auto v2 = Hilbert<1>::Offset(2, {+1});  // v1 is {2}.
//     auto v3 = Hilbert<1>::Offset(2, {+3});  // v1 is {3}.
template <std::size_t N, typename Int = int>
constexpr std::array<Int, N> OffsetV(std::size_t K,
                                     const std::array<Int, N>& center_v);

// The inverse operation of Offset() described above.  Example:
//     // Center the points of the 2nd iteration of a 1D Hilbert
//     // curve.
//     auto v0 = Hilbert<1>::Center(2, {0});  // v0 is {-3}.
//     auto v1 = Hilbert<1>::Center(2, {1});  // v1 is {-1}.
//     auto v2 = Hilbert<1>::Center(2, {2});  // v1 is {+1}.
//     auto v3 = Hilbert<1>::Center(2, {3});  // v1 is {+3}.
template <std::size_t N, typename Int = int>
constexpr std::array<Int, N> CenterV(std::size_t K,
                                     const std::array<Int, N>& offset_v);

template <typename Int>
constexpr /* static inline */ void Curve(std::size_t N,
                                         std::size_t K,
                                         Int* vs,
                                         const Int* pvs) {
  if (K == 0) {
    for (std::size_t j = 0; j < N; j++) {
      vs[j] = 0;
    }
    return;
  }

  size_t current = 0;
  for (std::size_t i = 0; i < (1U << N); i++) {
    for (const Int* p = pvs; p != pvs + N * (1 << N * (K - 1)); p += N) {
      std::size_t d = N - 1;
      if (i != 0 && i != (1U << N) - 1) {
        std::size_t j = (i - 1) >> 1;
        for (std::size_t bits = ~j & (j + 1); bits != 0; bits >>= 1) {
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
        v2[j] = p[order] * (sign ? 1 : -1) + offset;
      }
    }
  }
}

template <std::size_t N, std::size_t K, typename Int>
constexpr std::array<std::array<Int, N>, 1 << N * K> Curve() {
  std::array<std::array<Int, N>, 1 << N * K> curve{};
  Curve(K, curve.data());
  return curve;
}

template <std::size_t N, typename Int>
std::unique_ptr<std::array<Int, N>[]> Curve(std::size_t K) {
  std::unique_ptr<std::array<Int, N>[]> curve(
      new std::array<Int, N>[1 << N * K]);
  if (K == 0) {
    Curve<Int>(N, K, curve.get()[0].data(), nullptr);
    return curve;
  }
  auto prev = Curve<N, Int>(K - 1);
  Curve<Int>(N, K, curve.get()[0].data(), prev.get()[0].data());
  return curve;
}

template <std::size_t N, typename Int>
constexpr void Curve(std::size_t K, std::array<Int, N>* vs) {
  std::array<Int, N> v{};
  Curve<Int>(N, K, vs[0].data(), v.data());
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

template <std::size_t N, typename Int>
constexpr std::array<Int, N> IToV(std::size_t K, std::size_t i) {
  std::array<Int, N> v{};
  std::array<Int, N> orthant_v{};
  IToV<Int>(N, K, i, v.data(), orthant_v.data());
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

template <std::size_t N, typename Int>
constexpr std::size_t VToI(std::size_t K, const std::array<Int, N>& v) {
  std::array<Int, N> vec = v;
  std::array<Int, N> transformed{};
  std::array<Int, N> orthant_v{};
  return VToI<Int>(N, K, vec.data(), transformed.data(), orthant_v.data());
}

template <std::size_t N, typename Int>
constexpr std::array<Int, N> OffsetV(std::size_t K,
                                     const std::array<Int, N>& center_v) {
  std::array<Int, N> offset_v{};
  for (std::size_t i = 0; i < N; i++) {
    offset_v[i] = (center_v[i] + (1 << K)) >> 1;
  }
  return offset_v;
}

template <std::size_t N, typename Int>
constexpr std::array<Int, N> CenterV(std::size_t K,
                                     const std::array<Int, N>& offset_v) {
  std::array<Int, N> center_v{};
  for (std::size_t i = 0; i < N; i++) {
    center_v[i] = offset_v[i] * 2 - (1 << K) + 1;
  }
  return center_v;
}

}  // namespace hilbert

#endif  // HILBERT_HPP

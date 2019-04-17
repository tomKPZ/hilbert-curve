#pragma once

#include <array>
#include <cstdint>
#include <memory>
#include <type_traits>

template <std::size_t N, typename Int = int> class Hilbert {
public:
  using Vec = std::array<Int, N>;

  // Computes the K'th iteration of the Hilbert curve.  Returns an
  // array of 2^(N*K) Vec's.  Use this version of Curve() when K is a
  // small constant known at compile time.  Example:
  //     // Compute the 3rd iteration of a 2D Hilbert curve.
  //     constexpr auto curve = Hilbert<2>::Curve<3>();
  //     for (const auto& v : curve) { ... }
  template <std::size_t K> static constexpr std::array<Vec, 1 << N * K> Curve();

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
  static constexpr void Curve(Vec *vs, std::size_t K);

  // Returns the i'th vector of Curve(K).  Example:
  //     // Compute the 5th vector in the 3rd iteration of a 2D
  //     // Hilbert curve.
  //     auto v = Hilbert<2>::IToV(5, 3);  // v is {-7, -1}.
  static constexpr Vec IToV(std::size_t i, std::size_t K);

  // Returns the index that v would have in Curve(K).  Example:
  //     // Compute the index of the vector {-7, -1} in the 3rd
  //     // iteration of a 2D Hilbert curve.
  //     auto i = Hilbert<2>::VToI({-7, -1}, 3);  // i is 5.
  static constexpr std::size_t VToI(const Vec &v, std::size_t K);

  // Curve(), IToV(), and VToI() all operate on hilbert curves
  // centered at the origin with points separated a distance of 2.
  // For example, the 2nd iteration of a 1D hilbert curve would have
  // points at [{-3}, {-1}, {1}, {3}].  Sometimes this data is more
  // useful based at 0 with a distance 1 between points.  Example:
  //     // Offset the points of the 2nd iteration of a 1D Hilbert
  //     // curve.
  //     auto v0 = Hilbert<1>::Offset({-3}, 2);  // v0 is {0}.
  //     auto v1 = Hilbert<1>::Offset({-1}, 2);  // v1 is {1}.
  //     auto v2 = Hilbert<1>::Offset({+1}, 2);  // v1 is {2}.
  //     auto v3 = Hilbert<1>::Offset({+3}, 2);  // v1 is {3}.

  static constexpr Vec OffsetV(const Vec &center_v, std::size_t K);

  // The inverse operation of Offset() described above.  Example:
  //     // Center the points of the 2nd iteration of a 1D Hilbert
  //     // curve.
  //     auto v0 = Hilbert<1>::Center({0}, 2);  // v0 is {-3}.
  //     auto v1 = Hilbert<1>::Center({1}, 2);  // v1 is {-1}.
  //     auto v2 = Hilbert<1>::Center({2}, 2);  // v1 is {+1}.
  //     auto v3 = Hilbert<1>::Center({3}, 2);  // v1 is {+3}.
  static constexpr Vec CenterV(const Vec &offset_v, std::size_t K);

private:
  static_assert(std::is_signed_v<Int>);

  Hilbert() = delete;

  friend class Hilbert<N + 1, Int>;

  struct CompressedPermutationMatrix {
    std::array<std::size_t, N> order;
    std::array<bool, N> signs;
  };

  static constexpr bool IToVMap(std::size_t orthant, std::size_t j) {
    return (orthant + (1 << j)) & (1 << (j + 1));
  }

  static constexpr auto IToVTransform(std::size_t i) {
    CompressedPermutationMatrix ret{};

    std::size_t d = N - 1;
    if (i != 0 && i != (1 << N) - 1) {
      std::size_t j = (i - 1) >> 1;
      j = ~j & (j + 1);
      while (j != 0) {
        j >>= 1;
        d--;
      }
    }

    for (std::size_t j = 0; j < N; j++) {
      ret.order[j] = d;
      d = d == N - 1 ? 0 : d + 1;
    }

    for (std::size_t j = 0; j < N; j++) {
      std::size_t c = i + (j == 0 ? 1 : (1 << (j + 2)) - (1 << j) - 1);
      ret.signs[j] = c & (1 << (j + 1));
    }
    if constexpr (N > 0) {
      if (i == 0) {
        ret.signs[0] = 1;
      }
    }

    return ret;
  }

  static constexpr auto IToVTransforms() {
    std::array<CompressedPermutationMatrix, (1 << N)> ret{};
    for (std::size_t i = 0; i < (1 << N); i++) {
      ret[i] = IToVTransform(i);
    }
    return ret;
  }

  static constexpr auto VToITransforms() {
    std::array<CompressedPermutationMatrix, (1 << N)> ret{};
    for (std::size_t i = 0; i < (1 << N); i++) {
      for (std::size_t j = 0; j < N; j++) {
        ret[i].order[i_to_v_transforms[i].order[j]] = j;
        ret[i].signs[i_to_v_transforms[i].order[j]] =
            i_to_v_transforms[i].signs[j];
      }
    }
    return ret;
  }

  static constexpr auto i_to_v_transforms = IToVTransforms();
  static constexpr auto v_to_i_transforms = VToITransforms();
};

// static
template <std::size_t N, typename Int>
template <std::size_t K>
constexpr std::array<typename Hilbert<N, Int>::Vec, 1 << N * K>
Hilbert<N, Int>::Curve() {
  std::array<Vec, 1 << N * K> ret{};
  Curve(&ret[0], K);
  return ret;
}

// static
template <std::size_t N, typename Int>
std::unique_ptr<typename Hilbert<N, Int>::Vec[]>
Hilbert<N, Int>::Curve(std::size_t K) {
  std::unique_ptr<Hilbert<N, Int>::Vec[]> ret(
      new Hilbert<N, Int>::Vec[1 << (N * K)]);
  Curve(ret.get(), K);
  return ret;
}

// static
template <std::size_t N, typename Int>
constexpr void Hilbert<N, Int>::Curve(Vec *vs, std::size_t K) {
  if (K == 0) {
    vs[0] = Vec{};
    return;
  }

  Vec *prev_end = vs + (1 << N * (K - 1));
  Curve(vs, K - 1);
  size_t current = 1 << N * (K - 1);
  for (std::size_t i = 1; i < 1 << N; i++) {
    const auto &m = i_to_v_transforms[i];
    for (const Vec *p = vs; p != prev_end; p++) {
      Vec &v2 = vs[current++];
      for (std::size_t j = 0; j < N; j++) {
        v2[j] = (*p)[m.order[j]] * (m.signs[j] ? 1 : -1) +
                ((IToVMap(i, j) ? 1 : -1) << (K - 1));
      }
    }
  }

  current = 0;
  const auto &m = i_to_v_transforms[0];
  for (const Vec *p = vs; p != prev_end; p++) {
    const Vec v = *p;
    Vec &v2 = vs[current++];
    for (std::size_t j = 0; j < N; j++) {
      v2[j] = v[m.order[j]] * (m.signs[j] ? 1 : -1) +
              ((IToVMap(0, j) ? 1 : -1) << (K - 1));
    }
  }
}

// static
template <std::size_t N, typename Int>
constexpr typename Hilbert<N, Int>::Vec Hilbert<N, Int>::IToV(std::size_t i,
                                                              std::size_t K) {
  if (K == 0) {
    return {};
  }

  std::size_t orthant = i >> (N * (K - 1));
  std::size_t orthant_i = i & ((1 << N * (K - 1)) - 1);
  Vec orthant_v = IToV(orthant_i, K - 1);

  Vec v;
  const auto &m = i_to_v_transforms[orthant];
  for (std::size_t j = 0; j < N; j++) {
    v[j] = orthant_v[m.order[j]] * (m.signs[j] ? 1 : -1) +
           ((IToVMap(orthant, j) ? 1 : -1) << (K - 1));
  }
  return v;
}

// static
template <std::size_t N, typename Int>
constexpr std::size_t Hilbert<N, Int>::VToI(const Vec &v, std::size_t K) {
  if (K == 0) {
    return 0;
  }

  Vec transformed{};
  std::size_t orthant = 0;
  std::size_t parity = 0;
  for (std::size_t j = 0; j < N; j++) {
    parity ^= v[N - j - 1] > 0;
    orthant |= parity << (N - j - 1);
    transformed[j] = v[j] + ((v[j] > 0 ? -1 : 1) << (K - 1));
  }

  Vec orthant_v{};
  const auto &m = v_to_i_transforms[orthant];
  for (std::size_t j = 0; j < N; j++) {
    orthant_v[j] = transformed[m.order[j]] * (m.signs[j] ? 1 : -1);
  }

  std::size_t orthant_i = VToI(orthant_v, K - 1);
  return orthant_i + orthant * (1 << N * (K - 1));
}

// static
template <std::size_t N, typename Int>
constexpr typename Hilbert<N, Int>::Vec
Hilbert<N, Int>::OffsetV(const Vec &center_v, std::size_t K) {
  Vec offset_v{};
  for (std::size_t i = 0; i < N; i++) {
    offset_v[i] = (center_v[i] + (1 << K)) >> 1;
  }
  return offset_v;
}

// static
template <std::size_t N, typename Int>
constexpr typename Hilbert<N, Int>::Vec
Hilbert<N, Int>::CenterV(const Vec &offset_v, std::size_t K) {
  Vec center_v{};
  for (std::size_t i = 0; i < N; i++) {
    center_v[i] = offset_v[i] * 2 - (1 << K) + 1;
  }
  return center_v;
}

#pragma once

#include <array>
#include <cstdint>
#include <memory>
#include <type_traits>

template <std::size_t N, typename Int = int> class Hilbert {
public:
  using Vec = std::array<Int, N>;

  template <std::size_t K> static constexpr std::array<Vec, 1 << N * K> Curve();

  static std::unique_ptr<Vec[]> Curve(std::size_t K);

  static constexpr void Curve(Vec *vs, std::size_t K);

  static constexpr Vec IToV(std::size_t i, std::size_t K);

  static constexpr std::size_t VToI(const Vec &v, std::size_t K);

  static constexpr Vec Offset(const Vec &v, std::size_t K);

  static constexpr Vec Unoffset(const Vec &v, std::size_t K);

private:
  static_assert(std::is_signed_v<Int>);

  Hilbert() = delete;

  friend class Hilbert<N + 1, Int>;

  struct CompressedPermutationMatrix {
    std::array<std::size_t, N> order;
    std::array<bool, N> signs;
  };

  template <class InputIt, class OutputIt>
  static constexpr OutputIt Copy(InputIt first, InputIt last,
                                 OutputIt d_first) {
    while (first != last) {
      *d_first++ = *first++;
    }
    return d_first;
  }

  static constexpr std::array<std::size_t, 1 << N> BaseShape() {
    constexpr std::size_t TwoPowN = 1 << N;
    std::array<std::size_t, TwoPowN> ret{};
    if constexpr (N > 0) {
      auto np = Hilbert<N - 1, Int>::BaseShape();
      for (std::size_t i : {0, 1}) {
        for (std::size_t j = 0; j < TwoPowN / 2; j++) {
          auto &new_vec = ret[i * TwoPowN / 2 + j];
          new_vec = np[i ? TwoPowN / 2 - 1 - j : j];
          if (i) {
            new_vec |= 1 << (N - 1);
          }
        }
      }
    }
    return ret;
  }

  static constexpr std::array<std::size_t, (1 << N)>
  InvertPermutation(const std::array<std::size_t, (1 << N)> &arr) {
    std::array<std::size_t, 1 << N> ret{};
    for (std::size_t i = 0; i < 1 << N; i++) {
      ret[arr[i]] = i;
    }
    return ret;
  }

  static constexpr std::array<Vec, (1 << N) + 1> Transitions() {
    constexpr std::size_t TwoPowN = 1 << N;
    std::array<Vec, TwoPowN + 1> ret{};
    if constexpr (N > 0) {
      auto np = Hilbert<N - 1, Int>::Transitions();
      for (std::size_t j = 0; j < TwoPowN / 2; j++) {
        Vec &new_vec = ret[j];
        auto &old_vec = np[j];
        Copy(std::begin(old_vec), std::end(old_vec), std::begin(new_vec));
        new_vec[N - 1] = 0;
      }
      ret[TwoPowN / 2] = ret[TwoPowN / 2 - 1];
      ret[TwoPowN / 2][N - 1] = 1;
      for (std::size_t j = 0; j < TwoPowN / 2; j++) {
        Vec &new_vec = ret[TwoPowN / 2 + 1 + j];
        auto &old_vec = np[TwoPowN / 2 - 1 - j];
        Copy(std::begin(old_vec), std::end(old_vec), std::begin(new_vec));
        new_vec[N - 1] = 2;
      }
    }
    return ret;
  }

  static constexpr std::array<CompressedPermutationMatrix, (1 << N)>
  Transformations() {
    std::array<CompressedPermutationMatrix, (1 << N)> ret{};

    constexpr auto transitions = Transitions();
    for (std::size_t i = 0; i < (1 << N); i++) {
      std::size_t d = N;
      for (std::size_t j = 0; j < N; j++) {
        if (transitions[i + 1][j] - transitions[i][j] != 0) {
          d = N - j - 1;
          break;
        }
      }
      if constexpr (N > 0) {
        for (std::size_t j = 0; j < N; j++) {
          ret[i].order[j] = d;
          d = (d + 1) % N;
        }
      }
    }

    for (std::size_t j = 0; j < N; j++) {
      std::size_t c = j == 0 ? 1 : (1 << (j + 2)) - (1 << j) - 1;
      for (std::size_t i = 0; i < (1 << N); i++) {
        ret[i].signs[j] = c & (1 << (j + 1));
        c++;
      }
    }
    if constexpr (N > 0) {
      ret[0].signs[0] = 1;
    }

    return ret;
  }

  static constexpr std::array<CompressedPermutationMatrix, (1 << N)>
  InvertMatrices(const std::array<CompressedPermutationMatrix, (1 << N)> &ms) {
    std::array<CompressedPermutationMatrix, (1 << N)> ret{};
    for (std::size_t i = 0; i < (1 << N); i++) {
      for (std::size_t j = 0; j < N; j++) {
        ret[i].order[ms[i].order[j]] = j;
        ret[i].signs[ms[i].order[j]] = ms[i].signs[j];
      }
    }
    return ret;
  }

  static constexpr auto base_shape = BaseShape();
  static constexpr auto base_shape_inverse = InvertPermutation(base_shape);
  static constexpr auto transformations = Transformations();
  static constexpr auto transformations_inverse =
      InvertMatrices(transformations);
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
    const CompressedPermutationMatrix &m = transformations[i];
    for (const Vec *p = vs; p != prev_end; p++) {
      Vec &v2 = vs[current++];
      for (std::size_t j = 0; j < N; j++) {
        v2[j] = (*p)[m.order[j]] * (m.signs[j] ? 1 : -1) +
                ((base_shape[i] & (1 << j)) ? 1 : -1) * (1 << (K - 1));
      }
    }
  }

  current = 0;
  const CompressedPermutationMatrix &m = transformations[0];
  for (const Vec *p = vs; p != prev_end; p++) {
    const Vec v = *p;
    Vec &v2 = vs[current++];
    for (std::size_t j = 0; j < N; j++) {
      v2[j] = v[m.order[j]] * (m.signs[j] ? 1 : -1) +
              ((base_shape[0] & (1 << j)) ? 1 : -1) * (1 << (K - 1));
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

  std::size_t section = i / (1 << N * (K - 1));
  std::size_t section_i = i % (1 << N * (K - 1));
  Vec section_v = IToV(section_i, K - 1);

  Vec v;
  const CompressedPermutationMatrix &m = transformations[section];
  for (std::size_t j = 0; j < N; j++) {
    v[j] = section_v[m.order[j]] * (m.signs[j] ? 1 : -1) +
           ((base_shape[section] & (1 << j)) ? 1 : -1) * (1 << (K - 1));
  }
  return v;
}

// static
template <std::size_t N, typename Int>
constexpr std::size_t Hilbert<N, Int>::VToI(const Vec &v, std::size_t K) {
  if (K == 0) {
    return 0;
  }

  std::size_t coords = 0;
  Vec section_v{};
  for (std::size_t i = 0; i < N; i++) {
    if (v[i] > 0) {
      coords |= 1 << i;
    }
    if (v[i] > 0) {
      section_v[i] = v[i] - (1 << (K - 1));
    } else {
      section_v[i] = v[i] + (1 << (K - 1));
    }
  }

  std::size_t section = base_shape_inverse[coords];

  Vec transformed{};
  const CompressedPermutationMatrix &m = transformations_inverse[section];
  for (std::size_t i = 0; i < N; i++) {
    transformed[i] = section_v[m.order[i]] * (m.signs[i] ? 1 : -1);
  }

  std::size_t section_i = VToI(transformed, K - 1);
  return section_i + section * (1 << N * (K - 1));
}

// static
template <std::size_t N, typename Int>
constexpr typename Hilbert<N, Int>::Vec Hilbert<N, Int>::Offset(const Vec &v,
                                                                std::size_t K) {
  Vec ret{};
  for (std::size_t i = 0; i < N; i++) {
    ret[i] = (v[i] + (1 << K)) / 2;
  }
  return ret;
}

// static
template <std::size_t N, typename Int>
constexpr typename Hilbert<N, Int>::Vec
Hilbert<N, Int>::Unoffset(const Vec &v, std::size_t K) {
  Vec ret{};
  for (std::size_t i = 0; i < N; i++) {
    ret[i] = v[i] * 2 - (1 << K);
  }
  return ret;
}

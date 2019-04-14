#pragma once

#include <array>
#include <cstdint>
#include <vector>

template <std::size_t N> using Vector = std::array<int, N>;

template <std::size_t N> struct CompressedRotationMatrix {
  std::array<std::size_t, N> order;
  std::array<bool, N> signs;
};

template <class InputIt, class OutputIt>
constexpr OutputIt Copy(InputIt first, InputIt last, OutputIt d_first) {
  while (first != last) {
    *d_first++ = *first++;
  }
  return d_first;
}

template <class ForwardIt, class T>
constexpr void Fill(ForwardIt first, ForwardIt last, const T &value) {
  for (; first != last; ++first) {
    *first = value;
  }
}

template <std::size_t N>
constexpr Vector<N> operator+(const Vector<N> &v1, const Vector<N> &v2) {
  Vector<N> ret{};
  for (std::size_t i = 0; i < N; i++) {
    ret[i] = v1[i] + v2[i];
  }
  return ret;
}

template <std::size_t N>
constexpr Vector<N> operator-(const Vector<N> &v1, const Vector<N> &v2) {
  Vector<N> ret{};
  for (std::size_t i = 0; i < N; i++) {
    ret[i] = v1[i] - v2[i];
  }
  return ret;
}

template <std::size_t N>
constexpr Vector<N> operator*(const Vector<N> &v, int x) {
  Vector<N> ret{};
  for (std::size_t i = 0; i < N; i++) {
    ret[i] = v[i] * x;
  }
  return ret;
}

template <std::size_t N>
constexpr Vector<N> operator/(const Vector<N> &v, int x) {
  Vector<N> ret{};
  for (std::size_t i = 0; i < N; i++) {
    ret[i] = v[i] / x;
  }
  return ret;
}

template <std::size_t N>
constexpr Vector<N> operator*(const CompressedRotationMatrix<N> &m,
                              const Vector<N> &v) {
  Vector<N> ret{};
  for (std::size_t i = 0; i < N; i++) {
    ret[i] = (m.signs[i] ? 1 : -1) * v[m.order[i]];
  }
  return ret;
}

template <std::size_t N> constexpr Vector<N> Ones() {
  Vector<N> ret{};
  Fill(std::begin(ret), std::end(ret), 1);
  return ret;
}

constexpr std::size_t Pow(std::size_t b, std::size_t e) {
  std::size_t p = 1;
  for (std::size_t i = 0; i < e; i++) {
    p *= b;
  }
  return p;
}

template <std::size_t N> constexpr std::array<Vector<N>, (1 << N)> BaseShape() {
  if constexpr (N == 0) {
    return {{}};
  } else {
    constexpr std::size_t TwoPowN = 1 << N;
    std::array<Vector<N>, TwoPowN> ret{};
    std::array<Vector<N - 1>, TwoPowN / 2> np = BaseShape<N - 1>();
    for (int i : {0, 1}) {
      for (std::size_t j = 0; j < TwoPowN / 2; j++) {
        Vector<N> &new_vec = ret[i * TwoPowN / 2 + j];
        Vector<N - 1> &old_vec = np[i ? TwoPowN / 2 - 1 - j : j];
        Copy(std::begin(old_vec), std::end(old_vec), std::begin(new_vec));
        new_vec[N - 1] = i;
      }
    }
    return ret;
  }
}

template <std::size_t N>
constexpr std::array<Vector<N>, (1 << N) + 1> Transitions() {
  if constexpr (N == 0) {
    return {{}};
  } else {
    constexpr std::size_t TwoPowN = 1 << N;
    std::array<Vector<N>, TwoPowN + 1> ret{};
    std::array<Vector<N - 1>, TwoPowN / 2 + 1> np = Transitions<N - 1>();
    for (std::size_t j = 0; j < TwoPowN / 2; j++) {
      Vector<N> &new_vec = ret[j];
      Vector<N - 1> &old_vec = np[j];
      Copy(std::begin(old_vec), std::end(old_vec), std::begin(new_vec));
      new_vec[N - 1] = 0;
    }
    ret[TwoPowN / 2] = ret[TwoPowN / 2 - 1];
    ret[TwoPowN / 2][N - 1] = 1;
    for (std::size_t j = 0; j < TwoPowN / 2; j++) {
      Vector<N> &new_vec = ret[TwoPowN / 2 + 1 + j];
      Vector<N - 1> &old_vec = np[TwoPowN / 2 - 1 - j];
      Copy(std::begin(old_vec), std::end(old_vec), std::begin(new_vec));
      new_vec[N - 1] = 2;
    }
    return ret;
  }
}

template <std::size_t N>
constexpr std::array<CompressedRotationMatrix<N>, (1 << N)> Rotations() {
  std::array<CompressedRotationMatrix<N>, (1 << N)> ret{};

  constexpr auto transitions = Transitions<N>();
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

template <std::size_t N> struct HilbertConstants {
  constexpr static auto base_shape = BaseShape<N>();
  constexpr static auto rotations = Rotations<N>();
};

template <std::size_t N, std::size_t K>
constexpr std::array<Vector<N>, Pow(1 << N, K)> Hilbert() {
  std::array<Vector<N>, Pow(1 << N, K)> ret{};
  if constexpr (K > 0) {
    auto prev = Hilbert<N, K - 1>();
    std::size_t current = 0;
    for (std::size_t i = 0; i < (1 << N); i++) {
      for (const auto &v : prev) {
        constexpr Vector<N> offset = Ones<N>() * ((1 << (K - 1)) - 1);
        auto v2 = v * 2 - offset;
        v2 = HilbertConstants<N>::rotations[i] * v2;
        v2 = (v2 + offset) / 2;
        ret[current++] =
            v2 + HilbertConstants<N>::base_shape[i] * (1 << (K - 1));
      }
    }
  }
  return ret;
}

template <std::size_t N, std::size_t K> std::vector<Vector<N>> HilbertVector() {
  std::vector<Vector<N>> ret;
  if constexpr (K == 0) {
    ret.push_back(Vector<N>{});
  } else {
    auto prev = HilbertVector<N, K - 1>();
    for (std::size_t i = 0; i < (1 << N); i++) {
      for (const auto &v : prev) {
        constexpr Vector<N> offset = Ones<N>() * ((1 << (K - 1)) - 1);
        auto v2 = v * 2 - offset;
        v2 = HilbertConstants<N>::rotations[i] * v2;
        v2 = (v2 + offset) / 2;
        ret.push_back(v2 + HilbertConstants<N>::base_shape[i] * (1 << (K - 1)));
      }
    }
  }
  return ret;
}

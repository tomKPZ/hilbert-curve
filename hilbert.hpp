#pragma once

#include <array>
#include <cstdint>
#include <vector>

template <std::size_t N> class Hilbert {
public:
  using VecN = std::array<int, N>;

  template <std::size_t K>
  static constexpr std::array<VecN, 1 << N*K> Curve();

  template <std::size_t K> static std::vector<VecN> CurveVector();

private:
  Hilbert() = delete;

  friend class Hilbert<N + 1>;

  template <std::size_t M> using Vec = std::array<int, M>;

  struct CompressedRotationMatrix {
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

  static constexpr std::array<VecN, (1 << N)> BaseShape() {
    constexpr std::size_t TwoPowN = 1 << N;
    std::array<VecN, TwoPowN> ret{};
    if constexpr (N > 0) {
      std::array<Vec<N - 1>, TwoPowN / 2> np = Hilbert<N - 1>::BaseShape();
      for (int i : {0, 1}) {
        for (std::size_t j = 0; j < TwoPowN / 2; j++) {
          VecN &new_vec = ret[i * TwoPowN / 2 + j];
          Vec<N - 1> &old_vec = np[i ? TwoPowN / 2 - 1 - j : j];
          Copy(std::begin(old_vec), std::end(old_vec), std::begin(new_vec));
          new_vec[N - 1] = i;
        }
      }
    }
    return ret;
  }

  static constexpr std::array<VecN, (1 << N) + 1> Transitions() {
    constexpr std::size_t TwoPowN = 1 << N;
    std::array<VecN, TwoPowN + 1> ret{};
    if constexpr (N > 0) {
      std::array<Vec<N - 1>, TwoPowN / 2 + 1> np =
          Hilbert<N - 1>::Transitions();
      for (std::size_t j = 0; j < TwoPowN / 2; j++) {
        VecN &new_vec = ret[j];
        Vec<N - 1> &old_vec = np[j];
        Copy(std::begin(old_vec), std::end(old_vec), std::begin(new_vec));
        new_vec[N - 1] = 0;
      }
      ret[TwoPowN / 2] = ret[TwoPowN / 2 - 1];
      ret[TwoPowN / 2][N - 1] = 1;
      for (std::size_t j = 0; j < TwoPowN / 2; j++) {
        VecN &new_vec = ret[TwoPowN / 2 + 1 + j];
        Vec<N - 1> &old_vec = np[TwoPowN / 2 - 1 - j];
        Copy(std::begin(old_vec), std::end(old_vec), std::begin(new_vec));
        new_vec[N - 1] = 2;
      }
    }
    return ret;
  }

  static constexpr std::array<CompressedRotationMatrix, (1 << N)> Rotations() {
    std::array<CompressedRotationMatrix, (1 << N)> ret{};

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

  static constexpr auto base_shape = BaseShape();
  static constexpr auto rotations = Rotations();
};

template <std::size_t N>
template <std::size_t K>
constexpr std::array<typename Hilbert<N>::VecN, 1 << N*K>
Hilbert<N>::Curve() {
  std::array<VecN, 1 << N*K> ret{};
  if constexpr (K > 0) {
    auto prev = Curve<K - 1>();
    std::size_t current = 0;
    for (std::size_t i = 0; i < (1 << N); i++) {
      for (const auto &v : prev) {
	VecN& v2 = ret[current++];
	for (std::size_t j = 0; j < N; j++) {
	  constexpr int offset = (1 << (K - 1)) - 1;
	  const CompressedRotationMatrix& m = rotations[i];
	  int v2j = 2*v[m.order[j]] - offset;
	  v2j = (m.signs[j] ? 1 : -1) * v2j;
	  v2j = (v2j + offset)/2;
	  v2[j] = v2j + base_shape[i][j] * (1 << (K - 1));
	}
      }
    }
  }
  return ret;
}

template <std::size_t N>
template <std::size_t K>
std::vector<typename Hilbert<N>::VecN> Hilbert<N>::CurveVector() {
  std::vector<VecN> ret;
  if constexpr (K == 0) {
    ret.push_back(VecN{});
  } else {
    auto prev = CurveVector<K - 1>();
    for (std::size_t i = 0; i < (1 << N); i++) {
      for (const auto &v : prev) {
	VecN v2;
	for (std::size_t j = 0; j < N; j++) {
	  constexpr int offset = (1 << (K - 1)) - 1;
	  const CompressedRotationMatrix& m = rotations[i];
	  int v2j = 2*v[m.order[j]] - offset;
	  v2j = (m.signs[j] ? 1 : -1) * v2j;
	  v2j = (v2j + offset)/2;
	  v2[j] = v2j + base_shape[i][j] * (1 << (K - 1));
	}
	ret.push_back(v2);
      }
    }
  }
  return ret;
}

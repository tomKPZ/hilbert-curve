#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <vector>

template <std::size_t N> using Vector = std::array<int, N>;

template <class InputIt, class OutputIt>
constexpr OutputIt ConstexprCopy(InputIt first, InputIt last,
                                 OutputIt d_first) {
  while (first != last) {
    *d_first++ = *first++;
  }
  return d_first;
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
        ConstexprCopy(std::begin(old_vec), std::end(old_vec),
                      std::begin(new_vec));
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
      ConstexprCopy(std::begin(old_vec), std::end(old_vec),
                    std::begin(new_vec));
      new_vec[N - 1] = 0;
    }
    ret[TwoPowN / 2] = ret[TwoPowN / 2 - 1];
    ret[TwoPowN / 2][N - 1] = 1;
    for (std::size_t j = 0; j < TwoPowN / 2; j++) {
      Vector<N> &new_vec = ret[TwoPowN / 2 + 1 + j];
      Vector<N - 1> &old_vec = np[TwoPowN / 2 - 1 - j];
      ConstexprCopy(std::begin(old_vec), std::end(old_vec),
                    std::begin(new_vec));
      new_vec[N - 1] = 2;
    }
    return ret;
  }
}

int main() {
  std::cout << "BaseShape" << std::endl;
  constexpr auto vecs = BaseShape<3>();
  for (const auto &vec : vecs) {
    for (int i : vec) {
      std::cout << i << '\t';
    }
    std::cout << std::endl;
  }
  std::cout << "Transitions" << std::endl;
  constexpr auto transitions = Transitions<3>();
  for (const auto &vec : transitions) {
    for (int i : vec) {
      std::cout << i << '\t';
    }
    std::cout << std::endl;
  }
  return 0;
}

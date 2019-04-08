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

template <std::size_t N>
constexpr std::array<std::array<int, N>, 1 << N> BaseShape() {
  if constexpr (N == 0) {
    return {{}};
  } else {
    constexpr std::size_t RetSize = 1 << N;
    std::array<Vector<N>, RetSize> ret{};
    std::array<Vector<N - 1>, RetSize / 2> np = BaseShape<N - 1>();
    for (int i : {0, 1}) {
      for (std::size_t j = 0; j < RetSize / 2; j++) {
        Vector<N> &new_vec = ret[i * RetSize / 2 + j];
        Vector<N - 1> &old_vec = np[i ? RetSize / 2 - j - 1 : j];
        new_vec[0] = i;
        ConstexprCopy(std::begin(old_vec), std::end(old_vec),
                      std::begin(new_vec) + 1);
      }
    }
    return ret;
  }
}

int main() {
  constexpr auto vecs = BaseShape<3>();
  for (const auto &vec : vecs) {
    for (int i : vec) {
      std::cout << i << '\t';
    }
    std::cout << std::endl;
  }
  return 0;
}

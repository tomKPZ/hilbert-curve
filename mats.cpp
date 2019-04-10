#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <set>
#include <vector>

template <std::size_t N> using Vector = std::array<int, N>;
template <std::size_t N> using SquareMatrix = std::array<std::array<int, N>, N>;

template <std::size_t N>
constexpr SquareMatrix<N> operator*(const SquareMatrix<N> &l,
                                    const SquareMatrix<N> &r) {
  SquareMatrix<N> ret;
  for (std::size_t i = 0; i < N; i++) {
    for (std::size_t j = 0; j < N; j++) {
      int dot = 0;
      for (std::size_t k = 0; k < N; k++) {
        dot += l[i][k] * r[k][j];
      }
      ret[i][j] = dot;
    }
  }
  return ret;
}

template <class InputIt, class OutputIt>
constexpr OutputIt Copy(InputIt first, InputIt last, OutputIt d_first) {
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

template <std::size_t N> constexpr SquareMatrix<N> Identity() {
  SquareMatrix<N> ret{};
  for (std::size_t i = 0; i < N; i++) {
    ret[i][i] = 1;
  }
  return ret;
}

// O(N^2), but won't overflow.
template <std::size_t N, std::size_t K>
constexpr std::size_t BinomialCoefficient() {
  static_assert(N >= K);
  if constexpr (N == K || K == 0) {
    return 1;
  } else {
    return BinomialCoefficient<N - 1, K - 1>() +
           BinomialCoefficient<N - 1, K>();
  }
}

template <class ForwardIt, class T>
constexpr void Iota(ForwardIt first, ForwardIt last, T value) {
  while (first != last) {
    *first++ = value;
    ++value;
  }
}

template <std::size_t N, std::size_t K>
constexpr std::array<std::array<std::size_t, K>, BinomialCoefficient<N, K>()>
Choose() {
  static_assert(N >= K);
  constexpr std::size_t NChooseK = BinomialCoefficient<N, K>();
  std::array<std::array<std::size_t, K>, NChooseK> ret{};
  if constexpr (N == 0 || K == 0) {
    return ret;
  } else if constexpr (N == K) {
    Iota(std::begin(ret[0]), std::end(ret[0]), 0);
    return ret;
  } else {
    if constexpr (N - 1 >= K) {
      auto solutions = Choose<N - 1, K>();
      Copy(std::begin(solutions), std::end(solutions), std::begin(ret));
    }
    if constexpr (N - 1 >= K - 1) {
      auto solutions = Choose<N - 1, K - 1>();
      constexpr std::size_t offset = BinomialCoefficient<N - 1, K>();
      for (std::size_t i = 0; i < NChooseK - offset; i++) {
        Copy(std::begin(solutions[i]), std::end(solutions[i]),
             std::begin(ret[offset + i]));
        ret[offset + i][K - 1] = N - 1;
      }
    }
  }
  return ret;
}

template <std::size_t N>
constexpr std::array<SquareMatrix<N>, BinomialCoefficient<N, 2>()>
BaseRotationMatrices() {
  constexpr std::size_t RetSize = BinomialCoefficient<N, 2>();
  std::array<SquareMatrix<N>, RetSize> ret{};
  constexpr auto combinations = Choose<N, 2>();
  for (std::size_t i = 0; i < RetSize; i++) {
    const auto &combination = combinations[i];
    auto &matrix = ret[i];
    matrix = Identity<N>();
    matrix[combination[0]][combination[0]] = 0;
    matrix[combination[1]][combination[1]] = 0;
    // int d = ((N % 2) ^ (combination[0] % 2) ^ (combination[1] % 2)) ? 1 : -1;
    int d = 1;
    matrix[combination[0]][combination[1]] = d;
    matrix[combination[1]][combination[0]] = -d;
  }
  return ret;
}

void TestShapeAndTransitions() {
  constexpr std::size_t N = 5;
  std::cout << "BaseShape" << std::endl;
  constexpr auto vecs = BaseShape<N>();
  for (const auto &vec : vecs) {
    for (int i : vec) {
      std::cout << i << '\t';
    }
    std::cout << std::endl;
  }
  std::cout << "Transitions (absolute)" << std::endl;
  constexpr auto transitions = Transitions<N>();
  for (const auto &vec : transitions) {
    for (int i : vec) {
      std::cout << i << '\t';
    }
    std::cout << std::endl;
  }
  std::cout << "Transitions (relative)" << std::endl;
  for (std::size_t i = 0; i < std::size(transitions) - 1; i++) {
    for (std::size_t j = 0; j < N; j++) {
      std::cout << transitions[i + 1][j] - transitions[i][j] << '\t';
    }
    std::cout << std::endl;
  }
}

template <std::size_t B, std::size_t E> constexpr std::size_t Pow() {
  if constexpr (E == 0) {
    return 1;
  } else {
    return B * Pow<B, E - 1>();
  }
}

template <std::size_t N>
constexpr SquareMatrix<N> Pow(const SquareMatrix<N> &matrix,
                              std::size_t exponent) {
  SquareMatrix<N> ret = Identity<N>();
  for (std::size_t i = 0; i < exponent; i++) {
    ret = ret * matrix;
  }
  return ret;
}

template <std::size_t N, std::size_t M>
constexpr std::array<std::array<std::size_t, N>, Pow<M, N>()> Multipliers() {
  std::array<std::array<std::size_t, N>, Pow<M, N>()> ret{};
  if constexpr (N != 0 && M != 0) {
    auto solutions = Multipliers<N - 1, M>();
    std::size_t filled = 0;
    for (std::size_t i = 0; i < M; i++) {
      for (const auto &solution : solutions) {
        Copy(std::begin(solution), std::end(solution), std::begin(ret[filled]));
        ret[filled++][N - 1] = i;
      }
    }
  }
  return ret;
}

void TestChoose() {
  constexpr auto solutions = Choose<4, 2>();
  for (const auto &solution : solutions) {
    for (std::size_t x : solution) {
      std::cout << x << '\t';
    }
    std::cout << std::endl;
  }
}

template <std::size_t N> void PrintMatrix(const SquareMatrix<N> &matrix) {
  for (const auto &row : matrix) {
    for (int x : row) {
      std::cout << x << '\t';
    }
    std::cout << std::endl;
  }
}

void TestBaseRotationMatrices() {
  constexpr auto matrices = BaseRotationMatrices<4>();
  for (const auto &matrix : matrices) {
    PrintMatrix(matrix);
    std::cout << std::endl;
  }
}

template <std::size_t N> struct CompressedRotationMatrix {
  std::array<std::size_t, N> order;
  std::array<bool, N> signs;

  bool operator<(const CompressedRotationMatrix &other) const {
    if (order < other.order) {
      return true;
    }
    if (order > other.order) {
      return false;
    }
    return signs < other.signs;
  }

  bool operator==(const CompressedRotationMatrix &other) const {
    return order == other.order && signs == other.signs;
  }
};

template <std::size_t N>
void PrintCompressedRotationMatrix(const CompressedRotationMatrix<N> &matrix) {
  for (std::size_t i : matrix.order) {
    std::cout << i << '\t';
  }
  std::cout << std::endl;
  for (bool sign : matrix.signs) {
    std::cout << (sign ? '+' : '-') << '\t';
  }
  std::cout << std::endl;
}

template <std::size_t N>
constexpr CompressedRotationMatrix<N>
CompressRotationMatrix(const SquareMatrix<N> &matrix) {
  CompressedRotationMatrix<N> ret;

  for (std::size_t i = 0; i < N; i++) {
    for (std::size_t j = 0; j < N; j++) {
      if (matrix[i][j] != 0) {
        ret.order[j] = i;
        ret.signs[j] = matrix[i][j] > 0;
      }
    }
  }

  return ret;
}

template <typename T> constexpr void Swap(T &a, T &b) {
  T temp = a;
  a = b;
  b = temp;
}

// O(N^2).  Could make this O(N) if this is generalized to
// non-adjacent swaps.
template <std::size_t N>
constexpr std::size_t AdjacentSwaps(const std::array<std::size_t, N> &arr) {
  std::array<std::size_t, N> arr_{arr};
  std::size_t swaps = 0;
  for (std::size_t i = 0; i < N - 1; i++) {
    for (std::size_t j = 0; j < N - i - 1; j++) {
      if (arr_[j] > arr_[j + 1]) {
        Swap(arr_[j], arr_[j + 1]);
        swaps++;
      }
    }
  }
  return swaps;
}

template <std::size_t N>
constexpr bool Checksum(const CompressedRotationMatrix<N> &matrix) {
  bool checksum = (N % 2) ^ (AdjacentSwaps(matrix.order) % 2);
  for (std::size_t i = 0; i < N; i++) {
    checksum ^= matrix.signs[i];
  }
  return checksum;
}

template <typename BidirIt> constexpr bool Next(BidirIt first, BidirIt last) {
  if (first == last) {
    return false;
  }
  last--;
  *last = !*last;
  if (*last) {
    return true;
  }
  *last = false;
  return Next(first, last);
}

template <typename T, typename BidirIt>
constexpr bool Next(BidirIt first, BidirIt last, T max) {
  if (first == last) {
    return false;
  }
  last--;
  (*last)++;
  if (*last < max) {
    return true;
  }
  *last = 0;
  return Next(first, last, max);
}

// Slow reference implementation.
template <std::size_t N>
std::set<CompressedRotationMatrix<N>> SlowCompressedRotationMatrices() {
  constexpr auto base_matrices = BaseRotationMatrices<N>();
  std::set<CompressedRotationMatrix<N>> all_matrices;
  std::array<std::size_t, std::size(base_matrices)> multipliers{};
  do {
    SquareMatrix<N> matrix = Identity<N>();
    for (std::size_t i = 0; i < std::size(multipliers); i++) {
      matrix = matrix * Pow(base_matrices[i], multipliers[i]);
    }
    all_matrices.insert(CompressRotationMatrix(matrix));
  } while (Next(std::begin(multipliers), std::end(multipliers),
                static_cast<std::size_t>(4)));
  return all_matrices;
}

void TestMultipliers() {
  constexpr auto multipliers = Multipliers<4, 3>();
  for (const auto &mults : multipliers) {
    for (std::size_t x : mults) {
      std::cout << x << '\t';
    }
    std::cout << std::endl;
  }
}

template <class ForwardIt1, class ForwardIt2>
constexpr void IterSwap(ForwardIt1 a, ForwardIt2 b) {
  Swap(*a, *b);
}

template <class BidirIt> constexpr void Reverse(BidirIt first, BidirIt last) {
  while ((first != last) && (first != --last)) {
    IterSwap(first++, last);
  }
}

// Constexpr version copied from STL.
template <class BidirIt>
constexpr bool NextPermutation(BidirIt first, BidirIt last) {
  if (first == last) {
    return false;
  }
  BidirIt i = last;
  if (first == --i) {
    return false;
  }
  while (true) {
    BidirIt i1 = i;
    if (*--i < *i1) {
      BidirIt i2 = last;
      while (!(*i < *--i2)) {
      }
      IterSwap(i, i2);
      Reverse(i1, last);
      return true;
    }
    if (i == first) {
      Reverse(first, last);
      return false;
    }
  }
}

constexpr std::size_t Factorial(std::size_t x) {
  return x == 0 ? 1 : x * Factorial(x - 1);
}

constexpr std::size_t NumRotationMatrices(std::size_t N) {
  return Factorial(N) * (1 << (N - 1));
}

template <std::size_t N>
constexpr std::array<CompressedRotationMatrix<N>, NumRotationMatrices(N)>
CompressedRotationMatrices() {
  std::size_t i = 0;
  std::array<CompressedRotationMatrix<N>, NumRotationMatrices(N)> ret{};

  std::array<std::size_t, N> order{};
  std::array<bool, N - 1> signs{};
  Iota(std::begin(order), std::end(order), 0);
  do {
    do {
      ret[i].order = order;
      Copy(std::begin(signs), std::end(signs), std::begin(ret[i].signs));
      if (Checksum(ret[i])) {
        ret[i].signs[N - 1] = true;
      }
      i++;
    } while (Next(std::begin(signs), std::end(signs)));
  } while (NextPermutation(std::begin(order), std::end(order)));

  return ret;
}

void TestNext() {
  std::array<bool, 3> arr{};
  do {
    for (int j : arr) {
      std::cout << j << '\t';
    }
    std::cout << std::endl;
  } while (Next(std::begin(arr), std::end(arr)));
}

void TestCompressedRotationMatrices() {
  auto matrices = CompressedRotationMatrices<4>();
  for (const auto &matrix : matrices) {
    PrintCompressedRotationMatrix(matrix);
  }
}

void TestSlowCompressedRotationMatrices() {
  for (const auto &matrix : SlowCompressedRotationMatrices<4>()) {
    PrintCompressedRotationMatrix(matrix);
  }
}

void VerifyCompressedRotationMatrices() {
  constexpr std::size_t N = 5;
  const auto fast = CompressedRotationMatrices<N>();
  const auto slow = SlowCompressedRotationMatrices<N>();
  assert(std::equal(std::begin(fast), std::end(fast), std::begin(slow)));
}

int main() { return 0; }

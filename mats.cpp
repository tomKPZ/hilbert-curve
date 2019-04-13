#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <set>
#include <sstream>
#include <vector>

template <std::size_t N> using Vector = std::array<int, N>;
template <std::size_t N> using SquareMatrix = std::array<std::array<int, N>, N>;

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

template <class ForwardIt, class T>
constexpr void Fill(ForwardIt first, ForwardIt last, const T &value) {
  for (; first != last; ++first) {
    *first = value;
  }
}

template <std::size_t N> constexpr Vector<N> Ones() {
  Vector<N> ret{};
  Fill(std::begin(ret), std::end(ret), 1);
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

template <std::size_t N> constexpr Vector<N> operator-(const Vector<N> &v1) {
  Vector<N> ret{};
  for (std::size_t i = 0; i < N; i++) {
    ret[i] = -v1[i];
  }
  return ret;
}

template <std::size_t N>
constexpr SquareMatrix<N> operator*(const SquareMatrix<N> &l,
                                    const SquareMatrix<N> &r) {
  SquareMatrix<N> ret{};
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

template <std::size_t N>
constexpr Vector<N> operator*(const SquareMatrix<N> &m, const Vector<N> &v) {
  Vector<N> ret{};
  for (std::size_t i = 0; i < N; i++) {
    int dot = 0;
    for (std::size_t j = 0; j < N; j++) {
      dot += m[i][j] * v[j];
    }
    ret[i] = dot;
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
    matrix[combination[0]][combination[1]] = 1;
    matrix[combination[1]][combination[0]] = -1;
  }
  return ret;
}

template <std::size_t B, std::size_t E> constexpr std::size_t Pow() {
  std::size_t ret = 1;
  for (std::size_t i = 0; i < E; i++) {
    ret *= B;
  }
  return ret;
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

template <std::size_t N>
constexpr SquareMatrix<N>
UncompressRotationMatrix(const CompressedRotationMatrix<N> &matrix) {
  SquareMatrix<N> ret{};

  for (std::size_t i = 0; i < N; i++) {
    ret[i][matrix.order[i]] = matrix.signs[i] ? 1 : -1;
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

constexpr std::size_t NumRotationAndMirrorMatrices(std::size_t N) {
  return Factorial(N) * (1 << N);
}

template <std::size_t N>
constexpr std::array<CompressedRotationMatrix<N>,
                     NumRotationAndMirrorMatrices(N)>
RotationAndMirrorMatrices() {
  std::size_t i = 0;
  std::array<CompressedRotationMatrix<N>, NumRotationAndMirrorMatrices(N)>
      ret{};

  std::array<std::size_t, N> order{};
  std::array<bool, N> signs{};
  Iota(std::begin(order), std::end(order), 0);
  do {
    do {
      ret[i].order = order;
      ret[i].signs = signs;
      i++;
    } while (Next(std::begin(signs), std::end(signs)));
  } while (NextPermutation(std::begin(order), std::end(order)));

  return ret;
}

template <std::size_t N> void PrintNumberInBinary(std::size_t x) {
  for (std::size_t i = (1 << (N - 1)); i != 0; i >>= 1) {
    std::cout << ((i & x) ? 1 : 0) << '\t';
  }
  std::cout << std::endl;
}

template <std::size_t N> void PrintVector(const Vector<N> &vector) {
  for (int x : vector) {
    std::cout << x << '\t';
  }
  std::cout << std::endl;
}

template <std::size_t N> void PrintMatrix(const SquareMatrix<N> &matrix) {
  for (const auto &row : matrix) {
    for (int x : row) {
      std::cout << x << '\t';
    }
    std::cout << std::endl;
  }
}

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

void TestShapeAndTransitions() {
  constexpr std::size_t N = 3;
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

void TestChoose() {
  constexpr auto solutions = Choose<4, 2>();
  for (const auto &solution : solutions) {
    for (std::size_t x : solution) {
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

void TestMultipliers() {
  constexpr auto multipliers = Multipliers<4, 3>();
  for (const auto &mults : multipliers) {
    for (std::size_t x : mults) {
      std::cout << x << '\t';
    }
    std::cout << std::endl;
  }
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

void TestMatrixVectorMultiply() {
  SquareMatrix<3> m{{{7, 3, 19}, {5, 13, 17}, {23, 31, 29}}};
  Vector<3> v{{71, 61, 67}};
  Vector<3> product = m * v;
  PrintVector(product);
}

void VerifyMatrixVectorMultiply() {
  constexpr std::size_t N = 5;
  Vector<N> v{};
  for (std::size_t i = 0; i < N; i++) {
    v[i] = (i + 1) * (i + 1);
  }
  for (const auto &m : CompressedRotationMatrices<N>()) {
    assert(m * v == UncompressRotationMatrix(m) * v);
  }
}

void BruteForceRotations() {
  constexpr std::size_t N = 3;
  constexpr auto vecs = BaseShape<N>();
  constexpr auto transitions = Transitions<N>();
  // constexpr auto rotations = RotationAndMirrorMatrices<N>();
  constexpr auto rotations = CompressedRotationMatrices<N>();
  Vector<N> v1 = -Ones<N>();
  Vector<N> v2 = -Ones<N>();
  v2[N - 1] = 1;
  for (std::size_t i = 0; i < std::size(vecs); i++) {
    Vector<N> v1p = (transitions[i] - vecs[i]) * 2 - Ones<N>();
    Vector<N> v2p = (transitions[i + 1] - vecs[i]) * 2 - Ones<N>();
    PrintNumberInBinary<N>(i);
    PrintVector(vecs[i]);
    for (const auto &rotation : rotations) {
      if (rotation * v1p == v1 && rotation * v2p == v2) {
        PrintCompressedRotationMatrix(rotation);
      }
    }
    std::cout << std::endl;
  }
}

std::vector<std::vector<int>> Parse(const std::string &filename) {
  std::vector<std::vector<int>> ret;
  std::ifstream file(filename);
  std::string line;
  while (std::getline(file, line)) {
    std::vector<int> v;
    std::istringstream ss(line);
    std::istream s(ss.rdbuf());
    std::string num;
    while (std::getline(s, num, ',')) {
      v.push_back(std::atoi(num.c_str()));
    }
    ret.push_back(v);
  }
  return ret;
}

template <std::size_t N> void VerifyBaseShapes() {
  if constexpr (N == 0) {
    return;
  } else {
    VerifyBaseShapes<N - 1>();
    auto vs = Parse("hb/hb_1_" + std::to_string(N) + ".txt");
    auto base_shape = BaseShape<N>();
    assert(vs.size() == base_shape.size());
    for (std::size_t i = 0; i < vs.size(); i++) {
      assert(vs[i].size() == std::size(base_shape[i]));
    }
  }
}

template <std::size_t N>
std::array<CompressedRotationMatrix<N>, (1 << N)>
BruteForceRotationsFromDataFiles() {
  constexpr auto rotations = RotationAndMirrorMatrices<N>();
  constexpr auto base_shape = BaseShape<N>();
  auto base_shape_offset = base_shape;
  std::transform(std::begin(base_shape_offset), std::end(base_shape_offset),
                 std::begin(base_shape_offset),
                 [](Vector<N> &v) { return v * 2 - Ones<N>(); });
  auto vs = Parse("hb/hb_2_" + std::to_string(N) + ".txt");
  if constexpr (N == 2) {
    for (auto &v : vs) {
      std::reverse(v.begin(), v.end());
    }
  }
  std::array<CompressedRotationMatrix<N>, (1 << N)> ret{};
  for (std::size_t i = 0; i < std::size(base_shape); i++) {
    std::array<Vector<N>, base_shape.size()> rotated;
    for (std::size_t j = 0; j < std::size(base_shape); j++) {
      const auto &v = vs[std::size(base_shape) * i + j];
      for (std::size_t k = 0; k < N; k++) {
        rotated[j][k] = (v[k] - 2 * base_shape[i][k]) * 2 - 1;
      }
    }

    for (const auto &rotation : rotations) {
      bool satisfies = true;
      for (std::size_t j = 0; j < std::size(base_shape); j++) {
        if (rotation * base_shape_offset[j] != rotated[j]) {
          satisfies = false;
          break;
        }
      }
      if (satisfies) {
        ret[i] = rotation;
      }
    }
  }
  return ret;
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
    for (std::size_t j = 0; j < N; j++) {
      ret[i].order[j] = d;
      d = (d + 1) % N;
    }
  }

  for (std::size_t j = 0; j < N; j++) {
    std::size_t c = j == 0 ? 1 : (1 << (j + 2)) - (1 << j) - 1;
    for (std::size_t i = 0; i < (1 << N); i++) {
      ret[i].signs[j] = c & (1 << (j + 1));
      c++;
    }
  }
  ret[0].signs[0] = 1;

  return ret;
}

void TestRotations() {
  constexpr auto rotations = Rotations<4>();
  for (const auto &rotation : rotations) {
    PrintCompressedRotationMatrix(rotation);
  }
}

void VerifyRotations() {
  constexpr std::size_t N = 2;
  assert(BruteForceRotationsFromDataFiles<N>() == Rotations<N>());
  for (const auto &rotation : BruteForceRotationsFromDataFiles<N>()) {
    PrintCompressedRotationMatrix(rotation);
  }
  std::cout << std::endl;
  for (const auto &rotation : Rotations<N>()) {
    PrintCompressedRotationMatrix(rotation);
  }
}

constexpr std::size_t Pow(std::size_t b, std::size_t e) {
  std::size_t p = 1;
  for (std::size_t i = 0; i < e; i++) {
    p *= b;
  }
  return p;
}

template <std::size_t N, std::size_t K>
constexpr std::array<Vector<N>, Pow(1 << N, K)> Hilbert() {
  std::array<Vector<N>, Pow(1 << N, K)> ret{};
  if constexpr (K > 0) {
    constexpr auto base_shape = BaseShape<N>();
    constexpr auto rotations = Rotations<N>();
    auto prev = Hilbert<N, K - 1>();
    std::size_t current = 0;
    for (std::size_t i = 0; i < std::size(base_shape); i++) {
      for (const auto &v : prev) {
        constexpr Vector<N> offset = Ones<N>() * ((1 << (K - 1)) - 1);
        auto v2 = v * 2 - offset;
        v2 = rotations[i] * v2;
        v2 = (v2 + offset) / 2;
        ret[current++] = v2 + base_shape[i] * (1 << (K - 1));
      }
    }
  }
  return ret;
}

void TestHilbert() {
  constexpr auto hilbert = Hilbert<3, 3>();
  for (const auto &v : hilbert) {
    PrintVector(v);
  }
}

template <std::size_t N, std::size_t K> void VerifyHilbert() {
  constexpr auto hilbert = Hilbert<N, K>();
  auto vs =
      Parse("hb/hb_" + std::to_string(K) + "_" + std::to_string(N) + ".txt");
  if constexpr (N == 2) {
    for (auto &v : vs) {
      std::reverse(v.begin(), v.end());
    }
  }
  assert(hilbert.size() == vs.size());
  for (std::size_t i = 0; i < hilbert.size(); i++) {
    assert(vs[i].size() == hilbert[i].size());
    for (std::size_t j = 0; j < hilbert[i].size(); j++) {
      assert(vs[i][j] == hilbert[i][j]);
    }
  }
}

void VerifyHilbert() {
  // VerifyHilbert<1, 1>();
  // VerifyHilbert<1, 2>();
  // VerifyHilbert<1, 3>();
  // VerifyHilbert<1, 4>();
  // VerifyHilbert<1, 5>();
  // VerifyHilbert<1, 6>();
  // VerifyHilbert<1, 7>();
  // VerifyHilbert<1, 8>();
  // VerifyHilbert<1, 9>();
  // VerifyHilbert<1, 10>();

  // VerifyHilbert<2, 1>();
  // VerifyHilbert<2, 2>();
  // VerifyHilbert<2, 3>();
  // VerifyHilbert<2, 4>();
  // VerifyHilbert<2, 5>();

  // VerifyHilbert<3, 1>();
  // VerifyHilbert<3, 2>();
  // VerifyHilbert<3, 3>();

  // VerifyHilbert<4, 1>();
  // VerifyHilbert<4, 2>();

  // VerifyHilbert<5, 1>();
  // VerifyHilbert<5, 2>();

  // VerifyHilbert<6, 1>();

  // VerifyHilbert<7, 1>();

  // VerifyHilbert<8, 1>();

  // VerifyHilbert<9, 1>();
}

template <std::size_t N> struct HilbertConstants {
  constexpr static auto base_shape = BaseShape<N>();
  constexpr static auto rotations = Rotations<N>();
};

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

template <std::size_t N, std::size_t K> void VerifyHilbertVector() {
  std::cout << "Verifying " << N << "," << K << std::endl;
  auto hilbert = HilbertVector<N, K>();
  auto vs =
      Parse("hb/hb_" + std::to_string(K) + "_" + std::to_string(N) + ".txt");
  if constexpr (N == 2) {
    for (auto &v : vs) {
      std::reverse(v.begin(), v.end());
    }
  }
  assert(hilbert.size() == vs.size());
  for (std::size_t i = 0; i < hilbert.size(); i++) {
    assert(vs[i].size() == hilbert[i].size());
    for (std::size_t j = 0; j < hilbert[i].size(); j++) {
      assert(vs[i][j] == hilbert[i][j]);
    }
  }
}

void VerifyHilbertVector() {
  // VerifyHilbertVector<1, 1>();
  // VerifyHilbertVector<1, 10>();
  // VerifyHilbertVector<1, 2>();
  // VerifyHilbertVector<1, 3>();
  // VerifyHilbertVector<1, 4>();
  // VerifyHilbertVector<1, 5>();
  // VerifyHilbertVector<1, 6>();
  // VerifyHilbertVector<1, 7>();
  // VerifyHilbertVector<1, 8>();
  // VerifyHilbertVector<1, 9>();
  // VerifyHilbertVector<2, 1>();
  // VerifyHilbertVector<2, 10>();
  // VerifyHilbertVector<2, 2>();
  // VerifyHilbertVector<2, 3>();
  // VerifyHilbertVector<2, 4>();
  // VerifyHilbertVector<2, 5>();
  // VerifyHilbertVector<2, 6>();
  // VerifyHilbertVector<2, 7>();
  // VerifyHilbertVector<2, 8>();
  // VerifyHilbertVector<2, 9>();
  // VerifyHilbertVector<3, 1>();
  // VerifyHilbertVector<3, 2>();
  // VerifyHilbertVector<3, 3>();
  // VerifyHilbertVector<3, 4>();
  // VerifyHilbertVector<3, 5>();
  // VerifyHilbertVector<3, 6>();
  // VerifyHilbertVector<4, 1>();
  // VerifyHilbertVector<4, 2>();
  // VerifyHilbertVector<4, 3>();
  // VerifyHilbertVector<4, 4>();
  // VerifyHilbertVector<4, 5>();
  // VerifyHilbertVector<5, 1>();
  // VerifyHilbertVector<5, 2>();
  // VerifyHilbertVector<5, 3>();
  // VerifyHilbertVector<5, 4>();
  // VerifyHilbertVector<6, 1>();
  // VerifyHilbertVector<6, 2>();
  // VerifyHilbertVector<6, 3>();
  // VerifyHilbertVector<7, 1>();
  // VerifyHilbertVector<7, 2>();
  // VerifyHilbertVector<7, 3>();
  // VerifyHilbertVector<8, 1>();
  // VerifyHilbertVector<8, 2>();
  // VerifyHilbertVector<8, 3>();
  // VerifyHilbertVector<9, 1>();
  // VerifyHilbertVector<9, 2>();
  // VerifyHilbertVector<10, 1>();
  // VerifyHilbertVector<10, 2>();
  // VerifyHilbertVector<11, 1>();
  // VerifyHilbertVector<11, 2>();
}

int main() { return 0; }

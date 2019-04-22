#include "hilbert.hpp"

#include <memory>
#include <iostream>

// Dimension of curve used in examples.
constexpr int N = 2;

// Iteration of curve used in examples.
constexpr int K = 2;

// Compute and print the points of the curve.
void BasicCurveExample() {
  int curve[1 << N * K][N];
  Hilbert<>::Curve<N, K>(curve[0]);
  for (int i = 0; i < 1 << N * K; ++i) {
    for (int j = 0; j < N; ++j) {
      std::cout << curve[i][j] << '\t';
    }
    std::cout << std::endl;
  }
}

// Compute and point the points of the curve one-at-a-time and in
// reverse.
void IToVExample() {
  for (int i = (1 << N * K) - 1; i >0; --i) {
    int v[N];
    Hilbert<>::IToV<N, K>(i, v);
    for (int j = 0; j < N; ++j) {
      std::cout << v[j] << '\t';
    }
    std::cout << std::endl;
  }
}

// The same as BasicCurveExample(), except offsets the points before
// printing them.
void OffsetVExample() {
  int curve[1 << N * K][N];
  Hilbert<>::Curve<N, K>(curve[0]);
  for (int i = 0; i < 1 << N * K; ++i) {
    Hilbert<>::OffsetV<N, K>(curve[i], curve[i]);
    for (int j = 0; j < N; ++j) {
      std::cout << curve[i][j] << '\t';
    }
    std::cout << std::endl;
  }
}

// Compute the indices that points would have given their coordinates.
// Example limited to 2D for now.
void VToIAndCenterVExample() {
  for (int v0 = 0; v0 < 1 << K; v0++) {
    for (int v1 = 0; v1 < 1 << K; v1++) {
      int v[2] = {v0, v1};
      Hilbert<>::CenterV<2, K>(v, v);
      auto i = Hilbert<>::VToI<2, K>(v);
      std::cout << '(' << v0 << ",\t" << v1 << "):\t" << i << std::endl;
    }
  }
}

// Helper function for ConstexprCurveExample().
constexpr std::array<int, N << N * K> ConstexprCurveImpl() {
  std::array<int, N << N * K> curve{};
  Hilbert<>::Curve<N, K>(curve.data());
  return curve;
}

// Compute the points of the curve at compile-time and print them.
void ConstexprCurveExample() {
  constexpr auto curve = ConstexprCurveImpl();
  for (int i = 0; i < 1 << N * K; ++i) {
    for (int j = 0; j < N; ++j) {
      std::cout << curve[i * N + j] << '\t';
    }
    std::cout << std::endl;
  }
}

int main(void) {
  std::cout << "BasicCurveExample" << std::endl;
  BasicCurveExample();

  std::cout << "IToVExample" << std::endl;
  IToVExample();

  std::cout << "OffsetVExample" << std::endl;
  OffsetVExample();

  std::cout << "VToIAndCenterVExample" << std::endl;
  VToIAndCenterVExample();

  std::cout << "ConstexprCurveExample" << std::endl;
  ConstexprCurveExample();
  return 0;
}

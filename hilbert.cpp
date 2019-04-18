#include "hilbert.hpp"

#include <memory>

int main(void) {
  std::unique_ptr<int[]> curve{new int[2 * (1 << 2 * 14)]};
  Hilbert<>::Curve<2>(14, curve.get());
  return 0;
}

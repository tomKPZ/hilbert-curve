#include "hilbert.hpp"

#include <memory>

int main(void) {
  auto curve = std::make_unique<std::array<int, 2>[]>(2 << 2 * 14);
  Hilbert<>::Curve<2, 14>(curve[0].data());
  return 0;
}

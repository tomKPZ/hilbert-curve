#include "hilbert.hpp"

int main(void) {
  // Compute the 5th vector in the 3rd iteration of a 2D Hilbert curve.
  auto v = Hilbert<2>::IToV(5, 3);
}

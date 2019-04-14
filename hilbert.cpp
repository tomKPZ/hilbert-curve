#include "hilbert.hpp"

int main(void) {
  Hilbert<3>::VecN *buf = new Hilbert<3>::VecN[1 << (3 * 8)];
  Hilbert<3>::Curve(buf, 8);
  delete[] buf;
}

#include "hilbert.hpp"

int main(void) {
  Hilbert<2, short>::Curve(10);
  // constexpr auto vs = Hilbert<3>::Curve<3>();
}

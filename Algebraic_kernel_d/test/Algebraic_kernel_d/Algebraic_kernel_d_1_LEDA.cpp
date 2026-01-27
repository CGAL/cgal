// Test of Algebraic_kernel

#define CGAL_TEST_ALL_AK_VARIANTS 1

#include "Algebraic_kernel_d_1.h"

int main() {
#ifdef CGAL_HAS_LEDA_ARITHMETIC_KERNEL
  std::cout << " TEST AK1 USING LEDA " << std::endl;
  test_algebraic_kernel< CGAL::LEDA_arithmetic_kernel >();
#else
  std::cout << " NOTHING TESTED " << std::endl;
#endif

  return 0;
}

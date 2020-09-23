#include <iostream>

#include <CGAL/GMPXX_arithmetic_kernel.h>

#ifdef CGAL_USE_GMPXX

#include <CGAL/Test/_test_arithmetic_kernel.h>

int main() {
  std::cout << "TEST GMPXX_arithmetic_kernel" << std::endl;
  typedef CGAL::GMPXX_arithmetic_kernel AK;
  CGAL::test_arithmetic_kernel<AK>();
  return 0;
}

#else
int main() { return 0; }
#endif

#include <iostream>
#include <CGAL/basic.h>
#include <CGAL/GMP_arithmetic_kernel.h>

#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL

#include <CGAL/Test/_test_arithmetic_kernel.h>

int main() {
  std::cout << "TEST GMP_arithmetic_kernel" << std::endl; 
  typedef CGAL::GMP_arithmetic_kernel AK;
  CGAL::test_arithmetic_kernel<AK>();
  return 0;
}

#else
int main() { return 0; }
#endif

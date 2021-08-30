#include <iostream>

#include <CGAL/Arithmetic_kernel.h>

#if defined(CGAL_HAS_DEFAULT_ARITHMETIC_KERNEL)

#include <CGAL/Test/_test_arithmetic_kernel.h>

int main() {

  typedef CGAL::Arithmetic_kernel AK;
  CGAL::test_arithmetic_kernel<AK>();
  return 0;
}

#else
#warning CGAL has no default CGAL::Arithmetic kernel
int main() { return 0; }
#endif

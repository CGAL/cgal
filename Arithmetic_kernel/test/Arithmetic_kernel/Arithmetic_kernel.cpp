#include <iostream>
#include <CGAL/basic.h>
#include <CGAL/Arithmetic_kernel.h>

#if defined(CGAL_HAVE_DEFAULT_ARITHMETIC_KERNEL)

#include <CGAL/Test/_test_arithmetic_kernel.h>

int main() {

  typedef CGAL::Arithmetic_kernel AK;
  CGAL::test_arithmetic_kernel<AK>();
  return 0;
}

#else
int main() { return 0; }
#endif

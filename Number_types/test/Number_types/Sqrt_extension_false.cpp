
#include "Sqrt_extension.h"

int main(){ 
#ifdef CGAL_HAS_LEDA_ARITHMETIC_KERNEL
  sqrt_extension_test<CGAL::LEDA_arithmetic_kernel,CGAL::Tag_false>();
#endif // CGAL_HAS_LEDA_ARITHMETIC_KERNEL

#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL
  sqrt_extension_test<CGAL::CORE_arithmetic_kernel,CGAL::Tag_false>();
#endif // CGAL_HAS_CORE_ARITHMETIC_KERNEL
  test_nt_converter();
  return 0;
}



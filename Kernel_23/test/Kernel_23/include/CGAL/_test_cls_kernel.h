
#ifndef CGAL__TEST_KERNEL_C
#define CGAL__TEST_KERNEL_C

#include <CGAL/use.h>

template <class R>
bool
_test_kernel(const R& r)
{
  CGAL_USE(r);
  CGAL_USE_TYPE(typename R::RT);
  return true;
}

#endif


#ifndef CGAL__TEST_KERNEL_C
#define CGAL__TEST_KERNEL_C

#include <CGAL/use.h>

template <class R>
bool
_test_kernel(const R& r)
{
  typedef typename  R::RT RT;
  CGAL_USE_TYPE(RT);
  return true;
}

#endif

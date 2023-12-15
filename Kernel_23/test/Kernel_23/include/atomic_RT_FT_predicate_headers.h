#ifndef CGAL_KERNEL_23_TEST_ATOMIC_HEADERS_H
#define CGAL_KERNEL_23_TEST_ATOMIC_HEADERS_H

#define CGAL_NO_MPZF_DIVISION_OPERATOR

// These includes are there because this header is precompiled

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Mpzf.h>
#include <CGAL/tags.h>

#include <iostream>
#include <type_traits>

namespace CGAL {
namespace Kernel_23_tests {

struct Any
{
  template <class T>
  operator T();
};

} // namespace Kernel_23_tests
} // namespace CGAL

#endif // CGAL_KERNEL_23_TEST_ATOMIC_HEADERS_H

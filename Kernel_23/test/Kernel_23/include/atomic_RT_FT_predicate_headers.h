#ifndef CGAL_KERNEL_23_TEST_ATOMIC_HEADERS_H
#define CGAL_KERNEL_23_TEST_ATOMIC_HEADERS_H

#define CGAL_NO_MPZF_DIVISION_OPERATOR

#include <CGAL/Simple_cartesian.h>
#include <CGAL/tags.h>
#include <CGAL/Mpzf.h>

#include <iostream>

namespace CGAL {
namespace Kernel_23_tests {

struct Any {

  template <class T,
            typename = typename std::enable_if<!std::is_same<T, CGAL::RT_sufficient>::value>::type>
  operator T();
};

} // namespace Kernel_23_tests
} // namespace CGAL

#endif // CGAL_KERNEL_23_TEST_ATOMIC_HEADERS_H

#define CGAL_NO_DEPRECATION_WARNINGS 1 // because CGAL::copy_n is deprecated
#include <boost/config.hpp>

#if defined(BOOST_MSVC)
// avoid: warning C4996: 'CGAL::copy_n': was declared deprecated
#  pragma warning(disable:4996)
#endif
#include <CGAL/array.h>
#include <CGAL/tuple.h>
#include <CGAL/algorithm.h>

int main()
{
  CGAL::cpp0x::array<int, 3> arr;
  CGAL::cpp11::array<int, 3> arr2;

  CGAL::cpp0x::tuple<double, int> tuple;
  CGAL::cpp11::tuple<double, int> tuple2;

#ifndef CGAL_NO_DEPRECATED_CODE
  CGAL::copy_n(arr.begin(), 3, arr2.begin());
#endif // not CGAL_NO_DEPRECATED_CODE
  CGAL::cpp0x::copy_n(arr.begin(), 3, arr2.begin());
  CGAL::cpp11::copy_n(arr.begin(), 3, arr2.begin());
  
  CGAL::cpp0x::prev(arr.end());
  CGAL::cpp11::prev(arr.end());
  CGAL::cpp0x::next(arr.begin());
  CGAL::cpp11::next(arr.begin());
  return 0;
}

#include <CGAL/internal/disable_deprecation_warnings_and_errors.h> // because CGAL::copy_n is deprecated

#include <boost/config.hpp>

#if defined(BOOST_MSVC)
// avoid: warning C4996: 'CGAL::copy_n': was declared deprecated
#  pragma warning(disable:4996)
#endif
#include <CGAL/array.h>
#include <CGAL/tuple.h>
#include <CGAL/algorithm.h>
#include <CGAL/use.h>

int main()
{
  CGAL::cpp0x::array<int, 3> arr;
  std::array<int, 3> arr2;

  CGAL::cpp0x::tuple<double, int> tuple;
  std::tuple<double, int> tuple2;

  CGAL_USE(tuple);
  CGAL_USE(tuple2);

  CGAL::copy_n(arr.begin(), 3, arr2.begin());

  CGAL::cpp0x::copy_n(arr.begin(), 3, arr2.begin());
  std::copy_n(arr.begin(), 3, arr2.begin());
  
  CGAL::cpp0x::prev(arr.end());
  std::prev(arr.end());
  CGAL::cpp0x::next(arr.begin());
  std::next(arr.begin());
  return 0;
}

#include <CGAL/array.h>
#include <CGAL/tuple.h>
#include <CGAL/algorithm.h>

int main()
{
  CGAL::cpp0x::array<3, int> arr;
  CGAL::cpp11::array<3, int> arr2;

  CGAL::cpp0x::tuple<double, int> tuple;
  CGAL::cpp11::tuple<double, int> tuple2;

  CGAL::copy_n(arr.begin(), arr.end(), arr2.begin());
  
  return 0;
}



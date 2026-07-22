#include <CGAL/Exact_integer.h>
#include <CGAL/Simple_homogeneous.h>
#include <CGAL/Min_circle_2.h>
#include <CGAL/Min_circle_2_traits_2.h>

#include <array>
#include <iostream>

// typedefs
typedef  CGAL::Exact_integer             RT;
typedef  CGAL::Simple_homogeneous<RT>    K;
typedef  CGAL::Min_circle_2_traits_2<K>  Traits;
typedef  CGAL::Min_circle_2<Traits>      Min_circle;
typedef  K::Point_2                      Point;

int
main( int, char**)
{
  const int n = 100;
  std::array<Point, n> P;

  for ( int i = 0; i < n; ++i){
    P.at(i) = Point( (i%2 == 0 ? i : -i), 0, 1);
    // (0,0), (-1,0), (2,0), (-3,0), ...
  }

  Min_circle  mc1( P.begin(), P.end(), false);    // very slow
  Min_circle  mc2( P.begin(), P.end(), true);     // fast

  CGAL::IO::set_pretty_mode( std::cout);
  std::cout << mc2;

  return 0;
}

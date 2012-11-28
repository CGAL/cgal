#include <iostream>
#include <iterator>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_yz_3.h>
#include <CGAL/convex_hull_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K3;
typedef CGAL::Projection_traits_yz_3<K3> K;
typedef K::Point_2 Point_2;

int main()
{
  std::istream_iterator< Point_2 >  input_begin( std::cin );
  std::istream_iterator< Point_2 >  input_end;
  std::ostream_iterator< Point_2 >  output( std::cout, "\n" );
  CGAL::convex_hull_2( input_begin, input_end, output, K() );
  return 0;
}

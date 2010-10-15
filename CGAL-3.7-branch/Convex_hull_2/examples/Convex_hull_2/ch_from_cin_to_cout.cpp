#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/ch_graham_andrew.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;

int main()
{
  CGAL::set_ascii_mode(std::cin);
  CGAL::set_ascii_mode(std::cout);
  std::istream_iterator< Point_2 >  in_start( std::cin );
  std::istream_iterator< Point_2 >  in_end;
  std::ostream_iterator< Point_2 >  out( std::cout, "\n" );
  CGAL::ch_graham_andrew( in_start, in_end, out );
  return 0;
}

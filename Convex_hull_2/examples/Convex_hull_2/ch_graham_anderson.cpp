#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/ch_graham_andrew.h>
#include <vector>
#include <algorithm>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2      Point_2;

template <class InputIterator, class OutputIterator, class Traits>
OutputIterator
ch_graham_anderson( InputIterator  first, InputIterator  beyond,
                    OutputIterator result, const Traits&  ch_traits)
{
  typedef typename Traits::Point_2            Point_2;
  typedef typename Traits::Less_xy_2          Less_xy_2;
  typedef typename Traits::Less_rotate_ccw_2  Less_rotate_ccw_2;

  if (first == beyond) return result;
  std::vector< Point_2 >  V (first, beyond);
  typename std::vector< Point_2 >::iterator it =
               std::min_element(V.begin(), V.end(), Less_xy_2());
  const Point_2 p = *it;
  std::sort( V.begin(), V.end(), [&p](const Point_2& p1, const Point_2& p2){return Less_rotate_ccw_2()(p, p1, p2);} );
  if ( *(V.begin()) != *(V.rbegin()) )
  {
    result = CGAL::ch_graham_andrew_scan( V.begin(), V.end(), result, ch_traits);
  }
  // add the last point of the sequence that is
  // not added by ch_graham_andrew_scan
  *result++ = *(V.rbegin());
  return result;

}

int main()
{
  CGAL::IO::set_ascii_mode(std::cin);
  CGAL::IO::set_ascii_mode(std::cout);
  std::istream_iterator< Point_2 >  in_start( std::cin );
  std::istream_iterator< Point_2 >  in_end;
  std::ostream_iterator< Point_2 >  out( std::cout, "\n" );
  ch_graham_anderson(in_start, in_end, out, K());
  return 0;
}

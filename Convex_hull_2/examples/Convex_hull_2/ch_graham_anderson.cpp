#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/ch_graham_andrew.h>
#include <vector>
#include <algorithm>
#include <boost/bind.hpp>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2      Point_2;

template <class InputIterator, class OutputIterator, class Traits>
OutputIterator
ch_graham_anderson( InputIterator  first, InputIterator  beyond,
                    OutputIterator result, const Traits&  ch_traits)
{
  using namespace boost;

  typedef typename Traits::Point_2            Point_2;
  typedef typename Traits::Less_xy_2          Less_xy_2;
  typedef typename Traits::Less_rotate_ccw_2  Less_rotate_ccw_2;

  if (first == beyond) return result;
  std::vector< Point_2 >  V (first, beyond);
  typename std::vector< Point_2 >::iterator it =
               std::min_element(V.begin(), V.end(), Less_xy_2());
  std::sort( V.begin(), V.end(), boost::bind(Less_rotate_ccw_2(), *it, _1, _2) );
  if ( *(V.begin()) == *(V.rbegin()) )
  {
      *result = *(V.begin());  ++result;
      return result;
  }
  return CGAL::ch_graham_andrew_scan( V.begin(), V.end(), result, ch_traits);
}

int main()
{
  CGAL::set_ascii_mode(std::cin);
  CGAL::set_ascii_mode(std::cout);
  std::istream_iterator< Point_2 >  in_start( std::cin );
  std::istream_iterator< Point_2 >  in_end;
  std::ostream_iterator< Point_2 >  out( std::cout, "\n" );
  ch_graham_anderson(in_start, in_end, out, K());
  return 0;
}

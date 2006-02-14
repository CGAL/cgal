// file: examples/Convex_hull_2/ch_graham_anderson.C

#include <CGAL/Cartesian.h>
#include <CGAL/ch_graham_andrew.h>
#include <CGAL/functional.h>
#include <vector>
#include <algorithm>

typedef   CGAL::Point_2<CGAL::Cartesian<double> >        Point_2;

template <class InputIterator, class OutputIterator, class Traits>
OutputIterator
ch_graham_anderson( InputIterator  first, InputIterator  beyond,
                    OutputIterator result, const Traits&  ch_traits)
{
  typedef typename Traits::Less_xy_2          Less_xy_2;
  typedef typename Traits::Point_2            Point_2;
  typedef typename Traits::Less_rotate_ccw_2  Less_rotate_ccw_2;

  if (first == beyond) return result;
  std::vector< Point_2 >  V;
  std::copy( first, beyond, std::back_inserter(V) );
  typename std::vector< Point_2 >::iterator it =
               std::min_element(V.begin(), V.end(), Less_xy_2());
  std::sort( V.begin(), V.end(), CGAL::bind_1(Less_rotate_ccw_2(), *it) );
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
  ch_graham_anderson(in_start, in_end, out, CGAL::Cartesian<double>());
  return 0;
}

/*  
There must be a bug to fix in the sweepline code if this package.
The following data set hangs with the vertices in this particular order.
As the first thing the random_polygon_2 functions does
is a random_shuffle, we have put it in an #ifndef  for this testcase
 */

#define CGAL_DONT_SHUFFLE_IN_RANDOM_POLYGON_2 1

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/algorithm.h>



#include <fstream>
#include <list>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                                          Point_2;
typedef std::list<Point_2>                                  Container;
typedef CGAL::Polygon_2<K, Container>                       Polygon_2;


int main( )
{
   Polygon_2            polygon;
   std::list<Point_2>   point_set;

   std::string input = "-5 -7 46 59 55 -50 27 -93 46 -62 55 17 -71 0 -27 -41 -86 -73 7 -85 -71 93 41 -47 78 0 -13 52 -73 -25 -17 -11 56 -61 71 -69 -44 59 -9 5 -64 -96 -66 -45 39 54 -93 -22 5 -21 -32 32 -26 77 -99 -2 62 -73 -35 7 -88 60 39 56 2 3 -5 -50 3 68 82 16 -64 -38 -2 -16";


   std::stringstream input_stream(input);

   std::istream_iterator< Point_2 > in(input_stream);
   CGAL::copy_n_unique(in, 38,
                       std::back_inserter(point_set));

   std::ostream_iterator< Point_2 >  out( std::cout, "\n " );
  
   CGAL::random_polygon_2(point_set.size(), std::back_inserter(polygon),
                          point_set.begin());

   std::cout << "done" << std::endl;
   return 0;
}

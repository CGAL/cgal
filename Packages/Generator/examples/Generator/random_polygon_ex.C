//
// file : examples/Generator/random_polygon_ex.C
//
#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <fstream>

#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/Random.h>
#include <CGAL/copy_n.h>
#include <list>

typedef double                                            NT;
typedef CGAL::Cartesian<NT>                               K;
typedef CGAL::Polygon_traits_2<K>                         Traits;
typedef Traits::Point_2                                   Point_2;
typedef std::list<Point_2>                                Container;
typedef CGAL::Polygon_2<Traits, Container>                Polygon_2;
typedef CGAL::Creator_uniform_2<int, Point_2>             Creator;
typedef CGAL::Random_points_in_square_2<Point_2, Creator> Point_generator;

const double RADIUS = 100;
const int MAX_POLY_SIZE = 100;

int main( )
{
   Polygon_2            polygon;
   std::list<Point_2>   point_set;
   CGAL::Random         rand;

   int size = rand.get_int(4, MAX_POLY_SIZE);

   // copy size points from the generator, eliminating duplicates, so the
   // polygon will have <= size vertices
   CGAL::copy_n_unique(Point_generator(RADIUS), size,
                       std::back_inserter(point_set));
   CGAL::random_polygon_2(point_set.size(), std::back_inserter(polygon), 
                          point_set.begin());
   cout << polygon << endl;
   return 0;
}


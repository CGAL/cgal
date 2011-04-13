#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/all_furthest_neighbors_2.h>
#include <CGAL/IO/Ostream_iterator.h>
#include <iostream>
#include <vector>

using namespace std;
using CGAL::random_convex_set_2;
using CGAL::all_furthest_neighbors_2;

typedef double                                   FT;
typedef CGAL::Cartesian< FT >                    R;
typedef CGAL::Point_2< R >                       Point;
typedef CGAL::Polygon_traits_2< R >              P_traits;
typedef vector< Point >                          Point_cont;
typedef CGAL::Polygon_2< P_traits, Point_cont >  Polygon;
typedef CGAL::Creator_uniform_2< FT, Point >     Creator;
typedef CGAL::Random_points_in_square_2< Point, Creator >
  Point_generator;
typedef CGAL::Ostream_iterator< int, ostream >   Oiterator;

int
main()
{
  // generate random convex polygon:
  Polygon p;
  random_convex_set_2( 10, back_inserter( p), Point_generator( 1));

  // compute all furthest neighbors:
  all_furthest_neighbors_2(
    p.vertices_begin(),
    p.vertices_end(),
    Oiterator( cout));
  cout << endl;

  return 0;
} 

#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/all_furthest_neighbors_2.h>
#include <CGAL/IO/Ostream_iterator.h>
#include <iostream.h>
#include <vector.h>

typedef double                                 FT;
typedef CGAL_Cartesian< FT >                   R;
typedef CGAL_Point_2< R >                      Point_2;
typedef CGAL_Polygon_traits_2< R >             P_traits;
typedef vector< Point_2 >                      Point_cont;
typedef CGAL_Polygon_2< P_traits, Point_cont > Polygon_2;
typedef CGAL_Random_points_in_square_2<
  Point_2,
  CGAL_Creator_uniform_2< FT, Point_2 > >
Point_generator;
typedef CGAL_Ostream_iterator< int, ostream >  Ostream_iterator;

int
main()
{
  // generate random convex polygon:
  Polygon_2 p;
  CGAL_random_convex_set_2( 10,
                            back_inserter( p),
                            Point_generator( 1));

  // compute all furthest neighbors:
  CGAL_all_furthest_neighbors(
    p.vertices_begin(),
    p.vertices_end(),
    Ostream_iterator( cout));
  cout << endl;

} 

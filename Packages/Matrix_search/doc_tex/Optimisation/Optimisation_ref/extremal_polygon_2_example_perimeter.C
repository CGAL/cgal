#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/extremal_polygon_2.h>
#include <iostream.h>
#include <vector.h>

int main() {

  typedef double                            FT;
  typedef CGAL_Cartesian< FT >              R;
  typedef CGAL_Point_2< R >                 Point_2;
  typedef CGAL_Polygon_traits_2< R >        P_traits;
  typedef vector< Point_2 >                 Cont;
  typedef CGAL_Polygon_2< P_traits, Cont >  Polygon_2;
  typedef CGAL_Random_points_in_square_2<
    Point_2,
    CGAL_Creator_uniform_2< FT, Point_2 > >
  Point_generator;

  Polygon_2 p;
  int number_of_points( 10);
  int k( 5);

  CGAL_random_convex_set_2( number_of_points,
                            back_inserter( p),
                            Point_generator( 1));
  cout << "Generated Polygon:\n" << p << endl;

  Polygon_2 k_gon;
  CGAL_maximum_perimeter_inscribed_k_gon(
    p.vertices_begin(),
    p.vertices_end(),
    k,
    back_inserter( k_gon));
  cout << "Maximum perimeter " << k << "-gon:\n" << k_gon << endl;

} 

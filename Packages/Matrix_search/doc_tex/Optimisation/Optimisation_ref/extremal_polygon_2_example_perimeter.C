#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/extremal_polygon_2.h>
#include <iostream>
#include <vector>

using namespace std;
using CGAL::random_convex_set_2;
using CGAL::maximum_perimeter_inscribed_k_gon_2;

typedef double                                FT;
typedef CGAL::Cartesian< FT >                 R;
typedef CGAL::Point_2< R >                    Point;
typedef CGAL::Polygon_traits_2< R >           P_traits;
typedef vector< Point >                       Cont;
typedef CGAL::Polygon_2< P_traits, Cont >     Polygon;
typedef CGAL::Creator_uniform_2< FT, Point >  Creator;
typedef CGAL::Random_points_in_square_2< Point, Creator >
  Point_generator;

int main() {

  Polygon p;
  int number_of_points( 10);
  int k( 5);

  random_convex_set_2( number_of_points,
                       back_inserter( p),
                       Point_generator( 1));
  cout << "Generated Polygon:\n" << p << endl;

  Polygon k_gon;
  maximum_perimeter_inscribed_k_gon_2(
    p.vertices_begin(),
    p.vertices_end(),
    k,
    back_inserter( k_gon));
  cout << "Maximum perimeter " << k << "-gon:\n" << k_gon << endl;

  return 0;
} 

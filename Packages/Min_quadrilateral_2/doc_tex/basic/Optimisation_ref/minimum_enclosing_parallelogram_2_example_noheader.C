#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/min_quadrilateral_2.h>
#include <vector>
#include <iostream>

using CGAL::Random_points_in_square_2;
using CGAL::random_convex_set_2;
using CGAL::min_parallelogram_2;
using std::back_inserter;
using std::cout;
using std::endl;

typedef CGAL::Cartesian< double >                      R;
typedef R::Point_2                                     Point_2;
typedef R::Line_2                                      Line_2;
typedef CGAL::Polygon_traits_2< R >                    P_traits;
typedef std::vector< Point_2 >                         Cont;
typedef CGAL::Polygon_2< P_traits, Cont >              Polygon_2;
typedef CGAL::Creator_uniform_2< double, Point_2 >     Creator;
typedef Random_points_in_square_2< Point_2, Creator >  Point_generator;

int main()
{
  // build a random convex 20-gon p
  Polygon_2 p;
  random_convex_set_2(20, back_inserter(p), Point_generator(1.0));
  cout << p << endl;

  // compute the minimal enclosing parallelogram p_m of p
  Polygon_2 p_m;
  min_parallelogram_2(
    p.vertices_begin(), p.vertices_end(), back_inserter(p_m));
  cout << p_m << endl;

  return 0;
} 

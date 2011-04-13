// file: demo/Generator/random_poly_manual_demo.C
// ----------------------------------------------
// program generting a random simple polygon 

#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/IO/Window_stream.h>

typedef CGAL::Cartesian< double >                  K;
typedef K::Point_2                                 Point_2;
typedef std::list<Point_2>                         Container;
typedef CGAL::Polygon_traits_2<K>                  Traits;
typedef CGAL::Polygon_2<Traits, Container>         Polygon_2;
typedef CGAL::Random_points_in_square_2< Point_2 > Point_generator;

int main() {

  Polygon_2 polygon;

  // create 50-gon and write it into a window:
  CGAL::Window_stream W;
  W.init(-0.55, 0.55, -0.55);
  W.display();
  CGAL::cgalize(W);
  CGAL::random_polygon_2(50, std::back_inserter(polygon), Point_generator(0.5));
  W << polygon;

  W.read_mouse();         // wait for mouse-click
  return 0;
}

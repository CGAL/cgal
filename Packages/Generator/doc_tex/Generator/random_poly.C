#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/IO/Window_stream.h>
#include <CGAL/IO/Ostream_iterator.h>

typedef CGAL::Cartesian< double >                           R;
typedef CGAL::Point_2< R >                                  Point_2;
typedef std::list<Point_2>                                  Container;
typedef CGAL::Polygon_traits_2<R>                           Traits;
typedef CGAL::Polygon_2<Traits, Container>                  Polygon_2;
typedef CGAL::Creator_uniform_2< double, Point_2 >          Creator;
typedef CGAL::Random_points_in_square_2< Point_2, Creator > Point_generator;

int main() {

  Polygon_2 polygon;

  // create 50-gon and write it into a window:
  CGAL::Window_stream W;
  W.init(-0.55, 0.55, -0.55);
  W.display();
  CGAL::cgalize(W);
  CGAL::random_polygon_2(50, std::back_inserter(polygon), Point_generator(0.5));
  W << polygon;

  // wait for mouse-click:
  W.read_mouse();
}

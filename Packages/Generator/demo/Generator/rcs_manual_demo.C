//
// file: demo/Generator/rcs_manual_demo.C
// --------------------------------------
// generator of random convex set
//
#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/IO/Window_stream.h>
#include <CGAL/IO/Ostream_iterator.h>

typedef CGAL::Cartesian< double >   K;
typedef K::Point_2                  Point_2;
typedef CGAL::Random_points_in_square_2< 
     Point_2,
     CGAL::Creator_uniform_2< double, Point_2 > >
  Point_generator;

int main() {

  // create 500-gon and write it into a window:
  CGAL::Window_stream W;
  W.init( -0.55, 0.55, -0.55);
  W.display();
  CGAL::cgalize(W);
  CGAL::random_convex_set_2(
            500, 
            CGAL::Ostream_iterator< Point_2, CGAL::Window_stream >( W),
            Point_generator( 0.5));

  W.read_mouse();      // wait for mouse-click:
}

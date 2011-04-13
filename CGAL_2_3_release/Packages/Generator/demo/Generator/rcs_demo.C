// file : demo/Generator/rcs_demo.C
// --------------------------------
// program that generates a random convex set

#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <CGAL/IO/Window_stream.h>

typedef CGAL::Cartesian< double >                          K;
typedef K::Point_2                                         Point;
typedef CGAL::Polygon_traits_2< K >                        P_traits;
typedef std::vector< Point >                               Cont;
typedef CGAL::Polygon_2< P_traits, Cont >                  Polygon_2;
typedef CGAL::Creator_uniform_2< double, Point >           Creator;
typedef CGAL::Random_points_in_square_2< Point, Creator >  Point_generator;

int main(int argc, char* argv[])
{
  // this is not initialized on MIPSPRO
  CGAL::set_pretty_mode( std::cout);
  CGAL::set_pretty_mode( std::cerr);

  // take #points from command line:
  int n;
  if ( argc < 2 || (n = atoi( argv[1])) < 3) {
    std::cerr << "usage: " << argv[0] << " \"#points\" (>= 3)" << std::endl;
    return 1;
  }

  std::cout << "Test random_convex_set_2:\n" << std::endl;

  // build random n-gon:
  std::cout << "constructing random " << n << "-gon ..." << std::flush;
  Polygon_2 p;
  CGAL::random_convex_set_2( n, std::back_inserter( p), Point_generator( 1));
  std::cout << " done." << std::endl;

  // output polygon:
  std::cout << "\nHere is the result:" << std::endl;

  CGAL::Window_stream W;
  W.init( -1.05, 1.05, -1.05);
  W.display();
  CGAL::cgalize( W);
  W << p;

  W.read_mouse();        // wait for mouse-click

  // check convexity:
  if ( ! p.is_convex()) {
    std::cerr << "ERROR: polygon is not convex." << std::endl;
    return 1;
  }

  std::cout << "done." << std::endl;
  return 0;
} 

// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------


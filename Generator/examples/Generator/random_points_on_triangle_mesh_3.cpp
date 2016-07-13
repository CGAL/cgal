#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/point_generators_3.h>

#include <iostream>
#include <fstream>
using namespace CGAL;

typedef Simple_cartesian<double>                           K;
typedef K::Point_3                                         Point;
typedef K::FT                                              FT;
typedef Surface_mesh<Point>                                Surface_mesh;


int main()
{
 // Generated points are in that vector
  std::vector<Point> points;
  //Construct a Surface_mesh from an OFF file
  ::Surface_mesh sm;
  std::ifstream in("./data/star.off");
  in >> sm;
  CGAL_assertion(in && !sm.is_empty());

  // Create the generator, input is the Surface_mesh sm
  Random_points_on_triangle_mesh_3<Point, ::Surface_mesh>
      g(sm);

  // Get 100 random points in cdt
  CGAL::cpp11::copy_n( g, 100, std::back_inserter(points));

  // Check that we have really created 100 points.
  assert( points.size() == 100);

  // print the first point that was generated
  std::cout << points[0] << std::endl;

  return 0;
}


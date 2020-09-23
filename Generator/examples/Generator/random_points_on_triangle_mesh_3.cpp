#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/point_generators_3.h>

#include <iostream>
#include <fstream>
using namespace CGAL;
typedef Simple_cartesian<double>                           K;
typedef CGAL::Polyhedron_3<K>                              Polyhedron;
typedef K::Point_3                                         Point;
typedef K::FT                                              FT;


int main()
{
 // Generated points are in that vector
  std::vector<Point> points;
  // Create input polyhedron
  Polyhedron polyhedron;
  polyhedron.make_tetrahedron(Point(-1,0,0), Point(0,1,0), Point(1,0,0), Point(0,0,-1));

  // Create the generator, input is the Polyhedron polyhedron
  Random_points_in_triangle_mesh_3<Polyhedron>
      g(polyhedron);

  // Get 100 random points in cdt
  std::copy_n(g, 100, std::back_inserter(points));

  // Check that we have really created 100 points.
  assert( points.size() == 100);

  // print the first point that was generated
  std::cout << points[0] << std::endl;

  return 0;
}


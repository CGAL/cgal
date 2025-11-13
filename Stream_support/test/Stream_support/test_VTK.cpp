#include <CGAL/Simple_cartesian.h>

#include <CGAL/IO/VTK.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <list>

typedef CGAL::Simple_cartesian<double>                Kernel;
typedef Kernel::Point_3                               Point;
typedef std::vector<long unsigned int>                Face;

int main()
{
  const char* vtk_file = "test_VTK.vtp";

  std::vector<Point> points =  { Point(0,0,0), Point(1,0,0), Point(0,1,0) };
  std::vector<Face> polygons = { Face{0,1,2} };

  bool ok = CGAL::IO::write_VTP(vtk_file, points, polygons, CGAL::parameters::verbose(true).use_binary_mode(false));
  assert(ok);

  points.clear();
  polygons.clear();
  ok = CGAL::IO::read_VTP(vtk_file, points, polygons, CGAL::parameters::verbose(true));
  assert(ok);
  std::cout << points.size() << " points and " << polygons.size() << " polygons" << std::endl;

  std::cout << "Done!" << std::endl;
  return EXIT_SUCCESS;
}

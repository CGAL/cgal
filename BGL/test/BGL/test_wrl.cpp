#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/io.h>
#include <iostream>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Surface_mesh<Point>                            Mesh;

int main()
{
  Mesh sm;

  CGAL::make_tetrahedron(Point(0,0,0), Point(1,0,0), Point(1,1,0), Point(0,0,1), sm);
  CGAL::write_wrl(std::cout, sm);
  return 0;
}

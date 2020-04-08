#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/io.h>

#include <iostream>
#include <fstream>

#include <boost/graph/connected_components.hpp>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Surface_mesh<Point>                            Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

int main()
{
  Mesh sm;
  CGAL::make_quad(Point(0,0,0), Point(1,0,0),Point(1,1,0),Point(0,1,0), sm);
  CGAL::make_quad(Point(0,0,1), Point(1,0,1),Point(1,1,1),Point(0,1,1), sm);

  std::ofstream out("out.inp");
  out.precision(17);
  CGAL::write_inp(out, sm, "out.inp", "S4R");
  return 0;
}

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <iostream>
#include <fstream>


typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Surface_mesh<Point>                            Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

int main(int argc, char* argv[]) 
{
  Mesh sm;
  return 0;
}

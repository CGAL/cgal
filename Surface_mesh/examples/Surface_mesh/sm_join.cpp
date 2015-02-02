#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <iostream>
#include <fstream>

#include <boost/foreach.hpp>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Surface_mesh<Point>                            Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

int main(int /* argc */, char* argv[]) 
{
  Mesh sm1, sm2;
  std::ifstream in1(argv[1]);
  in1 >> sm1;
  std::ifstream in2(argv[2]);
  in2 >> sm2;
  sm1 += sm2;

  assert(sm1.is_valid());

  std::cout << sm1 << std::endl;
}

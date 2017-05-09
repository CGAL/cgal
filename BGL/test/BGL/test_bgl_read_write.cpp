#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/io.h>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Surface_mesh<Point>                            Mesh;


int main(int argc, char* argv[]) 
{
  Mesh sm;
  std::ifstream in((argc>1)?argv[1]:"data/prim.off");

  CGAL::read_off(in,sm);

  CGAL::write_off(std::cout, sm);

  return 0;
}

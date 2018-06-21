#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/Euler_operations.h>

#include <iostream>
#include <fstream>

#include <boost/foreach.hpp>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
//typedef CGAL::Surface_mesh<Point>                            Mesh;
typedef CGAL::Polyhedron_3<Kernel> Mesh;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;

int main()
{
  Point p(0,0,0), q(1,0,0), r(1,1,0), s(0,0,1);
  Mesh sm;
  CGAL::make_tetrahedron(p, q, r, s, sm);

  halfedge_descriptor hd1 = *(halfedges(sm).first);
  halfedge_descriptor hd2 = next(hd1,sm);

  CGAL::Euler::flip_edge(hd1,sm);
  std::cout << "after flip_edge" << std::endl;
  std::cout << sm << std::endl;
  CGAL::Euler::flip_edge(hd2,sm);
  std::cout << "after flip_edge" << std::endl;
  std::cout << sm << std::endl;

  std::cout << "done";

  return 0;
}

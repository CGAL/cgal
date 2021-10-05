#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <iostream>
#include <fstream>

#include <boost/graph/connected_components.hpp>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Surface_mesh<Point>                            Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/prim.off";

  Mesh sm;
  if(!CGAL::IO::read_polygon_mesh(filename, sm))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  Mesh::Property_map<vertex_descriptor,int> ccmap;
  ccmap = sm.add_property_map<vertex_descriptor,int>("v:CC").first;

  int num = connected_components(sm, ccmap);
  std::cout  << num << " connected components" << std::endl;
  for(vertex_descriptor v : vertices(sm)){
    std::cout  << v << " is in component " << ccmap[v] << std::endl;
  }

  return 0;
}

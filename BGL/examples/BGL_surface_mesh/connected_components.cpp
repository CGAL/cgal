#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <iostream>
#include <fstream>

#include <boost/graph/connected_components.hpp>
#include <boost/foreach.hpp>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Surface_mesh<Point>                            Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

int main(int argc, char* argv[]) 
{
  Mesh sm;
  std::ifstream in((argc>1)?argv[1]:"data/prim.off");
  in >> sm;

  Mesh::Property_map<vertex_descriptor,int> ccmap;
  ccmap = sm.add_property_map<vertex_descriptor,int>("v:CC").first;

  int num = connected_components(sm, ccmap);
  std::cout  << num << " connected components" << std::endl;
  BOOST_FOREACH(vertex_descriptor v, vertices(sm)){
    std::cout  << v << " is in component " << ccmap[v] << std::endl;
  }
  
  return 0;
}

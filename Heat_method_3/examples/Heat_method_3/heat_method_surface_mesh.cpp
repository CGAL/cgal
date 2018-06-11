#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Heat_method_3/Heat_method_3.h>

#include <iostream>
#include <fstream>

#include <boost/foreach.hpp>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Surface_mesh<Point>                            Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef Mesh::Property_map<vertex_descriptor, double> Vertex_distance_map; 
typedef CGAL::Heat_method_3::Heat_method_3<Mesh,Kernel,Vertex_distance_map> Heat_method;



int main(int argc, char* argv[]) 
{
  Mesh sm;
  std::ifstream in(argv[1]);
  in >> sm;
  Vertex_distance_map heat_intensity =
    sm.add_property_map<vertex_descriptor, double>("v:heat_intensity", 0).first; 
  Heat_method hm(sm, heat_intensity);

  vertex_descriptor source = *(vertices(sm).first);
  hm.add_source(source);

  BOOST_FOREACH(vertex_descriptor vd , vertices(sm)){
    std::cout << vd << "  is at distance " << hm.distance(vd) << " from " << source << std::endl;
  }

  std::cout << "done" << std::endl;
  return 0;
}

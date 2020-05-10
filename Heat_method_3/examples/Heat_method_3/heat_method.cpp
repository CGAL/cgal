#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>

#include <iostream>
#include <fstream>


typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point_3;
typedef CGAL::Surface_mesh<Point_3>                          Triangle_mesh;

typedef boost::graph_traits<Triangle_mesh>::vertex_descriptor vertex_descriptor;
typedef Triangle_mesh::Property_map<vertex_descriptor,double> Vertex_distance_map;

int main(int argc, char* argv[])
{
  Triangle_mesh tm;
  const char* filename = (argc > 1) ? argv[1] : "./data/elephant.off";
  std::ifstream input(filename);
  if (!input || !(input >> tm) || tm.is_empty()) {
    std::cerr << "Not a valid off file." << std::endl;
    return 1;
  }

  //property map for the distance values to the source set
  Vertex_distance_map vertex_distance = tm.add_property_map<vertex_descriptor, double>("v:distance", 0).first;

  vertex_descriptor source = *(vertices(tm).first);

  CGAL::Heat_method_3::estimate_geodesic_distances(tm, vertex_distance, source) ;

  std::cout << "Source vertex " << source << " at: " << tm.point(source) << std::endl;
  for(vertex_descriptor vd : vertices(tm)){
    std::cout << vd << " ("<< tm.point(vd) << ")"
              <<  " is at distance " << get(vertex_distance, vd) << std::endl;
  }

  return 0;
}

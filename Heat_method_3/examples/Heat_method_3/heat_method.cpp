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
  const char* filename = (argc > 1) ? argv[1] : "./data/elephant.off";

  Triangle_mesh tm;
  if(!CGAL::IO::read_polygon_mesh(filename, tm) ||
     CGAL::is_empty(tm) || !CGAL::is_triangle_mesh(tm))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }

  //property map for the distance values to the source set
  Vertex_distance_map vertex_distance = tm.add_property_map<vertex_descriptor, double>("v:distance", 0).first;

  vertex_descriptor source = *(vertices(tm).first);

  CGAL::Heat_method_3::estimate_geodesic_distances(tm, vertex_distance, source) ;

  std::cout << "Source vertex " << source << " at: " << tm.point(source) << std::endl;
  for(vertex_descriptor vd : vertices(tm))
  {
    std::cout << vd << " ("<< tm.point(vd) << ")"
              <<  " is at distance " << get(vertex_distance, vd) << std::endl;
  }

  return 0;
}

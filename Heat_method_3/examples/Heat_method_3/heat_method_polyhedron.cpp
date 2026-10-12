#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>

#include <unordered_map>
#include <fstream>
#include <iostream>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point_3;
typedef CGAL::Polyhedron_3<Kernel>                           Triangle_mesh;

typedef boost::graph_traits<Triangle_mesh>::vertex_descriptor vertex_descriptor;

int main(int argc, char* argv[])
{
  Triangle_mesh tm;
  const std::filesystem::path filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/elephant.off");
  std::ifstream input(filename);
  if (!input || !(input >> tm) || tm.is_empty()) {
    std::cerr << "Not a valid off file." << std::endl;
    return 1;
  }
  // map for the distance values to the source set
  std::unordered_map<vertex_descriptor, double> vertex_distance;

  vertex_descriptor source = *(vertices(tm).first);

  CGAL::Heat_method_3::estimate_geodesic_distances(tm,
                                                   boost::make_assoc_property_map(vertex_distance),
                                                   source) ;

  std::cout << "Source vertex at: " << source->point() << std::endl;
  double max = 0;
  vertex_descriptor vfar = source;
  for(vertex_descriptor vd : vertices(tm)){
    double d = vertex_distance[vd];
    if(d>max){
      max = d;
      vfar = vd;
    }
  }
  std::cout << "Furthest point at: " << vfar->point() << std::endl;

  return 0;
}

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>

#include <fstream>
#include <iostream>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point_3;
typedef CGAL::Surface_mesh<Point_3>                          Triangle_mesh;

typedef boost::graph_traits<Triangle_mesh>::vertex_descriptor vertex_descriptor;
typedef Triangle_mesh::Property_map<vertex_descriptor,double> Vertex_distance_map;
typedef CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Triangle_mesh> Heat_method;

int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/larger_sphere.off");

  Triangle_mesh tm;
  if(!CGAL::IO::read_polygon_mesh(filename, tm) ||
     CGAL::is_empty(tm) || !CGAL::is_triangle_mesh(tm))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }

  //property map for the distance values to the source set
  Vertex_distance_map vertex_distance = tm.add_property_map<vertex_descriptor, double>("v:distance", 0).first;

  Heat_method hm(tm);

  //add the first vertex as the source set
  vertex_descriptor source = *(vertices(tm).first);
  hm.add_source(source);
  hm.estimate_geodesic_distances(vertex_distance);

  Point_3 sp = tm.point(source);

  std::cout << "source: " << sp  << " " << source << std::endl;
  vertex_descriptor vfar;
  double sdistance = 0;

  for(vertex_descriptor vd : vertices(tm)){
    std::cout << vd << "  is at distance " << get(vertex_distance, vd) << " to " << source << std::endl;
    if(get(vertex_distance, vd) > sdistance){
      vfar = vd;
      sdistance = get(vertex_distance, vd);
    }
  }

  std::cout << "vfar: " << tm.point(vfar) << " " << vfar << std::endl;

  hm.add_source(vfar);
  hm.estimate_geodesic_distances(vertex_distance);

  for(vertex_descriptor vd : vertices(tm)){
    std::cout << vd << "  is at distance " << get(vertex_distance, vd) << "to the set of two sources" << std::endl;
  }

  std::cout << "done" << std::endl;
  return 0;
}

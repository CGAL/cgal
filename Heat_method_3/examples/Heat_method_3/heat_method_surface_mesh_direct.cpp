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

typedef CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Triangle_mesh, CGAL::Heat_method_3::Direct> Heat_method_idt;


int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/elephant.off");

  Triangle_mesh tm;
  if(!CGAL::IO::read_polygon_mesh(filename, tm) ||
     CGAL::is_empty(tm) || !CGAL::is_triangle_mesh(tm))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }

   //property map for the distance values to the source set
  Vertex_distance_map vertex_distance = tm.add_property_map<vertex_descriptor,double>("v:distance",0).first;

  //pass in the idt object and its vertex_distance_map
  Heat_method_idt hm_idt(tm);

  //add the first vertex as the source set
  vertex_descriptor source = *(vertices(tm).first);
  hm_idt.add_source(source);
  hm_idt.estimate_geodesic_distances(vertex_distance);

  for(vertex_descriptor vd : vertices(tm)){
    std::cout << vd << "  is at distance " << get(vertex_distance, vd) << " from " << source << std::endl;
  }

  return 0;
}

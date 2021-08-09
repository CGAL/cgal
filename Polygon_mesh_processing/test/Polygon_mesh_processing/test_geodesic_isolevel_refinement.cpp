#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>

#include <CGAL/Polygon_mesh_processing/internal/refine_mesh_at_isolevel.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>

#include <fstream>
#include <iostream>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point_3;
typedef CGAL::Surface_mesh<Point_3>                          Triangle_mesh;

typedef boost::graph_traits<Triangle_mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Triangle_mesh>::edge_descriptor edge_descriptor;

typedef Triangle_mesh::Property_map<vertex_descriptor,double> Vertex_distance_map;
typedef CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Triangle_mesh> Heat_method;

int main(int argc, char* argv[])
{
  const char* filename = argv[1];

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
  vertex_descriptor s = *(vertices(tm).first);
  hm.add_source(s);
  hm.estimate_geodesic_distances(vertex_distance);

  //property map for the constrained status of edges
  auto ecm = tm.add_property_map<edge_descriptor, bool>("e:is_constrained", 0).first;


  for (int i=2; i<argc; ++i)
    CGAL::Polygon_mesh_processing::experimental::refine_mesh_at_isolevel(tm, vertex_distance, atof(argv[i]), CGAL::parameters::edge_is_constrained_map(ecm));

  std::vector<Triangle_mesh> splitted;

  CGAL::Polygon_mesh_processing::split_connected_components(tm, splitted, CGAL::parameters::edge_is_constrained_map(ecm));

#ifdef CGAL_TEST_SUITE
  assert(splitted.size() == 22);
#else
  for(std::size_t i=0; i<splitted.size(); ++i)
    std::ofstream("out_"+std::to_string(i)+".off") << std::setprecision(17) << splitted[i];
#endif

  return 0;
}

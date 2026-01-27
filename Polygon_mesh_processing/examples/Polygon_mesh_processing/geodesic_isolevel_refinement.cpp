#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>

#include <CGAL/Polygon_mesh_processing/refine_mesh_at_isolevel.h>
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
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/elephant.off");

  Triangle_mesh tm;
  if(!CGAL::IO::read_polygon_mesh(filename, tm) ||
     CGAL::is_empty(tm) || !CGAL::is_triangle_mesh(tm))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }

  // default isovalues for cutting the mesh
  std::vector<double> isovalues = {0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 10};

  if (argc>2)
  {
    isovalues.clear();
    for (int i=2; i<argc; ++i)
      isovalues.push_back(atof(argv[i]));
  }

  //property map for the distance values to the source set
  Vertex_distance_map vertex_distance = tm.add_property_map<vertex_descriptor, double>("v:distance", 0).first;

  Heat_method hm(tm);

  //use heat method to compute approximated geodesic distances to the source vertex `s`
  vertex_descriptor s = *(vertices(tm).first);
  hm.add_source(s);
  hm.estimate_geodesic_distances(vertex_distance);

  // property map to flag new cut edge added in the mesh
  auto ecm = tm.add_property_map<edge_descriptor, bool>("e:is_constrained", 0).first;

  // refine the mesh along isovalues
  for (double isovalue : isovalues)
    CGAL::Polygon_mesh_processing::refine_mesh_at_isolevel(tm, vertex_distance, isovalue, CGAL::parameters::edge_is_constrained_map(ecm));

  // split the mesh in connected components bounded by the isocurves
  std::vector<Triangle_mesh> edges_split;
  CGAL::Polygon_mesh_processing::split_connected_components(tm, edges_split, CGAL::parameters::edge_is_constrained_map(ecm));

  assert(argc!=1 || edges_split.size() == 22);

  // export each submesh in a file
  for(std::size_t i=0; i<edges_split.size(); ++i)
    std::ofstream("out_"+std::to_string(i)+".off") << std::setprecision(17) << edges_split[i];

  return 0;
}

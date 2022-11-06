#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/Curvatures/interpolated_corrected_curvature_measures.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/property_map.h>

#include <boost/graph/graph_traits.hpp>

#include <iostream>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic_Kernel;
typedef CGAL::Surface_mesh<Epic_Kernel::Point_3> Surface_Mesh;
typedef boost::graph_traits<Surface_Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Surface_Mesh>::vertex_descriptor vertex_descriptor;

int main(int argc, char* argv[])
{
  Surface_Mesh g1;
  const std::string filename = (argc>1) ?
      argv[1] :
      CGAL::data_file_path("meshes/small_bunny.obj");

  if(!CGAL::IO::read_polygon_mesh(filename, g1))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }

  // creating and tying surface mesh property maps for curvatures (with defaults = 0)
  bool created = false;
  Surface_Mesh::Property_map<vertex_descriptor, Epic_Kernel::FT> mean_curvature_map, gaussian_curvature_map;
  boost::tie(mean_curvature_map, created) = g1.add_property_map<vertex_descriptor, Epic_Kernel::FT>("v:mean_curvature_map", 0);
  assert(created);

  boost::tie(gaussian_curvature_map, created) = g1.add_property_map<vertex_descriptor, Epic_Kernel::FT>("v:gaussian_curvature_map", 0);
  assert(created);

  // we use a tuple of 2 scalar values and 2 vectors for principal curvatures and directions
  Surface_Mesh::Property_map<vertex_descriptor, std::tuple<
      Epic_Kernel::FT,
      Epic_Kernel::FT,
      Epic_Kernel::Vector_3,
      Epic_Kernel::Vector_3
      >> principal_curvature_map;

  boost::tie(principal_curvature_map, created) = g1.add_property_map<vertex_descriptor, std::tuple<
      Epic_Kernel::FT,
      Epic_Kernel::FT,
      Epic_Kernel::Vector_3,
      Epic_Kernel::Vector_3
      >>("v:principal_curvature_map", { 0, 0,
          Epic_Kernel::Vector_3 (0,0,0),
          Epic_Kernel::Vector_3 (0,0,0)});
  assert(created);

  PMP::interpolated_corrected_mean_curvature(
      g1,
      mean_curvature_map
  );

  // uncomment this to compute a curvature while specifying named parameters
  // Example: an expansion ball radius of 0.5 and a vertex normals map (does not have to depend on positions)

  /*Surface_Mesh::Property_map<vertex_descriptor, Epic_Kernel::Vector_3> vnm;
  boost::tie(vnm, created) = g1.add_property_map<vertex_descriptor, Epic_Kernel::Vector_3>(
      "v:vnm", Epic_Kernel::Vector_3(0, 0, 0)
  );

  assert(created);

  PMP::interpolated_corrected_mean_curvature(
      g1,
      mean_curvature_map,
      CGAL::parameters::ball_radius(0.5).vertex_normal_map(vnm)
  );*/

  PMP::interpolated_corrected_gaussian_curvature(
      g1,
      gaussian_curvature_map
  );
  PMP::interpolated_corrected_principal_curvatures(
      g1,
      principal_curvature_map
  );

  for (vertex_descriptor v : vertices(g1))
  {
    auto PC = principal_curvature_map[v];
      std::cout << v.idx() << ": HC = " << mean_curvature_map[v]
                           << ", GC = " << gaussian_curvature_map[v] << "\n"
                           << ", PC = [ " << std::get<0>(PC) << " , " << std::get<1>(PC) << " ]\n";
  }
}

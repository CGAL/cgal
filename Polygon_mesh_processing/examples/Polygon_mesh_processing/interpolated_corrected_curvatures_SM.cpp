#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/interpolated_corrected_curvatures.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/property_map.h>

#include <boost/graph/graph_traits.hpp>

#include <iostream>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic_kernel;
typedef CGAL::Surface_mesh<Epic_kernel::Point_3> Surface_Mesh;
typedef boost::graph_traits<Surface_Mesh>::vertex_descriptor vertex_descriptor;

int main(int argc, char* argv[])
{
  Surface_Mesh smesh;
  const std::string filename = (argc > 1) ?
    argv[1] :
    CGAL::data_file_path("meshes/sphere.off");

  if (!CGAL::IO::read_polygon_mesh(filename, smesh))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }

  // creating and tying surface mesh property maps for curvatures (with defaults = 0)
  bool created = false;
  Surface_Mesh::Property_map<vertex_descriptor, Epic_kernel::FT>
    mean_curvature_map, Gaussian_curvature_map;

  boost::tie(mean_curvature_map, created) =
    smesh.add_property_map<vertex_descriptor, Epic_kernel::FT>("v:mean_curvature_map", 0);
  assert(created);

  boost::tie(Gaussian_curvature_map, created) =
    smesh.add_property_map<vertex_descriptor, Epic_kernel::FT>("v:Gaussian_curvature_map", 0);
  assert(created);

  // we use a tuple of 2 scalar values and 2 vectors for principal curvatures and directions
  Surface_Mesh::Property_map<vertex_descriptor, PMP::Principal_curvatures_and_directions<Epic_kernel>>
    principal_curvatures_and_directions_map;

  boost::tie(principal_curvatures_and_directions_map, created) =
    smesh.add_property_map<vertex_descriptor, PMP::Principal_curvatures_and_directions<Epic_kernel>>
    ("v:principal_curvatures_and_directions_map", { 0, 0,
        Epic_kernel::Vector_3(0,0,0),
        Epic_kernel::Vector_3(0,0,0) });
  assert(created);

  // user can call these fucntions to compute a specfic curvature type on each vertex.
  // (Note: if no ball radius is specified, the measure expansion of each vertex happens by
  // summing measures on faces adjacent to each vertex.)
  PMP::interpolated_corrected_mean_curvature(smesh, mean_curvature_map);

  PMP::interpolated_corrected_Gaussian_curvature(smesh, Gaussian_curvature_map);

  PMP::interpolated_corrected_principal_curvatures_and_directions(smesh, principal_curvatures_and_directions_map);

  // uncomment this to compute a curvature while specifying named parameters
  // Example: an expansion ball radius of 0.5 and a vertex normals map (does not have to depend on positions)

  /*Surface_Mesh::Property_map<vertex_descriptor, Epic_kernel::Vector_3> vnm;
  boost::tie(vnm, created) = smesh.add_property_map<vertex_descriptor, Epic_kernel::Vector_3>(
      "v:vnm", Epic_kernel::Vector_3(0, 0, 0)
  );

  assert(created);

  PMP::interpolated_corrected_mean_curvature(
      smesh,
      mean_curvature_map,
      CGAL::parameters::ball_radius(0.5).vertex_normal_map(vnm)
  );*/

  // This function can be used to compute multiple curvature types by specifiying them as named parameters
  // This is more efficient than computing each one separately (shared computations).
  PMP::interpolated_corrected_curvatures(
    smesh,
    CGAL::parameters::vertex_mean_curvature_map(mean_curvature_map)
    .vertex_principal_curvatures_and_directions_map(principal_curvatures_and_directions_map)
  );

  for (vertex_descriptor v : vertices(smesh))
  {
    auto PC = principal_curvatures_and_directions_map[v];
    std::cout << v.idx() << ": HC = " << mean_curvature_map[v]
      << ", GC = " << Gaussian_curvature_map[v] << "\n"
      << ", PC = [ " << PC.min_curvature << " , " << PC.max_curvature << " ]\n";
  }
}

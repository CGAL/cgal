#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/interpolated_corrected_curvatures.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/property_map.h>

#include <boost/graph/graph_traits.hpp>

#include <iostream>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic_kernel;
typedef CGAL::Surface_mesh<Epic_kernel::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

int main(int argc, char* argv[])
{
  Mesh smesh;
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
  Mesh::Property_map<vertex_descriptor, Epic_kernel::FT>
    mean_curvature_map, Gaussian_curvature_map;

  std::tie(mean_curvature_map, created) =
    smesh.add_property_map<vertex_descriptor, Epic_kernel::FT>("v:mean_curvature_map", 0);
  assert(created);

  std::tie(Gaussian_curvature_map, created) =
    smesh.add_property_map<vertex_descriptor, Epic_kernel::FT>("v:Gaussian_curvature_map", 0);
  assert(created);

  // we use a tuple of 2 scalar values and 2 vectors for principal curvatures and directions
  Mesh::Property_map<vertex_descriptor, PMP::Principal_curvatures_and_directions<Epic_kernel>>
    principal_curvatures_and_directions_map;

  std::tie(principal_curvatures_and_directions_map, created) =
    smesh.add_property_map<vertex_descriptor, PMP::Principal_curvatures_and_directions<Epic_kernel>>
    ("v:principal_curvatures_and_directions_map", { 0, 0,
        Epic_kernel::Vector_3(0,0,0),
        Epic_kernel::Vector_3(0,0,0) });
  assert(created);

  PMP::interpolated_corrected_curvatures(smesh,
    CGAL::parameters::vertex_mean_curvature_map(mean_curvature_map)
                     .vertex_Gaussian_curvature_map(Gaussian_curvature_map)
                     .vertex_principal_curvatures_and_directions_map(principal_curvatures_and_directions_map)
  // uncomment to use an expansion ball radius of 0.5 to estimate the curvatures
  //                 .ball_radius(0.5)
  );


  for (vertex_descriptor v : vertices(smesh))
  {
    auto PC = principal_curvatures_and_directions_map[v];
    std::cout << v.idx() << ": HC = " << mean_curvature_map[v]
      << ", GC = " << Gaussian_curvature_map[v] << "\n"
      << ", PC = [ " << PC.min_curvature << " , " << PC.max_curvature << " ]\n";
  }
  return 0;
}

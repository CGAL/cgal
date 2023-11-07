#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/acvd/acvd.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/property_map.h>
#include <CGAL/Polygon_mesh_processing/interpolated_corrected_curvatures.h>

#include <boost/graph/graph_traits.hpp>

#include <iostream>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic_kernel;
typedef CGAL::Surface_mesh<Epic_kernel::Point_3> Surface_Mesh;
typedef boost::graph_traits<Surface_Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::property_map<Surface_Mesh, CGAL::dynamic_vertex_property_t<CGAL::IO::Color> >::type VertexColorMap;


int main(int argc, char* argv[])
{
  Surface_Mesh smesh;
  const std::string filename = (argc > 1) ?
    CGAL::data_file_path(argv[1]) :
    CGAL::data_file_path("meshes/cactus.off");

  const int nb_clusters = (argc > 2) ? atoi(argv[2]) : 100;

  if (!CGAL::IO::read_polygon_mesh(filename, smesh))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }

  /// Uniform Isotropic ACVD Remeshing

  auto acvd_mesh = PMP::acvd_isotropic_simplification(smesh, nb_clusters);
  CGAL::IO::write_OFF("acvd_mesh.off", acvd_mesh);

  /// Adaptive Isotropic ACVD Remeshing

  bool created = false;
  Surface_Mesh::Property_map<vertex_descriptor, PMP::Principal_curvatures_and_directions<Epic_kernel>>
    principal_curvatures_and_directions_map;

  boost::tie(principal_curvatures_and_directions_map, created) =
    smesh.add_property_map<vertex_descriptor, PMP::Principal_curvatures_and_directions<Epic_kernel>>
    ("v:principal_curvatures_and_directions_map", { 0, 0,
        Epic_kernel::Vector_3(0,0,0),
        Epic_kernel::Vector_3(0,0,0) });
  assert(created);

  PMP::interpolated_corrected_principal_curvatures_and_directions(smesh, principal_curvatures_and_directions_map);

  const double gradation_factor = (argc > 3) ? atof(argv[3]) : 0;

  auto adaptive_acvd_mesh =
    PMP::acvd_isotropic_simplification(
      smesh,
      nb_clusters,
      CGAL::parameters::vertex_principal_curvatures_and_directions_map(principal_curvatures_and_directions_map)
        .gradation_factor(gradation_factor)
    );

  CGAL::IO::write_OFF("acvd_mesh_adaptive.off", adaptive_acvd_mesh);

  return 0;
}


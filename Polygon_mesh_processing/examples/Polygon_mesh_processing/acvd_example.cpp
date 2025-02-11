#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/acvd/acvd.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/property_map.h>
#include <CGAL/Polygon_mesh_processing/interpolated_corrected_curvatures.h>
#include <CGAL/Real_timer.h>

#include <boost/graph/graph_traits.hpp>

#include <iostream>
#include <filesystem>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic_kernel;
typedef CGAL::Surface_mesh<Epic_kernel::Point_3> Surface_Mesh;
typedef boost::graph_traits<Surface_Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::property_map<Surface_Mesh, CGAL::dynamic_vertex_property_t<CGAL::IO::Color> >::type VertexColorMap;

namespace params = CGAL::parameters;

int main(int argc, char* argv[])
{
  Surface_Mesh smesh;
  //CGAL::get_default_random() = CGAL::Random( 1739197120 ); //connexity constraint issue + infinite loop with cheese
  //CGAL::get_default_random() = CGAL::Random( 1739199762 ); //one small edge
  //CGAL::get_default_random() = CGAL::Random( 1739264620 ); //one very small edge
  CGAL::get_default_random() = CGAL::Random( 1739293586 ); // seed elephant non-manifold with 300 clusters

  std::cout << "Seed : " << CGAL::get_default_random().get_seed() << std::endl;
  const std::string filename = (argc > 1) ?
    argv[1] :
    CGAL::data_file_path("meshes/fandisk.off");

  const std::string stem = std::filesystem::path(filename).stem().string();
  const std::string extension = std::filesystem::path(filename).extension().string();

  const int nb_clusters = (argc > 2) ? atoi(argv[2]) : 3000;
  const std::string nbc = std::to_string(nb_clusters);

  if (!CGAL::IO::read_polygon_mesh(filename, smesh))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }

  ///// Uniform Isotropic ACVD
#if 0
  std::cout << "Uniform Isotropic ACVD ...." << std::endl;
  auto acvd_mesh = PMP::acvd_isotropic_remeshing(smesh, nb_clusters);
  CGAL::IO::write_polygon_mesh(stem+"_acvd_"+nbc+extension, acvd_mesh);

  std::cout << "Completed" << std::endl;
#endif

  //// With Post-Processing QEM Optimization
#if 0
  std::cout << "Uniform Isotropic ACVD with QEM optimization ...." << std::endl;

  auto acvd_mesh_qem_pp = PMP::acvd_isotropic_remeshing(smesh, nb_clusters, CGAL::parameters::postprocessing_qem(true));
  CGAL::IO::write_polygon_mesh(stem+"_acvd_qem-pp_"+nbc+extension, acvd_mesh_qem_pp);

  std::cout << "Completed" << std::endl;
#endif

  // With QEM Energy Minimization

  std::cout << "Uniform QEM ACVD ...." << std::endl;
CGAL::Real_timer timer; timer.start();
  auto acvd_mesh_qem = PMP::acvd_isotropic_remeshing(smesh, nb_clusters, params::use_qem_metric(true));
timer.stop();
  CGAL::IO::write_polygon_mesh( stem +"_acvd_qem_"+ std::to_string(nb_clusters) + extension, acvd_mesh_qem);

  std::cout << "Completed in " << timer.time() << "s\n";

  /// Adaptive Isotropic ACVD

  /*std::cout << "Adaptive Isotropic ACVD ...." << std::endl;

  const double gradation_factor = (argc > 3) ? atof(argv[3]) : 2;

  bool created = false;
  Surface_Mesh::Property_map<vertex_descriptor, PMP::Principal_curvatures_and_directions<Epic_kernel>>
   principal_curvatures_and_directions_map;

  boost::tie(principal_curvatures_and_directions_map, created) =
   smesh.add_property_map<vertex_descriptor, PMP::Principal_curvatures_and_directions<Epic_kernel>>
   ("v:principal_curvatures_and_directions_map", { 0, 0,
       Epic_kernel::Vector_3(0,0,0),
       Epic_kernel::Vector_3(0,0,0) });
  assert(created);

  PMP::interpolated_corrected_curvatures(smesh, CGAL::parameters::vertex_principal_curvatures_and_directions_map(principal_curvatures_and_directions_map));

  auto adaptive_acvd_mesh =
   PMP::acvd_isotropic_remeshing(
     smesh,
     nb_clusters,
     CGAL::parameters::vertex_principal_curvatures_and_directions_map(principal_curvatures_and_directions_map)
       .gradation_factor(gradation_factor)
   );

  CGAL::IO::write_OFF("fandisk_acvd_adaptive_3000.off", adaptive_acvd_mesh);

  std::cout << "Completed" << std::endl;*/

  return 0;
}


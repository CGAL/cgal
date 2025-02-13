#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/approximated_centroidal_Voronoi_diagram_remeshing.h>
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
typedef CGAL::Surface_mesh<Epic_kernel::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

namespace params = CGAL::parameters;

int main(int argc, char* argv[])
{
  Mesh smesh;
  //CGAL::get_default_random() = CGAL::Random( 1739197120 ); //connexity constraint issue + infinite loop with cheese
  //CGAL::get_default_random() = CGAL::Random( 1739199762 ); //one small edge
  //CGAL::get_default_random() = CGAL::Random( 1739264620 ); //one very small edge
  //CGAL::get_default_random() = CGAL::Random( 1739293586 ); // seed elephant non-manifold with 300 clusters + regression for fandisk 3000 with qem

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
  Mesh acvd_mesh = smesh;
  PMP::approximated_centroidal_Voronoi_diagram_remeshing(acvd_mesh, nb_clusters);
  CGAL::IO::write_polygon_mesh(stem+"_acvd_"+nbc+extension, acvd_mesh);

  std::cout << "Completed" << std::endl;
#endif

  //// With Post-Processing QEM Optimization
#if 0
  std::cout << "Uniform Isotropic ACVD with QEM optimization ...." << std::endl;

  Mesh acvd_mesh_qem_pp = smesh;
  PMP::approximated_centroidal_Voronoi_diagram_remeshing(acvd_mesh_qem_pp, nb_clusters, params::use_postprocessing_qem(true));
  CGAL::IO::write_polygon_mesh( stem +"_acvd_qem-pp_"+nbc+extension, acvd_mesh_qem_pp);

  std::cout << "Completed" << std::endl;
#endif

#if 1
  // With QEM Energy Minimization
  std::cout << "Uniform QEM ACVD ...." << std::endl;
  auto acvd_mesh_qem = smesh;
  PMP::approximated_centroidal_Voronoi_diagram_remeshing(acvd_mesh_qem, nb_clusters, params::use_qem_based_energy(true));
  CGAL::IO::write_polygon_mesh( stem +"_acvd_qem_"+ std::to_string(nb_clusters) + extension, acvd_mesh_qem);
#endif

#if 0
  /// Adaptive Isotropic ACVD
  std::cout << "Adaptive Isotropic ACVD ...." << std::endl;
  const double gradation_factor = 2;
  Mesh adaptive_acvd_mesh = smesh;
  PMP::approximated_centroidal_Voronoi_diagram_remeshing(adaptive_acvd_mesh, nb_clusters, params::gradation_factor(gradation_factor));
  CGAL::IO::write_OFF("fandisk_acvd_adaptive_3000.off", adaptive_acvd_mesh);
#endif

  std::cout << "Completed" << std::endl;

  return 0;
}


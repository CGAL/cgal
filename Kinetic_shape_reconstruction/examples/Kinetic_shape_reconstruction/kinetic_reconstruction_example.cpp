#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kinetic_shape_reconstruction_3.h>
#include <CGAL/IO/PLY_reader.h>
#include <CGAL/IO/PLY_writer.h>
#include <CGAL/Surface_mesh.h>

using SCF   = CGAL::Simple_cartesian<float>;
using SCD   = CGAL::Simple_cartesian<double>;
using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

using Kernel  = EPICK;
using Point_3 = typename Kernel::Point_3;

using Surface_mesh = CGAL::Surface_mesh<Point_3>;
using KSR = CGAL::Kinetic_shape_reconstruction_3<Kernel>;

int main(const int argc, const char** argv) {

  // Input.
  const auto kernel_name = boost::typeindex::type_id<Kernel>().pretty_name();
  std::string input_filename = (argc > 1 ? argv[1] : "data/reconstruction-test/syntetic-building.ply");
  std::ifstream input_file(input_filename);

  // TODO: reading data
  // TODO: getting point map
  // TODO: getting normal map
  // TODO: getting label map

  std::cout << std::endl;
  std::cout << "--- INPUT STATS: " << std::endl;
  std::cout << "* used kernel: "      << kernel_name        << std::endl;
  // std::cout << "* number of points: " << input_range.size() << std::endl;

  // Parameters.
  // const bool verbose = true;
  // const bool debug   = true;

  // Algorithm.
  // KSR ksr(verbose, debug);
  // const bool is_success = ksr.reconstruct(
  //   input_range, point_map, normal_map, label_map,
  //   CGAL::parameters::all_default());
  // assert(is_success);

  // Output.
  // ksr.output_reconstructed_model();

  std::cout << std::endl;
  std::cout << "--- OUTPUT STATS: " << std::endl;
  // std::cout << "* number of model vertices: " << num_vertices << std::endl;
  // std::cout << "* number of model edges: "    << num_edges    << std::endl;
  // std::cout << "* number of model faces: "    << num_faces    << std::endl;

  // Export.
  std::cout << std::endl;
  std::cout << "--- EXPORT: " << std::endl;
  // std::cout << "* model exported successfully" << std::endl;

  std::cout << std::endl << "3D KINETIC RECONSTRUCTION DONE!" << std::endl << std::endl;
  return EXIT_SUCCESS;
}

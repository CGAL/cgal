#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kinetic_shape_reconstruction_3.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/IO/PLY_writer.h>

using SCF   = CGAL::Simple_cartesian<float>;
using SCD   = CGAL::Simple_cartesian<double>;
using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

using Kernel  = EPICK;
using Point_3 = typename Kernel::Point_3;

using Point_set    = CGAL::Point_set_3<Point_3>;
using Point_map    = typename Point_set::Point_map;
using Vector_map   = typename Point_set::Vector_map;
using Label_map    = typename Point_set:: template Property_map<int>;
using Semantic_map = CGAL::KSR::Semantic_from_label_map<Label_map>;

using KSR = CGAL::Kinetic_shape_reconstruction_3<Kernel>;

int main(const int argc, const char** argv) {

  // Input.
  const auto kernel_name = boost::typeindex::type_id<Kernel>().pretty_name();

  const bool with_normals = true;
  Point_set point_set(with_normals);
  std::string input_filename = (argc > 1 ? argv[1] : "data/reconstruction-test/syntetic-building.ply");
  std::ifstream input_file(input_filename, std::ios_base::binary);
  input_file >> point_set;
  input_file.close();

  std::cout << std::endl;
  std::cout << "--- INPUT STATS: " << std::endl;
  std::cout << "* used kernel: "      << kernel_name      << std::endl;
  std::cout << "* number of points: " << point_set.size() << std::endl;

  // Parameters.
  const bool verbose = true;
  const bool debug   = true;

  // Define a map from a user-defined label to the semantic label.
  const Label_map label_map = point_set. template property_map<int>("label").first;
  const Semantic_map semantic_map(label_map, "0", "1", "2", "3", verbose);

  // Algorithm.
  KSR ksr(verbose, debug);
  const bool is_success = ksr.reconstruct(
    point_set,
    point_set.point_map(),
    point_set.normal_map(),
    semantic_map,
    CGAL::parameters::all_default());
  assert(is_success);

  // Output.
  std::vector<Point_3> output_vertices;
  std::vector< std::vector<std::size_t> > output_faces;
  ksr.output_reconstructed_model(
    std::back_inserter(output_vertices),
    std::back_inserter(output_faces));
  const std::size_t num_vertices = output_vertices.size();
  const std::size_t num_faces    = output_faces.size();

  std::cout << std::endl;
  std::cout << "--- OUTPUT STATS: " << std::endl;
  std::cout << "* number of vertices: " << num_vertices << std::endl;
  std::cout << "* number of faces: "    << num_faces    << std::endl;

  // Export.
  std::cout << std::endl;
  std::cout << "--- EXPORT: " << std::endl;

  // Model.
  const std::string output_filename = "reconstructed-model.ply";
  std::ofstream output_file_model(output_filename);
  output_file_model.precision(20);
  if (!CGAL::write_PLY(output_file_model, output_vertices, output_faces)) {
    std::cerr << "ERROR: can't write to the file " << output_filename << "!" << std::endl;
    return EXIT_FAILURE;
  }
  output_file_model.close();
  std::cout << "* the reconstructed model exported successfully" << std::endl;

  std::cout << std::endl << "3D KINETIC RECONSTRUCTION DONE!" << std::endl << std::endl;
  return EXIT_SUCCESS;
}

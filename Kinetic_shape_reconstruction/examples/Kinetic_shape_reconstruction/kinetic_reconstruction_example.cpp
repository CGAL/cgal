#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kinetic_shape_reconstruction_3.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/Real_timer.h>
#include <CGAL/IO/PLY.h>

#include "include/Parameters.h"
#include "include/Terminal_parser.h"

using SCF   = CGAL::Simple_cartesian<float>;
using SCD   = CGAL::Simple_cartesian<double>;
using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

using Kernel    = EPICK;
using FT        = typename Kernel::FT;
using Point_3   = typename Kernel::Point_3;
using Segment_3 = typename Kernel::Segment_3;

using Point_set    = CGAL::Point_set_3<Point_3>;
using Point_map    = typename Point_set::Point_map;
using Vector_map   = typename Point_set::Vector_map;
using Label_map    = typename Point_set:: template Property_map<int>;
using Semantic_map = CGAL::KSR::Semantic_from_label_map<Label_map>;

using KSR = CGAL::Kinetic_shape_reconstruction_3<Kernel>;

using Parameters      = CGAL::KSR::All_parameters<FT>;
using Terminal_parser = CGAL::KSR::Terminal_parser<FT>;
using Timer           = CGAL::Real_timer;

void parse_terminal(Terminal_parser& parser, Parameters& parameters) {
  // Set all parameters that can be loaded from the terminal.
  // add_str_parameter  - adds a string-type parameter
  // add_val_parameter  - adds a scalar-type parameter
  // add_bool_parameter - adds a boolean parameter

  std::cout << std::endl;
  std::cout << "--- INPUT PARAMETERS: " << std::endl;

  // Required parameters.
  parser.add_str_parameter("-data", parameters.data);

  // Label indices.
  parser.add_str_parameter("-gi", parameters.gi);
  parser.add_str_parameter("-bi", parameters.bi);
  parser.add_str_parameter("-ii", parameters.ii);
  parser.add_str_parameter("-vi", parameters.vi);

  // Main parameters.
  parser.add_val_parameter("-scale", parameters.scale);
  parser.add_val_parameter("-noise", parameters.noise);

  // Update.
  parameters.update_dependent();

  // Shape detection.
  parser.add_val_parameter("-kn"   , parameters.k_neighbors);
  parser.add_val_parameter("-dist" , parameters.distance_threshold);
  parser.add_val_parameter("-angle", parameters.angle_threshold);
  parser.add_val_parameter("-minp" , parameters.min_region_size);

  // Shape regularization.
  parser.add_bool_parameter("-regularize", parameters.regularize);

  // Partitioning.
  parser.add_val_parameter("-k", parameters.k_intersections);

  // Reconstruction.
  parser.add_val_parameter("-beta", parameters.graphcut_beta);
}

int main(const int argc, const char** argv) {

  // Parameters.
  std::cout.precision(20);
  std::cout << std::endl;
  std::cout << "--- PARSING INPUT: " << std::endl;
  const auto kernel_name = boost::typeindex::type_id<Kernel>().pretty_name();
  std::cout << "* used kernel: " << kernel_name << std::endl;
  const std::string path_to_save = "/Users/monet/Documents/gf/kinetic/logs/";
  Terminal_parser parser(argc, argv, path_to_save);

  Parameters parameters;
  parse_terminal(parser, parameters);

  // Input.
  Point_set point_set(parameters.with_normals);
  std::ifstream input_file(parameters.data, std::ios_base::binary);
  input_file >> point_set;
  input_file.close();

  std::cout << std::endl;
  std::cout << "--- INPUT STATS: " << std::endl;
  std::cout << "* number of points: " << point_set.size() << std::endl;

  // Define a map from a user-defined label to the semantic label.
  const Label_map label_map = point_set. template property_map<int>("label").first;
  const bool is_defined = point_set. template property_map<int>("label").second;
  const Semantic_map semantic_map(
    label_map,
    is_defined,
    parameters.gi,
    parameters.bi,
    parameters.ii,
    parameters.vi,
    parameters.verbose);

  // Algorithm.
  KSR ksr(parameters.verbose, parameters.debug);

  Timer timer;
  timer.start();
  const bool is_ksr_success = ksr.reconstruct(
    point_set,
    point_set.point_map(),
    point_set.normal_map(),
    semantic_map,
    CGAL::parameters::
    k_neighbors(parameters.k_neighbors).
    distance_threshold(parameters.distance_threshold).
    angle_threshold(parameters.angle_threshold).
    min_region_size(parameters.min_region_size).
    regularize(parameters.regularize).
    k_intersections(parameters.k_intersections).
    graphcut_beta(parameters.graphcut_beta));
  assert(is_ksr_success);
  const std::string success = is_ksr_success ? "SUCCESS" : "FAILED";
  timer.stop();
  const FT time = static_cast<FT>(timer.time());

  // Output.

  // Vertices.
  std::vector<Point_3> all_vertices;
  ksr.output_partition_vertices(
    std::back_inserter(all_vertices), -1);

  // Edges.
  std::vector<Segment_3> all_edges;
  ksr.output_partition_edges(
    std::back_inserter(all_edges), -1);

  // Faces.
  std::vector< std::vector<std::size_t> > all_faces;
  ksr.output_partition_faces(
    std::back_inserter(all_faces), -1, 6);

  // Model.
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

  // Edges.
  std::string output_filename = "partition-edges.polylines.txt";
  std::ofstream output_file_edges(output_filename);
  output_file_edges.precision(20);
  for (const auto& output_edge : all_edges)
    output_file_edges << "2 " << output_edge << std::endl;
  output_file_edges.close();
  std::cout << "* partition edges exported successfully" << std::endl;

  // Faces.
  output_filename = "partition-faces.ply";
  std::ofstream output_file_faces(output_filename);
  output_file_faces.precision(20);
  if (!CGAL::IO::write_PLY(output_file_faces, all_vertices, all_faces)) {
    std::cerr << "ERROR: can't write to the file " << output_filename << "!" << std::endl;
    return EXIT_FAILURE;
  }
  output_file_faces.close();
  std::cout << "* partition faces exported successfully" << std::endl;

  // Model.
  output_filename = "reconstructed-model.ply";
  std::ofstream output_file_model(output_filename);
  output_file_model.precision(20);
  if (!CGAL::IO::write_PLY(output_file_model, output_vertices, output_faces)) {
    std::cerr << "ERROR: can't write to the file " << output_filename << "!" << std::endl;
    return EXIT_FAILURE;
  }
  output_file_model.close();
  std::cout << "* the reconstructed model exported successfully" << std::endl;

  std::cout << std::endl << "3D KINETIC RECONSTRUCTION " << success <<
  " in " << time << " seconds!" << std::endl << std::endl;
  return EXIT_SUCCESS;
}

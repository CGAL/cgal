#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kinetic_shape_reconstruction_3.h>
#include <CGAL/Kinetic_shape_partitioning_Traits.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/Real_timer.h>
#include <CGAL/IO/PLY.h>
#include <sstream>

#ifdef DISABLED

#include "include/Parameters.h"
#include "include/Terminal_parser.h"

using Kernel    = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;
using FT        = typename Kernel::FT;
using Point_3   = typename Kernel::Point_3;
using Vector_3  = typename Kernel::Vector_3;
using Segment_3 = typename Kernel::Segment_3;

using Point_set    = CGAL::Point_set_3<Point_3>;
using Point_map    = typename Point_set::Point_map;
using Normal_map   = typename Point_set::Vector_map;
using Label_map    = typename Point_set:: template Property_map<int>;
using Semantic_map = CGAL::KSR::Semantic_from_label_map<Label_map>;
using Region_map = typename Point_set:: template Property_map<int>;

using Traits = typename CGAL::Kinetic_shape_partitioning_traits_3<Kernel, EPECK, Point_set, Point_map>;

using KSR = CGAL::Kinetic_shape_reconstruction_3<Traits, Normal_map>;

using Parameters      = CGAL::KSR::All_parameters<FT>;
using Terminal_parser = CGAL::KSR::Terminal_parser<FT>;
using Timer = CGAL::Real_timer;

template <typename T>
std::string to_stringp(const T a_value, const int n = 6)
{
  std::ostringstream out;
  out.precision(n);
  out << std::fixed << a_value;
  return out.str();
}

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

  // Debug.
  parser.add_bool_parameter("-debug", parameters.debug);

  // Verbose.
  parser.add_bool_parameter("-verbose", parameters.verbose);
}

int main(const int argc, const char** argv) {
  // Parameters.
  std::cout.precision(20);
  std::cout << std::endl;
  std::cout << "--- PARSING INPUT: " << std::endl;
  const auto kernel_name = boost::typeindex::type_id<Kernel>().pretty_name();
  std::cout << "* used kernel: " << kernel_name << std::endl;
  const std::string path_to_save = "";
  Terminal_parser parser(argc, argv, path_to_save);

  Parameters parameters;
  parse_terminal(parser, parameters);

  // Check if segmented point cloud already exists.
  std::string filename = parameters.data.substr(parameters.data.find_last_of("/\\") + 1);
  std::string base = filename.substr(0, filename.find_last_of("."));
  base = base + "_" + to_stringp(parameters.distance_threshold, 2) + "_" + to_stringp(parameters.angle_threshold, 2) + "_" + std::to_string(parameters.min_region_size) + ".ply";

  // Input.
  Point_set point_set(parameters.with_normals);
  std::ifstream segmented_file(base);
  if (segmented_file.is_open()) {
    segmented_file >> point_set;
    segmented_file.close();
  }
  else {
    std::ifstream input_file(parameters.data, std::ios_base::binary);
    input_file >> point_set;
    input_file.close();
  }

  for (std::size_t i = 0; i < point_set.size(); i++) {
    Vector_3 n = point_set.normal(i);
    if (abs(n * n) < 0.05)
      std::cout << "point " << i << " does not have a proper normal" << std::endl;
  }

  std::cout << std::endl;
  std::cout << "--- INPUT STATS: " << std::endl;
  std::cout << "* number of points: " << point_set.size() << std::endl;

  std::cout << "verbose " << parameters.verbose << std::endl;
  std::cout << "debug " << parameters.debug << std::endl;

  // Algorithm.
  KSR ksr(point_set, parameters.verbose, parameters.debug);

  const Region_map region_map = point_set. template property_map<int>("region").first;
  const bool is_segmented = point_set. template property_map<int>("region").second;

  Timer timer;
  timer.start();
  std::size_t num_shapes = ksr.detect_planar_shapes(point_set,
    CGAL::parameters::distance_threshold(parameters.distance_threshold)
    .angle_threshold(parameters.angle_threshold)
    .k_neighbors(parameters.k_neighbors)
    .min_region_size(parameters.min_region_size));

  std::cout << num_shapes << " detected planar shapes" << std::endl;

  num_shapes = ksr.regularize_shapes(CGAL::parameters::regularize_parallelism(true).regularize_coplanarity(true).regularize_axis_symmetry(false).regularize_orthogonality(false));

  std::cout << num_shapes << " detected planar shapes after regularization" << std::endl;

  bool is_ksr_success = ksr.initialize_partitioning();

  if (!is_ksr_success) {
    std::cout << "Initializing kinetic partitioning failed!" << std::endl;
    return 1;
  }

  is_ksr_success = ksr.partition(parameters.k_intersections);

  if (!is_ksr_success) {
    std::cout << "Initializing kinetic partitioning failed!" << std::endl;
    return 2;
  }

  ksr.setup_energyterms();

  ksr.reconstruct(parameters.graphcut_beta);
/*
  if (is_segmented)
    is_ksr_success = ksr.reconstruct(
      region_map,
      CGAL::parameters::
      k_neighbors(parameters.k_neighbors).
      distance_threshold(parameters.distance_threshold).
      angle_threshold(parameters.angle_threshold).
      min_region_size(parameters.min_region_size).
      regularize(parameters.regularize).
      k_intersections(parameters.k_intersections).
      graphcut_beta(parameters.graphcut_beta));
  else
    is_ksr_success = ksr.reconstruct(
    base,
    CGAL::parameters::
    k_neighbors(parameters.k_neighbors).
    distance_threshold(parameters.distance_threshold).
    angle_threshold(parameters.angle_threshold).
    min_region_size(parameters.min_region_size).
    regularize(parameters.regularize).
    k_intersections(parameters.k_intersections).
    graphcut_beta(parameters.graphcut_beta));*/
  assert(is_ksr_success);
  const std::string success = is_ksr_success ? "SUCCESS" : "FAILED";
  timer.stop();
  const FT time = static_cast<FT>(timer.time());

  const KSR::KSP& ksp = ksr.partitioning();

  // Output.
  CGAL::Linear_cell_complex_for_combinatorial_map<3, 3> lcc;
  ksp.get_linear_cell_complex(lcc);
/*

  // Vertices.
  std::vector<Point_3> all_vertices;
  ksp.output_partition_vertices(
    std::back_inserter(all_vertices), -1);

  // Edges.
  std::vector<Segment_3> all_edges;
  ksp.output_partition_edges(
    std::back_inserter(all_edges), -1);

  // Faces.
  std::vector< std::vector<std::size_t> > all_faces;
  ksp.output_partition_faces(
    std::back_inserter(all_faces), -1, 6);

  // Model.
  std::vector<Point_3> output_vertices;
  std::vector< std::vector<std::size_t> > output_faces;
  ksp.output_reconstructed_model(
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
  std::cout << "* the reconstructed model exported successfully" << std::endl;*/

  std::cout << std::endl << "3D KINETIC RECONSTRUCTION " << success <<
  " in " << time << " seconds!" << std::endl << std::endl;
  return EXIT_SUCCESS;
}

#endif

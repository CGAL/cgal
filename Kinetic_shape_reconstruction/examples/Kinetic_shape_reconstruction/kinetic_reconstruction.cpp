#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kinetic_shape_reconstruction_3.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/Real_timer.h>
#include <CGAL/IO/PLY.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <sstream>

#include "include/Parameters.h"
#include "include/Terminal_parser.h"

using Kernel    = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK     = CGAL::Exact_predicates_exact_constructions_kernel;
using FT        = typename Kernel::FT;
using Point_3   = typename Kernel::Point_3;
using Vector_3  = typename Kernel::Vector_3;
using Segment_3 = typename Kernel::Segment_3;

using Point_set    = CGAL::Point_set_3<Point_3>;
using Point_map    = typename Point_set::Point_map;
using Normal_map   = typename Point_set::Vector_map;
using Label_map    = typename Point_set:: template Property_map<int>;
using Region_map   = typename Point_set:: template Property_map<int>;


using KSR = CGAL::Kinetic_shape_reconstruction_3<Kernel, Point_set, Point_map, Normal_map>;

using Parameters      = CGAL::KSR::All_parameters<FT>;
using Terminal_parser = CGAL::KSR::Terminal_parser<FT>;
using Timer = CGAL::Real_timer;


double add_polys = 0, intersections = 0, iedges = 0, ifaces = 0, mapping = 0;

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

  if (!point_set.has_normal_map()) {
    point_set.add_normal_map();
    CGAL::pca_estimate_normals<CGAL::Parallel_if_available_tag>(point_set, 9);
    CGAL::mst_orient_normals(point_set, 9);
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

  auto param = CGAL::parameters::maximum_distance(parameters.distance_threshold)
    .maximum_angle(parameters.angle_threshold)
    .k_neighbors(parameters.k_neighbors)
    .minimum_region_size(parameters.min_region_size)
    .distance_tolerance(parameters.distance_threshold * 0.025)
    .debug(parameters.debug)
    .verbose(parameters.verbose)
    .regularize_parallelism(true)
    .regularize_coplanarity(true)
    .regularize_orthogonality(false)
    .regularize_axis_symmetry(false)
    .angle_tolerance(10)
    .maximum_offset(0.01);

  // Algorithm.
  KSR ksr(point_set, param);

/*
  auto rm = point_set. template property_map<int>("region");

  const Region_map region_map = point_set. template property_map<int>("region").value();
  const bool is_segmented = point_set. template property_map<int>("region").second;*/

  Timer timer;
  timer.start();
  std::size_t num_shapes = ksr.detect_planar_shapes(param);

  std::cout << num_shapes << " detected planar shapes" << std::endl;

  FT after_shape_detection = timer.time();

  //num_shapes = ksr.regularize_shapes(CGAL::parameters::regularize_parallelism(true).regularize_coplanarity(true).regularize_axis_symmetry(false).regularize_orthogonality(false));

  //std::cout << num_shapes << " detected planar shapes after regularization" << std::endl;

  ksr.initialize_partition(param);

  std::cout << add_polys << " add polys" << std::endl;
  std::cout << intersections << " intersections" << std::endl;
  std::cout << iedges << " iedges" << std::endl;
  std::cout << ifaces << " ifaces" << std::endl;
  std::cout << mapping << " mapping" << std::endl;

  FT after_init = timer.time();

  FT partition_time, finalization_time, conformal_time;

  ksr.partition(parameters.k_intersections, partition_time, finalization_time, conformal_time);

  FT after_partition = timer.time();

  ksr.setup_energyterms();

  FT after_energyterms = timer.time();

  ksr.reconstruct(parameters.graphcut_beta);
  FT after_reconstruction = timer.time();

  timer.stop();
  const FT time = static_cast<FT>(timer.time());

  std::vector<Point_3> vtx;
  std::vector<std::vector<std::size_t> > polylist;
  ksr.reconstructed_model_polylist(std::back_inserter(vtx), std::back_inserter(polylist));

  CGAL::KSR_3::dump_indexed_polygons(vtx, polylist, "polylist");

  ksr.reconstruct(0.3);

  vtx.clear();
  polylist.clear();
  ksr.reconstructed_model_polylist(std::back_inserter(vtx), std::back_inserter(polylist));

  CGAL::KSR_3::dump_indexed_polygons(vtx, polylist, "polylist_b0.3");
  ksr.reconstruct(0.5);

  vtx.clear();
  polylist.clear();
  ksr.reconstructed_model_polylist(std::back_inserter(vtx), std::back_inserter(polylist));

  CGAL::KSR_3::dump_indexed_polygons(vtx, polylist, "polylist_b0.5");
  ksr.reconstruct(0.7);

  vtx.clear();
  polylist.clear();
  ksr.reconstructed_model_polylist(std::back_inserter(vtx), std::back_inserter(polylist));

  CGAL::KSR_3::dump_indexed_polygons(vtx, polylist, "polylist_b0.7");
  ksr.reconstruct(0.95);

  vtx.clear();
  polylist.clear();
  ksr.reconstructed_model_polylist(std::back_inserter(vtx), std::back_inserter(polylist));

  CGAL::KSR_3::dump_indexed_polygons(vtx, polylist, "polylist_b0.95");

  // Output.
  //ksr.partition().get_linear_cell_complex();

  std::cout << "Shape detection:        " << after_shape_detection << " seconds!" << std::endl;
  std::cout << "Kinetic partition:      " << (after_partition - after_shape_detection) << " seconds!" << std::endl;
  std::cout << " initialization:        " << (after_init - after_shape_detection) << " seconds!" << std::endl;
  std::cout << " partition:             " << (partition_time) << " seconds!" << std::endl;
  std::cout << " finalization:          " << (finalization_time) << " seconds!" << std::endl;
  std::cout << " making conformal:      " << (conformal_time) << " seconds!" << std::endl;
  std::cout << "Kinetic reconstruction: " << (after_reconstruction - after_partition) << " seconds!" << std::endl;
  std::cout << "Total time:             " << time << " seconds!" << std::endl << std::endl;

  return EXIT_SUCCESS;
}

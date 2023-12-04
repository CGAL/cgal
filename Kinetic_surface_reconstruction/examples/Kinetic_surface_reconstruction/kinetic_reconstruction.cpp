#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kinetic_shape_reconstruction_3.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/Real_timer.h>
#include <CGAL/IO/PLY.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/bounding_box.h>
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

  parser.add_str_parameter("-data", parameters.data);

  // Shape detection.
  parser.add_val_parameter("-kn"   , parameters.k_neighbors);
  parser.add_val_parameter("-dist" , parameters.distance_threshold);
  parser.add_val_parameter("-angle", parameters.angle_threshold);
  parser.add_val_parameter("-minp" , parameters.min_region_size);

  parser.add_val_parameter("-odepth", parameters.max_octree_depth);
  parser.add_val_parameter("-osize", parameters.max_octree_node_size);

  // Shape regularization.
  parser.add_bool_parameter("-regularize", parameters.regularize);

  // Shape regularization.
  parser.add_bool_parameter("-reorient", parameters.reorient);

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

  // If no input data is provided, use input from data directory.
  if (parameters.data.empty()) {
    parameters.data = CGAL::data_file_path("points_3/building.ply");
  }

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

  if (parameters.distance_threshold == 0) {
    CGAL::Bbox_3 bbox = CGAL::bbox_3(CGAL::make_transform_iterator_from_property_map(point_set.begin(), point_set.point_map()),
      CGAL::make_transform_iterator_from_property_map(point_set.end(), point_set.point_map()));

    FT d = CGAL::approximate_sqrt
    ((bbox.xmax() - bbox.xmin()) * (bbox.xmax() - bbox.xmin())
      + (bbox.ymax() - bbox.ymin()) * (bbox.ymax() - bbox.ymin())
      + (bbox.zmax() - bbox.zmin()) * (bbox.zmax() - bbox.zmin()));

    parameters.distance_threshold = d * 0.03;
  }

  if (parameters.min_region_size == 0)
    parameters.min_region_size = point_set.size() * 0.01;

  std::cout << std::endl;
  std::cout << "--- INPUT STATS: " << std::endl;
  std::cout << "* number of points: " << point_set.size() << std::endl;

  std::cout << "verbose " << parameters.verbose << std::endl;
  std::cout << "distance_threshold " << parameters.distance_threshold << std::endl;
  std::cout << "angle_threshold " << parameters.angle_threshold << std::endl;
  std::cout << "min_region_size " << parameters.min_region_size << std::endl;
  std::cout << "k " << parameters.k_intersections << std::endl;
  std::cout << "graphcut_beta " << parameters.graphcut_beta << std::endl;

  auto param = CGAL::parameters::maximum_distance(parameters.distance_threshold)
    .maximum_angle(parameters.angle_threshold)
    .k_neighbors(parameters.k_neighbors)
    .minimum_region_size(parameters.min_region_size)
    .debug(parameters.debug)
    .verbose(parameters.verbose)
    .max_octree_depth(parameters.max_octree_depth)
    .max_octree_node_size(parameters.max_octree_node_size)
    .reorient_bbox(false)
    .regularize_parallelism(false)
    .regularize_coplanarity(false)
    .regularize_orthogonality(false)
    .regularize_axis_symmetry(false)
    .angle_tolerance(5)
    .maximum_offset(0.02);

  // Algorithm.
  KSR ksr(point_set, param);

  const Region_map region_map = point_set. template property_map<int>("region").first;
  const bool is_segmented = point_set. template property_map<int>("region").second;

  Timer timer;
  timer.start();
  std::size_t num_shapes = ksr.detect_planar_shapes(false, param);

  std::cout << num_shapes << " detected planar shapes" << std::endl;

  FT after_shape_detection = timer.time();

  ksr.initialize_partition(param);

  FT after_init = timer.time();

  ksr.partition(parameters.k_intersections);

  FT after_partition = timer.time();

  std::vector<Point_3> vtx;
  std::vector<std::vector<std::size_t> > polylist;

  ksr.reconstruct(parameters.graphcut_beta, std::back_inserter(vtx), std::back_inserter(polylist));
  FT after_reconstruction = timer.time();

  if (polylist.size() > 0)
    CGAL::IO::write_polygon_soup("polylist_" + std::to_string(parameters.graphcut_beta), vtx, polylist);

  timer.stop();
  const FT time = static_cast<FT>(timer.time());

  std::vector<FT> betas{0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99};

  for (FT b : betas) {
    if (b == parameters.graphcut_beta)
      continue;

    vtx.clear();
    polylist.clear();

    ksr.reconstruct(b, std::back_inserter(vtx), std::back_inserter(polylist));

    if (polylist.size() > 0)
      CGAL::IO::write_polygon_soup("polylist_" + std::to_string(parameters.graphcut_beta), vtx, polylist);
  }

  std::cout << "Shape detection:        " << after_shape_detection << " seconds!" << std::endl;
  std::cout << "Kinetic partition:      " << (after_partition - after_shape_detection) << " seconds!" << std::endl;
  std::cout << " initialization:        " << (after_init - after_shape_detection) << " seconds!" << std::endl;
  std::cout << " partition:             " << (after_partition - after_init) << " seconds!" << std::endl;
  std::cout << "Kinetic reconstruction: " << (after_reconstruction - after_partition) << " seconds!" << std::endl;
  std::cout << "Total time:             " << time << " seconds!" << std::endl << std::endl;

  return EXIT_SUCCESS;
}

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Kinetic_surface_reconstruction_3.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/Real_timer.h>
#include <CGAL/IO/PLY.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/bounding_box.h>
#include <sstream>
#include <filesystem>
#include <CGAL/boost/graph/IO/OFF.h>
#include <CGAL/boost/graph/graph_traits_Linear_cell_complex_for_combinatorial_map.h>

#include "include/Parameters.h"
#include "include/Terminal_parser.h"

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT = typename Kernel::FT;
using Point_3 = typename Kernel::Point_3;
using Vector_3 = typename Kernel::Vector_3;
using Segment_3 = typename Kernel::Segment_3;
using Plane_3 = typename Kernel::Plane_3;

using Point_set = CGAL::Point_set_3<Point_3>;
using Point_map = typename Point_set::Point_map;
using Normal_map = typename Point_set::Vector_map;

using KSR = CGAL::Kinetic_surface_reconstruction_3<Kernel, Point_set, Point_map, Normal_map>;

using Parameters = CGAL::KSR::All_parameters<FT>;
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
  parser.add_val_parameter("-kn", parameters.k_neighbors);
  parser.add_val_parameter("-dist", parameters.maximum_distance);
  parser.add_val_parameter("-angle", parameters.maximum_angle);
  parser.add_val_parameter("-minp", parameters.min_region_size);


  // Shape regularization.
  parser.add_bool_parameter("-regparallel", parameters.regparallel);
  parser.add_bool_parameter("-regcoplanar", parameters.regcoplanar);
  parser.add_bool_parameter("-regorthogonal", parameters.regorthogonal);
  parser.add_bool_parameter("-regsymmetric", parameters.regsymmetric);

  parser.add_val_parameter("-regoff", parameters.maximum_offset);
  parser.add_val_parameter("-regangle", parameters.angle_tolerance);

  // Shape regularization.
  parser.add_bool_parameter("-reorient", parameters.reorient);

  // Partitioning.
  parser.add_val_parameter("-k", parameters.k_intersections);
  parser.add_val_parameter("-odepth", parameters.max_octree_depth);
  parser.add_val_parameter("-osize", parameters.max_octree_node_size);

  // Reconstruction.
  parser.add_val_parameter("-lambda", parameters.graphcut_lambda);
  parser.add_val_parameter("-ground", parameters.use_ground);

  // Debug.
  parser.add_bool_parameter("-debug", parameters.debug);

  // Verbose.
  parser.add_bool_parameter("-verbose", parameters.verbose);
}

void trim(std::string &str) {
  str.erase(str.find_last_not_of(" ") + 1);
  str.erase(0, str.find_first_not_of(" "));
}

std::size_t load_planes(const std::string &filename, Point_set &ps, std::vector<std::pair<Plane_3, std::vector<typename Point_set::Index>>> &regions) {
  std::ifstream in(filename);
  if (!in.is_open())
    return 0;

  ps.clear();
  regions.clear();

  // scan points
  int lineNumber = 0; // line counter
  std::string line; // line buffer
  std::istringstream iss;

  while (line.empty()) {
    if (!getline(in, line))
      return 0;
    trim(line);
  }

  std::size_t num_points = 0;
  iss.str(line);

  if (!(iss >> std::string("num_points:") >> num_points))
    return 0;

  ps.reserve(num_points);
  lineNumber = 0;
  while (getline(in, line)) {
    trim(line);
    if (line.length() == 0 || line[0] == '#')
      continue;

    FT x, y, z;
    iss.clear();
    iss.str(line);
    if (iss >> CGAL::IO::iformat(x) >> CGAL::IO::iformat(y) >> CGAL::IO::iformat(z)) {
      Point_3 point(x, y, z);
      ps.insert(point);
      lineNumber++;
      if (lineNumber == num_points)
        break;
    }
  }

  line.clear();
  while (line.empty()) {
    if (!getline(in, line))
      return 0;
    trim(line);
  }

  iss.str(line);
  std::size_t num_colors = 0;

  // Check for colors
  auto it = line.find("num_colors:");
  if (it != std::string::npos) {
    lineNumber = 0;
    while (getline(in, line)) {
      trim(line);
      if (line.length() == 0 || line[0] == '#')
        continue;
      lineNumber++;
      if (lineNumber == num_points)
        break;
    }

    line.clear();
    while (line.empty()) {
      if (!getline(in, line))
        return 0;
      trim(line);
    }
  }

  iss.clear();
  iss.str(line);

  std::size_t num_normals = 0;
  if (!(iss >> std::string("num_normals:") >> num_normals))
    return 0;

  lineNumber = 0;
  while (getline(in, line)) {
    trim(line);
    if (line.length() == 0 || line[0] == '#')
      continue;

    FT nx, ny, nz;
    iss.clear();
    iss.str(line);
    if (iss >> CGAL::IO::iformat(nx) >> CGAL::IO::iformat(ny) >> CGAL::IO::iformat(nz)) {
      Vector_3 normal(nx, ny, nz);
      ps.normal(lineNumber) = normal;
      lineNumber++;
      if (lineNumber == num_points)
        break;
    }
  }

  std::size_t num_groups = 0;

  line.clear();
  while (line.empty()) {
    if (!getline(in, line))
      return 0;
    trim(line);
  }

  if (line.empty())
    return 0;

  iss.clear();
  iss.str(line);
  if (!(iss >> std::string("num_groups:") >> num_groups))
    return 0;

  bool check = true;

  std::size_t group_counter = 0;
  regions.reserve(num_groups);
  while (group_counter < num_groups) {
    FT a = 0, b = 0, c = 0, d = 0;
    Plane_3 plane;
    int red = 0, green = 0, blue = 0;
    std::size_t group_num_points = 0;
    std::string s;

    while (getline(in, line)) {
      trim(line);
      if (line.length() == 0 || line[0] == '#')
        continue;

      iss.clear();
      iss.str(line);

      std::vector<std::string> token;
      for (std::string t; iss >> t; ) token.push_back(t);

      if (token[0] == "group_parameters:" && token.size() == 5) {
        a = std::stod(token[1]);
        b = std::stod(token[2]);
        c = std::stod(token[3]);
        d = std::stod(token[4]);
        if (c < 0) {
          a = -a;
          b = -b;
          c = -c;
          d = -d;
        }
        FT l = a * a + b * b + c * c;
        if (abs(l) < 0.9 || abs(l) > 1.1)
          std::cout << "warning: plane normal of region " << group_counter << " is not normalized " << l << std::endl;
        plane = Plane_3(a, b, c, d);
        continue;
      }

      if (token[0] == "group_color:" && token.size() == 4) {
        red = std::stod(token[1]);
        green = std::stod(token[2]);
        blue = std::stod(token[3]);
        continue;
      }

      if (token[0] == "group_num_points:" && token.size() == 2) {
        group_num_points = std::stod(token[1]);
        continue;
      }

      if (line[0] >= '0' && line[0] <= '9') {
        // Indices are coming, check if al necessary parameters are set.
        if (group_num_points == 0 || (a * a + b * b + c * c) == 0) {
          std::cout << "skipping region " << group_counter << " due to missing parameters" << std::endl;
          break;
        }

        regions.push_back(std::make_pair(Plane_3(a, b, c, d), std::vector<typename Point_set::Index>()));
        std::size_t idx = 0;

        regions.back().second.reserve(group_num_points);

        for (const std::string& s : token)
          regions.back().second.push_back(static_cast<typename Point_set::Index>(std::stoi(s)));

        std::vector<Point_3> pts;
        for (auto idx : regions.back().second)
          pts.push_back(ps.point(idx));

        Plane_3 pl;
        CGAL::linear_least_squares_fitting_3(pts.begin(), pts.end(), pl, CGAL::Dimension_tag<0>());
        if (pl.c() < 0) {
          pl = Plane_3(-pl.a(), -pl.b(), -pl.c(), -pl.d());
        }
        regions.back().first = pl;

        group_counter++;

        if (regions.back().second.size() != group_num_points)
          std::cout << "warning: region " << group_counter << " has inconsistent number of points " << group_num_points << " " << regions.back().second.size() << std::endl;
        break;
      }
    }
  }

  return group_counter;
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

  //If no input data is provided, use input from data directory.
  if (parameters.data.empty())
    parameters.data = CGAL::data_file_path("points_3/building.ply");

  Point_set point_set;
  std::vector<std::pair<Plane_3, std::vector<typename Point_set::Index>>> regions;
  //load_planes("../GF-test/saddle_point_cloud.vg", point_set, regions);
  std::filesystem::path path = parameters.data;

  //CGAL::IO::write_point_set("saddle_point.ply", point_set, CGAL::parameters::stream_precision(17).use_binary_mode(false));
  std::string vg_file = path.filename().generic_string() + "_" + to_stringp(parameters.maximum_distance) + "_" + to_stringp(parameters.maximum_angle) + "_" + std::to_string(parameters.min_region_size) + ".vg";
//   if (std::filesystem::exists(vg_file))
//     load_planes(vg_file, point_set, regions);

  // Input.
  if (!regions.empty())
    std::cout << regions.size() << " planar shapes loaded" << std::endl;

  if (regions.empty() && !parameters.data.empty())
    CGAL::IO::read_point_set(parameters.data, point_set);

  if (point_set.size() == 0) {
    std::cout << "input file not found or empty!" << std::endl;
    return EXIT_FAILURE;
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

  if (parameters.maximum_distance == 0) {
    CGAL::Bbox_3 bbox = CGAL::bbox_3(CGAL::make_transform_iterator_from_property_map(point_set.begin(), point_set.point_map()),
      CGAL::make_transform_iterator_from_property_map(point_set.end(), point_set.point_map()));

    FT d = CGAL::approximate_sqrt
    ((bbox.xmax() - bbox.xmin()) * (bbox.xmax() - bbox.xmin())
      + (bbox.ymax() - bbox.ymin()) * (bbox.ymax() - bbox.ymin())
      + (bbox.zmax() - bbox.zmin()) * (bbox.zmax() - bbox.zmin()));

    parameters.maximum_distance = d * 0.03;
  }

  if (parameters.min_region_size == 0)
    parameters.min_region_size = static_cast<std::atomic_size_t>(point_set.size() * 0.01);

  std::cout << std::endl;
  std::cout << "--- INPUT STATS: " << std::endl;
  std::cout << "* number of points: " << point_set.size() << std::endl;

  std::cout << "verbose " << parameters.verbose << std::endl;
  std::cout << "maximum_distance " << parameters.maximum_distance << std::endl;
  std::cout << "maximum_angle " << parameters.maximum_angle << std::endl;
  std::cout << "min_region_size " << parameters.min_region_size << std::endl;
  std::cout << "k " << parameters.k_intersections << std::endl;
  std::cout << "graphcut_lambda " << parameters.graphcut_lambda << std::endl;

  auto param = CGAL::parameters::maximum_distance(parameters.maximum_distance)
    .maximum_angle(parameters.maximum_angle)
    .k_neighbors(parameters.k_neighbors)
    .minimum_region_size(parameters.min_region_size)
    .debug(parameters.debug)
    .verbose(parameters.verbose)
    .max_octree_depth(parameters.max_octree_depth)
    .max_octree_node_size(parameters.max_octree_node_size)
    .reorient_bbox(parameters.reorient)
    .regularize_parallelism(parameters.regparallel)
    .regularize_coplanarity(parameters.regcoplanar)
    .regularize_orthogonality(parameters.regorthogonal)
    .regularize_axis_symmetry(parameters.regsymmetric)
    .angle_tolerance(parameters.angle_tolerance)
    .maximum_offset(parameters.maximum_offset)
    .bbox_dilation_ratio(1.1);

  // Algorithm.
  KSR ksr(point_set, param);

  typename KSR::LCC lcc_input;
  std::string lcc_file = path.filename().generic_string() + "_" + to_stringp(parameters.maximum_distance) + "_" + to_stringp(parameters.maximum_angle) + "_" + std::to_string(parameters.min_region_size) + ".lcc";
/*
  if (!regions.empty() && std::filesystem::exists(lcc_file)) {
    std::ifstream lccfile(lcc_file);
    if (lccfile.is_open()) {
      lccfile >> lcc_input;
      lccfile.close();
    }
  }*/

  Timer timer;
  timer.start();

  if (regions.empty())
    ksr.detect_planar_shapes(param);
  else if (lcc_input.template one_dart_per_cell<3>().empty() && lcc_input.template one_dart_per_cell<3>().empty())
      ksr.insert_planar_shapes(regions);

  if (lcc_input.template one_dart_per_cell<3>().empty()) {
    std::cout << ksr.planar_shapes().size() << " planar shapes regularized into ";
    ksr.regularize_planar_shapes(param);
    std::cout << ksr.planar_shapes().size() << std::endl;
  }

  FT after_shape_detection = 0;
  FT after_partition = 0;

  if (!lcc_input.template one_dart_per_cell<3>().empty())
    ksr.insert_planar_shapes_and_linear_cell_complex(regions, lcc_input);
  else {
    after_shape_detection = timer.time();

    ksr.partition(parameters.k_intersections);

    after_partition = timer.time();

/*
    const typename KSR::LCC& lcc = ksr.get_linear_cell_complex();
    std::ofstream file(lcc_file);
    file << lcc;
    file.close();*/
  }

  std::vector<Point_3> vtx;
  std::vector<std::vector<std::size_t> > polylist;

  std::map<typename KSR::KSP::Face_support, bool> external_nodes;

  if (parameters.use_ground) {
    external_nodes[KSR::KSP::Face_support::ZMIN] = false;
    ksr.reconstruct_with_ground(parameters.graphcut_lambda, std::back_inserter(vtx), std::back_inserter(polylist));
  }
  else
    ksr.reconstruct(parameters.graphcut_lambda, external_nodes, std::back_inserter(vtx), std::back_inserter(polylist));

  FT after_reconstruction = timer.time();

  std::cout << polylist.size() << " polygons, " << vtx.size() << " vertices" << std::endl;

  if (polylist.size() > 0)
    CGAL::IO::write_polygon_soup("polylist_" + std::to_string(parameters.graphcut_lambda) + (parameters.use_ground ? "_g" : "_") + ".off", vtx, polylist);

  timer.stop();
  const FT time = static_cast<FT>(timer.time());

  std::vector<FT> lambdas{ 0.3, 0.4, 0.5, 0.6, 0.7, 0.73, 0.75, 0.77, 0.8, 0.9, 0.95, 0.99 };

  bool non_empty = false;

  for (FT l : lambdas) {
    if (l == parameters.graphcut_lambda)
      continue;

    vtx.clear();
    polylist.clear();

    if (parameters.use_ground)
      ksr.reconstruct_with_ground(l, std::back_inserter(vtx), std::back_inserter(polylist));
    else
      ksr.reconstruct(l, external_nodes, std::back_inserter(vtx), std::back_inserter(polylist));

    if (polylist.size() > 0) {
      non_empty = true;
      CGAL::IO::write_polygon_soup("polylist_" + std::to_string(l) + (parameters.use_ground ? "_g" : "_") + ".off", vtx, polylist);
    }
  }

  std::cout << "Shape detection and initialization\nof kinetic partition:     " << after_shape_detection << " seconds!" << std::endl;
  std::cout << "Kinetic partition:        " << (after_partition - after_shape_detection) << " seconds!" << std::endl;
  std::cout << "Kinetic reconstruction:   " << (after_reconstruction - after_partition) << " seconds!" << std::endl;
  std::cout << "Total time:               " << after_reconstruction << " seconds!" << std::endl << std::endl;

  return (non_empty) ? EXIT_SUCCESS : EXIT_FAILURE;
}

#include "CGAL/_test_pipeline.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/property_map.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <list>
#include <vector>

namespace SS3 = CGAL::Straight_skeletons_3;

using K = CGAL::Exact_predicates_exact_constructions_kernel;
using FT = K::FT;
using Point_3 = K::Point_3;

using Mesh = CGAL::Surface_mesh<Point_3>;
using face_descriptor = typename boost::graph_traits<Mesh>::face_descriptor;

int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  // Argument parsing
  char* mesh_filename;
  char* weights_filename;
  std::filesystem::path save_path = std::filesystem::current_path();
  std::vector<FT> save_times;

  if (!SS3::utils::parse_args(argc, argv, mesh_filename, weights_filename, save_path, save_times)) {
    return EXIT_FAILURE;
  }

  // do nothing if commented
  if (std::string_view(mesh_filename) == "#" || std::string_view(mesh_filename) == "//") {
    return EXIT_SUCCESS;
  }

  std::string mesh_file_stem = std::filesystem::path(mesh_filename).stem().string();
  std::string log_path = "logs/" + mesh_file_stem + ".log";
  std::cout << "Writing to " << log_path << std::endl;
  std::ofstream log(log_path);
  log.precision(17);
  CGAL::Straight_skeletons_3::internal::set_log_stream(log);

  Mesh sm;
  if (!SS3::utils::preprocess_input(mesh_filename, sm)) {
    return EXIT_FAILURE;
  }

  CGAL::IO::write_polygon_mesh("input.obj", sm, CGAL::parameters::stream_precision(17));

  // assign weights
  auto res = sm.add_property_map<face_descriptor, FT>("f:weight");
  if(!res.second) {
    std::cerr << "Error: failed to add property map?" << std::endl;
    return EXIT_FAILURE;
  }

  auto fwm = res.first;
  for (face_descriptor f : faces(sm)) {
    put (fwm, f, 1);
  }

  bool result = SS3::utils::run(sm, fwm, save_times, save_path);
  std::cout << "result: " << result << std::endl;

  std::cout << "Done" << std::endl;
  return result ? EXIT_SUCCESS : EXIT_FAILURE;
}

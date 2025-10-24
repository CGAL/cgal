#include "CGAL/_cardinal_weights.h"
#include "CGAL/_color_input.h"
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

// For the weight file format, check _cardinal_weights.h
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

  // ---

  Mesh sm;
  if (!SS3::utils::preprocess_input(mesh_filename, sm)) {
    return EXIT_FAILURE;
  }

  // assign weights
  auto res = sm.add_property_map<face_descriptor, double>("f:weight");
  if(!res.second) {
    std::cerr << "Error: failed to add property map?" << std::endl;
    return EXIT_FAILURE;
  }

  auto fwm = res.first;

  if (!SS3::utils::assign_cardinal_weights(weights_filename, sm, fwm)) {
    std::cerr << "Error: failed to assign weights " << (weights_filename ? weights_filename : "(default)") << std::endl;
    return EXIT_FAILURE;
  }

#if 0
  // @tmp hardcode weights
  for (face_descriptor f : faces(sm)) {
    bool bottom = true;
    for (auto v : vertices_around_face(halfedge(f, sm), sm)) {
      if (sm.point(v).z() != -2) {
        bottom = false;
        break;
      }
    }

    put (fwm, f, bottom ? 10.0 : 1.0);
  }
  // @tmp
#endif

  bool result = SS3::utils::run(sm, fwm, save_times, save_path);

  std::cout << "Done" << std::endl;
  return result ? EXIT_SUCCESS : EXIT_FAILURE;
}

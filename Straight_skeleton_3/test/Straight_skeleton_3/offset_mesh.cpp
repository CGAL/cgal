#include "CGAL/_cardinal_weights.h"
#include "CGAL/_color_input.h"

#include <CGAL/Straight_skeleton_3/face_offset.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/property_map.h>
#include <CGAL/Real_timer.h>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <list>
#include <vector>

namespace PMP = CGAL::Polygon_mesh_processing;
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

  if (argc == 1) {
    std::cout << "Usage: " << argv[0] << " <input filename> [--weight-path <weights filename>] [--save-path <output directory>] [--save-times <time_0> <time_1> ...]" << std::endl;
    std::cout << "  <input filename>                   Path to input mesh file (required)" << std::endl;
    std::cout << "  --weight-path <weights filename>   Path to weights file (optional)" << std::endl;
    std::cout << "  --save-path <output directory>     Directory to save results (optional, defaults to current directory)" << std::endl;
    std::cout << "  --save-times <times>           List of times (required, all must be non-zero and same sign)" << std::endl;
    std::exit(0);
  }

  // Argument parsing
  const char* mesh_filename = argv[1];
  const char* weights_filename = nullptr;
  std::filesystem::path save_path = std::filesystem::current_path();
  std::vector<std::string> time_strings;

  for (int i = 2; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--weight-path" && i + 1 < argc) {
      std::filesystem::path candidate(argv[i + 1]);
      if (std::filesystem::exists(candidate) && std::filesystem::is_regular_file(candidate)) {
        weights_filename = argv[i + 1];
      } else {
        weights_filename = nullptr;
        std::cerr << "Warning: '" << argv[i + 1] << "' is not a valid file. Using default weights." << std::endl;
      }
      ++i;
    } else if (arg == "--save-path" && i + 1 < argc) {
      std::filesystem::path candidate(argv[i + 1]);
      if (std::filesystem::is_directory(candidate)) {
        save_path = candidate;
      } else {
        std::cerr << "Warning: '" << argv[i + 1] << "' is not a valid directory. Using current directory instead." << std::endl;
        save_path = std::filesystem::current_path();
      }
      ++i;
    } else if (arg == "--save-times") {
      // Read all following values until next flag or end
      int j = i + 1;
      while (j < argc && std::string(argv[j]).find("--") != 0) {
        time_strings.push_back(argv[j]);
        ++j;
      }
      i = j - 1;
    }
  }

  // times at which a result should be produced
  std::set<FT> unique_save_times;
  for (const auto& s : time_strings) {
    try {
      double val = std::stod(s);
      unique_save_times.insert(val);
    } catch (const std::exception& e) {
      std::cerr << "Error: invalid time value '" << s << "'" << std::endl;
      return EXIT_FAILURE;
    }
  }

  // Check that all times are of the same sign (ignoring zeros)
  std::vector<FT> save_times(std::cbegin(unique_save_times), std::cend(unique_save_times));

  bool has_positive_times = false, has_negative_times = false;
  for (auto v : save_times) {
    if (v == 0) {
      std::cerr << "Error: time should be non-zero" << std::endl;
      return EXIT_FAILURE;
    } else if (v > 0) {
      has_positive_times = true;
    } else if (v < 0) {
      has_negative_times = true;
    }
  }

  if (has_positive_times && has_negative_times) {
    std::cerr << "Error: times must all be positive or all negative." << std::endl;
    return EXIT_FAILURE;
  }

  std::sort(save_times.begin(), save_times.end(), [](const FT& a, const FT& b) {
    return CGAL::abs(a) < CGAL::abs(b);
  });

  // ---

  Mesh sm;
  if(!CGAL::IO::read_polygon_mesh(mesh_filename, sm)) {
    std::cerr << "Error: failed to read input " << mesh_filename << std::endl;
    return EXIT_FAILURE;
  }

  if(CGAL::is_empty(sm) || !is_valid_face_graph(sm)) {
    std::cerr << "Error: invalid input " << mesh_filename << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Input mesh: " << vertices(sm).size() << " NV " << faces(sm).size() << " NF" << std::endl;

  // let's tolerate untriangulated inputs
  if (!CGAL::is_triangle_mesh(sm)) {
    PMP::triangulate_faces(sm);
  }

  // some sanity checks
  if (CGAL::is_triangle_mesh(sm) && PMP::does_self_intersect(sm)) {
    std::cerr << "Error: input has self intersections" << std::endl;
    return EXIT_FAILURE;
  }

#if 0 // not an issue for the algorithm, but sometimes useful to analyze the input
  auto vol_id_map = sm.add_property_map<face_descriptor, std::size_t>().first;
  std::size_t vccn = PMP::volume_connected_components(sm, vol_id_map,
                                                      CGAL::parameters::do_orientation_tests(false));
  std::size_t ccn = PMP::internal::number_of_connected_components(sm);
  if (vccn != ccn) {
    std::cerr << "Error: input has nested connected components" << std::endl;
    return EXIT_FAILURE;
  }
#endif

  // assign weights
  auto res = sm.add_property_map<face_descriptor, double>("f:weight");
  if(!res.second) {
    std::cerr << "Error: failed to add property map?" << std::endl;
    return EXIT_FAILURE;
  }

  auto fwm = res.first;

  if (!CGAL::utils::assign_cardinal_weights(weights_filename, sm, fwm)) {
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

  std::vector<Mesh> results;
  results.reserve(save_times.size());

  // Main call
  CGAL::Real_timer timer;
  timer.start();

  bool success = SS3::face_offset(sm, save_times,
                                  results,
                                  CGAL::parameters::face_weight_map(fwm));

  timer.stop();
  std::cout << "Elapsed: " << timer.time() << std::endl;

  // Check the sanity of the results, and save them
  for (std::size_t i=0; i<results.size(); ++i) {
    const FT save_time = save_times[i];
    const Mesh& sm = results[i];

    if (!sm.is_valid() || !is_valid_face_graph(sm)) {
      std::cerr << "Error: broken output" << std::endl;
      return EXIT_FAILURE;
    }
    if (!CGAL::is_closed(sm)) {
      std::cerr << "open output" << std::endl;
      return EXIT_FAILURE;
    }
    if (!CGAL::is_triangle_mesh(sm)) {
      std::cerr << "Error: non-triangulated output" << std::endl;
      return EXIT_FAILURE;
    }
    if (PMP::has_degenerate_faces(sm)) {
      std::cerr << "Error: degenerate faces in output" << std::endl;
      return EXIT_FAILURE;
    }
    if (PMP::does_self_intersect(sm)) {
      std::cerr << "Error: output self-intersects" << std::endl;
      return EXIT_FAILURE;
    }

    std::stringstream out_ss;
    out_ss << save_path.string() << "/result_" << save_time << ".obj";

    // this writes 'double', not FT (aka EPECK::FT), but this is what we want
    if (!CGAL::IO::write_polygon_mesh(out_ss.str(), sm,
                                      CGAL::parameters::stream_precision(17))) {
      std::cerr << "Error: failed to write result " << out_ss.str() << std::endl;
      return EXIT_FAILURE;
    }
  }

  std::cout << "Done" << std::endl;
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}

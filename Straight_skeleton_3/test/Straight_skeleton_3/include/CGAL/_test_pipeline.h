// Copyright (c) 2024-2025 GeometryFactory (France)
// This file is part of CGAL (www.cgal.org)
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_STRAIGHT_SKELETON_3_TEST_OFFSET_PIPELINE_H
#define CGAL_STRAIGHT_SKELETON_3_TEST_OFFSET_PIPELINE_H

#include <CGAL/Straight_skeleton_3/face_offset.h>

#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Real_timer.h>

#include <filesystem>
#include <set>
#include <vector>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace utils {

template <typename FT>
bool parse_args(int argc, char** argv,
                char*& mesh_filename,
                char*& weights_filename,
                std::filesystem::path& save_path,
                std::vector<FT>& save_times)
{
  if (argc == 1) {
    std::cout << "Usage: " << argv[0] << " <input filename> [--weight-path <weights filename>] [--save-path <output directory>] [--save-times <time_0> <time_1> ...]" << std::endl;
    std::cout << "  <input filename>                   Path to input mesh file (required)" << std::endl;
    std::cout << "  --weight-path <weights filename>   Path to weights file (optional)" << std::endl;
    std::cout << "  --save-path <output directory>     Directory to save results (optional, defaults to current directory)" << std::endl;
    std::cout << "  --save-times <times>               List of times (optional, all must be non-zero and same sign)" << std::endl;
    return false;
  }

  mesh_filename = argv[1];

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
      FT val = std::stod(s);
      unique_save_times.insert(val);
    } catch (const std::exception& e) {
      std::cerr << "Error: invalid time value '" << s << "'" << std::endl;
      return false;
    }
  }

  // Check that all times are of the same sign (ignoring zeros)
  save_times = std::vector<FT>(std::cbegin(unique_save_times), std::cend(unique_save_times));

  bool has_positive_times = false, has_negative_times = false;
  for (auto v : save_times) {
    if (v == 0) {
      std::cerr << "Error: time should be non-zero" << std::endl;
      return false;
    } else if (v > 0) {
      has_positive_times = true;
    } else if (v < 0) {
      has_negative_times = true;
    }
  }

  if (has_positive_times && has_negative_times) {
    std::cerr << "Error: times must all be positive or all negative." << std::endl;
    return false;
  }

  std::sort(save_times.begin(), save_times.end(), [](const FT& a, const FT& b) {
    return CGAL::abs(a) < CGAL::abs(b);
  });

  return true;
}

template <typename Mesh>
bool preprocess_input(const char* mesh_filename,
                      Mesh& sm)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  if(!CGAL::IO::read_polygon_mesh(mesh_filename, sm)) {
    std::cerr << "Error: failed to read input " << mesh_filename << std::endl;
    return false;
  }

  if(CGAL::is_empty(sm) || !is_valid_face_graph(sm)) {
    std::cerr << "Error: invalid input " << mesh_filename << std::endl;
    return false;
  }

  std::cout << "Input mesh: " << vertices(sm).size() << " NV " << faces(sm).size() << " NF" << std::endl;

  // let's tolerate untriangulated inputs
  if (!CGAL::is_triangle_mesh(sm)) {
    PMP::triangulate_faces(sm);
  }

  // some sanity checks
  if (CGAL::is_triangle_mesh(sm) && PMP::does_self_intersect(sm)) {
    std::cerr << "Error: input has self intersections" << std::endl;
    return false;
  }

#if 0 // not an issue for the algorithm, but sometimes useful to analyze the input
  auto vol_id_map = sm.add_property_map<face_descriptor, std::size_t>().first;
  std::size_t vccn = PMP::volume_connected_components(sm, vol_id_map,
                                                      CGAL::parameters::do_orientation_tests(false));
  std::size_t ccn = PMP::internal::number_of_connected_components(sm);
  if (vccn != ccn) {
    std::cerr << "Error: input has nested connected components" << std::endl;
    return false;
  }
#endif

  return true;
}

template <typename Mesh,
          typename FaceWeightMap,
          typename FT>
bool test_datum(const Mesh& sm,
                const FaceWeightMap fwm,
                const std::vector<FT>& save_times,
                const std::filesystem::path& save_path)
{
  namespace PMP = CGAL::Polygon_mesh_processing;
  namespace SS3 = CGAL::Straight_skeletons_3;

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

  if (!success) {
    std::cerr << "Error: face_offset() returned 'false'" << std::endl;
    return false;
  }

  // Check the sanity of the results, and save them
  for (std::size_t i=0; i<results.size(); ++i) {
    const FT save_time = save_times[i];
    const Mesh& sm = results[i];

    if (!sm.is_valid() || !is_valid_face_graph(sm)) {
      std::cerr << "Error: broken output" << std::endl;
      return false;
    }
    if (!CGAL::is_closed(sm)) {
      std::cerr << "open output" << std::endl;
      return false;
    }
    if (!CGAL::is_triangle_mesh(sm)) {
      std::cerr << "Error: non-triangulated output" << std::endl;
      return false;
    }
    if (PMP::has_degenerate_faces(sm)) {
      std::cerr << "Error: degenerate faces in output" << std::endl;
      return false;
    }
    if (PMP::does_self_intersect(sm)) {
      std::cerr << "Error: output self-intersects" << std::endl;
      return false;
    }

    std::stringstream out_ss;
    out_ss << save_path.string() << "/result_" << save_time << ".obj";

    // this writes 'double', not FT (aka EPECK::FT), but this is what we want
    if (!CGAL::IO::write_polygon_mesh(out_ss.str(), sm,
                                      CGAL::parameters::stream_precision(17))) {
      std::cerr << "Error: failed to write result " << out_ss.str() << std::endl;
      return false;
    }
  }

  return true;
}

} // namespace utils
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif // CGAL_STRAIGHT_SKELETON_3_TEST_OFFSET_PIPELINE_H

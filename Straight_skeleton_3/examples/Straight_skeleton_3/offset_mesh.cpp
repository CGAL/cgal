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

namespace CGAL {
namespace utils {

// Simple helper function to draw a mesh whose faces are colored according to the weights (speeds).
template<typename PolygonMesh, typename Values>
void save_colored_mesh(const PolygonMesh& pmesh,
                       const Values& values,
                       const std::string fullpath)
{
  using Color = CGAL::IO::Color;

  using face_descriptor = typename boost::graph_traits<PolygonMesh>::face_descriptor;

  using value_type = typename CGAL::cpp20::remove_cvref<decltype(values[face_descriptor()])>::type;

  std::cout << "Saving colored mesh to " << fullpath << std::endl;

  // get a unique vector of values
  std::vector<value_type> unique_values;
  for (auto f : faces(pmesh)) {
    unique_values.push_back(values[f]);
  }

  std::sort(unique_values.begin(), unique_values.end());
  unique_values.erase(std::unique(unique_values.begin(), unique_values.end()), unique_values.end());

  srand(static_cast<unsigned int>(time(NULL)));

  std::map<std::size_t, CGAL::Color> colors;
  for (const auto& value : unique_values) {
    colors[value] = Color(static_cast<unsigned char>(rand() % 256),
                          static_cast<unsigned char>(rand() % 256),
                          static_cast<unsigned char>(rand() % 256));
    std::cout << " value " << value << " has color " << colors[value] << std::endl;
  }

  auto& nc_pmesh = const_cast<PolygonMesh&>(pmesh);
  auto face_color = nc_pmesh.template add_property_map<face_descriptor, Color>("f:color").first;

  for (auto f : faces(pmesh)) {
    std::cout << "facet " << f << " with value " << values[f]
              << " gets color " << colors[values[f]] << std::endl;
    put(face_color, f, colors[values[f]]);
  }

  std::ofstream out(fullpath);
  CGAL::IO::write_PLY(out, pmesh, CGAL::parameters::face_color_map(face_color));
}

// This is an helper function that computes the weight associated to a face
// using its normal and the provided weights in the x, y and z directions.
// Weights are written in the property map `fwm`.
//
// \tparam PolygonMesh must be a model of `FaceListGraph`
// \tparam FaceWeightMap must be a model of `ReadWritePropertyMap`
//
// \param weights_filename the name of the file containing the weights
//                         in the format:
//                         x1: <value> x2: <value> y1: <value> y2: <value>
//                         [bottom: <value> top: <value>]
// \param pmesh the polygon mesh to assign weights to
// \param fwm the face weight map to write the weights to
//
// \return `true` if all weights could be read, are valid, and were successfully assigned; `false` otherwise.
template <typename PolygonMesh, typename FaceWeightMap>
bool assign_weights(const char* weights_filename,
                    const PolygonMesh& pmesh,
                    FaceWeightMap fwm)
{
  using face_descriptor = typename boost::graph_traits<PolygonMesh>::face_descriptor;

  using GT = typename CGAL::GetGeomTraits<PolygonMesh>::type;
  using FT = typename GT::FT;
  using Point_3 = typename GT::Point_3;
  using Vector_3 = typename GT::Vector_3;

  std::ifstream weights_in(weights_filename);

  if (!weights_filename || !weights_in) {
    CGAL_SS3_TRACE_V(1, "Warning: no input weights provided; all weights are set to '1'.");
    for (face_descriptor f : faces(pmesh))
      put(fwm, f, 1.);

    return true;
  }

  std::string x1_str, x2_str, y1_str, y2_str, bot_str, top_str;
  FT x1_val, x2_val, y1_val, y2_val, bot_val, top_val;
  x1_val = x2_val = y1_val = y2_val = bot_val = top_val = 0;

  if(!(weights_in >> x1_str >> x1_val
                  >> x2_str >> x2_val
                  >> y1_str >> y1_val
                  >> y2_str >> y2_val)) {
    std::cerr << "Error: failed to read weights" << std::endl;
    return false;
  }

  CGAL_assertion(x1_str == "x1:" && x2_str == "x2:" && y1_str == "y1:" && y2_str == "y2:");

  if(weights_in >> bot_str >> bot_val
                >> top_str >> top_val) {
    CGAL_SS3_TRACE_V(8, "bottom & top weight info detected");
    CGAL_assertion(bot_str == "bottom:" && top_str == "top:");
  } else {
    if(x2_val != y2_val) {
      std::cerr << "Warning: unknown z-speeds, and x-speed and y-speed differ..." << std::endl;
      // arbitrary choice
      top_val = y2_val;
    } else {
      // assign the uniform speed to the top
      top_val = x2_val;
    }
  }

  if(x1_val < 0 || x2_val < 0 || y1_val < 0 || y2_val < 0 || bot_val < 0 || top_val < 0) {
    std::cerr << "Error: negative weights are not allowed" << std::endl;
    return false;
  }

  CGAL_SS3_TRACE_V(8, "x1_val = " << x1_val);
  CGAL_SS3_TRACE_V(8, "x2_val = " << x2_val);
  CGAL_SS3_TRACE_V(8, "y1_val = " << y1_val);
  CGAL_SS3_TRACE_V(8, "y2_val = " << y2_val);
  CGAL_SS3_TRACE_V(8, "bot_val = " << bot_val);
  CGAL_SS3_TRACE_V(8, "top_val = " << top_val);

  if(is_zero(x1_val) && is_zero(x2_val) &&
     is_zero(y1_val) && is_zero(y2_val) &&
     is_zero(bot_val) && is_zero(top_val)) {
    std::cerr << "Error: all weights are zero" << std::endl;
    return false;
  }

  FT eps_weight = std::numeric_limits<double>::max(); // 'double' on purpose
  if(x1_val > 0) eps_weight = (std::min)(eps_weight, x1_val);
  if(x2_val > 0) eps_weight = (std::min)(eps_weight, x2_val);
  if(y1_val > 0) eps_weight = (std::min)(eps_weight, y1_val);
  if(y2_val > 0) eps_weight = (std::min)(eps_weight, y2_val);
  if(bot_val > 0) eps_weight = (std::min)(eps_weight, bot_val);
  if(top_val > 0) eps_weight = (std::min)(eps_weight, top_val);

  CGAL_SS3_TRACE_V(8, "min_weight = " << eps_weight);

  // @todo handle true zero
  eps_weight = 1e-10 * eps_weight;

  if(x1_val == 0) { CGAL_SS3_TRACE_V(16, "x1_val to eps weight " << eps_weight); x1_val = eps_weight; }
  if(x2_val == 0) { CGAL_SS3_TRACE_V(16, "x2_val to eps weight " << eps_weight); x2_val = eps_weight; }
  if(y1_val == 0) { CGAL_SS3_TRACE_V(16, "y1_val to eps weight " << eps_weight); y1_val = eps_weight; }
  if(y2_val == 0) { CGAL_SS3_TRACE_V(16, "y2_val to eps weight " << eps_weight); y2_val = eps_weight; }
  if(bot_val == 0) { CGAL_SS3_TRACE_V(16, "bot_val to eps weight " << eps_weight); bot_val = eps_weight; }
  if(top_val == 0) { CGAL_SS3_TRACE_V(16, "top_val to eps weight " << eps_weight); top_val = eps_weight; }

  for (face_descriptor f : faces(pmesh))
  {
    // internal stuff, we don't need to normalize and introduce inexactness
    Vector_3 v = CGAL::NULL_VECTOR;
    PMP::internal::sum_normals<Point_3>(pmesh, f, get(CGAL::vertex_point, pmesh), v, GT());
    FT sq_n = v.squared_length();

    CGAL_SS3_TRACE_V(16, "facet: " << f << " normal: " << v);

#if 1
    FT sq_cos_x = CGAL::square(v.x()) / sq_n;
    FT sq_cos_y = CGAL::square(v.y()) / sq_n;
    FT sq_cos_z = CGAL::square(v.z()) / sq_n;

    // The weight is a weighted sum of all speed contributions, where the weights are
    // the squared cosines of the angles between the normal and the axes
    FT weight_x = (v.x() >= 0) ? x1_val : x2_val;
    FT weight_y = (v.y() >= 0) ? y1_val : y2_val;
    FT weight_z = (v.z() >= 0) ? top_val : bot_val;
    FT weight = weight_x*sq_cos_x + weight_y*sq_cos_y + weight_z*sq_cos_z;
#else
    if(v.x() == 0 && v.y() == 0) {
      if(v.z() > 0)
        weight = vz2;
      else
        weight = vz1;
    } else {
      if(v.x() >= 0) {
        const FT sq_cos = CGAL::square(CGAL::scalar_product(v, west)) / sq_n;
        if(v.y() >= 0) {
          // north east quadrant
          weight = vy2 * (1 - sq_cos) + vx2 * sq_cos;
        } else {
          // south east quadrant
          weight = vy1 * (1 - sq_cos) + vx2 * sq_cos;
        }
      } else { // x < 0
        const FT sq_cos = CGAL::square(CGAL::scalar_product(v, east)) / sq_n;
        if(v.y() >= 0) {
          // north west quadrant
          weight = vy2 * (1 - sq_cos) + vx1 * sq_cos;
        } else {
          // south west quadrant
          weight = vy1 * (1 - sq_cos) + vx1 * sq_cos;
        }
      }
    }
#endif

    CGAL_SS3_TRACE_V(16, "facet: " << f << " weight: " << weight);

    // @todo currently 'double', but only because the run_and_compare.sh pipeline
    // with multiple .cpp needs to save to a file and the f:weight pmap will not
    // be detected if its value type is, e.g., EPECK::FT
    // Could be FT if there were no intermediary saving
    put(fwm, f, CGAL::to_double(weight));
    CGAL_postcondition(get(fwm, f) > 0);
  }

  CGAL_SS3_TRACE_V(8, "E-W-S-N weights: " << x1_val << " " << x2_val << " " << y1_val << " " << y2_val);

  utils::save_colored_mesh(pmesh, fwm, "results/weighted.ply");

  return true;
}

} // namespace utils
} // namespace CGAL

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
    std::cout << "Usage: " << argv[0] << " <input filename> [--weight-path <weights filename>] [--save-path <output directory>] [--save-offsets <offset_0> <offset_1> ...]" << std::endl;
    std::cout << "  <input filename>                   Path to input mesh file (required)" << std::endl;
    std::cout << "  --weight-path <weights filename>   Path to weights file (optional)" << std::endl;
    std::cout << "  --save-path <output directory>     Directory to save results (optional, defaults to current directory)" << std::endl;
    std::cout << "  --save-offsets <offsets>           List of offsets (required, all must be non-zero and same sign)" << std::endl;
    std::exit(0);
  }

  // Argument parsing
  const char* mesh_filename = argv[1];
  const char* weights_filename = nullptr;
  std::filesystem::path save_path = std::filesystem::current_path();
  std::vector<std::string> offset_strings;

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
    } else if (arg == "--save-offsets") {
      // Read all following values until next flag or end
      int j = i + 1;
      while (j < argc && std::string(argv[j]).find("--") != 0) {
        offset_strings.push_back(argv[j]);
        ++j;
      }
      i = j - 1;
    }
  }

  // offsets at which a result should be produced
  std::set<FT> unique_save_offsets;
  for (const auto& s : offset_strings) {
    try {
      double val = std::stod(s);
      unique_save_offsets.insert(val);
    } catch (const std::exception& e) {
      std::cerr << "Error: invalid offset value '" << s << "'" << std::endl;
      return EXIT_FAILURE;
    }
  }

  // Check that all offsets are of the same sign (ignoring zeros)
  std::vector<FT> save_offsets(std::cbegin(unique_save_offsets), std::cend(unique_save_offsets));

  bool has_positive_offsets = false, has_negative_offsets = false;
  for (auto v : save_offsets) {
    if (v == 0) {
      std::cerr << "Error: offset should be non-zero" << std::endl;
      return EXIT_FAILURE;
    } else if (v > 0) {
      has_positive_offsets = true;
    } else if (v < 0) {
      has_negative_offsets = true;
    }
  }

  if (has_positive_offsets && has_negative_offsets) {
    std::cerr << "Error: offsets must all be positive or all negative." << std::endl;
    return EXIT_FAILURE;
  }

  std::sort(save_offsets.begin(), save_offsets.end(), [](const FT& a, const FT& b) {
    return CGAL::abs(a) < CGAL::abs(b);
  });

  // ---

  Mesh sm;
  if(!CGAL::IO::read_polygon_mesh(mesh_filename, sm) || CGAL::is_empty(sm) || !is_valid_face_graph(sm)) {
    std::cerr << "Error: failed to read input " << mesh_filename << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Input mesh: " << vertices(sm).size() << " NV " << faces(sm).size() << " NF" << std::endl;

  // let's tolerate untriangulated inputs
  if (!CGAL::is_triangle_mesh(sm)) {
    PMP::triangulate_faces(sm);
  }

  // some sanity checks
  if (PMP::does_self_intersect(sm)) {
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

  if (!CGAL::utils::assign_weights(weights_filename, sm, fwm)) {
    std::cerr << "Error: failed to assign weights " << (weights_filename ? weights_filename : "(default)") << std::endl;
    return EXIT_FAILURE;
  }

  std::vector<Mesh> results;
  results.reserve(save_offsets.size());

  // Main call
  CGAL::Real_timer timer;
  timer.start();

  bool success = SS3::face_offset(sm, save_offsets,
                                  results,
                                  CGAL::parameters::face_weight_map(fwm));

  timer.stop();
  std::cout << "Elapsed: " << timer.time() << std::endl;

  // Check the sanity of the results, and save them
  for (std::size_t i=0; i<results.size(); ++i) {
    const FT save_offset = save_offsets[i];
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
    out_ss << save_path.string() << "/result_" << save_offset << ".obj";

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

#include "util/Configuration.h"
#include "cgal_kernel.h"
#include "algo/3d/SimpleStraightSkel.h"
#include "algo/3d/PolyhedronTransformation.h"
#include "algo/3d/SelfIntersection.h"
#include "data/3d/Polyhedron.h"
#include "db/3d/Surface_meshIO.h"

#include <CGAL/Surface_mesh.h>
#include <CGAL/property_map.h>

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

using Mesh = CGAL::Surface_mesh<CGAL::Point3>;
using face_descriptor = typename boost::graph_traits<Mesh>::face_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

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

  std::ifstream weights_in(weights_filename);

  if(!weights_filename || !weights_in) {
    CGAL_SS3_TRACE("Warning: no input weights provided; all weights are set to '1'.");
    for(face_descriptor f : faces(pmesh))
      put(fwm, f, 1.);

    return true;
  }

  std::string x1_str, x2_str, y1_str, y2_str, bot_str, top_str;
  CGAL::FT x1_val, x2_val, y1_val, y2_val, bot_val, top_val;
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
    CGAL_SS3_TRACE("bottom & top weight info detected");
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

  CGAL_SS3_TRACE("x1_val = " << x1_val);
  CGAL_SS3_TRACE("x2_val = " << x2_val);
  CGAL_SS3_TRACE("y1_val = " << y1_val);
  CGAL_SS3_TRACE("y2_val = " << y2_val);
  CGAL_SS3_TRACE("bot_val = " << bot_val);
  CGAL_SS3_TRACE("top_val = " << top_val);

  if(is_zero(x1_val) && is_zero(x2_val) &&
     is_zero(y1_val) && is_zero(y2_val) &&
     is_zero(bot_val) && is_zero(top_val)) {
    std::cerr << "Error: all weights are zero" << std::endl;
    return false;
  }

  CGAL::FT eps_weight = std::numeric_limits<double>::max(); // 'double' on purpose
  if(x1_val > 0) eps_weight = (std::min)(eps_weight, x1_val);
  if(x2_val > 0) eps_weight = (std::min)(eps_weight, x2_val);
  if(y1_val > 0) eps_weight = (std::min)(eps_weight, y1_val);
  if(y2_val > 0) eps_weight = (std::min)(eps_weight, y2_val);
  if(bot_val > 0) eps_weight = (std::min)(eps_weight, bot_val);
  if(top_val > 0) eps_weight = (std::min)(eps_weight, top_val);

  CGAL_SS3_TRACE("min_weight = " << eps_weight);

  // @todo handle true zero
  eps_weight = 1e-10 * eps_weight;

  if(x1_val == 0) { CGAL_SS3_TRACE("x1_val to eps weight " << eps_weight); x1_val = eps_weight; }
  if(x2_val == 0) { CGAL_SS3_TRACE("x2_val to eps weight " << eps_weight); x2_val = eps_weight; }
  if(y1_val == 0) { CGAL_SS3_TRACE("y1_val to eps weight " << eps_weight); y1_val = eps_weight; }
  if(y2_val == 0) { CGAL_SS3_TRACE("y2_val to eps weight " << eps_weight); y2_val = eps_weight; }
  if(bot_val == 0) { CGAL_SS3_TRACE("bot_val to eps weight " << eps_weight); bot_val = eps_weight; }
  if(top_val == 0) { CGAL_SS3_TRACE("top_val to eps weight " << eps_weight); top_val = eps_weight; }

  for(face_descriptor f : faces(pmesh))
  {
    // internal stuff, we don't need to normalize and introduce inexactness
    Vector3 v = CGAL::NULL_VECTOR;
    CGAL::Polygon_mesh_processing::internal::sum_normals<Point3>(pmesh, f, get(CGAL::vertex_point, pmesh), v, CGAL::K());
    CGAL::FT sq_n = v.squared_length();

    CGAL_SS3_TRACE("facet: " << f << " normal: " << v);

#if 1
    CGAL::FT sq_cos_x = CGAL::square(v.x()) / sq_n;
    CGAL::FT sq_cos_y = CGAL::square(v.y()) / sq_n;
    CGAL::FT sq_cos_z = CGAL::square(v.z()) / sq_n;

    // The weight is a weighted sum of all speed contributions, where the weights are
    // the squared cosines of the angles between the normal and the axes
    CGAL::FT weight_x = (v.x() >= 0) ? x1_val : x2_val;
    CGAL::FT weight_y = (v.y() >= 0) ? y1_val : y2_val;
    CGAL::FT weight_z = (v.z() >= 0) ? top_val : bot_val;
    CGAL::FT weight = weight_x*sq_cos_x + weight_y*sq_cos_y + weight_z*sq_cos_z;
#else
    if(v.x() == 0 && v.y() == 0) {
      if(v.z() > 0)
        weight = vz2;
      else
        weight = vz1;
    } else {
      if(v.x() >= 0) {
        const CGAL::FT sq_cos = CGAL::square(CGAL::scalar_product(v, west)) / sq_n;
        if(v.y() >= 0) {
          // north east quadrant
          weight = vy2 * (1 - sq_cos) + vx2 * sq_cos;
        } else {
          // south east quadrant
          weight = vy1 * (1 - sq_cos) + vx2 * sq_cos;
        }
      } else { // x < 0
        const CGAL::FT sq_cos = CGAL::square(CGAL::scalar_product(v, east)) / sq_n;
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

    CGAL_SS3_TRACE("facet: " << f << " weight: " << weight);

    // @todo currently 'double', but only because the run_and_compare.sh pipeline
    // with multiple .cpp needs to save to a file and the f:weight pmap will not
    // be detected if its value type is, e.g., EPECK::FT
    // Could be CGAL::FT if there were no intermediary saving
    put(fwm, f, CGAL::to_double(weight));
    CGAL_postcondition(get(fwm, f) != 0);
  }

  CGAL_SS3_TRACE("E-W-S-N weights: " << x1_val << " " << x2_val << " " << y1_val << " " << y2_val);

  utils::save_colored_mesh(pmesh, fwm, "results/weighted.ply");

  return true;
}

} // namespace utils

// Given a range of values, this function generates face-offsets of an input mesh.
// The mesh is constructed by translating the faces of the input mesh by a distance
// that depends on the specified offset and their respective weights (speeds); during
// offsetting, elements of the mesh interact with each other (merging, splitting, etc.).
//
// Positive offset values correspond to outward offsets, while negative values correspond to inward offsets.
//
// \warning A tiny perturbation is always applied to the input mesh's geometry as to avoid so-called
// degenerate positions, meaning a configuration where any three supporting planes of the facets
// of the input mesh do not intersect in a single point.
//
// \tparam TriangleMesh must be a model of `MutableFaceGraph`, `FaceListGraph`, `HalfedgeListGraph`
// \tparam NamedParametersIn must be a model of `NamedParameters`
// \tparam NamedParametersOut must be a model of `NamedParameters`
//
// \param tmesh the input mesh to offset
// \param save_offsets the offsets to apply to the mesh
// \param results the output vector of meshes that will contain the offset meshes
// \param np_in the input named parameters
// \param np_out the output named parameters
//
// \pre offsets value should be non-zero and all of the same sign.
//
// \return `true` if offset meshes were successfully constructed; `false` otherwise.
template <typename TriangleMesh,
          typename NamedParametersIn = parameters::Default_named_parameters,
          typename NamedParametersOut = parameters::Default_named_parameters>
bool face_offset(TriangleMesh& tmesh,
                 std::vector<CGAL::FT> save_offsets, // intentional copy
                 std::vector<TriangleMesh>& results,
                 const NamedParametersIn& np_in = parameters::default_values(),
                 const NamedParametersOut& np_out = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  const std::filesystem::path save_path = choose_parameter(get_parameter(np_in, internal_np::io_path),
                                                           std::filesystem::current_path());

  CGAL_precondition(!CGAL::is_empty(tmesh));
  CGAL_precondition(CGAL::is_closed(tmesh));
  CGAL_precondition(CGAL::is_triangle_mesh(tmesh));
  CGAL_precondition(!PMP::does_self_intersect(tmesh));

  const bool outwards = (!save_offsets.empty() && CGAL::is_positive(save_offsets.front()));

  // The underlying implementation always shrinks the mesh, so for outwards offsets,
  // we reverse the face orientations, whether the mesh was outward oriented or not.
  if (outwards) {
    CGAL_SS3_TRACE("Reversing face orientations...");
    PMP::reverse_face_orientations(tmesh);

    // also switch to negative offsets
    for (CGAL::FT& offset : save_offsets) {
      offset = -offset;
    }
  }

  // Convert the suface mesh into the SLS3-specific data structure that allows faces with multiple
  // borders and disconnected facet connected components
  algo::_3d::PolyhedronSPtr p = algo::_3d::PolyhedronTransformation::convert(tmesh, np_in);

  CGAL_SS3_TRACE("Post merge: " << p->vertices().size() << " NV " << p->facets().size() << " NF");

  // algo::_3d::PolyhedronTransformation::truncatePrecision(p);

  // @todo In safe mode, this should be:
  // if (bad) -> reduce amplitude of perturbation and try again"
  algo::_3d::PolyhedronTransformation::normalizeFacetPlanes(p);
  algo::_3d::PolyhedronTransformation::randTiltPlanesv3(p);

  CGAL_SS3_TRACE("Post randomization: " << p->vertices().size() << " NV " << p->facets().size() << " NF");

  // visitor to collect the results
  struct Mesh_collector_visitor
      : public  algo::_3d::utils::Default_mesh_offset_visitor
  {
      Mesh_collector_visitor(std::vector<algo::_3d::PolyhedronSPtr>& results) : results_(results) { }

      void on_save_offset_event(algo::_3d::PolyhedronSPtr polyhedron, CGAL::FT) override {
          algo::_3d::PolyhedronSPtr other_polyhedron = polyhedron->clone();
          results_.push_back(other_polyhedron);
      }

  private:
      std::vector<algo::_3d::PolyhedronSPtr>& results_;
  };

  // run the skeleton code
  algo::ControllerSPtr controller = { };
  algo::_3d::SimpleStraightSkelSPtr algoskel3d = algo::_3d::SimpleStraightSkel::create(p, controller, save_offsets, save_path);

  std::vector<algo::_3d::PolyhedronSPtr> results_p;
  Mesh_collector_visitor visitor(results_p);
  algoskel3d->setVisitor(&visitor);

  bool success = algoskel3d->run();
  if (!success) {
      return false;
  }

  if (results_p.size() != save_offsets.size()) {
    return false;
  }

  if (outwards) {
    for (CGAL::FT& offset : save_offsets) {
      offset = -offset;
    }
  }

  for (std::size_t i=0; i<results_p.size(); ++i) {
    const CGAL::FT save_offset = save_offsets[i];

    CGAL_SS3_TRACE("Post processing of result @ " << save_offset)

    algo::_3d::PolyhedronSPtr result_p = results_p[i];
    CGAL_assertion(result_p && result_p->isConsistent());

    CGAL_SS3_TRACE("At offset: " << save_offset << ", result with " << result_p->vertices().size()
                     << " vertices and " << result_p->facets().size() << " faces");

    // Convert back to Surface_mesh structure to save the results.
    // This could be avoided if a polygonal output was prefered, but then we
    // need some specific conversion code.
    TriangleMesh result_t;
    bool success = db::_3d::Surface_meshIO::save(result_p, result_t,
                                                 CGAL::parameters::do_not_triangulate_faces(false) /*triangulate*/);
    if (!success) {
      std::cerr << "Error: failed to convert back to Surface_mesh" << std::endl;
      return false;
    }

    CGAL_postcondition(result_t.is_valid());
    CGAL_postcondition(is_valid_face_graph(result_t));
    CGAL_postcondition(CGAL::is_closed(result_t));
    CGAL_postcondition(CGAL::is_triangle_mesh(result_t));
    CGAL_postcondition(!PMP::has_degenerate_faces(result_t));
    CGAL_postcondition(!PMP::does_self_intersect(result_t));

    CGAL_SS3_TRACE("At offset: " << save_offset << ", result with " << num_vertices(result_t)
                     << " vertices and " << num_faces(result_t) << " faces");

    if (outwards) {
      CGAL_SS3_TRACE("Reversing face orientations...");
      PMP::reverse_face_orientations(result_t);
    }

    results.push_back(result_t);
  }

  CGAL_SS3_TRACE("Offset mesh(es) constructed");
  return true;
}

} // namespace CGAL

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

  // SLS3 config file
  util::ConfigurationSPtr config = util::Configuration::getInstance();
  std::string str_conf_file = config->findDefaultFilename();
  if (!config->load(str_conf_file)) {
      std::cerr << "Error: Config file '" << str_conf_file << "' not found." << std::endl;
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
  std::set<CGAL::FT> unique_save_offsets;
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
  std::vector<CGAL::FT> save_offsets(std::cbegin(unique_save_offsets), std::cend(unique_save_offsets));

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

  std::sort(save_offsets.begin(), save_offsets.end(), [](const CGAL::FT& a, const CGAL::FT& b) {
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

#if 0 // not an issue, but sometimes useful to debug
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
  bool success = CGAL::face_offset(sm, save_offsets,
                                   results,
                                   CGAL::parameters::face_weight_map(fwm));

  // Check the sanity of the results, and save them
  for (std::size_t i=0; i<results.size(); ++i) {
    const CGAL::FT save_offset = save_offsets[i];
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

    // this writes 'double', not CGAL::FT (aka EPECK::FT), but this is what we want
    if (!CGAL::IO::write_polygon_mesh(out_ss.str(), sm,
                                      CGAL::parameters::stream_precision(17))) {
      std::cerr << "Error: failed to write result " << out_ss.str() << std::endl;
      return EXIT_FAILURE;
    }
  }

  std::cout << "Done" << std::endl;
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}

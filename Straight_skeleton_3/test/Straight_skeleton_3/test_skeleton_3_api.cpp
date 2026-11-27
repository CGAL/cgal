#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/create_offset_polyhedra.h>
#include <CGAL/create_straight_skeleton_3.h>
#include <CGAL/Straight_skeleton_3.h>
#include <CGAL/Straight_skeleton_3/IO/OBJ.h>

#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/property_map.h>

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace PMP = CGAL::Polygon_mesh_processing;

using K = CGAL::Exact_predicates_exact_constructions_kernel;
using FT = K::FT;
using Point_3 = K::Point_3;

using Mesh = CGAL::Surface_mesh<Point_3>;
using Polyhedron = CGAL::Polyhedron_3<K>;

using Straight_skeleton_3 = CGAL::Straight_skeleton_3<K>;
using NodeSPtr = std::shared_ptr<Straight_skeleton_3::Node>;
using SpokeSPtr = std::shared_ptr<Straight_skeleton_3::Spoke>;

bool check(bool b, const std::string& msg)
{
  if (!b)
    std::cerr << "Error: " << msg << std::endl;
  return b;
}

template <typename PolygonMesh>
PolygonMesh read_triangle_mesh(const std::string& filename)
{
  PolygonMesh pmesh;
  if (!CGAL::IO::read_polygon_mesh(filename, pmesh) || CGAL::is_empty(pmesh)) {
    std::cerr << "Error: failed to read " << filename << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (!CGAL::is_triangle_mesh(pmesh))
    PMP::triangulate_faces(pmesh);
  return pmesh;
}

bool test_API(const Mesh& sm)
{
  std::cout << "== API" << std::endl;

  // faces move inward if weights are negative, and outward if weights are positive
  CGAL::Constant_property_map<Mesh::Face_index, FT> inward_weights(-2);
  CGAL::Constant_property_map<Mesh::Face_index, FT> outward_weights(2);
  bool ok = true;

  ok &= check(CGAL::create_straight_skeleton_3(sm) != nullptr, "(tmesh)");
  ok &= check(CGAL::create_straight_skeleton_3(sm, CGAL::parameters::face_weight_map(inward_weights)) != nullptr, "(tmesh, np)");
  ok &= check(CGAL::create_straight_skeleton_3(sm, FT(0.5)) != nullptr, "(tmesh, FT)");
  ok &= check(CGAL::create_straight_skeleton_3(sm, FT(0.5), CGAL::parameters::face_weight_map(outward_weights)) != nullptr, "(tmesh, FT, np)");
  ok &= check(CGAL::create_straight_skeleton_3(sm, 0.5) != nullptr, "(tmesh, double)");
  ok &= check(CGAL::create_straight_skeleton_3(sm, 0.5, CGAL::parameters::face_weight_map(outward_weights)) != nullptr, "(tmesh, double, np)");

  std::vector<double> times = { 0.25, 0.5 };
  std::vector<Mesh> results;
  ok &= check(CGAL::create_straight_skeleton_and_offset_polyhedra_3(sm, times, results), "offsets with double times");
  ok &= check(results.size() == times.size(), "number of offsets with double times");

  return ok;
}

// Inputs that are not index-based (here, Polyhedron_3) must be supported, as input and as output
bool test_polyhedron_3(const std::string& filename)
{
  std::cout << "== test polyhedron_3" << std::endl;

  const Mesh sm = read_triangle_mesh<Mesh>(filename);
  const Polyhedron poly = read_triangle_mesh<Polyhedron>(filename);
  bool ok = true;

  auto ss_sm = CGAL::create_straight_skeleton_3(sm);
  auto ss_poly = CGAL::create_straight_skeleton_3(poly);
  ok &= check(ss_poly != nullptr, "skeleton of a Polyhedron_3");
  if (ss_sm && ss_poly)
    ok &= check(ss_sm->nodes().size() == ss_poly->nodes().size(), "same skeleton for Surface_mesh and Polyhedron_3");

  CGAL::Constant_property_map<boost::graph_traits<Polyhedron>::face_descriptor, FT> weights(-1);
  std::vector<FT> times = { 0.25 };

  std::vector<Mesh> results_sm;
  ok &= check(CGAL::create_straight_skeleton_and_offset_polyhedra_3(sm, times, results_sm), "offset of a Surface_mesh");

  std::vector<Polyhedron> results_poly;
  std::map<boost::graph_traits<Polyhedron>::face_descriptor,
           boost::graph_traits<Polyhedron>::face_descriptor> f2f;
  ok &= check(CGAL::create_straight_skeleton_and_offset_polyhedra_3(poly, times, results_poly,
                                                                    CGAL::parameters::face_weight_map(weights),
                                                                    CGAL::parameters::face_to_face_map(boost::make_assoc_property_map(f2f))),
              "offset of a Polyhedron_3");

  if (!results_sm.empty() && !results_poly.empty()) {
    ok &= check(PMP::volume(results_sm[0]) == PMP::volume(results_poly[0]), "same offset for Surface_mesh and Polyhedron_3");
    ok &= check(f2f.size() == num_faces(results_poly[0]), "face_to_face_map filled for every output face");
  }

  return ok;
}

bool is_close(const FT& a, const FT& b)
{
  return CGAL::abs(a - b) < FT(1e-6);
}

FT max_node_time(const std::shared_ptr<Straight_skeleton_3>& ss)
{
  FT max_time = 0;
  for (const NodeSPtr& node : ss->nodes())
    max_time = (std::max)(max_time, node->time());
  return max_time;
}

// all node times are in [0; max_time], and time increases from the source to the target of spokes
bool check_times(const std::shared_ptr<Straight_skeleton_3>& ss, const FT& max_time)
{
  for (const NodeSPtr& node : ss->nodes()) {
    if (node->time() < 0 || node->time() > max_time + FT(1e-6))
      return false;
  }
  for (const SpokeSPtr& spoke : ss->spokes()) {
    if (spoke->has_target() && spoke->source()->time() > spoke->target()->time())
      return false;
  }
  return true;
}

// Time is positive, and the sign of the weights gives the direction. 'sm' is the cube [-1;1]^3.
bool test_time_and_weights(const Mesh& sm)
{
  std::cout << "== time and weights" << std::endl;

  using face_descriptor = Mesh::Face_index;
  CGAL::Constant_property_map<face_descriptor, FT> inward_weights(-2);
  CGAL::Constant_property_map<face_descriptor, FT> outward_weights(2);
  CGAL::Constant_property_map<face_descriptor, FT> zero_weights(0);

  std::map<face_descriptor, FT> mixed_weights;
  for (face_descriptor f : faces(sm))
    mixed_weights[f] = (f.idx() % 2 == 0) ? 1 : -1;

  bool ok = true;

  // interior skeletons: the cube collapses at t = 1 / |weight|
  auto ss = CGAL::create_straight_skeleton_3(sm);
  ok &= check(ss && check_times(ss, 1) && is_close(max_node_time(ss), 1), "interior skeleton times in [0;1]");
  ss = CGAL::create_straight_skeleton_3(sm, CGAL::parameters::face_weight_map(inward_weights));
  ok &= check(ss && check_times(ss, 0.5) && is_close(max_node_time(ss), 0.5), "interior skeleton with weight -2 collapses at 0.5");
  ss = CGAL::create_straight_skeleton_3(sm, FT(0.5));
  ok &= check(ss && check_times(ss, 0.5), "partial interior skeleton times in [0;0.5]");

  // exterior skeleton
  ss = CGAL::create_straight_skeleton_3(sm, FT(1), CGAL::parameters::face_weight_map(outward_weights));
  ok &= check(ss && check_times(ss, 1), "exterior skeleton times in [0;1]");
  if (ss)
    ok &= check(CGAL::Straight_skeletons_3::IO::write_OBJ("exterior_skeleton.obj", ss, FT(1)), "write exterior skeleton");

  // offsets: the cube has side 2 - 2*0.25 inward (weight -1), and 2 + 2*2*0.25 outward (weight 2)
  std::vector<FT> times = { 0.25 };
  std::vector<Mesh> results;
  ok &= check(CGAL::create_straight_skeleton_and_offset_polyhedra_3(sm, times, results) &&
              results.size() == 1 && is_close(PMP::volume(results[0]), FT(3.375)) &&
              PMP::is_outward_oriented(results[0]), "inward offset");
  results.clear();
  ok &= check(CGAL::create_straight_skeleton_and_offset_polyhedra_3(sm, times, results,
                                                                    CGAL::parameters::face_weight_map(outward_weights)) &&
              results.size() == 1 && is_close(PMP::volume(results[0]), FT(27)) &&
              PMP::is_outward_oriented(results[0]), "outward offset");
  results.clear();

  // invalid times and weights
  ok &= check(CGAL::create_straight_skeleton_3(sm, FT(-1)) == nullptr, "negative time is rejected");
  ok &= check(CGAL::create_straight_skeleton_3(sm, FT(0)) == nullptr, "zero time is rejected");
  ok &= check(CGAL::create_straight_skeleton_3(sm, CGAL::parameters::face_weight_map(outward_weights)) == nullptr,
              "positive weights without a maximum time are rejected");
  ok &= check(CGAL::create_straight_skeleton_3(sm, CGAL::parameters::face_weight_map(zero_weights)) == nullptr,
              "zero weights are rejected");
  ok &= check(CGAL::create_straight_skeleton_3(sm, FT(1), CGAL::parameters::face_weight_map(
                                                             boost::make_assoc_property_map(mixed_weights))) == nullptr,
              "weights of different signs are rejected");
  ok &= check(!CGAL::create_straight_skeleton_and_offset_polyhedra_3(sm, std::vector<FT>{ }, results,
                                                                     CGAL::parameters::face_weight_map(outward_weights)),
              "outward offsets without times are rejected");
  ok &= check(!CGAL::create_straight_skeleton_and_offset_polyhedra_3(sm, std::vector<FT>{ -0.25 }, results),
              "negative offset time is rejected");

  return ok;
}

std::size_t number_of_connected_components(Mesh sm)
{
  auto fcc = sm.add_property_map<Mesh::Face_index, std::size_t>("f:cc", 0).first;
  return PMP::connected_components(sm, fcc);
}

// The input is made of two convex polyhedra that collide as they grow, which creates an event
// between elements of different connected components (e.g., a pierce or an edge split event).
// An outward offset is then the union of the outward offsets of the two polyhedra.
bool test_outward_collision(const std::string& filename)
{
  std::cout << "== outward collision: " << filename << std::endl;

  const Mesh sm = read_triangle_mesh<Mesh>(filename);

  std::vector<Mesh> ccs;
  PMP::split_connected_components(sm, ccs);
  if (!check(ccs.size() == 2, "two connected components"))
    return false;

  CGAL::Constant_property_map<Mesh::Face_index, FT> outward_weights(1);
  const auto np = CGAL::parameters::face_weight_map(outward_weights);

  // before, and after the collision
  const std::vector<FT> times = { 0.1, 0.3 };
  std::vector<Mesh> results, results_0, results_1;
  if (!check(CGAL::create_straight_skeleton_and_offset_polyhedra_3(sm, times, results, np) &&
             CGAL::create_straight_skeleton_and_offset_polyhedra_3(ccs[0], times, results_0, np) &&
             CGAL::create_straight_skeleton_and_offset_polyhedra_3(ccs[1], times, results_1, np) &&
             results.size() == 2 && results_0.size() == 2 && results_1.size() == 2, "outward offsets"))
    return false;

  bool ok = true;
  ok &= check(number_of_connected_components(results[0]) == 2, "two connected components before the collision");
  ok &= check(number_of_connected_components(results[1]) == 1, "one connected component after the collision");

  for (std::size_t i=0; i<times.size(); ++i) {
    Mesh union_mesh;
    ok &= check(PMP::corefine_and_compute_union(results_0[i], results_1[i], union_mesh), "union of the offsets");
    ok &= check(PMP::is_outward_oriented(results[i]), "outward-oriented offset");
    ok &= check(is_close(PMP::volume(results[i]), PMP::volume(union_mesh)), "offset is the union of the offsets");
  }

  return ok;
}

int main(int argc, char** argv)
{
  const std::string cube_filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/cube.off");
  const std::string pyramid_filename = (argc > 2) ? argv[2] : "data/pyramid.off";

  const Mesh cube = read_triangle_mesh<Mesh>(cube_filename);
  const Mesh pyramid = read_triangle_mesh<Mesh>(pyramid_filename);

  bool ok = true;
  ok &= test_API(cube);
  ok &= test_polyhedron_3(pyramid_filename);
  ok &= test_time_and_weights(cube);
  ok &= test_outward_collision("data/events/EdgeSplitEvent_outward.obj");
  ok &= test_outward_collision("data/events/PierceEvent_outward.obj");

  std::cout << (ok ? "Done" : "FAILED") << std::endl;
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}

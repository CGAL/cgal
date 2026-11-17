#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/create_offset_polyhedra.h>
#include <CGAL/create_straight_skeleton_3.h>
#include <CGAL/Straight_skeleton_3.h>

#include <CGAL/Polygon_mesh_processing/measure.h>
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

  CGAL::Constant_property_map<Mesh::Face_index, FT> weights(2);
  bool ok = true;

  ok &= check(CGAL::create_straight_skeleton_3(sm) != nullptr, "(tmesh)");
  ok &= check(CGAL::create_straight_skeleton_3(sm, CGAL::parameters::face_weight_map(weights)) != nullptr, "(tmesh, np)");
  ok &= check(CGAL::create_straight_skeleton_3(sm, FT(-0.5)) != nullptr, "(tmesh, FT)");
  ok &= check(CGAL::create_straight_skeleton_3(sm, FT(0.5), CGAL::parameters::face_weight_map(weights)) != nullptr, "(tmesh, FT, np)");
  ok &= check(CGAL::create_straight_skeleton_3(sm, -0.5) != nullptr, "(tmesh, double)");
  ok &= check(CGAL::create_straight_skeleton_3(sm, 0.5, CGAL::parameters::face_weight_map(weights)) != nullptr, "(tmesh, double, np)");

  std::vector<double> times = { -0.25, -0.5 };
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

  CGAL::Constant_property_map<boost::graph_traits<Polyhedron>::face_descriptor, FT> weights(1);
  std::vector<FT> times = { -0.25 };

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

int main(int argc, char** argv)
{
  const std::string cube_filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/cube.off");
  const std::string pyramid_filename = (argc > 2) ? argv[2] : "data/pyramid.off";

  const Mesh cube = read_triangle_mesh<Mesh>(cube_filename);
  const Mesh pyramid = read_triangle_mesh<Mesh>(pyramid_filename);

  bool ok = true;
  ok &= test_API(cube);
  ok &= test_polyhedron_3(pyramid_filename);

  std::cout << (ok ? "Done" : "FAILED") << std::endl;
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}

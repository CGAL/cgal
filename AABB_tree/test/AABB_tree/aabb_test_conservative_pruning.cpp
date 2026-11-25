/*! \file
 *
 * Tests the conservative pruning of the traversal of an AABB tree for linear
 * queries (lines, rays, and segments):
 * 1. the conservative test of a box never rejects a box that the query
 *    intersects (compared with an exact kernel), and
 * 2. the queries of the tree return exactly the primitives that the same test
 *    of every primitive (brute force) returns, also on degenerate queries.
 */

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_tree/internal/AABB_conservative_bbox_test.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/Random.h>
#include <CGAL/intersection_3.h>

using EK = CGAL::Simple_cartesian<CGAL::Exact_rational>;
using Epick = CGAL::Exact_predicates_inexact_constructions_kernel;
using Epeck = CGAL::Exact_predicates_exact_constructions_kernel;

/*! tests the conservative test against the exact test, on random boxes and
 * queries; with small integer coordinates, most queries touch faces, edges, or
 * corners of the box, or are parallel to coordinate axes or planes.
 * \return `true` if the conservative test never rejects a box that the query
 *         intersects; `false` otherwise.
 */
bool test_conservative_bbox_test() {
  using Test = CGAL::internal::AABB_tree::Conservative_linear_bbox_test;
  CGAL::Random rng(0);
  std::size_t num_failures = 0, num_intersecting = 0;
  const int n = 300000;
  for (int i = 0; i < n; ++i) {
    const bool integer = (i % 4 != 0);
    auto coord = [&](int lo, int hi) { return integer ? double(rng.get_int(lo, hi + 1)) : rng.get_double(lo, hi); };
    double c[6];
    for (int j = 0; j < 3; ++j) {
      c[j] = coord(-3, 3);
      c[j + 3] = c[j] + std::abs(coord(0, 3));
    }
    const CGAL::Bbox_3 bbox(c[0], c[1], c[2], c[3], c[4], c[5]);
    const double p[3] = { coord(-5, 5), coord(-5, 5), coord(-5, 5) };
    const double v[3] = { coord(-2, 2), coord(-2, 2), coord(-2, 2) };
    if ((std::max)({ std::abs(v[0]), std::abs(v[1]), std::abs(v[2]) }) < 1e-3) continue;

    const EK::Point_3 ep(p[0], p[1], p[2]);
    const EK::Vector_3 ev(v[0], v[1], v[2]);
    CGAL::Interval_nt<> ip[3], iv[3];
    for (int j = 0; j < 3; ++j) {
      ip[j] = p[j];
      iv[j] = v[j];
    }
    // The part of the line in the box has a parameter below 11 / max|v_i| < 1e6.
    const bool kinds[3][2] = { { false, false }, { true, false }, { true, true } };
    for (const auto& kind : kinds) {
      bool expected;
      if (! kind[0]) expected = CGAL::do_intersect(EK::Segment_3(ep - 1e6 * ev, ep + 1e6 * ev), bbox);
      else if (! kind[1]) expected = CGAL::do_intersect(EK::Segment_3(ep, ep + 1e6 * ev), bbox);
      else expected = CGAL::do_intersect(EK::Segment_3(ep, ep + ev), bbox);
      if (! expected) continue;
      ++num_intersecting;
      if (! Test(ip, iv, kind[0], kind[1])(bbox)) {
        if (++num_failures <= 10)
          std::cerr << "ERROR: the conservative test rejects the box " << bbox << " for the query (" << p[0] << " "
                    << p[1] << " " << p[2] << ") + t (" << v[0] << " " << v[1] << " " << v[2] << ")"
                    << (kind[1] ? ", t in [0, 1]" : kind[0] ? ", t >= 0" : "") << std::endl;
      }
    }
  }
  std::cout << "Conservative box test: " << num_failures << " failure(s) in " << num_intersecting
            << " intersecting cases" << std::endl;
  return num_failures == 0;
}

/*! tests the queries of the tree against brute force, on queries through mesh
 * vertices (lines, rays, and segments; also parallel to coordinate axes).
 * \tparam K the kernel.
 * \param filename the name of the input mesh file.
 * \param kernel_name the name of the kernel (for the report).
 * \return `true` if the tree and brute force agree on all queries; `false`
 *         otherwise.
 */
template <typename K>
bool test_tree(const std::string& filename, const std::string& kernel_name) {
  using Point = typename K::Point_3;
  using Mesh = CGAL::Surface_mesh<Point>;
  using Primitive = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
  using Traits = CGAL::AABB_traits_3<K, Primitive>;
  using Tree = CGAL::AABB_tree<Traits>;
  using Face = typename boost::graph_traits<Mesh>::face_descriptor;

  Mesh mesh;
  if (! CGAL::IO::read_polygon_mesh(filename, mesh)) {
    std::cerr << "ERROR: cannot read " << filename << std::endl;
    return false;
  }
  Tree tree(faces(mesh).first, faces(mesh).second, mesh);
  std::vector<typename K::Triangle_3> triangles;
  std::vector<Face> face_of;
  for (Face f : faces(mesh)) {
    // The same triangle as the primitive (the same vertex order: with inexact
    // constructions, the computed intersection depends on it).
    auto h = halfedge(f, mesh);
    triangles.emplace_back(mesh.point(target(h, mesh)), mesh.point(target(next(h, mesh), mesh)),
                           mesh.point(source(h, mesh)));
    face_of.push_back(f);
  }

  CGAL::Random rng(1);
  const int nv = static_cast<int>(num_vertices(mesh));
  std::size_t num_failures = 0, num_queries = 0;

  auto check = [&](const auto& query) {
    ++num_queries;
    // The primitives that the query intersects (do_intersect()), and those for
    // which the intersection functor of the kernel (the one that the traits
    // use) returns an object. The two may differ: intersection(Ray_3,
    // Triangle_3) ignores the direction of a ray in the coplanar case, and with
    // inexact constructions the functor may return no object.
    std::vector<Face> expected, expected_objects;
    for (std::size_t k = 0; k < triangles.size(); ++k) {
      if (CGAL::do_intersect(query, triangles[k])) expected.push_back(face_of[k]);
      if (K().intersect_3_object()(query, triangles[k])) expected_objects.push_back(face_of[k]);
    }
    std::vector<Face> found;
    tree.all_intersected_primitives(query, std::back_inserter(found));
    std::sort(expected.begin(), expected.end());
    std::sort(found.begin(), found.end());
    bool ok = (found == expected) &&
              (tree.number_of_intersected_primitives(query) == expected.size()) &&
              (tree.do_intersect(query) == ! expected.empty());
    std::vector<typename Tree::template Intersection_and_primitive_id<std::decay_t<decltype(query)>>::Type> xs;
    tree.all_intersections(query, std::back_inserter(xs));
    std::vector<Face> found_objects;
    for (const auto& x : xs) found_objects.push_back(x.second);
    std::sort(found_objects.begin(), found_objects.end());
    // Every primitive that the query intersects and for which intersection()
    // returns an object is reported, and only objects that intersection()
    // returns. (The spurious objects of intersection(), see above, may lie in
    // boxes that the query does not intersect; the traversal reports those
    // that lie in visited leaves. With inexact constructions, intersection()
    // may return no object for a primitive that the query intersects.)
    std::vector<Face> required;
    std::set_intersection(expected.begin(), expected.end(), expected_objects.begin(), expected_objects.end(),
                          std::back_inserter(required));
    ok = ok && std::includes(found_objects.begin(), found_objects.end(), required.begin(), required.end()) &&
      std::includes(expected_objects.begin(), expected_objects.end(), found_objects.begin(), found_objects.end());
    auto any = tree.any_intersected_primitive(query);
    ok = ok && (bool(any) == ! expected.empty()) &&
      (! any || std::binary_search(expected.begin(), expected.end(), *any));
    if (! ok && ++num_failures <= 10)
      std::cerr << "ERROR: " << kernel_name << ": the tree and brute force differ for a query on " << filename
                << ": brute force " << expected.size() << ", all_intersected_primitives " << found.size()
                << ", number_of_intersected_primitives " << tree.number_of_intersected_primitives(query)
                << ", do_intersect " << tree.do_intersect(query) << ", all_intersections " << xs.size() << " (brute force "
                << expected_objects.size() << ")"
                << ", any_intersected_primitive " << bool(any) << std::endl;
  };

  for (int i = 0; i < 300; ++i) {
    const Point& a = mesh.point(typename Mesh::Vertex_index(rng.get_int(0, nv)));
    const Point& b = mesh.point(typename Mesh::Vertex_index(rng.get_int(0, nv)));
    if (a == b) continue;
    check(typename K::Line_3(a, b));
    check(typename K::Ray_3(a, b));
    check(typename K::Segment_3(a, b));
    // Parallel to a coordinate axis (zero direction components), through a vertex.
    const typename K::Vector_3 axes[3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
    const auto& axis = axes[i % 3];
    check(typename K::Line_3(a, axis));
    check(typename K::Ray_3(a, axis));
    check(typename K::Segment_3(a, a + axis));
  }
  std::cout << kernel_name << ", " << filename << ": " << num_failures << " failure(s) in " << num_queries
            << " queries" << std::endl;
  return num_failures == 0;
}

int main(int argc, char* argv[]) {
  bool ok = test_conservative_bbox_test();
  const std::vector<std::string> filenames = (argc > 1) ?
    std::vector<std::string>(argv + 1, argv + argc) :
    std::vector<std::string>{ CGAL::data_file_path("meshes/cube_meshed.off"), CGAL::data_file_path("meshes/cross.off"),
                              CGAL::data_file_path("meshes/joint.off") };
  for (const auto& filename : filenames) {
    ok &= test_tree<Epick>(filename, "Epick");
    ok &= test_tree<Epeck>(filename, "Epeck");
  }
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}

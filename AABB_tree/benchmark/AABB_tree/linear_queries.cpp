/*! \file
 *
 * Benchmark of the queries of an AABB tree with linear queries (lines, rays,
 * and segments) on a triangle mesh, with Epick and Epeck. Two query sets are
 * used: through two random points in the bounding box (general position), and
 * through two random vertices of the mesh (degenerate).
 *
 * For every kernel, query set, and query type, the program prints the time of
 * number_of_intersected_primitives(), all_intersections(), do_intersect(), and
 * any_intersected_primitive(), and a checksum of their results, which allows
 * to compare the results of two builds (e.g., before and after a change).
 *
 * Usage: linear_queries [mesh] [number of queries]
 */

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Random.h>
#include <CGAL/Real_timer.h>

/*! times the queries of `tree` for the queries that `make` constructs from
 * the pairs of points in `points`, and prints one line of results.
 */
template <typename Tree, typename Points, typename Make>
void time_queries(const Tree& tree, const Points& points, Make make, const std::string& label) {
  using Query = decltype(make(points.front().first, points.front().second));
  using Intersection = typename Tree::template Intersection_and_primitive_id<Query>::Type;

  std::size_t checksum = 0;
  double times[4];
  CGAL::Real_timer timer;

  timer.start();
  for (const auto& [a, b] : points) checksum += tree.number_of_intersected_primitives(make(a, b));
  timer.stop();
  times[0] = timer.time();

  std::vector<Intersection> intersections;
  timer.reset();
  timer.start();
  for (const auto& [a, b] : points) {
    intersections.clear();
    tree.all_intersections(make(a, b), std::back_inserter(intersections));
    checksum += intersections.size();
  }
  timer.stop();
  times[1] = timer.time();

  timer.reset();
  timer.start();
  for (const auto& [a, b] : points) checksum += tree.do_intersect(make(a, b));
  timer.stop();
  times[2] = timer.time();

  timer.reset();
  timer.start();
  for (const auto& [a, b] : points) checksum += bool(tree.any_intersected_primitive(make(a, b)));
  timer.stop();
  times[3] = timer.time();

  std::cout << std::left << std::setw(28) << label << std::right << std::fixed << std::setprecision(3);
  for (double t : times) std::cout << std::setw(12) << t;
  std::cout << std::setw(14) << checksum << std::endl;
}

/*! runs the benchmark with the kernel `K` on the mesh in `filename` with `n`
 * queries of every kind.
 */
template <typename K>
bool run(const std::string& filename, std::size_t n, const std::string& kernel_name) {
  using Point = typename K::Point_3;
  using Mesh = CGAL::Surface_mesh<Point>;
  using Primitive = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
  using Tree = CGAL::AABB_tree<CGAL::AABB_traits_3<K, Primitive>>;

  Mesh mesh;
  if (! CGAL::IO::read_polygon_mesh(filename, mesh) || is_empty(mesh)) {
    std::cerr << "Error: cannot read the mesh " << filename << std::endl;
    return false;
  }
  Tree tree(faces(mesh).first, faces(mesh).second, mesh);
  tree.build();

  const auto bb = CGAL::Polygon_mesh_processing::bbox(mesh);
  const int nv = static_cast<int>(num_vertices(mesh));
  for (const std::string type : { "random", "vertices" }) {
    CGAL::Random rng(0);
    std::vector<std::pair<Point, Point>> points;
    points.reserve(n);
    while (points.size() < n) {
      Point a, b;
      if (type == "random") {
        a = Point(rng.get_double(bb.xmin(), bb.xmax()), rng.get_double(bb.ymin(), bb.ymax()),
                  rng.get_double(bb.zmin(), bb.zmax()));
        b = Point(rng.get_double(bb.xmin(), bb.xmax()), rng.get_double(bb.ymin(), bb.ymax()),
                  rng.get_double(bb.zmin(), bb.zmax()));
      }
      else {
        a = mesh.point(typename Mesh::Vertex_index(rng.get_int(0, nv)));
        b = mesh.point(typename Mesh::Vertex_index(rng.get_int(0, nv)));
      }
      if (a != b) points.emplace_back(a, b);
    }
    const std::string label = kernel_name + " " + type + " ";
    time_queries(tree, points, [](const Point& a, const Point& b) { return typename K::Line_3(a, b); },
                 label + "line");
    time_queries(tree, points, [](const Point& a, const Point& b) { return typename K::Ray_3(a, b); },
                 label + "ray");
    time_queries(tree, points, [](const Point& a, const Point& b) { return typename K::Segment_3(a, b); },
                 label + "segment");
  }
  return true;
}

int main(int argc, char* argv[]) {
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/armadillo.off");
  const std::size_t n = (argc > 2) ? std::stoul(argv[2]) : 50000;

  std::cout << "Mesh: " << filename << ", queries: " << n << " of each kind\n"
            << std::left << std::setw(28) << "kernel, query set, type" << std::right << std::setw(12) << "count"
            << std::setw(12) << "all" << std::setw(12) << "do_inters." << std::setw(12) << "any"
            << std::setw(14) << "checksum" << std::endl;
  bool ok = run<CGAL::Exact_predicates_inexact_constructions_kernel>(filename, n, "Epick");
  ok = ok && run<CGAL::Exact_predicates_exact_constructions_kernel>(filename, n, "Epeck");
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}

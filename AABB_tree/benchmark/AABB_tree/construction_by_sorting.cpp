
#include <iostream>

#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_traits_construct_by_sorting.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/Timer.h>
#include <CGAL/point_generators_3.h>

#include <boost/core/demangle.hpp>

static std::size_t C = 100;
static std::size_t T = 10000;

template<typename K>
std::vector<CGAL::Segment_3<K>> generate_queries(std::size_t n) {
  typedef CGAL::Point_3<K> Point_3;
  typedef CGAL::Segment_3<K> Segment_3;

  // Generate some points
  CGAL::Random r(23);
  CGAL::Random_points_in_cube_3<Point_3, CGAL::Creator_uniform_3<typename K::FT, Point_3> > g(2.0, r);
  std::vector<Point_3> points;
  points.reserve(n * 2);
  std::copy_n(g, n * 2, std::back_inserter(points));

  // Combine those points into Segments
  std::vector<Segment_3> segments;
  segments.reserve(n);
  for (std::size_t i = 0, j = points.size() - 1; i < j; ++i, --j)
    segments.push_back(Segment_3(points[i], points[j]));

  return segments;
}

template<typename K>
double benchmark_recursive_partitioning_traversal(std::string input_path) {
  typedef CGAL::Polyhedron_3<K> Polyhedron;
  typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
  typedef CGAL::AABB_traits<K, Primitive> Traits;
  typedef CGAL::AABB_tree<Traits> Tree;

  std::ifstream in(input_path);
  Polyhedron polyhedron;
  in >> polyhedron;

  Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
  tree.build();

  auto queries = generate_queries<K>(T);
  typedef typename Tree::AABB_traits::template Intersection_and_primitive_id<decltype(queries[0])>::Type Result_type;
  std::vector<Result_type> v;
  v.reserve(queries.size());

  CGAL::Timer t;
  t.start();
  for (const auto &query : queries)
    tree.template all_intersections(query, std::back_inserter(v));
  t.stop();

  return t.time() / queries.size();
}

template<typename K>
double benchmark_sorting_traversal(std::string input_path) {
  typedef CGAL::Polyhedron_3<K> Polyhedron;
  typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
  typedef CGAL::AABB_traits_construct_by_sorting<K, Primitive> Traits;
  typedef CGAL::AABB_tree<Traits> Tree;

  std::ifstream in(input_path);
  Polyhedron polyhedron;
  in >> polyhedron;

  Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
  tree.build();

  auto queries = generate_queries<K>(T);
  typedef typename Tree::AABB_traits::template Intersection_and_primitive_id<decltype(queries[0])>::Type Result_type;
  std::vector<Result_type> v;
  v.reserve(queries.size());

  CGAL::Timer t;
  t.start();
  for (const auto &query : queries)
    tree.template all_intersections(query, std::back_inserter(v));
  t.stop();

  return t.time() / queries.size();
}

template<typename Traits>
double benchmark_construction(std::string input_path) {
  typedef typename Traits::Geom_traits K;
  typedef CGAL::AABB_tree<Traits> Tree;
  typedef CGAL::Polyhedron_3<K> Polyhedron;
  typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;

  std::ifstream in(input_path);
  Polyhedron polyhedron;
  in >> polyhedron;

  CGAL::Timer t;
  t.start();
  for (int i = 0; i < C; ++i) {
    Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
    tree.build();
  }
  t.stop();

  return t.time() / C;
}

template<typename K>
void benchmark(std::string input_path) {
  typedef CGAL::Polyhedron_3<K> Polyhedron;
  typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;

  typedef CGAL::AABB_traits<K, Primitive> Traits_construct_by_splitting;
  typedef CGAL::AABB_traits_construct_by_sorting<K, Primitive> Traits_construct_by_sorting;

  std::cout << "{| class=\"wikitable\"\n";
  std::cout << "! " << boost::core::demangle(typeid(K).name()) << " !! Recursive Partition !! Hilbert Sort"
            << "\n|-\n";
  std::cout << "| construction || " << std::flush
            << benchmark_construction<Traits_construct_by_splitting>(input_path) << " s || " << std::flush
            << benchmark_construction<Traits_construct_by_sorting>(input_path) << " s"
            << "\n|-\n";
  std::cout << "| traversal || " << std::flush
            << benchmark_recursive_partitioning_traversal<K>(input_path) << " s || " << std::flush
            << benchmark_sorting_traversal<K>(input_path) << " s"
            << "\n|-\n";
  std::cout << "|}\n\n";
}

int main(int argc, char **argv) {

  // Determine our data source, with a default if no path is provided
  std::string input_path = argc > 1 ? argv[1] : "data/handle.off";

  benchmark<CGAL::Simple_cartesian<float>>(input_path);
  benchmark<CGAL::Cartesian<float>>(input_path);
  benchmark<CGAL::Simple_cartesian<double>>(input_path);
  benchmark<CGAL::Cartesian<double>>(input_path);
  benchmark<CGAL::Exact_predicates_inexact_constructions_kernel>(input_path);

}

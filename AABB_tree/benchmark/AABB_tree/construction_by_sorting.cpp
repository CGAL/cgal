
#include <iostream>

#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh/IO/PLY.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_traits_construct_by_sorting.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/Real_timer.h>
#include <CGAL/point_generators_3.h>

#include <boost/core/demangle.hpp>

static std::size_t C = 10;
static std::size_t T = 100000;

template<typename K>
std::vector<CGAL::Ray_3<K>> generate_queries(std::size_t n) {
  typedef CGAL::Point_3<K> Point_3;
  typedef CGAL::Ray_3<K> Ray_3;

  // Generate some points
  CGAL::Random r(23);
  CGAL::Random_points_in_cube_3<Point_3, CGAL::Creator_uniform_3<typename K::FT, Point_3> > g(2.0, r);
  std::vector<Point_3> points;
  points.reserve(n * 2);
  std::copy_n(g, n * 2, std::back_inserter(points));

  // Combine those points into Rays
  std::vector<Ray_3> queries;
  queries.reserve(n);
  for (std::size_t i = 0, j = points.size() - 1; i < j; ++i, --j)
    queries.push_back(Ray_3(points[i], points[j]));

  return queries;
}

template<typename Traits, typename Data, typename Queries>
double benchmark_traversal(Data polyhedron, Queries queries) {
  typedef CGAL::AABB_tree<Traits> Tree;

  Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
  tree.build();

  typedef typename Tree::AABB_traits::template Intersection_and_primitive_id<decltype(queries[0])>::Type Result_type;
  std::vector<Result_type> v;
  v.reserve(queries.size());

  CGAL::Real_timer t;
  t.start();
  for (const auto &query : queries)
    tree.template all_intersections(query, std::back_inserter(v));
  t.stop();

  return t.time() / queries.size();
}

template<typename Traits, typename Data>
double benchmark_construction(Data polyhedron) {
  typedef CGAL::AABB_tree<Traits> Tree;

  CGAL::Real_timer t;
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
  typedef CGAL::Surface_mesh<CGAL::Point_3<K>> Polyhedron;
  typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;

  typedef CGAL::AABB_traits<K, Primitive, CGAL::Default> Traits_construct_by_splitting;
  typedef CGAL::AABB_traits_construct_by_sorting<K, Primitive, CGAL::Default, CGAL::Parallel_tag> Traits_construct_by_sorting;

  std::ifstream in(input_path);
  Polyhedron polyhedron;
//  in >> polyhedron;
  CGAL::IO::read_PLY(in, polyhedron);

  auto queries = generate_queries<K>(T);

  std::cout << "{| class=\"wikitable\"\n";
  std::cout << "! " << boost::core::demangle(typeid(K).name()) << " !! Recursive Partition !! Hilbert Sort"
            << "\n|-" << std::endl;
  std::cout << "| construction || " << std::flush
            << benchmark_construction<Traits_construct_by_splitting>(polyhedron) << " s || " << std::flush
            << benchmark_construction<Traits_construct_by_sorting>(polyhedron) << " s"
            << "\n|-" << std::endl;
  std::cout << "| traversal || " << std::flush
            << benchmark_traversal<Traits_construct_by_splitting>(polyhedron, queries) << " s || " << std::flush
            << benchmark_traversal<Traits_construct_by_sorting>(polyhedron, queries) << " s"
            << "\n|-" << std::endl;
  std::cout << "|}\n" << std::endl;
}

int main(int argc, char **argv) {

  // Determine our data source, with a default if no path is provided
  std::string input_path = argc > 1 ? argv[1] : "data/handle.off";

//  benchmark<CGAL::Simple_cartesian<float>>(input_path);
//  benchmark<CGAL::Simple_cartesian<double>>(input_path);
  benchmark<CGAL::Exact_predicates_inexact_constructions_kernel>(input_path);

}

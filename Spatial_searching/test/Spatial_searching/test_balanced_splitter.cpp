#include <cmath>
#include <vector>
#include <fstream>

#include <CGAL/Real_timer.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Euclidean_distance.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Kd_tree.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/point_generators_3.h>

using SCF   = CGAL::Simple_cartesian<float>;
using SCD   = CGAL::Simple_cartesian<double>;
using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;

using Kernel = SCD;

using FT         = typename Kernel::FT;
using Point_3    = typename Kernel::Point_3;
using Triangle_3 = typename Kernel::Triangle_3;
using STraits    = CGAL::Search_traits_3<Kernel>;

using Splitter = CGAL::Sliding_midpoint<STraits>;
using KD_tree  = CGAL::Kd_tree<STraits, Splitter>;
using Timer    = CGAL::Real_timer;

using Iterator  = typename std::vector<Triangle_3>::const_iterator;
using Primitive = CGAL::AABB_triangle_primitive<Kernel, Iterator>;
using ABtraits  = CGAL::AABB_traits<Kernel, Primitive>;
using AB_tree   = CGAL::AABB_tree<ABtraits>;

void test_balanced_tree(const std::vector<Point_3>& /* points */) {

  std::vector<Point_3> points = {
    Point_3(2,3,3), Point_3(5,4,2), Point_3(9,6,7), Point_3(4,7,9), Point_3(8,1,5),
    Point_3(7,2,6), Point_3(9,4,1), Point_3(8,4,2), Point_3(9,7,8), Point_3(6,3,1),
    Point_3(3,4,5), Point_3(1,6,8), Point_3(9,5,3), Point_3(2,1,3), Point_3(8,7,6),
    Point_3(5,4,2), Point_3(6,3,1), Point_3(8,7,6), Point_3(9,6,7), Point_3(2,1,3),
    Point_3(7,2,6), Point_3(4,7,9), Point_3(1,6,8), Point_3(3,4,5), Point_3(9,4,1)
  };
  assert(points.size() == 25);
  std::cout << "* num points: " << points.size() << std::endl;

  // using Splitter = CGAL::Sliding_midpoint<STraits>;

  using Splitter = CGAL::Balanced_splitter<STraits>;
  using Kd_tree  = CGAL::Kd_tree<STraits, Splitter>;

  const unsigned int bucket_size = 5; // TODO: better test bucket size!
  std::cout << "* bucket size: " << bucket_size << std::endl;
  Splitter splitter(bucket_size);

  std::cout << "* building the balanced tree ... " << std::endl;
  Kd_tree tree(points.begin(), points.end(), splitter);
  // tree.preprocess(0.00001); // smaller -> more duplicates
  tree.build();
  std::cout << "* building done" << std::endl << std::endl;
  tree.statistics(std::cout);
  tree.print(std::cout);

  // Search.
  // Test for default fuzzy search.
  std::cout << "* testing search: " << std::endl;
  using Fuzzy_sphere = CGAL::Fuzzy_sphere<STraits>;
  Fuzzy_sphere fuzzy_sphere(Point_3(2,3,3), 2);

  std::vector<Point_3> result;
  tree.search(std::back_inserter(result), fuzzy_sphere);
  std::cout << "- found points (fuzzy): " << result.size() << std::endl;
  for (const auto& point : result) {
    std::cout << point << std::endl;
  }

  // Orthogonal search.
  // Test if it works for orthogonal search when extended node is required.
  using Distance = CGAL::Euclidean_distance<STraits>;
  using Neighbor_search = CGAL::Orthogonal_k_neighbor_search<STraits, Distance, Splitter, Kd_tree>;
  const Point_3 query(2,3,3);
  Neighbor_search nsearch(tree, query, 2);
  std::cout << "- found points (orth k): " << std::distance(nsearch.begin(), nsearch.end()) << std::endl;
  for (auto it = nsearch.begin(); it != nsearch.end(); ++it) {
    std::cout << it->first << " with dist: " << CGAL::sqrt(it->second) << std::endl;
  }
}

void call_kd_tree(
  const std::vector<Point_3>& points) {

  const unsigned int bucket_size = 10;
  std::cout << "* bucket size: " << bucket_size << std::endl;
  Splitter splitter(bucket_size);

  std::cout << "* building the tree with preprocessing ... ";
  KD_tree tree(points.begin(), points.end(), splitter);
  // It is off for the moment to compare against kd-tree-single-thread!
  // tree.preprocess();
  tree.build();
  std::cout << "finished" << std::endl;
  tree.statistics(std::cout);
}

void test_runtime_kd_tree_build(
  const std::vector<Point_3>& points,
  const std::size_t num_iters) {

  Timer timer;
  std::cout << "* kd tree build tests, num iterations: " << num_iters << std::endl;

  double mean_time = 0.0;
  for (std::size_t k = 0; k < num_iters; ++k) {
    timer.reset();
    timer.start();
    KD_tree tree(points.begin(), points.end());
    tree.build();
    timer.stop();
    mean_time += timer.time();
  }
  mean_time /= static_cast<double>(num_iters);
  std::cout << "mean time without preprocessing: "
    << mean_time << " sec." << std::endl;

  mean_time = 0.0;
  for (std::size_t k = 0; k < num_iters; ++k) {
    timer.reset();
    timer.start();
    KD_tree tree(points.begin(), points.end());
    tree.preprocess();
    tree.build();
    timer.stop();
    mean_time += timer.time();
  }
  mean_time /= static_cast<double>(num_iters);
  std::cout << "mean time with preprocessing: "
    << mean_time << " sec." << std::endl;
}

void test_runtime_aabb_tree_build(
  const std::vector<Point_3>& points,
  const std::vector<Triangle_3>& triangles,
  const std::size_t num_iters) {

  Timer timer;
  std::cout << "* aabb tree build tests, num iterations: " << num_iters << std::endl;

  double mean_time = 0.0;
  for (std::size_t k = 0; k < num_iters; ++k) {
    timer.reset();
    timer.start();
    AB_tree tree(triangles.begin(), triangles.end());
    tree.build();
    tree.accelerate_distance_queries(points.begin(), points.end(), false /* preprocess */);
    timer.stop();
    mean_time += timer.time();
  }
  mean_time /= static_cast<double>(num_iters);
  std::cout << "mean time without preprocessing: "
    << mean_time << " sec." << std::endl;

  mean_time = 0.0;
  for (std::size_t k = 0; k < num_iters; ++k) {
    timer.reset();
    timer.start();
    AB_tree tree(triangles.begin(), triangles.end());
    tree.build();
    tree.accelerate_distance_queries(points.begin(), points.end(), true /* preprocess */);
    timer.stop();
    mean_time += timer.time();
  }
  mean_time /= static_cast<double>(num_iters);
  std::cout << "mean time with preprocessing: "
    << mean_time << " sec." << std::endl;
}

void test_runtime_aabb_tree_query(
  const std::vector<Point_3>& points,
  const std::vector<Triangle_3>& triangles,
  const std::size_t num_iters) {

  Timer timer;
  std::cout << "* aabb tree query tests, num iterations: " << num_iters << std::endl;

  int nb_points = 100;
  double size = 100.0;
  std::vector<Point_3> queries;
  queries.reserve(nb_points);
  CGAL::Random_points_in_cube_3<Point_3> generator(size);
  for (std::size_t i = 0; i < nb_points; ++i) {
    queries.push_back(*generator++);
  }

  AB_tree tree1(triangles.begin(), triangles.end());
  tree1.build();
  tree1.accelerate_distance_queries(points.begin(), points.end(), false /* preprocess */);

  double mean_time = 0.0;
  for (std::size_t k = 0; k < num_iters; ++k) {
    timer.reset();
    timer.start();
    for (const auto& query : queries) {
      tree1.closest_point(query);
    }
    timer.stop();
    mean_time += timer.time();
  }
  mean_time /= static_cast<double>(num_iters);
  std::cout << "mean time without preprocessing: "
    << mean_time << " sec." << std::endl;

  AB_tree tree2(triangles.begin(), triangles.end());
  tree2.build();
  tree2.accelerate_distance_queries(points.begin(), points.end(), true /* preprocess */);

  mean_time = 0.0;
  for (std::size_t k = 0; k < num_iters; ++k) {
    timer.reset();
    timer.start();
    for (const auto& query : queries) {
      tree2.closest_point(query);
    }
    timer.stop();
    mean_time += timer.time();
  }
  mean_time /= static_cast<double>(num_iters);
  std::cout << "mean time with preprocessing: "
    << mean_time << " sec." << std::endl;
}

int main(int argc, char* argv[]) {

  std::cout.precision(20);
  std::vector<Point_3> points;
  const std::string filename = (argc > 1 ? argv[1] : "data/failure-tiny.xyz"); // "data/clean-result.xyz"
  std::cout << "* input points: " << filename << std::endl;
  std::ifstream in(filename);
  CGAL::IO::set_ascii_mode(in);

  Point_3 p, q; std::size_t count = 0;
  while (in >> p >> q) {
    // Use it to compare against kd-tree-single-thread.
    // if (count == 320000) break;
    points.push_back(p);
    ++count;
  }

  test_balanced_tree(points);
  return EXIT_SUCCESS;

  std::cout << "* num points: " << points.size() << std::endl;

  // Default call.
  if (argc == 1) {
    call_kd_tree(points);
    return EXIT_SUCCESS;
  }

  // Runtime tests.
  call_kd_tree(points);

  const std::size_t num_iters = 10;
  test_runtime_kd_tree_build(points, num_iters);
  if (argc > 2) {
    const std::string meshfile = (argc > 2 ? argv[2] : "data/clean-mesh.off");

    std::vector<Point_3> vertices;
    std::vector< std::vector<std::size_t> > faces;
    std::cout << "* input mesh: " << meshfile << std::endl;
    CGAL::IO::read_polygon_soup(meshfile, vertices, faces);
    std::cout << "* num faces: " << faces.size() << std::endl;

    std::vector<Triangle_3> triangles;
    triangles.reserve(faces.size());
    for (const auto& face : faces) {
      assert(face.size() == 3);
      Triangle_3 triangle(
        vertices[face[0]], vertices[face[1]], vertices[face[2]]);
      triangles.push_back(triangle);
    }
    assert(triangles.size() == faces.size());

    test_runtime_aabb_tree_build(points, triangles, num_iters);
    test_runtime_aabb_tree_query(points, triangles, num_iters);
  }
  return EXIT_SUCCESS;
}

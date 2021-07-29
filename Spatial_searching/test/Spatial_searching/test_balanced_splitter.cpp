// #define KD_TREE_DEBUG // use this to output more debug info

#include <cmath>
#include <vector>
#include <fstream>

#include <CGAL/Real_timer.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Euclidean_distance.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Kd_tree.h>

using SCF   = CGAL::Simple_cartesian<float>;
using SCD   = CGAL::Simple_cartesian<double>;
using CDD   = CGAL::Cartesian_d<double>;
using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;
using Timer = CGAL::Real_timer;

template<typename Kernel>
void run_tests_2d() {

  std::cout << "* testing 2D ..." << std::endl;
  using Point_2  = typename Kernel::Point_2;
  using Traits   = CGAL::Search_traits_2<Kernel>;
  using Splitter = CGAL::Balanced_splitter<Traits>;
  using Kd_tree  = CGAL::Kd_tree<Traits, Splitter>;

  // No duplicates but it is usually the worst case for median splitters.
  std::vector<Point_2> points = {
    Point_2(0,1) , Point_2(0,2) , Point_2(0,3) , Point_2(0,4) , Point_2(0,5) ,
    Point_2(0,6) , Point_2(0,7) , Point_2(0,8) , Point_2(0,9) , Point_2(0,10),
    Point_2(0,11), Point_2(0,12), Point_2(0,13), Point_2(0,14), Point_2(0,15),
    Point_2(0,16), Point_2(0,17), Point_2(0,18), Point_2(0,19), Point_2(19,0)
  };
  assert(points.size() == 20);

  // Bucket size 5.
  Splitter splitter_bs_5(5);
  Kd_tree tree_0(points.begin(), points.end(), splitter_bs_5);
  tree_0.build();
  assert(tree_0.root()->depth() == 4);
  assert(tree_0.root()->num_items() == points.size());
  assert(tree_0.root()->num_nodes() == 5);

  // Default bucket size = 10.
  Kd_tree tree_1(points.begin(), points.end());
  tree_1.build();
  assert(tree_1.root()->depth() == 3);
  assert(tree_1.root()->num_items() == points.size());
  assert(tree_1.root()->num_nodes() == 3);

  // One line of collinear points with duplicates.
  points.clear();
  points = {
    Point_2(0, 0), Point_2(1, 0), Point_2(2, 0), Point_2(3, 0), Point_2(4, 0),
    Point_2(5, 0), Point_2(6, 0), Point_2(7, 0), Point_2(8, 0), Point_2(9, 0),
    Point_2(0, 0), Point_2(1, 0), Point_2(2, 0), Point_2(3, 0), Point_2(4, 0),
    Point_2(5, 0), Point_2(6, 0), Point_2(7, 0), Point_2(8, 0), Point_2(9, 0)
  };
  assert(points.size() == 20);

  Kd_tree tree_2(points.begin(), points.end(), splitter_bs_5);
  tree_2.build();
  assert(tree_2.root()->depth() == 3);
  assert(tree_2.root()->num_items() == points.size());
  assert(tree_2.root()->num_nodes() == 4);

  // Two lines of collinear points without duplicates.
  points.clear();
  points = {
    Point_2(0, 0), Point_2(1, 0), Point_2(2, 0), Point_2(3, 0), Point_2(4, 0),
    Point_2(5, 0), Point_2(6, 0), Point_2(7, 0), Point_2(8, 0), Point_2(9, 0),
    Point_2(0, 1), Point_2(1, 1), Point_2(2, 1), Point_2(3, 1), Point_2(4, 1),
    Point_2(5, 1), Point_2(6, 1), Point_2(7, 1), Point_2(8, 1), Point_2(9, 1)
  };
  assert(points.size() == 20);

  Kd_tree tree_3(points.begin(), points.end(), splitter_bs_5);
  tree_3.build();
  assert(tree_3.root()->depth() == 4);
  assert(tree_3.root()->num_items() == points.size());
  assert(tree_3.root()->num_nodes() == 5);
}

template<typename Kernel>
void run_tests_3d() {

  std::cout << "* testing 3D ..." << std::endl;
  using Point_3  = typename Kernel::Point_3;
  using Traits   = CGAL::Search_traits_3<Kernel>;
  using Splitter = CGAL::Balanced_splitter<Traits>;
  using Kd_tree  = CGAL::Kd_tree<Traits, Splitter>;

  // This data set has duplicates.
  std::vector<Point_3> points = {
    Point_3(2,3,3), Point_3(5,4,2), Point_3(9,6,7), Point_3(4,7,9), Point_3(8,1,5),
    Point_3(7,2,6), Point_3(9,4,1), Point_3(8,4,2), Point_3(9,7,8), Point_3(6,3,1),
    Point_3(3,4,5), Point_3(1,6,8), Point_3(9,5,3), Point_3(2,1,3), Point_3(8,7,6),
    Point_3(5,4,2), Point_3(6,3,1), Point_3(8,7,6), Point_3(9,6,7), Point_3(2,1,3),
    Point_3(7,2,6), Point_3(4,7,9), Point_3(1,6,8), Point_3(3,4,5), Point_3(9,4,1)
  };
  assert(points.size() == 25);

  // Decrease the bucket size but not too much.
  Splitter splitter_bs_5(5);
  Kd_tree tree_0(points.begin(), points.end(), splitter_bs_5);
  tree_0.build();
  assert(tree_0.root()->depth() == 3);
  assert(tree_0.root()->num_items() == points.size());
  assert(tree_0.root()->num_nodes() == 4);

  // Default bucket size = 10.
  Kd_tree tree_1(points.begin(), points.end());
  tree_1.build();
  assert(tree_1.root()->depth() == 3);
  assert(tree_1.root()->num_items() == points.size());
  assert(tree_1.root()->num_nodes() == 4);

  // Increase the bucket size but not too much.
  Splitter splitter_bs_15(15);
  Kd_tree tree_2(points.begin(), points.end(), splitter_bs_15);
  tree_2.build();
  assert(tree_2.root()->depth() == 2);
  assert(tree_2.root()->num_items() == points.size());
  assert(tree_2.root()->num_nodes() == 2);

  // The bucket size covers all points.
  Splitter splitter_bs_30(30);
  Kd_tree tree_3(points.begin(), points.end(), splitter_bs_30);
  tree_3.build();
  assert(tree_3.root()->depth() == 1);
  assert(tree_3.root()->num_items() == points.size());
  assert(tree_3.root()->num_nodes() == 1);
}

template<typename Kernel>
void run_tests_kd() {

  std::cout << "* testing KD ..." << std::endl;
  using Traits   = CDD;
  using Point_d  = CGAL::Point_d<Traits>;
  using Splitter = CGAL::Balanced_splitter<Traits>;
  using Kd_tree  = CGAL::Kd_tree<Traits, Splitter>;

  using Point_generator = CGAL::Random_points_in_cube_d<Point_d>;

  const int dim = 5;
  std::vector<Point_d> points;
  Point_generator generator(dim);
  std::copy_n(generator, 50, std::back_inserter(points));
  assert(points.size() == 50);

  Kd_tree tree(points.begin(), points.end());
  assert(tree.root()->depth() < 10);
  assert(tree.root()->num_items() == points.size());
  assert(tree.root()->num_nodes() > 0);
}

template<typename Kernel>
void run_all_tests() {

  run_tests_2d<Kernel>();
  run_tests_3d<Kernel>();
  run_tests_kd<Kernel>();
}

template<typename Kernel>
void test_balanced_tree(
  const std::string filename,
  const int ref_depth,
  const std::size_t ref_num_leaves,
  const bool verbose = false) {

  std::cout << std::endl;
  std::cout << "- testing " << filename << " ... " << std::endl;

  using Traits   = CGAL::Search_traits_3<Kernel>;
  using Splitter = CGAL::Balanced_splitter<Traits>;
  using Kd_tree  = CGAL::Kd_tree<Traits, Splitter>;

  using Distance        = CGAL::Euclidean_distance<Traits>;
  using Fuzzy_sphere    = CGAL::Fuzzy_sphere<Traits>;
  using Neighbor_search = CGAL::Orthogonal_k_neighbor_search<Traits, Distance, Splitter, Kd_tree>;

  using Point_3 = typename Kernel::Point_3;
  std::vector<Point_3> points;
  std::ifstream in(filename);
  CGAL::IO::set_ascii_mode(in);

  Point_3 p; std::size_t count = 0;
  while (in >> p) {
    points.push_back(p);
    ++count;
  }
  assert(points.size() > 0);

  std::cout << "* num points: " << points.size() << std::endl;
  std::cout << "* building the balanced tree ... " << std::endl;

  // Building a tree.
  Splitter splitter(5); // bucket size = 5
  Kd_tree tree(points.begin(), points.end(), splitter);

  Timer timer;
  timer.start();
  tree.build();
  timer.stop();

  std::cout << "* building done in " << timer.time() << " sec." << std::endl;
  if (verbose) {
    std::cout << std::endl;
    tree.statistics(std::cout);
    tree.print(std::cout);
  }

  if (ref_depth > 0) {
    assert(tree.root()->depth() == ref_depth);
  } else {
    assert(tree.root()->depth() < 10);
  }

  assert(tree.root()->num_items() == points.size());

  if (ref_num_leaves != std::size_t(-1)) {
    assert(tree.root()->num_nodes() == ref_num_leaves);
  } else {
    assert(tree.root()->num_nodes() > 0);
  }

  // Search.

  // Run fuzzy sphere search.
  std::cout << "* testing fuzzy search ... " << std::endl;

  const auto radius = 2;
  const Point_3 center(2,3,3);
  Fuzzy_sphere fuzzy_sphere(center, radius);

  std::vector<Point_3> result;
  tree.search(std::back_inserter(result), fuzzy_sphere);
  if (verbose) {
    std::cout << "- center: " << center << std::endl;
    std::cout << "- radius: " << radius << std::endl;
    std::cout << "- found points (fuzzy): " << result.size() << std::endl;
    for (const auto& point : result) {
      std::cout << point << std::endl;
    }
    std::cout << std::endl;
  }

  if (filename == "data/balanced.xyz") {
    assert(result.size() == 3);
    assert(result[0] == Point_3(2,3,3));
    assert(result[1] == Point_3(2,1,3));
    assert(result[2] == Point_3(2,1,3));
  }

  // Run orthogonal search when an extended node is required.
  std::cout << "* testing orthogonal search ... " << std::endl;
  result.clear();
  const Point_3 query(2,3,3);
  Neighbor_search nsearch(tree, query, 3);

  if (verbose) {
    std::cout << "- query point: " << query << std::endl;
    std::cout << "- found points (orth k): " << std::distance(nsearch.begin(), nsearch.end()) << std::endl;
  }
  for (auto it = nsearch.begin(); it != nsearch.end(); ++it) {
    if (verbose) std::cout << it->first << " with dist: " << CGAL::sqrt(it->second) << std::endl;
    result.push_back(it->first);
  }
  if (verbose) std::cout << std::endl;

  if (filename == "data/balanced.xyz") {
    assert(result.size() == 3);
    assert(result[0] == Point_3(2,3,3));
    assert(result[1] == Point_3(2,1,3));
    assert(result[2] == Point_3(2,1,3));
  }
}

int main(int argc, char* argv[]) {

  // Basic test on real input data.
  std::cout.precision(20);

  if (argc > 1) {
    test_balanced_tree<EPICK>(std::string(argv[1]), -1, -1, true);
    return EXIT_SUCCESS;
  } else {
    // This data set has duplicates.
    test_balanced_tree<EPICK>("data/balanced.xyz", 3, 4, false);
    // This data set has duplicates.
    test_balanced_tree<EPICK>("data/failure-tiny.xyz", 8, 82, false);
  }

  // Run tests.
  std::cout << std::endl;
  std::cout << "- testing SCFKR kernel ... " << std::endl;
  run_all_tests<SCF>();
  std::cout << "- done with success! " << std::endl;

  std::cout << "- testing SCDKR kernel ... " << std::endl;
  run_all_tests<SCD>();
  std::cout << "- done with success! " << std::endl;

  std::cout << "- testing EPICK kernel ... " << std::endl;
  run_all_tests<EPICK>();
  std::cout << "- done with success! " << std::endl;

  // TODO: Should we fix EPECK for running this test?
  // std::cout << "- testing EPECK kernel ... ";
  // run_all_tests<EPECK>(); // triggers a warning about tmp ref return and seg fault
  // std::cout << "done with success! " << std::endl;
  std::cout << std::endl;

  return EXIT_SUCCESS;
}

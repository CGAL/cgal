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

using SCF   = CGAL::Simple_cartesian<float>;
using SCD   = CGAL::Simple_cartesian<double>;
using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;

using Kernel = SCD;

using Point_3 = typename Kernel::Point_3;
using STraits = CGAL::Search_traits_3<Kernel>;
using Timer   = CGAL::Real_timer;

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

  // Other splitters.
  // using Splitter = CGAL::Sliding_midpoint<STraits>;

  using Splitter = CGAL::Balanced_splitter<STraits>;
  using Kd_tree  = CGAL::Kd_tree<STraits, Splitter>;

  const unsigned int bucket_size = 5;
  std::cout << "* bucket size: " << bucket_size << std::endl;
  Splitter splitter(bucket_size);

  std::cout << "* building the balanced tree ... " << std::endl;
  Kd_tree tree(points.begin(), points.end(), splitter);
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

int main(int argc, char* argv[]) {

  std::cout.precision(20);
  std::vector<Point_3> points;
  const std::string filename = (argc > 1 ? argv[1] : "data/failure-tiny.xyz"); // "data/clean-result.xyz"
  std::cout << "* input points: " << filename << std::endl;
  std::ifstream in(filename);
  CGAL::IO::set_ascii_mode(in);

  Point_3 p, q; std::size_t count = 0;
  while (in >> p >> q) {
    points.push_back(p);
    ++count;
  }

  test_balanced_tree(points);
  return EXIT_SUCCESS;
}

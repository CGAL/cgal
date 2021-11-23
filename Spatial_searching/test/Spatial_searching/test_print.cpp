#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>

#include <CGAL/Real_timer.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Splitters.h>
#include <CGAL/Kd_tree.h>

using Kernel  = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = typename Kernel::Point_3;

using Traits   = CGAL::Search_traits_3<Kernel>;
using Splitter = CGAL::Sliding_midpoint<Traits>;
using Kd_tree  = CGAL::Kd_tree<Traits, Splitter>;

using Timer = CGAL::Real_timer;

void test_print(const std::string filename) {

  std::cout << std::endl;
  std::cout << "- testing " << filename << " ... " << std::endl;

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
  std::cout << "* building the tree ... " << std::endl;

  // Building the tree.
  Splitter splitter(5); // bucket size = 5
  Kd_tree tree(points.begin(), points.end(), splitter);

  Timer timer;
  timer.start();
  tree.build();
  timer.stop();
  std::cout << "* building done in " << timer.time() << " sec." << std::endl;

  tree.statistics(std::cout);
  // Use this command to print in png:
  // dot -Tpng tree.graphviz > tree.png
  std::ofstream outfile("tree.graphviz");
  tree.write_graphviz(outfile);

  assert(tree.root()->num_items() == points.size());
  assert(tree.root()->num_nodes() == 8);
  assert(tree.root()->depth() == 5);
  std::cout << std::endl;
}

int main() {

  std::cout.precision(20);
  test_print("data/balanced.xyz");
  return EXIT_SUCCESS;
}

#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/MEDIT.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <deque>

typedef CGAL::Simple_cartesian<double>  Kernel;
typedef Kernel::Point_3                 Point_3;


int main() {
  const std::string filename = "./data/polyhedral_complex.mesh";
  std::ifstream input(filename);
  std::vector<Point_3> points;
  std::vector<std::array<int,4>> cells;
  std::vector<int> subdomains;
  std::vector<std::array<int,3>> edges;
  std::deque<CGAL::IO::internal::Corner_with_index<int>> corners;
  bool verbose = false;
  bool success = CGAL::IO::read_MEDIT(input, points, cells, CGAL::parameters::subdomains(std::ref(subdomains)).edges_with_indices(std::ref(edges)).corners_with_indices(std::ref(corners)).verbose(verbose));
  assert(success);

  for(auto c : corners)
    std ::cout << c.v << std::endl;

  std::ostringstream output;
  output.precision(17);
  CGAL::IO::write_MEDIT(output, points, cells, CGAL::parameters::subdomains(std::cref(subdomains)).edges_with_indices(std::cref(edges)).corners_with_indices(std::cref(corners)));
  std::istringstream input2(output.str());
  std::vector<Point_3> points2;
  std::vector<std::array<int,4>> cells2;
  std::vector<int> subdomains2;
  std::vector<std::array<int,3>> edges2;
  std::deque<CGAL::IO::internal::Corner_with_index<int>> corners2;
  success = CGAL::IO::read_MEDIT(input2, points2, cells2, CGAL::parameters::subdomains(std::ref(subdomains2)).edges_with_indices(std::ref(edges2)).corners_with_indices(std::ref(corners2)).verbose(verbose));

  assert(points == points2);
  assert(cells == cells2);
  assert(subdomains == subdomains2);
  assert(corners == corners2);
  assert(edges == edges2);
  assert(success);
  return 0;
}

#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/MEDIT.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

typedef CGAL::Simple_cartesian<double>  Kernel;
typedef Kernel::Point_3                 Point_3;

int main() {
  const std::string filename = CGAL::data_file_path("meshes/elephant.mesh");
  std::ifstream input(filename);
  std::vector<Point_3> points;
  std::vector<int> subdomains;
  std::vector<std::array<int,4>> cells;
  bool verbose = false;
  bool success = CGAL::IO::read_MEDIT(input, points, cells, subdomains, verbose);
  assert(success);

  std::ostringstream output;
  output.precision(17);
  CGAL::IO::write_MEDIT(output, points, cells, subdomains);
  std::istringstream input2(output.str());
  std::vector<Point_3> points2;
  std::vector<std::array<int,4>> cells2;
  std::vector<int> subdomains2;
  success = CGAL::IO::read_MEDIT(input2, points2, cells2, subdomains2, verbose);
  assert(points == points2);
  assert(cells == cells2);
  assert(subdomains == subdomains2);
  assert(success);
  return 0;
}

#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/MEDIT.h>
#include <iostream>
#include <fstream>
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
  return 0;
}

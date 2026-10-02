#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/MEDIT.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <deque>
#include <array>
#include <tuple>

typedef CGAL::Simple_cartesian<double>  Kernel;
typedef Kernel::Point_3                 Point_3;

int main() {
  const std::string filename = "./data/polyhedral_complex.mesh";
  std::ifstream input(filename);
  std::vector<Point_3> points;
  std::vector<std::array<int,4>> cells;
  std::vector<int> subdomains;
  std::vector<std::array<int,3>> edges;
  std::vector<std::tuple<int,int,int,int>> facets;
  std::deque<CGAL::IO::internal::Vertex_with_corner_index<int>> corners;
  bool verbose = false;
  bool success = CGAL::IO::read_MEDIT(input, points, cells, CGAL::parameters::subdomain_indices(std::ref(subdomains))
                                                                             .facets_with_patch_index(std::ref(facets))
                                                                             .edges_with_curve_index(std::ref(edges))
                                                                             .vertices_with_corner_index(std::ref(corners))
                                                                             .verbose(verbose));
  assert(success);

  std::ostringstream output;
  output.precision(17);
  CGAL::IO::write_MEDIT(output, points, cells, CGAL::parameters::subdomain_indices(std::cref(subdomains))
                                                                .facets_with_patch_index(std::cref(facets))
                                                                .edges_with_curve_index(std::cref(edges))
                                                                .vertices_with_corner_index(std::cref(corners)));

  std::istringstream input2(output.str());
  std::vector<Point_3> points2;
  std::vector<std::array<int,4>> cells2;
  std::vector<int> subdomains2;
  std::vector<std::array<int,3>> edges2;
  std::vector<std::tuple<int,int,int,int>> facets2;
  std::deque<CGAL::IO::internal::Vertex_with_corner_index<int>> corners2;
  success = CGAL::IO::read_MEDIT(input2, points2, cells2, CGAL::parameters::subdomain_indices(std::ref(subdomains2))
                                                                           .facets_with_patch_index(std::ref(facets2))
                                                                           .edges_with_curve_index(std::ref(edges2))
                                                                           .vertices_with_corner_index(std::ref(corners2))
                                                                           .verbose(verbose));

  assert(points == points2);
  assert(cells == cells2);
  assert(subdomains == subdomains2);
  assert(corners == corners2);
  assert(edges == edges2);
  assert(facets == facets2);
  assert(success);
  std::cout << "done" << std::endl;
  return 0;
}

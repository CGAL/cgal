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
  std::vector<std::pair<Point_3,int>> points_with_ref;
  std::vector<std::array<int,5>> cells_with_ref;
  std::vector<std::array<int,3>> edges_with_ref;
  std::vector<std::tuple<int,int,int,int>> facets_with_ref;
  std::vector<int> ridges;
  std::deque<int> corners;
  bool success = CGAL::IO::read_MEDIT(input, points_with_ref, CGAL::parameters::tetrahedra(std::ref(cells_with_ref))
                                                                               .triangles(std::ref(facets_with_ref))
                                                                               .edges(std::ref(edges_with_ref))
                                                                               .ridges(std::ref(ridges))
                                                                               .corners(std::ref(corners))
                                                                               .verbose(false));
  assert(success);

  std::ostringstream output;
  output.precision(17);
  CGAL::IO::write_MEDIT(output, points_with_ref, CGAL::parameters::tetrahedra(std::cref(cells_with_ref))
                                                                  .triangles(std::cref(facets_with_ref))
                                                                  .edges(std::cref(edges_with_ref))
                                                                  .ridges(std::cref(ridges))
                                                                  .corners(std::cref(corners)));

  std::istringstream input_bis(output.str());
  std::vector<std::pair<Point_3,int>> points_with_ref_bis;
  std::vector<std::array<int,5>> cells_with_ref_bis;
  std::vector<std::array<int,3>> edges_with_ref_bis;
  std::vector<std::tuple<int,int,int,int>> facets_with_ref_bis;
  std::vector<int> ridges_bis;
  std::deque<int> corners_bis;
  success = CGAL::IO::read_MEDIT(input_bis, points_with_ref_bis, CGAL::parameters::tetrahedra(std::ref(cells_with_ref_bis))
                                                                                   .triangles(std::ref(facets_with_ref_bis))
                                                                                   .edges(std::ref(edges_with_ref_bis))
                                                                                   .ridges(std::ref(ridges_bis))
                                                                                   .corners(std::ref(corners_bis))
                                                                                   .verbose(false));

  assert(points_with_ref== points_with_ref_bis);
  assert(cells_with_ref == cells_with_ref_bis);
  assert(facets_with_ref == facets_with_ref_bis);
  assert(edges_with_ref == edges_with_ref_bis);
  assert(corners == corners_bis);
  assert(ridges == ridges_bis);
  assert(success);
  std::cout << "done" << std::endl;
  return 0;
}

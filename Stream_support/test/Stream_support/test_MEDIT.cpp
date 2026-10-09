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
  std::vector<std::array<int,3>> facets;
  std::vector<std::array<int,2>> edges;
  std::vector<int> ridges;
  std::deque<int> corners;
  std::vector<int> point_refs, cell_refs, facet_refs, edge_refs;

  std::cout << "A  " << & cell_refs << std::endl;

  bool success = CGAL::IO::read_MEDIT(input, points, CGAL::parameters::tetrahedra(std::ref(cells))
                                                                               .triangles(std::ref(facets))
                                                                               .edges(std::ref(edges))
                                                                               .ridges(std::ref(ridges))
                                                                               .corners(std::ref(corners))
                                                                               .tetrahedra_ref(std::ref(cell_refs))
                                                                               .triangles_ref(std::ref(facet_refs))
                                                                               .edges_ref(std::ref(edge_refs))
                                                                               .vertices_ref(std::ref(point_refs))
                                                                               .verbose(false));
  assert(success);

  std::ostringstream output;
  output.precision(17);
  CGAL::IO::write_MEDIT(output, points, CGAL::parameters::tetrahedra(std::cref(cells))
                                                                  .triangles(std::cref(facets))
                                                                  .edges(std::cref(edges))
                                                                  .ridges(std::cref(ridges))
                                                                  .corners(std::cref(corners))
                                                                               .tetrahedra_ref(std::cref(cell_refs))
                                                                               .triangles_ref(std::cref(facet_refs))
                                                                               .edges_ref(std::cref(edge_refs))
                                                                               .vertices_ref(std::cref(point_refs)));

  std::istringstream input_bis(output.str());
  std::vector<Point_3> points_bis;
  std::vector<std::array<int,4>> cells_bis;
  std::vector<std::array<int,3>> facets_bis;
  std::vector<std::array<int,2>> edges_bis;
  std::vector<int> ridges_bis;
  std::deque<int> corners_bis;
  std::vector<int> point_refs_bis, cell_refs_bis, facet_refs_bis, edge_refs_bis;

  success = CGAL::IO::read_MEDIT(input_bis, points_bis, CGAL::parameters::tetrahedra(std::ref(cells_bis))
                                                                                   .triangles(std::ref(facets_bis))
                                                                                   .edges(std::ref(edges_bis))
                                                                                   .ridges(std::ref(ridges_bis))
                                                                                   .corners(std::ref(corners_bis))
                                                                               .tetrahedra_ref(std::ref(cell_refs_bis))
                                                                               .triangles_ref(std::ref(facet_refs_bis))
                                                                               .edges_ref(std::ref(edge_refs_bis))
                                                                               .vertices_ref(std::ref(point_refs_bis))
                                                                                   .verbose(false));

  assert(points== points_bis);
  assert(cells == cells_bis);
  assert(facets == facets_bis);
  assert(edges == edges_bis);
  assert(ridges == ridges_bis);
  assert(corners == corners_bis);
  assert(cell_refs == cell_refs_bis);
  assert(facet_refs == facet_refs_bis);
  assert(edge_refs == edge_refs_bis);
  assert(point_refs == point_refs_bis);
  assert(success);
  std::cout << "done" << std::endl;
  return 0;
}

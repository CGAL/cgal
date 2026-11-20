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

  bool success = CGAL::IO::read_MEDIT(input, points, CGAL::parameters::point_refs(std::ref(point_refs))
                                                                      .tetrahedra(std::ref(cells))
                                                                      .tetrahedron_refs(std::ref(cell_refs))
                                                                      .triangles(std::ref(facets))
                                                                      .triangle_refs(std::ref(facet_refs))
                                                                      .edges(std::ref(edges))
                                                                      .edge_refs(std::ref(edge_refs))
                                                                      .ridges(std::ref(ridges))
                                                                      .corners(std::ref(corners))
                                                                      .verbose(false));
  assert(success);

  std::ostringstream output;
  output.precision(17);
  CGAL::IO::write_MEDIT(output, points, CGAL::parameters::point_refs(std::cref(point_refs))
                                                         .tetrahedra(std::cref(cells))
                                                         .tetrahedron_refs(std::cref(cell_refs))
                                                         .triangles(std::cref(facets))
                                                         .triangle_refs(std::cref(facet_refs))
                                                         .edges(std::cref(edges))
                                                         .edge_refs(std::cref(edge_refs))
                                                         .ridges(std::cref(ridges))
                                                         .corners(std::cref(corners)));

  std::istringstream input_bis(output.str());
  std::vector<Point_3> points_bis;
  std::vector<std::array<int,4>> cells_bis;
  std::vector<std::array<int,3>> facets_bis;
  std::vector<std::array<int,2>> edges_bis;
  std::vector<int> ridges_bis;
  std::deque<int> corners_bis;
  std::vector<int> point_refs_bis, cell_refs_bis, facet_refs_bis, edge_refs_bis;

  success = CGAL::IO::read_MEDIT(input_bis, points_bis, CGAL::parameters::point_refs(std::ref(point_refs_bis))
                                                                         .tetrahedra(std::ref(cells_bis))
                                                                         .tetrahedron_refs(std::ref(cell_refs_bis))
                                                                         .triangles(std::ref(facets_bis))
                                                                         .triangle_refs(std::ref(facet_refs_bis))
                                                                         .edges(std::ref(edges_bis))
                                                                         .edge_refs(std::ref(edge_refs_bis))
                                                                         .ridges(std::ref(ridges_bis))
                                                                         .corners(std::ref(corners_bis))
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



  std::ostringstream output_bis;
  points.resize(4);
  cells.clear();
  std::array<int,4> cell = {0, 1, 2, 3};
  cells.emplace_back(cell);
  CGAL::IO::write_MEDIT(output_bis, points, CGAL::parameters::tetrahedra(std::cref(cells)));

  std::istringstream input_ter(output_bis.str());
  points_bis.clear();
  cells_bis.clear();
  facets_bis.clear();
  edges_bis.clear();
  ridges_bis.clear();
  corners_bis.clear();
  point_refs_bis.clear();
  cell_refs_bis.clear();
  facet_refs_bis.clear();
  edge_refs_bis.clear();

  success = CGAL::IO::read_MEDIT(input_ter, points_bis, CGAL::parameters::point_refs(std::ref(point_refs_bis))
                                                                         .tetrahedra(std::ref(cells_bis))
                                                                         .tetrahedron_refs(std::ref(cell_refs_bis))
                                                                         .triangles(std::ref(facets_bis))
                                                                         .triangle_refs(std::ref(facet_refs_bis))
                                                                         .edges(std::ref(edges_bis))
                                                                         .edge_refs(std::ref(edge_refs_bis))
                                                                         .ridges(std::ref(ridges_bis))
                                                                         .corners(std::ref(corners_bis))
                                                                         .verbose(false));
  assert(points_bis.size() == 4);
  assert(point_refs_bis.size() == 4);
  assert(cells_bis.size() ==  1);
  assert(cell_refs_bis.size() ==  1);
  assert(facets_bis.empty());
  assert(edges_bis.empty());
  assert(ridges_bis.empty());
  assert(corners_bis.empty());
  assert(facet_refs_bis.empty());
  assert(edge_refs_bis.empty());

  std::cout << std::endl << "done" << std::endl;
  return 0;
}

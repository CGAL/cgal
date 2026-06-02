// Reproduction case for issue #9327 (Basic Viewer edge/vertex sizing).
//
// initialize_vertices_and_edges_size() derives the edge width and vertex size
// from the bounding box diagonal. This hexahedron spans ~500 units per side, so
// the diagonal is ~813 and the edges/vertices are drawn hundreds of world units
// thick, obscuring the geometry. Compare with draw_lcc.cpp (same shape at ~5
// units), where the sizes look reasonable. The contrast demonstrates that the
// sizing is tied to object scale instead of a viewport-relative (pixel) metric.

#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/draw_linear_cell_complex.h>

using LCC=CGAL::Linear_cell_complex_for_combinatorial_map<3>;
using Point=LCC::Point;

int main()
{
  LCC lcc;
  lcc.make_hexahedron(Point(0,0,0), Point(500,0,0),
                      Point(500,500,0), Point(0,500,0),
                      Point(0,500,400), Point(0,0,400),
                      Point(500,0,400), Point(500,500,400));
  CGAL::draw(lcc);

  return EXIT_SUCCESS;
}
